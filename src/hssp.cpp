#include "hssp.h"

Pool::Pool(Problem *P_) : stop(false) {
	this->P = P_;

	for(int t = 0; t < P->cores; ++t)
		threads.push_back(std::thread(&Pool::work, this, t));
}

void Pool::join(){
	for(std::thread &t : threads)
		t.join();
}

void Pool::terminate(){
	stop = true;
	cv.notify_all();
}

int Pool::working(){
	assert_with_log(workers >= 0 && workers <= P->cores, "number of workers out of bounds");
	return workers;
}

Task Pool::snapshot(Point *S, Point *subset, double *C, double *Ce, double *Cr, double hv, double ubound1, int cur_pos, int subset_size, bool is_new){
	Task task {
		new Point[P->n], new Point[P->k],
		new double[P->n], new double[P->n], new double[P->n], hv, ubound1,
		cur_pos, subset_size,
		is_new
	};

	std::copy(S+cur_pos, S+P->n, task.S+cur_pos);
	std::copy(subset, subset+subset_size, task.subset);

	std::copy(C+cur_pos, C+P->n, task.C+cur_pos);
	std::copy(Ce+cur_pos, Ce+P->n, task.Ce+cur_pos);
	std::copy(Cr+cur_pos, Cr+P->n, task.Cr+cur_pos);

	// NOTE: the return value should be optimized by the compiler
	return task;
}

void Pool::schedule(Point *S, Point *subset, double *C, double *Ce, double *Cr, double hv, double ubound1, int cur_pos, int subset_size, bool is_new){
	const Task task = snapshot(S, subset, C, Ce, Cr, hv, ubound1, cur_pos, subset_size, is_new);

	std::unique_lock<std::mutex> ul(tasks_mtx);
	tasks.push(task);
	cv.notify_one();
	ul.unlock();
}

void Pool::work(int id){
	bool dismiss = true;
	Point aux1[P->n];
	double aux2[P->n];

	while(true){
		bool idle = true;
		std::unique_lock<std::mutex> ul(tasks_mtx);

		if(dismiss)
			dismiss = false;
		else {
			// if the queue has tasks do not decrease the number of workers
			if(!tasks.empty())
				idle = false;
			else if(--workers == 0 /* && tasks.empty() */ )
				terminate();
		}

		while(!stop && tasks.empty()){
			assert_with_log(idle == true, "");
			cv.wait(ul);
		}

		if(stop)
			return;

		// if the queue had tasks do not increase the number of workers
		workers += idle;

		Task task = tasks.top();
		tasks.pop();
		ul.unlock();

		hvc_s *hvcs = build(task.S, task.subset, aux1, task.subset_size, 0, P->n, P->dim, P->ref, false);
		hvc_s *full = build(task.S+task.cur_pos, task.subset, aux1, task.subset_size, P->n - task.cur_pos, P->n, P->dim, P->ref, true);

		branch(task.S+task.cur_pos, task.subset+task.subset_size, task.subset, aux1, task.S, task.C, task.Ce, task.Cr, aux2, task.is_new, task.hv, task.ubound1, id, this, hvcs, full, P);

		delete [] task.S; delete [] task.subset;
		delete [] task.C; delete [] task.Ce; delete [] task.Cr;

		if(P->dim == 3){
			dealloc(hvcs);
			dealloc(full);
		}

		P->ntasks++;
	}
}

double contribution(Point *q, Point *S, Point *aux1, int s, double hv, Problem *P){
	aux1[s] = *q;
	std::copy(S, S+s, aux1);
	return hypervolume(aux1, s+1, P->dim, P->ref) - hv;
}

#ifndef NDEBUG
// TODO: move test function to test.cpp
void test(Point *S, Point *cur, Point *subset, Point *aux1, double *C, double &hv, double &ubound1, int subset_size, int points_left, int cur_pos, bool is_new, Problem *P){
	if(!P->devmode)
		return;

	if(is_new){
		std::vector<int> uniq;

		for(Point *it = subset; it != subset+subset_size; ++it){
			const int pos = initial_position(it, P);
			uniq.push_back(pos);
		}

		// test branching
		std::sort(uniq.begin(), uniq.end());
		std::lock_guard<std::mutex> lg(P->subsets_mtx);
		auto ret = P->U.insert(uniq);
		assert_with_log(ret.second == true, "same subset tested twice in branch");
	}

	// test hypervolumes
	std::copy(subset, subset+subset_size, aux1);
	assert_with_log(eq__(hv, hypervolume(aux1, subset_size, P->dim, P->ref)), "wrong hypervolume indicator");

	const double cache = ubound1;
	std::copy(subset, subset+subset_size, aux1);
	std::copy(cur, cur+points_left, aux1+subset_size);
	ubound1 = hypervolume(aux1, subset_size+points_left, P->dim, P->ref);
	assert_with_log(eq__(cache, ubound1), "wrong hypervolume indicator (faulty all_contributions or problem in the incremental calculation of this bound)");

	// test subset
	for(Point *it = cur; it != S + P->n; ++it){
		Point *e = std::find_if(subset, subset+subset_size, [it, P](const Point &p) -> bool {
			return std::equal(it->values, it->values + P->dim, p.values, eq__);
		});

		assert_with_log(subset+subset_size == e, "subset contains unassigned points");
	}

	// test hypervolume contributions
	for(int i = cur_pos; i < P->n; ++i){
		const double c = contribution(S+i, subset, aux1, subset_size, hv, P);

		assert_with_log(eq__(C[i], c), "incorrect hypervolume contribution");
		assert_with_log(C[cur_pos] > C[i] || eq__(C[cur_pos], C[i]), "maximal contributer in wrong position");
	}

	// test exclusive hypervolume contributions
	// TODO: write test to verify the exclusive hypervolume contributions
}
#endif

// check if the hypervolume of the first set with the second set reaches best
double bound1(Point *cur, Point *subset, Point *aux1, double *Cr, double *Ce, double *aux2, double &ubound1, int subset_size, int points_left, int cur_pos, bool is_new, hvc_s *full, Problem *P){
	const int missing = P->k - subset_size;
	const int excess = points_left - missing;
	assert_with_log(excess >= 0, "invalid negative excess");

	// NOTE: disable the incremental calculation of this upper bound
	if(P->dim != 3 && P->dim != 4){
		std::copy(subset, subset+subset_size, aux1);
		std::copy(cur, cur+points_left, aux1+subset_size);
		ubound1 = hypervolume(aux1, subset_size+points_left, P->dim, P->ref);
	}

	if(ubound1 <= P->best)
		return ubound1;

	// NOTE: currently all_contributions algorithm is only available for d=3 and d=4
	if(P->dim == 3 || P->dim == 4){
		if(!is_new){
			const int all_size = subset_size + points_left;
			double all[all_size], *tmp = all+subset_size;

			if(P->dim == 3){
				double *contribs = getContributions(full);
				int *selected = getSelected(full);
				assert_with_log(all_size == getSize(full), "");

				for(int i = 0; i < all_size; ++i)
					if(!selected[i] /* if not in the subset */ )
						*tmp++ = contribs[i];
			}

			if(P->dim != 3){
				std::copy(subset, subset+subset_size, aux1);
				std::copy(cur, cur+points_left, aux1+subset_size);
				all_contributions(aux1, all_size, P->dim, P->ref, all);
			}

			std::copy(all+subset_size, all+all_size, Ce+cur_pos);
		}

		std::copy(Ce+cur_pos, Ce+P->n, aux2);
	}
	else {
		// for d>4 fallback for the exclusive hypervolume contributions at the root
		std::copy(Cr+cur_pos, Cr+P->n, aux2);
	}

	std::nth_element(aux2, aux2+excess, aux2+points_left, std::less<double>());
	const double exclusive_excess = std::accumulate(aux2, aux2+excess, 0.0);
	return ubound1 - exclusive_excess;
}

double bound2(double *C, double *aux2, double &hv, int subset_size, int points_left, int cur_pos, Problem *P){
	const int missing = P->k - subset_size;

	std::copy(C+cur_pos, C+P->n, aux2);
	std::nth_element(aux2, aux2+missing, aux2+points_left, std::greater<double>());
	return hv + std::accumulate(aux2, aux2+missing, 0.0);
}

double bound3(Point *S, Point *subset, Point *aux1, double *C, double *aux2, double hv, int subset_size, int cur_pos, int points_left, Problem *P){
	Point aux3[P->n];
	bool selected[P->n];
	std::fill(selected, selected + P->n, false);

	std::copy(C+cur_pos, C+P->n, aux2+cur_pos);
	std::copy(subset, subset+subset_size, aux3);

	for(int s = subset_size; s < P->k; ++s){
		const int p = cur_pos + argmax(aux2+cur_pos, points_left);
		selected[p] = true;
		aux3[s] = S[p];
		hv += aux2[p];

		for(int q = cur_pos; q < P->n; ++q)
			aux2[q] = !selected[q] ? contribution(S+q, aux3, aux1, s+1, hv, P) : 0;
	}

	const double approximation = 1 - (1.0 / std::exp(1.0)); // ~0.6321
	return hv / approximation;
}

double bound2_extension(Point *S, Point *subset, Point *aux1, double *C, double &hv, int subset_size, int points_left, int cur_pos, Problem *P){
	const int missing = P->k - subset_size;

	// NOTE: this is disabled for now. two problems exist: nth_subsetsum is too slow and it is also
	// very demanding to calculate this many hypervolumes
	std::vector<HeapEntry> solutions = nth_subsetsum(C+cur_pos, points_left, missing, points_left);

	if(solutions.size() && hv + solutions.back().sum >= P->best)
		return INF;

	for(HeapEntry &h : solutions){
		std::copy(subset, subset+subset_size, aux1);

		for(int i = 0; i < missing; ++i){
			const int pos = h.subset[i].data_pos;
			aux1[subset_size+i] = S[cur_pos+pos];
		}

		const double this_hv = hypervolume(aux1, P->k, P->dim, P->ref);

		// make sure that a new best solution is not lost
		if(this_hv > P->best)
			return INF;
		if(hv + h.sum <= P->best)
			return -INF;
	}

	return INF;
}

void update(Point *next, Point *subset, Point *aux1, double *C, double hv, int subset_size, int cur_pos, int next_pos, int points_left, hvc_s *hvcs, Problem *P){
	if(P->dim == 3)
		update_contributions(next, points_left, hvcs, C+next_pos);

	if(P->dim != 3){
		double *it = C+next_pos;

		for(Point *p = next; p != next + points_left; ++p)
			*it++ = contribution(p, subset, aux1, subset_size, hv + C[cur_pos], P);
	}
}

// restores the exclusive hypervolume contribution of the current point. this is needed for the
// incremental update of bound1. also we ignore the order for other points since we do not need
// to know what exclusive hypervolume contribution a point has, but the all values themselves
void restore_current(double *Ce, double *cur_point, int cur_pos, hvc_s *full, Problem *P){
	if(P->dim != 3)
		return;

	const int pos = findPointPos(full, cur_point);
	double *Ce_ = getContributions(full);
	assert_with_log(pos >= 0, "failed to find current point");

	// NOTE: finding the position of the contribution by value should be avoided
	for(int p = cur_pos; p < P->n; ++p){
		if(eq__(Ce[p], Ce_[pos])){
			std::iter_swap(Ce+cur_pos, Ce+p);
			return;
		}
	}

	assert_with_log(false, "failed to find the contribution");
}

void swap_next(Point *S, double *C, double *Ce, double *Cr, int pos1, int pos2, Problem *P){
	if(pos1 >= P->n || pos2 >= P->n)
		return;

	// swap points
	std::iter_swap(S+pos1, S+pos2);
	// swap hypervolume contributions
	std::iter_swap(C+pos1, C+pos2);
	// swap exclusive hypervolume contributions
	std::iter_swap(Ce+pos1, Ce+pos2);
	// swap root exclusive hypervolume contributions
	std::iter_swap(Cr+pos1, Cr+pos2);
}

// TODO: later use the Task struct as the node. the code will need to be fully refactored
void branch(Point *cur, Point *end, Point *subset, Point *aux1, Point *S, double *C, double *Ce, double *Cr, double *aux2, bool is_new, double hv, double ubound1, int id, Pool *pool, hvc_s *hvcs, hvc_s *full, Problem *P){
	const int subset_size = gap(subset, end);
	const int cur_pos     = gap(S, cur);
	const int points_left = P->n - cur_pos;
	const int next_pos    = cur_pos+1;
	P->nodes[id]++;
	ping(P, id);

#ifndef NDEBUG
	test(S, cur, subset, aux1, C, hv, ubound1, subset_size, points_left, cur_pos, is_new, P);
#endif

	if(subset_size == P->k){
		// guard invalidating new best solutions when two or more are found at the same time
		std::lock_guard<std::mutex> lg(P->solution_mtx);

		if(hv > P->best){
			P->best = hv;
			search_info("new", P, id);
			std::copy(subset, end, P->solution);
		}

		return;
	}

	if(!points_left)
		return;
	if(subset_size + points_left < P->k)
		return;

	if(bound2(C, aux2, hv, subset_size, points_left, cur_pos, P) <= P->best)
		return;
	if(bound1(cur, subset, aux1, Cr, Ce, aux2, ubound1, subset_size, points_left, cur_pos, is_new, full, P) <= P->best)
		return;
	//if(bound2_extension(S, subset, aux1, C, hv, subset_size, points_left, cur_pos, P) <= P->best)
	//	return;
	/*if(bound3(S, subset, aux1, C, aux2, hv, subset_size, cur_pos, points_left, P) <= P->best){
		logger::info("branch pruned by bound3");
		return;
	}*/

	// if(d <= 4 && subset_size + points_left - 1 == P->k)
	//	logger::info("base case found: remove the point with the least exclusive hypervolume contribution ");

	double Cb[P->n];
	std::copy(C+next_pos, C+P->n, Cb+next_pos);

	double *cur_point = cur->values;
	restore_current(Ce, cur_point, cur_pos, full, P);

	if(P->dim == 3 && hvcs)
		addPoint(hvcs, cur_point, 0);
	if(P->dim == 3 && full)
		enablePoint(full, cur_point);

	// accept the current point
	*end = *cur;

	// update the hypervolume contributions given the newly accepted point
	update(cur+1, subset, aux1, C, hv, subset_size+1, cur_pos, next_pos, points_left-1, hvcs, P);

	const int accept_pos = next_pos + argmax(C+next_pos, points_left - 1);
	const int ignore_pos = next_pos + argmax(Cb+next_pos, points_left - 1);

	if(P->cores > 1 && pool && pool->working() < P->cores && P->best > 0.0 && subset_size > 0){
		swap_next(S, C, Ce, Cr, next_pos, accept_pos, P);
		pool->schedule(S, subset, C, Ce, Cr, hv + C[cur_pos], ubound1, next_pos, subset_size+1, true);
		swap_next(S, C, Ce, Cr, next_pos, accept_pos, P);
	}
	else {
		swap_next(S, C, Ce, Cr, next_pos, accept_pos, P);
		branch(cur+1, end+1, subset, aux1, S, C, Ce, Cr, aux2, true, hv + C[cur_pos], ubound1, id, pool, hvcs, full, P);
		swap_next(S, C, Ce, Cr, next_pos, accept_pos, P);
	}

	if(P->dim == 3 && hvcs)
		removePointAt(hvcs, subset_size, 0); // removePoint(hvcs, cur_point, 0);
	if(P->dim == 3 && full)
		removePoint(full, cur_point, 0); // disablePoint(full, cur_point); is superfluous

	if(P->cores > 1 && pool && pool->working() < P->cores && P->best > 0.0 && subset_size > 0){
		swap_next(S, Cb, Ce, Cr, next_pos, ignore_pos, P);
		pool->schedule(S, subset, Cb, Ce, Cr, hv, ubound1 - Ce[cur_pos], next_pos, subset_size, false);
		swap_next(S, Cb, Ce, Cr, next_pos, ignore_pos, P);
	}
	else {
		swap_next(S, Cb, Ce, Cr, next_pos, ignore_pos, P);
		branch(cur+1, end+0, subset, aux1, S, Cb, Ce, Cr, aux2, false, hv, ubound1 - Ce[cur_pos], id, pool, hvcs, full, P);
		swap_next(S, Cb, Ce, Cr, next_pos, ignore_pos, P);
	}

	if(P->dim == 3 && full)
		addPoint(full, cur_point, 0);
}

// decremental greedy
void heuristic(Point *S, Point *aux1, double *aux2, double full, Problem *P){
	if(P->dim != 3 && P->dim != 4)
		return;

	std::copy(S, S + P->n, aux1);

	for(int s = P->n; s != P->k; --s){
		all_contributions(aux1, s, P->dim, P->ref, aux2);
		const int minp = argmin(aux2, s);
		std::iter_swap(aux1+minp, aux1+s-1);
		full -= aux2[minp];
	}

	if(full > P->best){
		P->best = full;
		std::copy(aux1, aux1 + P->k, P->solution);
		search_info("init", P, 0);
	}
}

void compute_root_contributions(Point *S, Point *subset, Point *aux1, double *C, double *Ce, double *Cr, Problem *P){
	if(P->dim == 3){
		get_contributions(S, subset, 0, P->n, P->dim, P->ref, C);
		all_contributions(S, P->n, P->dim, P->ref, Ce);
	}

	if(P->dim != 3){
		std::copy(S, S + P->n, aux1);
		const double hs = hypervolume(aux1, P->n, P->dim, P->ref);

		for(Point *p = S; p != S + P->n; ++p)
			*C++ = contribution(p, subset, aux1, 0, 0, P);

		double *it = Ce;
		for(int i = 0; i < P->n; i++){
			std::copy(S, S + P->n, aux1);
			std::iter_swap(aux1+0, aux1+i);
			*it++ = hs - hypervolume(aux1+1, P->n-1, P->dim, P->ref);
		}
	}

	std::copy(Ce, Ce+P->n, Cr);
}

double root(Problem *P){
	Point S[P->n], subset[P->k], aux1[P->n];
	double C[P->n], Ce[P->n], Cr[P->n], aux2[P->n];

	std::copy(P->X, P->X + P->n, S);
	compute_root_contributions(S, subset, aux1, C, Ce, Cr, P);

	if(P->k == 1){
		const int maxp = argmax(C, P->n);
		P->solution[0] = S[maxp];
		P->best = C[maxp];
	}
	else if(P->k == P->n){
		std::copy(S, S + P->n, P->solution);
		P->best = hypervolume(S, P->n, P->dim, P->ref);
	}
	else if(P->k == P->n - 1){
		const int minp = argmin(Ce, P->n);
		std::copy(S, S + P->n, P->solution);
		std::iter_swap(P->solution+minp, P->solution + P->n - 1);
		P->best = hypervolume(S, P->n, P->dim, P->ref) - Ce[minp];
	}
	else {
		std::copy(S, S + P->n, aux1);
		const int r = argmax(C, P->n);
		const double hs = hypervolume(aux1, P->n, P->dim, P->ref);

		heuristic(S, aux1, aux2, hs, P);
		swap_next(S, C, Ce, Cr, 0, r, P);

		if(P->cores == 1){
			hvc_s *hvcs = build(S, subset, aux1, 0, 0, P->n, P->dim, P->ref, false);
			hvc_s *full = build(S, subset, aux1, 0, P->n, P->n, P->dim, P->ref, true);
			branch(S, subset, subset, aux1, S, C, Ce, Cr, aux2, false, 0.0, hs, 0, nullptr, hvcs, full, P);

			if(P->dim == 3 && hvcs)
				dealloc(hvcs);
			if(P->dim == 3 && full)
				dealloc(full);
		}
		else {
			Pool pool(P);
			pool.schedule(S, subset, C, Ce, Cr, 0.0, hs, 0, 0, false);
			pool.join();
		}

		swap_next(S, C, Ce, Cr, 0, r, P);
	}

	search_info("end", P, 0);
	concurrency_info(P);
	return P->best;
}

double hssp(std::vector<std::string> args, std::vector<int> &sol, long int &nodes){
	Problem P {};
	std::ifstream file;
	std::vector<double> ref;

	P.t0 = high_resolution_clock::now();
	P.ping = P.t0;
	P.best = -INF;
	P.cores = 1u;

	if(!read_args(args, ref, file, &P))
		return -1;
	if(!read_input(file, ref, &P))
		return -1;
	if(P.printpointsonly){
		print_problem(&P);
		clear_memory(&P);
		return -1;
	}

	P.nodes.resize(P.cores);
	const double volume = root(&P);
	sol = solution_to_vector(&P);
	nodes = std::accumulate(P.nodes.begin(), P.nodes.end(), 0);
	clear_memory(&P);
	return volume;
}
