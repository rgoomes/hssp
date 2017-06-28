#include "hypervolume.h"

double hv2D(Point *S, const int n){
	std::sort(S, S+n, [](const Point &a, const Point &b){ return a.y() > b.y(); });

	double volume = 0.0;
	Point * const last_point = S+n-1;

	for(Point *it = S; it != last_point; ++it)
		volume += (it->y() - std::next(it)->y()) * it->x();

	return volume + last_point->x() * last_point->y();
}

double hv3D(Point *S, const int n){
	std::sort(S, S+n, [](const Point &a, const Point &b){ return a.z() > b.z(); });

	Tree T {{S->x() , S->y()}, {0, INF}};
	double area = S->x() * S->y(), volume = 0.0;
	Point * const end = S+n, * const last_point = end-1;

	for(Point *it = S+1; it != end; ++it){
		Point2D plane {it->x(), it->y()};
		Tree::iterator ub = T.upper_bound(plane), pred = std::prev(ub);
		double ub_y = ub != T.end() ? ub->y : 0.0;

		volume += area * (std::prev(it)->z() - it->z());

		while(it->y() >= pred->y){
			Tree::iterator dm = pred--;
			area -= (dm->x - pred->x) * (dm->y - ub_y);
			T.erase(dm);
		}

		if(it != last_point)
			T.emplace_hint(pred, plane);

		area += (it->x() - pred->x) * (it->y() - ub_y);
	}

	return volume + area * last_point->z();
}

double hvND(Point *S, const int n, const int dim, const double *ref){
	double data[n*dim];
	flatten_set(S, data, n, 1, dim);
	return fpli_hv(data, dim, n, ref);
}

double hypervolume(Point *S, int n, const int dim, const double *ref){
	if(n <= 0)
		return 0.0;

	switch(dim){
		// NOTE: hv2D and hv3D are only for maximization
		//case 2:
		//	return hv2D(S, n);
		//case 3:
		//	return hv3D(S, n);
		default:
			return hvND(S, n, dim, ref);
	}
}

void all_contributions(Point *S, const int n, const int dim, double *ref, double *contributions){
	double data[n*dim];
	flatten_set(S, data, n, 1, dim);
	hvc(data, dim, n, ref, contributions, 0);
}

void get_contributions(Point *S, Point *subset, int subset_size, int n, int dim, double *ref, double *contributions){
	double data[subset_size*dim];
	flatten_set(subset, data, subset_size, 1, dim);
	hvc_s *hvcs = init(data, dim, subset_size, n, ref);

	for(Point *p = S; p != S+n; ++p)
		*contributions++ = oneContribution(hvcs, p->values);

	dealloc(hvcs);
}

void update_contributions(Point *S, int n, hvc_s *hvcs, double *contributions){
	for(Point *p = S; p != S+n; ++p)
		*contributions++ = oneContribution(hvcs, p->values);
}

hvc_s *build(Point *S, Point *subset, Point *aux1, int subset_size, int points_left, int n, int dim, double *ref, bool enable_subset){
	if(dim != 3)
		return nullptr;

	const int all_size = subset_size + points_left;
	std::copy(subset, subset+subset_size, aux1);
	std::copy(S, S+points_left, aux1+subset_size);

	double data[all_size*dim];
	flatten_set(aux1, data, all_size, 1, dim);
	hvc_s *hvcs = init(data, dim, all_size, n, ref);

	if(!enable_subset)
		return hvcs;
	for(int p = 0; p < subset_size; ++p)
		enablePoint(hvcs, &data[p*dim]);

	return hvcs;
}
