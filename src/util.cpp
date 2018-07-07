#include "util.h"

const std::string tryhelp_str = "Try 'hssp --help' for more information.\n";
const std::string usage_str =
	"usage:  hssp [OPTION...] [FILE]\n\n"
	"      --help           display this help and exit\n"
	"  -a  --maximize       use reference point for maximization problem\n"
	"  -r, --reference R    use R as the reference point (default reference is the origin)\n"
	"  -k, --subsetsize K   select K points (default subsetsize is half the number of points)\n"
	"  -j, --concurrency T  use T parallel workers (default concurrency is 1)\n"
	"  -v, --verbose        displays extra information\n";

std::mutex print_mtx__;

namespace logger {
	void fail(std::string msg, bool to_console){
		std::lock_guard<std::mutex> lguard(print_mtx__);
		if(to_console)
			std::cout << "[ FAIL ] " << msg << std::endl;
		else
			std::cerr << "[ \033[1;31m" << "FAIL" << "\033[0m ] " << msg << std::endl;
	}

	void okay(std::string msg, bool to_console){
		std::lock_guard<std::mutex> lguard(print_mtx__);
		if(to_console)
			std::cout << "[  OK  ] " << msg << std::endl;
		else
			std::cerr << "[ \033[1;32m" << " OK " << "\033[0m ] " << msg << std::endl;
	}

	void warn(std::string msg, bool to_console){
		std::lock_guard<std::mutex> lguard(print_mtx__);
		if(to_console)
			std::cout << "[ WARN ] " << msg << std::endl;
		else
			std::cerr << "[ \033[1;33m" << "WARN" << "\033[0m ] " << msg << std::endl;
	}

	void info(std::string msg, bool to_console){
		std::lock_guard<std::mutex> lguard(print_mtx__);
		if(to_console)
			std::cout << "[ INFO ] " << msg << std::endl;
		else
			std::cerr << "[ \033[1;37m" << "INFO" << "\033[0m ] " << msg << std::endl;
	}
}

void search_info(std::string msg, Problem *P, int id){
	if(!P->verbose)
		return;

	const duration<float> fs = high_resolution_clock::now() - P->t0;
	const long int total_nodes = std::accumulate(P->nodes.begin(), P->nodes.end(), 0);

	std::ostringstream lbound;
	lbound << P->best;

	std::stringstream ss;
	ss << "lbound " << std::fixed << std::setprecision(PRECISION) << lbound.str() << " nodes " << total_nodes << " time " << std::fixed << std::setprecision(6) << fs.count() << "s" << " nps " << int(total_nodes / fs.count()) << " thread " << id << " type " << msg;
	logger::info(ss.str());

	if(P->pingmode){
		P->ping = high_resolution_clock::now();

		if(msg == "ping")
			P->ping_time *= 2.0;
		if(msg == "new")
			P->ping_time = 1.0;
	}
}

void ping(Problem *P, int id){
	if(!P->pingmode)
		return;

	const duration<float> fs = high_resolution_clock::now() - P->ping;
	if(fs.count() > P->ping_time){
		// lock multiple threads attempting to ping at the same time. this avoids duplicated ping logs
		// when concurrency is enabled
		std::unique_lock<std::mutex> ulock(P->ping_mtx, std::try_to_lock);
		if(ulock.owns_lock())
			search_info("ping", P, id);

	}
}

void concurrency_info(Problem *P){
	if(!P->verbose)
		return;

	for(int t = 0; t < P->cores; ++t)
		logger::info("thread " + std::to_string(t) + " nodes " + std::to_string(P->nodes[t]));

	logger::info("thread tasks " + std::to_string(P->ntasks));
}

Point *new_points(int size, int dim, bool alloc_coordinates){
	Point *S = new Point[size];

	for(Point *it = S; alloc_coordinates && it != S+size; ++it)
		it->values = new double[dim];

	return S;
}

bool eq__(double a, double b){
	return std::fabs(a - b) <= EPSILON;
}

// sets dest variable. returns true on failure
bool stoi__(int *dest, std::string num){
	try {
		*dest = stoi(num);
		return false;
	} catch(...){
		return true;
	}
}

void assert__(std::string expr_str, bool expr, std::string file, int line, std::string msg){
	if(!expr){
		std::stringstream ss;
		ss << "assertion '" << expr_str << "' failed in " << file << " line " << line << ": " << msg;
		logger::fail(ss.str());
		std::abort();
	}
}

void flatten_set(Point *S, double *to, int n, int direction, int dim){
	for(Point *it = S; it != S+n; ++it)
		for(int d = 0; d < dim; ++d)
			*to++ = direction * it->values[d];
}

int initial_position(Point *point, Problem *P){
	assert_with_log(P->initialized, "problem not yet initialized");

	for(Point *it = P->X; it != P->X + P->n; ++it)
		if(std::equal(it->values, it->values + P->dim, point->values, eq__))
			return gap(P->X, it);

	return -1;
}

bool is_comment(std::string str){
	size_t starts_with = str.find_first_not_of(" \t\r\n");
	return str[starts_with] == '#';
}

void print_problem(Problem *P){
	const std::string sep = " ";

	std::cout << "#" << std::endl;
	for(Point *it = P->X; it != P->X + P->n; ++it){
		for(int d = 0; d < P->dim; ++d)
			std::cout << std::fixed << std::setprecision(PRECISION) << it->values[d] << (d+1 < P->dim ? sep : "");

		std::cout << std::endl;
	}
	std::cout << "#" << std::endl;
}

std::vector<int> solution_to_vector(Problem *P){
	std::vector<int> solution(P->n, 0);

	for(Point *it = P->solution; it != P->solution + P->k; ++it){
		const int pos = initial_position(it, P);
		assert_with_log(pos >= 0, "point is not in the original set");
		solution[pos] = 1;
	}

	return solution;
}

void clear_memory(Problem *P){
	for(Point *it = P->X; it && it != P->X + P->n; ++it)
		if(it->values)
			delete [] it->values;

	if(P->X)
		delete [] P->X;
	if(P->ref)
		delete [] P->ref;
	if(P->solution)
		delete [] P->solution;
}

bool ERROR(Problem *P, std::string error_msg = ""){
	printf("hssp: %s", error_msg.c_str());
	clear_memory(P);
	return false;
}

bool read_args(std::vector<std::string> args, std::vector<double> &ref, std::ifstream &file, Problem *P){
	if(contains(args, std::string("--help") ))
		return ERROR(P, usage_str + "\n");

	for(std::vector<std::string>::iterator it = args.begin(); it != args.end(); ++it){
		if(*it == "-r" || *it == "--reference"){
			if(std::next(it) != args.end()){
				std::istringstream iss(*++it);
				ref.assign(std::istream_iterator<double>(iss), std::istream_iterator<double>());
			}
		}
		else if(*it == "-k" || *it == "subsetsize"){
			if(std::next(it) != args.end() && !stoi__(&P->k, *std::next(it))){
				P->kset = true;
				++it;
			}
		}
		else if(*it == "-j" || *it == "--concurrency"){
			if(std::next(it) != args.end() && !stoi__(&P->cores, *std::next(it)))
				++it;
		}
		else if(*it == "-a" || *it == "--maximize")
			P->maximize = true;
		else if(*it == "-v" || *it == "--verbose")
			P->verbose = true;
		else if(*it == "--ping")
			P->pingmode = true;
		else if(*it == "--points")
			P->printpointsonly = true;
		else if(*it == "--devmode")
			P->devmode = true;
		else {
			file.open(*it, std::ifstream::in);

			// is this an input file?
			if(!file.good() && *it != "--help")
				return ERROR(P, "unrecognized option '" + *it + "'\n" + tryhelp_str);
		}
	}

	if(!file.good())
		return ERROR(P, "no file input specified\n" + tryhelp_str);
	if(P->cores <= 0)
		return ERROR(P, "invalid concurrency value \'" + std::to_string(P->cores) + "\'\n" + tryhelp_str);

	return true; // no error
}

bool read_input(std::ifstream &file, std::vector<double> &ref, Problem *P){
	std::string line;
	std::vector<double> values;
	std::back_insert_iterator<std::vector<double> > it = std::back_inserter(values);

	while(std::getline(file, line)){
		if(is_comment(line))
			continue;

		std::istringstream iss(line);
		size_t vector_size = values.size();
		it = std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), it);
		size_t insert_dims = values.size() - vector_size;

		if(P->dim && insert_dims != (size_t) P->dim)
			return ERROR(P, "input data points must have the same dimension\n");
		if(ref.size() && insert_dims != ref.size())
			return ERROR(P, "invalid input data point or reference point dimension\n");

		P->dim = (int) insert_dims;
	}

	P->n = P->dim == 0 ? 0 : (int) values.size() / P->dim;

	if(P->n <= 0)
		return ERROR(P, "no input data given\n" + tryhelp_str);
	if(P->k > P->n || P->k < 0 || (P->k == 0 && P->kset))
		return ERROR(P, "invalid subset size \'" + std::to_string(P->k) + "\'\n" + tryhelp_str);

	if(P->k == 0 && !P->kset)
		P->k = P->n/2;

	ref.resize(P->dim, 0.0);
	P->ref = new double[P->dim]();
	P->X = new_points(P->n, P->dim, true);

	for(int i = 0; i < P->n; ++i){
		for(int j = 0; j < P->dim; ++j){
			const double value = values[i*P->dim+j];
			if((P->maximize && (ref[j] >= value)) || (!P->maximize && (ref[j] <= value)))
				return ERROR(P, "reference point must dominate input data points\n");

			// normalize data to the origin reference point and transform it to minimization problem
			P->X[i].values[j] = (value - ref[j]) * (!P->maximize - P->maximize);
		}
	}

	P->solution = new_points(P->n, P->dim, false);
	P->initialized = true;

	file.close();
	return true;
}
