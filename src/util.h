#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <float.h>

#include <string>    // string, getline, to_string
#include <sstream>   // istringstream, stringstream, ostringstream
#include <fstream>   // ifstream
#include <iterator>  // iterator, istream_iterator, back_inserter, back_insert_iterator, next
#include <vector>    // vector
#include <iomanip>   // setprecision
#include <iostream>  // fixed, cerr, cout, endl
#include <cmath>     // fabs
#include <algorithm> // find, equal, max_element, min_element, copy
#include <numeric>   // accumulate
#include <chrono>    // high_resolution_clock, time_point, duration
#include <mutex>     // mutex, unique_lock, lock_guard, try_to_lock
#include <set>       // set
#include <cstdlib>   // abort
#include <limits>    // numeric_limits

#define INF std::numeric_limits<double>::infinity()
#define PRECISION 15
#define EPSILON 1.0e-15

#ifndef NDEBUG
#define assert_with_log(expr, msg) assert__(#expr, expr, __FILE__, __LINE__, msg)
#else
#define assert_with_log(expr, msg) {;}
#endif

using namespace std::chrono;

namespace logger {
	void fail(std::string , bool = false );
	void okay(std::string , bool = false );
	void warn(std::string , bool = false );
	void info(std::string , bool = false );
}

struct Point {
	double *values;

	double x() const { return values[0]; }
	double y() const { return values[1]; }
	double z() const { return values[2]; }
};

struct Problem {
	std::set<std::vector<int> > U;
	std::mutex ping_mtx, subsets_mtx, solution_mtx;
	Point *X, *solution;
	double *ref, best, ping_time = 1.0;
	std::vector<long int> nodes;
	long int ntasks;
	int n, k, dim, maximize, cores;
	bool initialized, verbose, kset, pingmode, printpointsonly, devmode;
	high_resolution_clock::time_point t0, ping;
};

template <typename T> inline int gap(T *beg, T *end){ return end - beg; }
template <typename T> inline int argmax(T *beg, int s){ return gap(beg, std::max_element(beg, beg+s)); }
template <typename T> inline int argmin(T *beg, int s){ return gap(beg, std::min_element(beg, beg+s)); }
template <typename T> inline bool contains(std::vector<T> &A, T v){
	return std::find(A.begin(), A.end(), v) != A.end();
}

void search_info(std::string , Problem *, int );
void ping(Problem * , int );
void concurrency_info(Problem * );
Point *new_points(int , int , bool );
bool eq__(double , double );
bool stoi__(int * , std::string );
void assert__(std::string , bool ,  std::string , int , std::string );
void flatten_set(Point * , double * , int , int , int );
int initial_position(Point *, Problem * );
bool is_comment(std::string );
void print_problem(Problem * );
std::vector<int> solution_to_vector(Problem * );
void clear_memory(Problem * );
bool read_args(std::vector<std::string> , std::vector<double> & , std::ifstream & , Problem * );
bool read_input(std::ifstream & , std::vector<double> & , Problem * );

#endif
