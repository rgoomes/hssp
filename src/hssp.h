#ifndef HSSP_H
#define HSSP_H

#include "util.h"
#include "hypervolume.h"
#include "nth_subsetsum.h"

#include <algorithm>          // copy, nth_element, iter_swap, sort, find_if, equal
#include <functional>         // greater, less
#include <numeric>            // accumulate
#include <chrono>             // high_resolution_clock, now
#include <sstream>            // stringstream
#include <cmath>              // exp
#include <vector>             // vector
#include <string>             // string
#include <fstream>            // ifstream
#include <mutex>              // mutex, unique_lock, lock_guard
#include <condition_variable> // condition_variable
#include <thread>             // thread
#include <queue>              // queue, priority_queue
#include <atomic>             // atomic

extern "C"
{
#include "HVC/hvc-class.h"
}

using namespace std::chrono;

struct Task {
	Point *S, *subset;
	double *C, *Ce, *Cr, hv, ubound1;
	int cur_pos, subset_size;
	bool is_new;

	bool operator< (const Task &t) const {
		return hv < t.hv;
	}
};

class Pool {
	Problem *P;
	std::mutex tasks_mtx;
	std::condition_variable cv;
	std::priority_queue<Task> tasks;
	std::vector<std::thread> threads;

	int workers = 0;
	std::atomic<bool> stop;

public:
	Pool(Problem * );
	void join();
	void terminate();
	void work(int );
	int working();
	void schedule(Point *, Point *, double *, double *, double *, double , double , int , int , bool );
	Task snapshot(Point *, Point *, double *, double *, double *, double , double , int , int , bool );
};

double hssp(std::vector<std::string> , std::vector<int> & , long int & );
void branch(Point * , Point * , Point * , Point * , Point * , double * , double * , double * , double * , bool , double , double , int , Pool * , hvc_s * , hvc_s * , Problem * );

#endif
