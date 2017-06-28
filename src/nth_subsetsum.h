#ifndef NTH_SUBSETSUM_H
#define NTH_SUBSETSUM_H

#include <vector>    // vector
#include <queue>     // priority_queue
#include <algorithm> // sort

struct DataPoint {
	double val;
	int data_pos;
	int sorted_pos;

	bool operator< (const DataPoint &r) const {
		return val > r.val;
	}
};

struct HeapEntry {
	std::vector<DataPoint> subset;
	double sum = 0.0f; // ++

	bool operator< (const HeapEntry &t) const {
		return sum < t.sum;
	}
};

std::vector<HeapEntry> nth_subsetsum(double * , const int , const int , const unsigned int );

#endif
