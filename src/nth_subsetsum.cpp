#include "nth_subsetsum.h"

unsigned int choose(int n, int k){
	return k == 0 ? 1 : n*choose(n - 1, k - 1) / k;
}

bool ascending_index__(const DataPoint &a, const DataPoint &b){
	return a.sorted_pos < b.sorted_pos;
}

// nth_subsetsum: given an array of size n find the m greatest subset sums of size k
std::vector<HeapEntry> nth_subsetsum(double *data, const int n, const int k, const unsigned int m){
	if(n < k || k < 1)
		return {};
	if(m > choose(n, k))
		return {};

	HeapEntry best;
	DataPoint sorted_data[n];
	std::vector<HeapEntry> solutions;
	std::priority_queue<HeapEntry> heap;

	for(int i = 0; i < n; ++i)
		sorted_data[i] = {data[i], i, -1};

	std::sort(sorted_data, sorted_data+n);

	for(int idx = 0; idx < n; ++idx)
		sorted_data[idx].sorted_pos = idx;

	for(int i = 0; i < k; ++i){
		best.subset.push_back(sorted_data[i]);
		best.sum += sorted_data[i].val;
	}

	heap.push(best);

	while(solutions.size() != m){
		HeapEntry top = heap.top();
		solutions.push_back(top);
		heap.pop();

		// NOTE: sorting here can be easily removed. the improvement is very small for the use case
		// and therefeore we leave here the sort for readability
		std::sort(top.subset.begin(), top.subset.end(), ascending_index__);

		const int next_pos = top.subset[k-1].sorted_pos + 1;
		if(next_pos >= n)
			continue;

		int pos;
		for(pos = k; --pos ; )
			if(top.subset[pos].sorted_pos != top.subset[pos-1].sorted_pos + 1)
				break;

		for( ; pos < k; ++pos){
			HeapEntry succ = top;
			succ.sum += sorted_data[next_pos].val - succ.subset[pos].val;
			succ.subset[pos] = sorted_data[next_pos];
			heap.push(succ);
		}
	}

	return solutions;
}
