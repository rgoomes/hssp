#ifndef HYPERVOLUME_H
#define HYPERVOLUME_H

#include "util.h"

extern "C"
{
#include "hv-2.0rc2-src/hv.h"

#include "HVC/hvc.h"
#include "HVC/hvc-class.h"
#include "HVC/io.h"
}

#include <algorithm> // sort, copy
#include <iterator>  // prev, next
#include <set>       // multiset

struct Point2D {
	double x, y;

	bool operator < (const Point2D &p) const {
		return x < p.x;
	}
};

typedef std::multiset<Point2D> Tree;

double hypervolume(Point *, int , const int , const double * );
void all_contributions(Point *, const int , const int , double * , double * );
void get_contributions(Point * , Point * , int , int , int , double * , double * );
void update_contributions(Point * , int , hvc_s * , double * );
hvc_s *build(Point * , Point * , Point * , int , int , int , int , double * , bool );

#endif
