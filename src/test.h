#ifndef TEST_H
#define TEST_H

#include "util.h"
#include "hssp.h"

#include <string>    // string, getline, stod, to_string
#include <fstream>   // ifstream
#include <sstream>   // stringstream
#include <vector>    // vector
#include <chrono>    // high_resolution_clock, now, duration
#include <iomanip>   // setprecision
#include <iostream>  // fixed
#include <iterator>  // istream_iterator

using namespace std::chrono;

struct File {
	std::string name; // file name
	int size; // n of the problem
	int output_type; // 0 means no output
	std::string path; // path to file
	std::string ref; // ex "0 0 0"
	std::string problem_type; // "-i" for minimization
};

const std::vector<File> test_files {
	// convex
	{"convex.1s.3d.10.dat", 10, 1, "test/validation/", "1 1 1", ""},
	{"convex.1s.3d.20.dat", 20, 1, "test/validation/", "1 1 1", ""},
	{"convex.1s.3d.30.dat", 30, 1, "test/validation/", "1 1 1", ""},
	{"convex.1s.3d.40.dat", 40, 1, "test/validation/", "1 1 1", ""},
	{"convex.1s.3d.50.dat", 50, 0, "test/validation/", "1 1 1", ""},
	{"convex.1s.3d.60.dat", 60, 0, "test/validation/", "1 1 1", ""},
	{"convex.1s.3d.70.dat", 70, 0, "test/validation/", "1 1 1", ""},
	//{"convex.1s.3d.80.dat", 80, 0, "test/validation/", "1 1 1", ""},

	// cliff
	{"cliff.1s.3d.10.dat", 10, 1, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.20.dat", 20, 1, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.30.dat", 30, 1, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.40.dat", 40, 1, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.50.dat", 50, 1, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.60.dat", 60, 0, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.70.dat", 70, 0, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.80.dat", 80, 0, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.90.dat", 90, 0, "test/validation/", "1 1 1", ""},
	{"cliff.1s.3d.100.dat", 100, 0, "test/validation/", "1 1 1", ""},
	//{"cliff.1s.3d.200.dat", 200, 0, "test/validation/", "1 1 1", ""},

	// linear
	{"linear.1s.3d.10.dat", 10, 1, "test/validation/", "1 1 1", ""},
	{"linear.1s.3d.20.dat", 20, 1, "test/validation/", "1 1 1", ""},
	{"linear.1s.3d.30.dat", 30, 0, "test/validation/", "1 1 1", ""},
	{"linear.1s.3d.40.dat", 40, 0, "test/validation/", "1 1 1", ""},
	{"linear.1s.3d.50.dat", 50, 0, "test/validation/", "1 1 1", ""},
	{"linear.1s.3d.60.dat", 60, 0, "test/validation/", "1 1 1", ""},
	{"linear.1s.3d.70.dat", 70, 0, "test/validation/", "1 1 1", ""},
	//{"linear.1s.3d.80.dat", 80, 0, "test/validation/", "1 1 1", ""},

	// concave
	{"concave.1s.3d.10.dat", 10, 1, "test/validation/", "1 1 1", ""},
	{"concave.1s.3d.20.dat", 20, 1, "test/validation/", "1 1 1", ""},
	{"concave.1s.3d.30.dat", 30, 1, "test/validation/", "1 1 1", ""},
	{"concave.1s.3d.40.dat", 40, 1, "test/validation/", "1 1 1", ""},
	{"concave.1s.3d.50.dat", 50, 1, "test/validation/", "1 1 1", ""},
	{"concave.1s.3d.60.dat", 60, 0, "test/validation/", "1 1 1", ""},
	{"concave.1s.3d.70.dat", 70, 0, "test/validation/", "1 1 1", ""},
	//{"concave.1s.3d.80.dat", 80, 0, "test/validation/", "1 1 1", ""},
};

void run_validation_tests();

#endif
