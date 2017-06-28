#include "test.h"

const std::string results_folder_path__ = "results/";

std::vector<double> read_results_file(std::string path_to_results){
	std::string line;
	std::vector<double> values;
	std::ifstream file(path_to_results);

	while(std::getline(file, line)){
		if(is_comment(line))
			continue;

		values.push_back(std::stod(line));
	}

	return values;
}

std::vector<std::vector<int> > read_position_file(std::string path_to_positions){
	std::string line;
	std::vector<std::vector<int> > positions;
	std::ifstream file(path_to_positions);

	while(std::getline(file, line)){
		if(is_comment(line))
			continue;

		std::istringstream iss(line);
		std::vector<int> values((std::istream_iterator<int>(iss)), std::istream_iterator<int>());
		positions.push_back(values);
	}

	return positions;
}

std::vector<double> load_outputs(File file){
	std::string test_type, results_file, results_file_path;

	switch(file.output_type){
		case 0:
			return {};
		case 1:
			test_type = file.name.substr(0, file.name.find("."));
			results_file = "myDataSets.1.all." + test_type + "." + "ILP.3D" + "." + std::to_string(file.size) + ".out2";
			results_file_path = file.path + results_folder_path__ + results_file;
			return read_results_file(results_file_path);
		case 2:
			return {};
		default:
			return {};
	}
}

std::vector<std::vector<int> > load_positions(File file){
	std::string test_type, positions_file, positions_file_path;

	switch(file.output_type){
		case 0:
			return {};
		case 1:
			test_type = file.name.substr(0, file.name.find("."));
			positions_file = "myDataSets.1.all." + test_type + "." + "ILP.3D" + "." + std::to_string(file.size) + ".sub";
			positions_file_path = file.path + results_folder_path__ + positions_file;
			return read_position_file(positions_file_path);
		case 2:
			return {};
		default:
			return {};
	}
}

void run_validation_tests(){
	long int total_nodes = 0;
	const int max_repetitions = 1;

	for(File file : test_files){
		std::string path_to_file = file.path + file.name;
		std::vector<double> volumes = load_outputs(file);
		std::vector<std::vector<int> > positions = load_positions(file);

		// if(volumes.empty()) continue;

		for(int k = 1; k <= file.size; ++k){
			std::vector<std::string> args{"-k", std::to_string(k), path_to_file, "-r", file.ref, "-j", "4"};

			if(file.problem_type.size())
				args.push_back(file.problem_type);

			for(int repeat = 1; repeat <= max_repetitions; ++repeat){
				long int nodes;
				std::vector<int> solution;

				auto t1 = high_resolution_clock::now();
				double got_volume = hssp(args, solution, nodes);
				auto t2 = high_resolution_clock::now();
				const duration<float> elap = t2 - t1;

				int type = -1;
				double expected_volume = 0.0;

				if(!volumes.empty() && !positions.empty()){
					expected_volume = volumes[k-1];
					type = eq__(got_volume, expected_volume) && (solution == positions[k-1]);
				}

				if(repeat != max_repetitions)
					continue;

				std::stringstream ss;
				ss << std::fixed << std::setprecision(6) << "got " << got_volume << " expected " << expected_volume << " file " << file.name << " nodes " << nodes << " time " << std::fixed << std::setprecision(6) << elap.count() << "s" << " k " << k;

				if(type == -1)
					logger::warn(ss.str());
				else if(type == 0)
					logger::fail(ss.str());
				else if(type == 1)
					logger::okay(ss.str());

				total_nodes += nodes;
			}
		}
	}

	logger::info("total nodes " + std::to_string(total_nodes));
}
