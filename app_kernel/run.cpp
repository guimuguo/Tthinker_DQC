//########################################################################
//## Copyright 2019 Da Yan http://www.cs.uab.edu/yanda
//##
//## Licensed under the Apache License, Version 2.0 (the "License");
//## you may not use this file except in compliance with the License.
//## You may obtain a copy of the License at
//##
//## //http://www.apache.org/licenses/LICENSE-2.0
//##
//## Unless required by applicable law or agreed to in writing, software
//## distributed under the License is distributed on an "AS IS" BASIS,
//## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//## See the License for the specific language governing permissions and
//## limitations under the License.
//########################################################################

//########################################################################
//## Contributors
//##
//##
//########################################################################

#include "kernel_app.h"
#include <math.h> // used for ceil method, as it won't compile without it on Windows. 
// #include <chrono>
using namespace std::chrono;

int main(int argc, char **argv)
{
    auto start = steady_clock::now();

    if(argc != 9 && argc != 10){
		cout<<"arg1 = input path, arg2 = number of threads"
				<<", arg3 = out-going degree ratio, arg4 = incoming degree ratio, arg5 = min_size, arg6 = time delay threshold, "
						<<"arg7 = kernel file, arg8 = prime_k, arg9 = use trie check"<<endl;
		return -1;
	}
    char* input_path = argv[1];
    num_compers = atoi(argv[2]); //number of compers
	gdmin_deg_ratio_o = atof(argv[3]);
	gdmin_deg_ratio_i = atof(argv[4]);
	gnmin_size = atoi(argv[5]);
	gnmax_size = INT_MAX;
    gnmin_deg_o = ceil(gdmin_deg_ratio_o * (gnmin_size - 1));
    gnmin_deg_i = ceil(gdmin_deg_ratio_i * (gnmin_size - 1));
	TIME_THRESHOLD = atof(argv[6]);
	kernel_file = argv[7];
	prime_k = atoi(argv[8]);
	if(argc == 10)
		trie_check = atoi(argv[9]);


    QCWorker worker(num_compers);
    worker.load_data(input_path);
    // worker.initialize_tasks();
    worker.run();

    auto end = steady_clock::now();
    float duration = (float)duration_cast<milliseconds>(end - start).count() / 1000;
    cout << "Execution Time:" << duration << endl;

    cout << "trie count: " << trie->print_result() << endl;

    log("Done");
}

//./run ...
