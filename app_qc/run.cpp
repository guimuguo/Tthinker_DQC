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

#include "qc_app.h"
#include <math.h> // used for ceil method, as it won't compile without it on Windows.

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include "systemI.h"

using namespace std;

// #include <chrono>
using namespace std::chrono;

extern std::atomic<bool> global_end_label_t;

int main(int argc, char **argv)
{   
    char outputfileTime[1000];
	char outputfile[1000];
	char outputfilePeakMem[1000];
    
    int comper = atoi(argv[2]);
    sprintf(outputfile, "%s_%d_opsm_output.txt", string(argv[1]).c_str(),comper);
    sprintf(outputfileTime, "%s_%d_running_time_output.txt", string(argv[1]).c_str(),comper);
    sprintf(outputfilePeakMem, "%s_%d_maxmem_output.txt",string(argv[1]).c_str(),comper);

	
    ofstream fout(outputfile);
	ofstream foutPeakMem(outputfilePeakMem);
	ofstream foutTime(outputfileTime);

    GetCurrentPid();
    thread t = thread(info, GetCurrentPid(), ref(foutPeakMem));  // launch a thread to record

    auto start_t = steady_clock::now();
    cout << "Start testing memory and time..." << endl;

    
    if (argc < 7)
    {
        cout << "arg1 = input path, arg2 = number of threads"
             << ", arg3 = out degree ratio, arg4 = in degree ratio, arg5 = min_size, arg6 = time delay threshold" << endl;
        return -1;
    }
    char *input_path = argv[1];
    num_compers = atoi(argv[2]); //number of compers
    gdmin_deg_ratio_o = atof(argv[3]);
    gdmin_deg_ratio_i = atof(argv[4]);
    gnmin_size = atoi(argv[5]);
    gnmax_size = INT_MAX;
    gnmin_deg_o = ceil(gdmin_deg_ratio_o * (gnmin_size - 1));
    gnmin_deg_i = ceil(gdmin_deg_ratio_i * (gnmin_size - 1));
    TIME_THRESHOLD = atof(argv[6]); // tau_time
    MAX_RESULT_FOUND = atoi(argv[7]);
    
    if (argc > 7)
        BIGTASK_THRESHOLD = atof(argv[7]); // tau_split
    if (argc > 8)
        tasks_per_fetch_g = atoi(argv[8]);
    if (argc > 9)
        Qbig_capacity = atoi(argv[9]);
    if (argc > 10)
        Qreg_capacity = atoi(argv[10]);
    if (argc > 11)
        BT_TASKS_PER_FILE = atoi(argv[11]);
    if (argc > 12)
        MINI_BATCH_NUM = atoi(argv[12]);
    if (argc > 13)
        RT_TASKS_PER_FILE = atoi(argv[13]);
    if (argc > 14)
        BT_THRESHOLD_FOR_REFILL = atoi(argv[14]);
    if (argc > 15)
        RT_THRESHOLD_FOR_REFILL = atoi(argv[15]);

    cout << "input_path:" << argv[1] << endl;
    cout << "num_compers:" << num_compers << endl;
    cout << "out gdmin_deg_ratio:" << gdmin_deg_ratio_o << endl;
    cout << "in gdmin_deg_ratio:" << gdmin_deg_ratio_i << endl;
    cout << "gnmin_size:" << gnmin_size << endl;
    cout << "TIME_THRESHOLD:" << TIME_THRESHOLD << endl;
    cout << "BIGTASK_THRESHOLD:" << BIGTASK_THRESHOLD << endl;
    cout << "tasks_per_fetch_g:" << tasks_per_fetch_g << endl;
    cout << "Qbig_capacity:" << Qbig_capacity << endl;
    cout << "Qreg_capacity:" << Qreg_capacity << endl;
    cout << "BT_TASKS_PER_FILE:" << BT_TASKS_PER_FILE << endl;
    cout << "MINI_BATCH_NUM:" << MINI_BATCH_NUM << endl;
    cout << "RT_TASKS_PER_FILE:" << RT_TASKS_PER_FILE << endl;
    cout << "BT_THRESHOLD_FOR_REFILL:" << BT_THRESHOLD_FOR_REFILL << endl;
    cout << "RT_THRESHOLD_FOR_REFILL:" << RT_THRESHOLD_FOR_REFILL << endl;
    

    QCWorker worker(num_compers);

    auto start = steady_clock::now();
    worker.load_data(input_path);
    auto end = steady_clock::now();
    float duration = (float)duration_cast<milliseconds>(end - start).count() / 1000;
    cout << "load_data() execution time:" << duration << endl;
    worker.latest_elapsed_time = duration; // used to calculate step time inside the worker

    start = steady_clock::now();
    worker.run();
    end = steady_clock::now();
    duration = (float)duration_cast<milliseconds>(end - start).count() / 1000;
    cout << "run() execution Time:" << duration - worker.init_time << endl;
    // cout << "trie count: " << trie->print_result() << endl;
    // cout << "================================================" << endl;


    auto end_t = steady_clock::now();
    cout << "disk time:"<<disk_time<<endl;
    cout << "2-hop time:"<<hop2_time<<endl;
    cout << "Running time(without disk time): " << (float)duration_cast<milliseconds>(end_t - start_t).count() / 1000 - disk_time<< " s" << endl;
    cout << "Running time(without disk time and 2hop): " << (float)duration_cast<milliseconds>(end_t - start_t).count() / 1000 - disk_time - hop2_time<< " s" << endl;

    foutTime << (double)duration_cast<milliseconds>(end_t - start_t).count() / 1000 << " s" << endl;
    foutTime.close();

    global_end_label_t = false;
	sleep(1);
	t.join();
	foutPeakMem.close();
	
	fout.close();

    cout << "Testing ends..." << endl;

    log("Done");
}

//./run ...
