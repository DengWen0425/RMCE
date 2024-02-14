#include "parmce.hpp"
#include <tbb/tbb.h>
#include <iostream>
#include <string>


using namespace std;

int main(int argc, char *argv[])
{
    int num_thread = tbb::task_scheduler_init::default_num_threads();
    if (argc < 2) {
        cout << "Too few parameters!!!!" << endl;
        return 0;
    } 

    if (argc == 2) {
        cout << "Use default thread number: " << num_thread << endl;
    } else if (argc == 3) {
        cout << "Use thread number: " << argv[2] << endl;
        num_thread = atoi(argv[2]);
    }

    const string file_path = argv[1];

    parmce_main(file_path, num_thread);

    return 0;

}