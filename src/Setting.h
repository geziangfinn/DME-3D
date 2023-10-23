#pragma once

#include "global.h"
using namespace std;



class Setting {
public:
    // basic
    string input_file_name;
    string output_file_name;
    string preOrderfile;
    string inOrderfile;
    string layerfile;
    string benchmark_path = "../benchmarks";
    string topo_choice; // NS, CL, RGM
    int metric = 1; // 1=> L1, 2=>L2(euclidean)
    int refine_M = 6; // by default
    string get_case_name();
};
extern Setting setting;