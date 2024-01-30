#ifndef PARSER_H
#define PARSER_H
#include <iostream>
#include <boost/format.hpp>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include "global.h"
using namespace std;

struct parser_sink{
    double x;
    double y;
    int layer=1;
    double capacitance;
};
struct parser_blockage{
    double llx;//ll: lower left, ur: upper right
    double lly;
    double urx;
    double ury;
    int layer=1;
};

class ispd2009_parser{
    public:
    string filePath;
    vector<parser_sink> sinks;
    vector<parser_blockage> blockages;
    ispd2009_parser()
    {
        filePath="";
    }
    ispd2009_parser(string _filePath)
    {
        filePath=_filePath;
    }
    void setFilePath(string path)
    {
        filePath=path;
    }
    vector<parser_sink> get_sinks();
    vector<parser_blockage> get_blockages();
    void parse();
};
#endif