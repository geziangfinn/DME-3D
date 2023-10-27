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

struct sink{
    double x;
    double y;
    int layer;
    double capacitance;
};

vector<sink> parse(string filePath)
{
    ifstream file(filePath);
    string line, label;
    int numSink;
    vector<sink> sinks;
    while(getline(file, line)){
        if(line=="" || line[0]=='#') continue;
        stringstream ss(line);
        ss >> label;
        if(label == "num"){ // Netlist (Insts)
            ss >>label >>numSink;
            assert(label=="sink");
            int countInst = 0;
            while(countInst < numSink){
                getline(file, line);
                if(line=="") continue;
                stringstream ss2(line);
                // add inst
                sink tempsink;
                ss2 >> label >>label>> tempsink.x >> tempsink.y>>tempsink.layer>>tempsink.capacitance;
                sinks.emplace_back(tempsink);
                countInst++;
            }
            break;
        }
    }
    cout.setf(ios::fixed,ios::floatfield);
    for(sink sink:sinks)
    {
        cout<<sink.x<<" "<<sink.y<<" "<<sink.layer<<" "<<sink.capacitance<<endl;
    }
    return sinks;
}
#endif