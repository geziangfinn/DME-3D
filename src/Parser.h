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
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
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