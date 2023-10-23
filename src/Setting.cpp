#include "Setting.h"

string Setting::get_case_name()
{
    string caseName;
    string fileName;
    string inFileName =input_file_name;
    char sep = '/'; 
    size_t i = inFileName.rfind(sep, inFileName.length());
    if (i != string::npos) {
        fileName = inFileName.substr(i+1, inFileName.length() - i);
    } else {
        fileName = inFileName;
    }
    char sep2 = '.';
    i = fileName.rfind(sep2, fileName.length());
    if (i != string::npos) {
        caseName = fileName.substr(0, i);
    } else{
        caseName = fileName;
    }
    return caseName; 
}