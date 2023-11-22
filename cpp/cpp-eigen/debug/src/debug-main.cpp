#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <variant>
#include <string>
#include <cstddef>

#include <Eigen/Sparse>

using std::string;
using ParameterValue = std::variant<size_t, double, string>;

std::unordered_map<string, ParameterValue> parseInputFile(const string filename);

int main(){
    std::ofstream simlog;
    simlog.open("debug-log.txt");
    simlog << "Beginning simulation" << std::endl;

    
    // Parse input file
    simlog << "Parsing input file" << std::endl;
    std::unordered_map<string, ParameterValue> inputParameters = parseInputFile("../debug/debug.inp");
    size_t N = std::get<size_t>(inputParameters["N"]);
    size_t Nx = std::get<size_t>(inputParameters["Nx"]);
    size_t Nt = std::get<size_t>(inputParameters["Nt"]);
    size_t W = std::get<size_t>(inputParameters["W"]);
    double vprime = std::get<double>(inputParameters["vprime"]);


    // Initialize simulation - particle species, and grid

    // Initial step

    // PIC Loop

    simlog << "Closing log" << std::endl;
    simlog.close();
    return 0;
}

std::unordered_map<string, ParameterValue> parseInputFile(const string filename){
    std::unordered_map<string, ParameterValue> parameters;
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file: " << filename << std::endl;
        return parameters;
    }

    string line, paramName, paramValueStr;

    while(getline(inputFile, line)){
        size_t delimiterPos = line.find("=");
        if (delimiterPos != string::npos){
            string paramName = line.substr(0,delimiterPos);
            string paramValueStr = line.substr(delimiterPos + 1);
            if (paramName == "N" || paramName == "Nt" || paramName == "Nx" || paramName == "W"){
                size_t paramValue = std::stoul(paramName);
                parameters[paramName] = paramValue;
            } 
            else if (paramName == "vprime"){
                double paramValue = std::stod(paramName);
                parameters[paramName] = paramValue;
            }
        }
    }
    return parameters;
}