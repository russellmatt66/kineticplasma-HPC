#define _USE_MATH_DEFINES

#include<fstream>

int main(){
    std::ofstream simlog;
    simlog.open("debug.log");
    simlog << "Beginning simulation" << std::endl;
    simlog << "Parsing input file" << std::endl;
    // Parse input file
    
    // Initialize simulation - particle species, and grid

    // Initial step

    // PIC Loop

    simlog << "Closing log" << std::endl;
    simlog.close();
    return 0;
}