// #include "simulation_parameters.h"
// #include <iostream>
// #include <string>

void runner() {

    int repetitions = 1;

    std::ofstream outputFile("output_new.txt", std::ios::app);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        exit(0);
    }
    outputFile << "BKG PROBABILITY;BKG ERROR;SIGNAL PROBABILITY;SIGNAL ERROR;CORR;MAN_CORR" << std::endl;

    double counter = 0;
    for (size_t i = 0; i < repetitions; ++i) {
        gROOT->ProcessLine(".x tester.c");
        cout << counter/repetitions << endl;
        counter += 1;
    }
}