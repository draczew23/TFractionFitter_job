void runner() {

    int repetitions = 10000;

    std::ofstream outputFile("output_new.txt", std::ios::app);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        exit(0);
    }
    outputFile << "BKG PROBABILITY;BKG ERROR;SIGNAL PROBABILITY;SIGNAL ERROR;CORR;MAN_CORR" << std::endl;

    for (size_t i = 0; i < repetitions; ++i) {
        gROOT->ProcessLine(".x tester.c");
    }
}