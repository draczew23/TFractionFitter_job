#include <iostream>
#include <fstream>

void outfile() {
    // Open a file for appending
    std::ofstream outputFile("output.txt", std::ios::app);

    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
    }

    // Use a for loop to write numbers to the file
    for (int i = 1; i <= 10; ++i) {
        // Write the number to the file
        outputFile << i << ";" << endl;
    }

    // Close the file
    outputFile.close();

    std::cout << "Data appended to file successfully!" << std::endl;
}
