#include <iostream>
using std::ofstream;
using std::ifstream;
#include <fstream>
#include <sstream>
using std::string;
#include "data.h"
template<class T>
T fromString(const std::string& s) {
    std::istringstream stream (s);
    T t;
    stream >> t;
    return t;
}

void delim_to_array(float* data[], string filename, int length, int width, char delimiter) {
    string line;
    int row = 0, col = 0;
    ifstream file;
    file.open(filename);

    if (file.is_open()) {
        while (getline(file, line)) {
            if (row >= length) {
                std::cerr << "delim_to_array: Input file truncated (rows exceed " << length << ")" << std::endl;
                std::cin.get(); // Wait before exit
                exit(EXIT_FAILURE);
            }
            std::stringstream lineStream(line);
            std::string item;
            col = 0;
            while (getline(lineStream, item, delimiter)) {
                if (col >= width) {
                    std::cerr << "delim_to_array: Dimensions mismatch on row " << row << " (columns exceed " << width << ")" << std::endl;
                    std::cin.get(); // Wait before exit
                    exit(EXIT_FAILURE);
                }
                (*data)[row*width+col] = fromString<float>(item);
                col++;
            }
            row++;
        }

        file.close();
    } else {
        std::cerr << "delim_to_array: Failed to open " << filename << std::endl;
        //std::cerr << "delim_to_array: " << strerror(errno) << std::endl;
        std::cin.get(); // Wait before exit
        exit(EXIT_FAILURE);
    }
}

void csv_to_array(float* data[], string filename, int length, int width) {
    delim_to_array(data, filename, length, width, ',');
}

