
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allunder/underDIAGV3.h"
#include "alldef/allstruct.h"
#include "fstream"

using namespace std;

vector<DOUBLE> loadfile(string filename)
{
    ifstream inputFile(filename);
    if(!inputFile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(0);
    }

    vector<DOUBLE> result;
    DOUBLE value;
    while(inputFile >> value)
    {
        result.push_back(value);
    }

    inputFile.close();

    return result;
}

vector<DOUBLE> loadfile_apart(string filename,int apart)
{
    ifstream inputFile(filename);
    if(!inputFile.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(0);
    }

    vector<DOUBLE> result;
    result.reserve(4096);
    
    DOUBLE value;

    int count = 0;
    while (inputFile >> value) {
        
        if (count % apart == 0) {
            result.push_back(value);
        }
        count++;
    }

    inputFile.close();

    return result;
}

vector<DOUBLE> loadfile_2(string filename1,string filename2)
{
    std::ifstream file1(filename1);
    if (!file1.is_open()) {
        std::cerr << "Error opening file: filename1" << std::endl;
        exit(0);
    }

    std::vector<std::vector<DOUBLE>> cricket_graeme_1;
    DOUBLE value;
    std::string line;
    while (std::getline(file1, line)) {
        std::istringstream iss(line);
        std::vector<DOUBLE> row;
        while (iss >> value) {
            row.push_back(value);
        }
        cricket_graeme_1.push_back(row);
    }
    file1.close();

    std::ifstream file2(filename2);
    if (!file2.is_open()) {
        std::cerr << "Error opening file: filename2" << std::endl;
        exit(0);
    }

    std::vector<std::vector<DOUBLE>> cricket_sebastian_1;
    while (std::getline(file2, line)) {
        std::istringstream iss(line);
        std::vector<DOUBLE> row;
        while (iss >> value) {
            row.push_back(value);
        }
        cricket_sebastian_1.push_back(row);
    }
    file2.close();

    std::vector<std::vector<DOUBLE>> tmp = cricket_graeme_1;
    tmp.insert(tmp.end(), cricket_sebastian_1.begin(), cricket_sebastian_1.end());

    std::vector<DOUBLE> TS4;
    for (const auto& row : tmp) {
        if (!row.empty()) {
            TS4.push_back(row[0]);
        }
    }

    return  TS4;
}

vector<DOUBLE> loadfile_(string filename)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open file!" << std::endl;
        exit(0);
    }

    std::vector<std::vector<DOUBLE>> matrix;

    std::string line;
    while (std::getline(inputFile, line)) {
        std::vector<DOUBLE> row;
        std::istringstream iss(line);

        DOUBLE value;
        while (iss >> value) {
            row.push_back(value);
        }

        matrix.push_back(row);
    }

    inputFile.close();

    std::vector<DOUBLE> TS5;
    size_t l1 = matrix[0].size();
    size_t l2 = matrix.size();
    for (size_t row = 0; row < l2; ++row) {
        for (size_t col = 0; col < l1; ++col) {
            TS5.push_back(matrix[row][col]);
        }
    }

    return TS5;
}