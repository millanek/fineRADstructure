//
//  utils.h
//  RADTAGpainter
//
//  Created by Milan Malinsky on 29/01/2016.
//  Copyright (c) 2016 Milan Malinsky. All rights reserved.
//

#ifndef __RADTAGpainter__utils__
#define __RADTAGpainter__utils__

#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <vector>

#include <map>
#include <sstream>
#include <fstream>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <stdexcept>
#include <math.h>
#include <numeric>
#include <ctime>
#include <regex>

using std::string;

#define BIN "RADpainter"
#define BUGREPORT "mm812@cam.ac.uk"
#define GZIP_EXT ".gz"
#define THIS_AUTHOR "Milan Malinsky"
#define V "0.1 r100"


std::vector<std::string> split(const std::string &s, char delim);
std::string stripExtension(const std::string& filename);
void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_size);

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
}

// Print an arbitrary matrix (vector of vectors)
template <class T> void print_matrix(T matrix, std::ofstream& outFile) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (j == (matrix[i].size()-1))
            outFile << matrix[i][j] << std::endl;
            else
            outFile << matrix[i][j] << "\t";
        }
    }
}

// Print an arbitrary matrix (vector of vectors)
template <class T> void print_matrix_wNames(T matrix, std::ofstream& outFile, std::vector<std::string>& names) {
    for (int i = 0; i < matrix.size(); i++) {
        outFile << names[i] << "\t";
        for (int j = 0; j < matrix[i].size(); j++) {
            if (j == (matrix[i].size()-1))
            outFile << matrix[i][j] << std::endl;
            else
            outFile << matrix[i][j] << "\t";
        }
    }
}


// Print an arbitrary vector to a file
template <class T> void print_vector(T vector, std::ofstream& outFile, char delim = '\t') {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1))
        outFile << vector[i] << std::endl;
        else
        outFile << vector[i] << delim;
    }
}

// Print an arbitrary vector to an output stream
template <class T> void print_vector_stream(T vector, std::ostream& outStream, char delim = '\t') {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1))
        outStream << vector[i] << std::endl;
        else
        outStream << vector[i] << delim;
    }
}



#endif /* defined(__RADTAGpainter__utils__) */
