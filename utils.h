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
#define BUGREPORT "milan.malinsky@unibas.ch"
#define GZIP_EXT ".gz"
#define THIS_AUTHOR "Milan Malinsky"
#define V "0.2 r102"


std::vector<std::string> split(const std::string &s, char delim);
std::string stripExtension(const std::string& filename);
void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_size);


class BlockCoancestries
{
public:
    BlockCoancestries(int rows, int columns)
    : matrix (rows, std::vector < std::vector <double> > (columns, std::vector<double> (1))),
    p_ij_matrix (rows, std::vector < std::vector <double> > (columns, std::vector<double> (1))),
    c_prime_matrix (rows, std::vector < std::vector <double> > (columns, std::vector<double> (1))),
    matrixNewMissing (rows, std::vector < std::vector <double> > (columns, std::vector<double> (1)))
    {}
    
    std::vector <std::vector <std::vector <double > > > matrix;
    std::vector <std::vector <std::vector <double > > > p_ij_matrix;
    std::vector <std::vector <std::vector <double > > > c_prime_matrix;
    std::vector <std::vector <std::vector <double > > > matrixNewMissing;
};



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

// -------------------------------------    SOME BASIC MATH/STATS  ----------------------------------------

template <class T> double vector_average(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    double average = (double)sum / (double)vector.size();
    return average;
}

template <class T> double vector_sum(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    return sum;
}


inline void copy_except(int i, std::vector<double>& inVec, std::vector<double>& outVec) {
    std::copy(inVec.begin(), inVec.begin() + i, outVec.begin());
    std::copy(inVec.begin() + i + 1, inVec.end(), outVec.begin()+i);
    //std::cerr << "copying:" << i << " "; print_vector_stream(inVec, std::cerr);
    //std::cerr << "copied: " << i << " "; print_vector_stream(outVec, std::cerr);
}

// jackknive standard error
template <class T> double jackknive_std_err(T& vector) {
    std::vector<double> jackkniveAverages;
    std::vector<double> JregionDs; JregionDs.resize(vector.size()-1);
    for (std::vector<double>::size_type i = 0; i != vector.size(); i++) {
        // std::cerr << "copying " << i << std::endl;
        copy_except(i, vector, JregionDs);
        jackkniveAverages.push_back(vector_average(JregionDs));
        JregionDs.clear(); JregionDs.resize(vector.size()-1);
    }
    double jackkniveOverallMean = vector_average(jackkniveAverages);
    double sum = 0;
    for (int i = 0; i < jackkniveAverages.size(); i++) {
        sum += pow((jackkniveAverages[i] - jackkniveOverallMean), 2.0);
    }
    double var = ((double)(jackkniveAverages.size()-1)/(double)jackkniveAverages.size()) * sum;
    double Dstd_err = sqrt(var);
    return Dstd_err;
}

// jackknive standard error of sum
template <class T> double jackknive_std_err_sum(T& vector) {
    std::vector<double> jackkniveSums;
    std::vector<double> JregionDs; JregionDs.resize(vector.size()-1);
    for (std::vector<double>::size_type i = 0; i != vector.size(); i++) {
        // std::cerr << "copying " << i << std::endl;
        copy_except(i, vector, JregionDs);
        jackkniveSums.push_back(vector_sum(JregionDs));
        JregionDs.clear(); JregionDs.resize(vector.size()-1);
    }
    double jackkniveOverallMean = vector_average(jackkniveSums);
    double sum = 0;
    for (int i = 0; i < jackkniveSums.size(); i++) {
        sum += pow((jackkniveSums[i] - jackkniveOverallMean), 2.0);
    }
    double var = ((double)(jackkniveSums.size()-1)/(double)jackkniveSums.size()) * sum;
    double Dstd_err = sqrt(var);
    return Dstd_err;
}

#endif /* defined(__RADTAGpainter__utils__) */
