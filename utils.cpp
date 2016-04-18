//
//  utils.cpp
//  RADTAGpainter
//
//  Created by Milan Malinsky on 29/01/2016.
//  Copyright (c) 2016 Milan Malinsky. All rights reserved.
//

#include "utils.h"

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// Initialize a matrix
void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_size) {
    for (int i = 0; i < m_size; i++) {
        std::vector<double> v(m_size,0);
        m.push_back(v);
    }
}

// Initialize a matrix
void initialize_matrix_int(std::vector<std::vector<int> >& m, int m_size) {
    for (int i = 0; i < m_size; i++) {
        std::vector<int> v(m_size,0);
        m.push_back(v);
    }
}

// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
    return filename; // no suffix
    else
    return filename.substr(0, suffixPos);
}