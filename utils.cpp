//
//  utils.cpp
//  RADTAGpainter
//
//  Created by Milan Malinsky on 29/01/2016.
//  Copyright (c) 2016 Milan Malinsky. All rights reserved.
//

#include "utils.h"

double stringToDouble(std::string s) {
    double d;
    std::stringstream ss(s); //turn the string into a stream
    ss >> d; //convert
    return d;
}

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


std::string suffix(const std::string& seq, size_t len)
{
    assert(seq.length() >= len);
    return seq.substr(seq.length() - len);
}

// Returns true if the filename has an extension indicating it is compressed
bool isGzip(const std::string& filename)
{
    size_t suffix_length = sizeof(GZIP_EXT) - 1;
    
    // Assume files without an extension are not compressed
    if(filename.length() < suffix_length)
        return false;
    
    std::string extension = suffix(filename, suffix_length);
    return extension == GZIP_EXT;
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }
}

// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }
}


void assertGZOpen(gzstreambase& gh, const std::string& fn)
{
    if(!gh.good())
    {
        std::cerr << "Error: could not open " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
}


// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        igzstream* pGZ = new igzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
        assertFileOpen(*pReader, filename);
        return pReader;
    }
}


double calculateInbreedingCoefficient(std::vector<int>& numericGenotypes) {
    int naa = 0; int nAa = 0; int nAA = 0;
    int nSamples = (int)numericGenotypes.size();
    for (std::vector<int>::size_type i = 0; i != nSamples; i++) {
        if (numericGenotypes[i] == 0) naa++;
        if (numericGenotypes[i] == 1) nAa++;
        if (numericGenotypes[i] == 2) nAA++;
    }
    
    // Get the proportions of alt-hom and hets
    double pAA = (double)nAA / nSamples;
    double pAa = (double)nAa / nSamples;
    
    // Allele frequencies
    double p = pAA + (0.5 * pAa);
    double q = 1 - p;
    
    // Get the Hardy-Weinberg prediction for expected number of heterozygotes:
    double HWAa = 2*p*q;
    
    
    // Get the inbreeding coefficient
    double F = (HWAa - pAa) / HWAa;
    return F;
}

