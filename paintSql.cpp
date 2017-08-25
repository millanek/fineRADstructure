//
//  paintSql.cpp
//  RADTAGpainter
//
//  Created by Milan Malinsky on 29/01/2016.
//  Copyright (c) 2016 Milan Malinsky. All rights reserved.
//

#include "paintSql.h"
#include "utils.h"

#define SUBPROGRAM "paint"

static const int initMissing = 2;

static const char *PAINTSQL_USAGE_MESSAGE =
"Usage: " BIN " " SUBPROGRAM " [OPTIONS] INPUT.txt\n"
"Generate a co-ancestry matrix from RAD data\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -t, --tetraploid                        the tetraploid mode (see manual for input format)\n"
"       -c, --chr                               output per-chromosome coancestry matrices\n"
"       -n, --run-name                          run-name will be included in the output file name(s)\n\n"
"       -m, --missing2                          (deprecated) output a conancestry matrix with missing data treated\n"
"                                               as if any pair of individuals are equally distant\n"
"\n\n"
"\nReport bugs to " BUGREPORT "\n\n";


// Options
static const char* shortopts = "hn:mtc";
static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "chr",   no_argument, NULL, 'c' },
    { "tetraploid",   no_argument, NULL, 't' },
    { "help",   no_argument, NULL, 'h' },
    { "missing2", no_argument, NULL, 'm' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string runName = "";
    static string sqlFileName = "";
    static bool bTetraploid = false;
    static bool bOutputChr = false;
    static bool bMissing2 = false;
}

// Just recording where either donor or recipient are "missing data"
static std::vector<std::vector<double> > outChunksNoMissing;
static std::vector<std::vector<double> > missingnessMatrix;
static std::vector<std::vector<double> > localChunksNoMissing;
static std::vector<std::vector<double> > localMissingnessMatrix;


void incrementMissingessMatrixBothHaplotypes(std::vector<std::vector<double> >& missingnessMatrix, int i) {
    int Nindividuals = (int)missingnessMatrix[i].size();
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            missingnessMatrix[i][j] = missingnessMatrix[i][j] + 1.0;
        }
    }
}
void incrementMissingessMatrixOneHaplotype(std::vector<std::vector<double> >& missingnessMatrix, int i) {
    int Nindividuals = (int)missingnessMatrix[i].size();
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            missingnessMatrix[i][j] = missingnessMatrix[i][j] + 0.5;
        }
    }
}

// Helper functions
bool checkSimilaritiesSumToOne(const std::vector<double>& simVectorPerInd) {
    double sum = 0;
    for (int i = 0; i < simVectorPerInd.size(); i++) {
        sum = sum + simVectorPerInd[i];
    }
    // std::cout << "Sum: " << sum << std::endl;
    if (sum > 0.9999999999 && sum < 1.00000000001)
        return true;
    else
        return false;
}


// get the proportion of difference between two haplotypes
double compareSeqs(const std::string& seq1, const std::string& seq2) {
    if ((seq1.length() != seq2.length())) {
        std::cerr << "seq1: " << seq1 << "!" << seq1.length() << std::endl;
        std::cerr << "seq2: " << seq2 << "!" << seq2.length()<< std::endl;
    }
    assert(seq1.length() == seq2.length());
    int numDiff = 0; int numNs = 0;
    for (int i = 0; i < seq1.length(); i++) {
        if (seq1[i] != 'N' && seq2[i] != 'N') {
            if (seq1[i] != seq2[i]) {
                numDiff++;
            }
        } else {
            numNs++;
        }
    }
    double diffR = (double)numDiff/(seq1.length() - numNs);
    return diffR;
}

std::vector<int> findMinDiffIndices(std::vector<double>& diffs, double minDiff) {
    std::vector<int> indices;
    std::vector<double>::iterator iter = diffs.begin();
    while ((iter = std::find(iter, diffs.end(), minDiff)) != diffs.end())
    {
        indices.push_back((int)(iter - diffs.begin()));
        iter++;
    }
    //std::cerr << "indices.size(): " << indices.size() << std::endl;
    assert(indices.size() > 0);
    return indices;
}


double sumCoancestryReceivedSoFar(int nSamples, int i, std::vector<std::vector<double> >& coancestryM) {
    double coancestrySoFarSum = 0;
    for (int j = 0; j < nSamples; j++) {
        if (j != i) {
            coancestrySoFarSum = coancestrySoFarSum + coancestryM[i][j];
        }
    }
    return coancestrySoFarSum;
}

// For missing recipient assume that is is eqaully similar to all the possible donors
std::vector<double> getMissingRecipientSimVector(int nSamples, int i) {
    std::vector<double> recipientSimVector(nSamples, 0.0);
    for (int j = 0; j < nSamples; j++) {
        if (j != i) {
            recipientSimVector[j] = (1.0/(double)(nSamples-1));
        }
    }
    return recipientSimVector;
}

// For missing recipient assume that its similarity to the possible donors is given by the co-ancestry matrix
std::vector<double> getMissingRecipientSimVector(int nSamples, int i, std::vector<std::vector<double> >& coancestryM, double coancestrySoFarSum) {
    std::vector<double> recipientSimVector(nSamples, 0.0);
    for (int j = 0; j < nSamples; j++) {
        if (j != i) {
            recipientSimVector[j] = coancestryM[i][j]/coancestrySoFarSum;
        }
    }
    return recipientSimVector;
}



bool checkIfTagInformative(const std::vector<std::string>& fields) {
    int missing = 0;
    std::regex Ns("N+");
    std::regex NsHet("N+/N+");
    for (int i = 0; i < fields.size(); i++) {
        if (fields[i] == "" || fields[i] == " " || std::regex_match(fields[i], Ns) || std::regex_match(fields[i], NsHet)) {
            missing++;
        }
    }
    if (missing >= fields.size() - 1) {
        return false;
    } else {
        return true;
    }
}


// Divide similarity equally between samples that are the closest - to equal 1 in total
// For a missing donor, assign similarity 1/(N-1) (N the number of samples), and take a corresponding amount away from all the closest ones
std::vector<double> calculateSimilarity(const std::vector<std::string>& allHaps, const std::string& recipientHap, int thisIndI, int nSamples, std::vector<std::vector<double> >& coancestryM, int tagsSoFar, int numRecipientHaps) {
    double sumOfCoancestryReceivedSoFar = sumCoancestryReceivedSoFar(nSamples, thisIndI, coancestryM);

    if (std::regex_match(recipientHap, std::regex("N+"))) {
	//std::cerr << "recipientHap is N: " << recipientHap << std::endl;
        incrementMissingessMatrixOneHaplotype(missingnessMatrix, thisIndI);
        incrementMissingessMatrixOneHaplotype(localMissingnessMatrix, thisIndI);
        if (tagsSoFar < initMissing) return getMissingRecipientSimVector(nSamples, thisIndI);
        else return getMissingRecipientSimVector(nSamples, thisIndI, coancestryM, sumOfCoancestryReceivedSoFar);
    }
    std::vector<double> diffVector(nSamples*2, 0.0); // Proportion of difference between all haplotypes and the recipient
    std::vector<double> simVectorPerInd(nSamples, 0.0);
    std::vector<double> simVectorPerIndNoMissing(nSamples, 0.0);
    std::vector<double> thisMissing(nSamples, 0.0);
    double numMissing = 0;
    double totalToMissing = 0;
    double thisMissingBasic = 0;
    const double missingRandomDiff = 1-pow(0.25,recipientHap.length());
    assert(nSamples == allHaps.size());
    for (int i = 0; i < allHaps.size(); i++) {
        if (i != thisIndI) {
            std::vector<std::string> donorHaps = split(allHaps[i], '/');
            assert(donorHaps.size() < 3);
            if (donorHaps.size() == 2) {
                if (std::regex_match(donorHaps[0], std::regex("N+"))) {
                    numMissing = numMissing + 0.5;
                    thisMissing[i] += (1.0/numRecipientHaps)/2;
                    diffVector[2*i] = missingRandomDiff;
                    thisMissingBasic = thisMissingBasic + (1/(double)(nSamples-1))/2;
                    if (tagsSoFar < initMissing) {
                        simVectorPerInd[i] += (1/(double)(nSamples-1))/2;
                        totalToMissing += (1/(double)(nSamples-1))/2;
                    } else {
                        simVectorPerInd[i] += (coancestryM[thisIndI][i]/sumOfCoancestryReceivedSoFar)/2;
                        totalToMissing += (coancestryM[thisIndI][i]/sumOfCoancestryReceivedSoFar)/2;
                    }
                } else {
                    diffVector[2*i] = compareSeqs(recipientHap, donorHaps[0]);
                    // outChunksNoMissing[thisIndI]
                }
                if (std::regex_match(donorHaps[1], std::regex("N+"))) {
                    numMissing = numMissing + 0.5; thisMissing[i] += (1.0/numRecipientHaps)/2;
                    // std::cerr << "numMissing: " << numMissing << std::endl;
                    diffVector[(2*i)+1] = missingRandomDiff;
                    thisMissingBasic = thisMissingBasic + (1/(double)(nSamples-1))/2;
                    if (tagsSoFar < initMissing) {
                        simVectorPerInd[i] += (1/(double)(nSamples-1))/2;
                        totalToMissing += (1/(double)(nSamples-1))/2;
                    } else {
                        simVectorPerInd[i] += (coancestryM[thisIndI][i]/sumOfCoancestryReceivedSoFar)/2;
                        totalToMissing += (coancestryM[thisIndI][i]/sumOfCoancestryReceivedSoFar)/2;
                    }
                } else {
                    diffVector[(2*i)+1] = compareSeqs(recipientHap, donorHaps[1]);
                }
            } else {
                if (allHaps[i] == "" || allHaps[i] == " " || std::regex_match(allHaps[i], std::regex("N+"))) {
                    numMissing = numMissing + 1; thisMissing[i] += (1.0/numRecipientHaps);
                    diffVector[2*i] = missingRandomDiff; diffVector[(2*i)+1] = diffVector[2*i];
                    thisMissingBasic = thisMissingBasic + (1/(double)(nSamples-1));
                    if (tagsSoFar < initMissing) {
                        simVectorPerInd[i] += (1/(double)(nSamples-1));
                        totalToMissing += (1/(double)(nSamples-1));
                    } else {
                        simVectorPerInd[i] += (coancestryM[thisIndI][i]/sumOfCoancestryReceivedSoFar);
                        totalToMissing += (coancestryM[thisIndI][i]/sumOfCoancestryReceivedSoFar);
                    }
                } else {
                    diffVector[2*i] = compareSeqs(recipientHap, donorHaps[0]);
                    diffVector[(2*i)+1] = diffVector[2*i];
                }
            }
        } else {
            diffVector[2*i] = 1.1; // Comparing with itself - so assigning difference proportion above one as we are not interested in within-individual comparisons
            diffVector[(2*i)+1] = 1.1; // Comparing with itself - so assigning difference proportion above one as we are not interested in within-individual comparisons
        }
    }
    // Find the smallest diff:
    double minDiff = *std::min_element(diffVector.begin(),diffVector.end());
    int numClosest = (int)std::count(diffVector.begin(),diffVector.end(), minDiff);
    std::vector<int> indicesMin = findMinDiffIndices(diffVector, minDiff);
    assert(numClosest == indicesMin.size());
    double subtractFromClosestNoMissing;
    if (totalToMissing > 0) { assert(minDiff <= missingRandomDiff);}
    if (minDiff < missingRandomDiff) {
        for (int i = 0; i < nSamples; i++) {
            missingnessMatrix[thisIndI][i] += thisMissing[i];
            localMissingnessMatrix[thisIndI][i] += thisMissing[i];
        }
        subtractFromClosestNoMissing = totalToMissing/(double)numClosest; // Take this away from the closest (assigned to missing)
    } else {
        // the missing ones are actually the "closest"
        subtractFromClosestNoMissing = 0;
    }
    double subtractFromClosest = totalToMissing/(double)numClosest; // Take this away from the closest (assigned to missing)
    //std::cerr << "thisIndI: " << thisIndI << std::endl;
    //std::cerr << "numClosest: " << numClosest << std::endl;
    //std::cerr << "numMissing: " << numMissing << std::endl;
    
    // Fill in the per-individual co-ancestry vector
    //std::cerr << "totalToMissing: " << totalToMissing << std::endl;
    double add = (1/(double)numClosest) - subtractFromClosest;
    double addNoMissing = (1/(double)numClosest) - subtractFromClosestNoMissing;
    //std::cerr << "add: " << add << std::endl;
    for (int i = 0; i < simVectorPerInd.size(); i++) {
        if (std::find(indicesMin.begin(), indicesMin.end(), 2*i) != indicesMin.end()) {
            //if (thisMissing[i] != 1) {
            simVectorPerInd[i] += add;
            simVectorPerIndNoMissing[i] += addNoMissing;
            outChunksNoMissing[thisIndI][i] += addNoMissing;
            localChunksNoMissing[thisIndI][i] += addNoMissing;
            //}
        }
        if (std::find(indicesMin.begin(), indicesMin.end(), (2*i)+1) != indicesMin.end()) {
            simVectorPerInd[i] += add;
            simVectorPerIndNoMissing[i] += addNoMissing;
            outChunksNoMissing[thisIndI][i] += addNoMissing;
            localChunksNoMissing[thisIndI][i] += addNoMissing;
        }
    }
    
    if (minDiff == missingRandomDiff) {
        //std::cerr << "thisIndI: " << thisIndI << std::endl;
        if(checkSimilaritiesSumToOne(simVectorPerIndNoMissing) == false) {
            double sum = 0;
            for (int i = 0; i < simVectorPerIndNoMissing.size(); i++) {
                sum = sum;
            }
            std::cerr << "sum: " << sum << std::endl;
            print_vector_stream(diffVector, std::cerr);
            print_vector_stream(simVectorPerIndNoMissing, std::cerr);
            
        } assert(checkSimilaritiesSumToOne(simVectorPerIndNoMissing));
    }
    
    if(checkSimilaritiesSumToOne(simVectorPerInd) == false) {
        double sum = 0;
        for (int i = 0; i < simVectorPerInd.size(); i++) {
            sum = sum + simVectorPerInd[i];
        }
        std::cerr << "sum: " << sum << std::endl;
        print_vector_stream(diffVector, std::cerr);
        print_vector_stream(simVectorPerInd, std::cerr);
        
    } assert(checkSimilaritiesSumToOne(simVectorPerInd));
    // print_vector_stream(simVectorPerInd, std::cout);
    return simVectorPerInd;
}

// For missing recipient assume that is is eqaully similar to all the possible donors
void incrementMissingRecipient(std::vector<std::vector<double> >& outChunksMatrix, int i) {
    int Nindividuals = (int)outChunksMatrix[i].size();
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            outChunksMatrix[i][j] = outChunksMatrix[i][j] + (1.0/(double)(Nindividuals-1));
        }
    }
}




// For missing recipient assume that its similarity to possible donors is given by the co-ancestry matrix seen so far
void incrementMissingRecipientOnCoancestry(std::vector<std::vector<double> >& outChunksMatrix, int i) {
    int Nindividuals = (int)outChunksMatrix[i].size();
    double coancestrySoFarSum = 0;
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            coancestrySoFarSum = coancestrySoFarSum + outChunksMatrix[i][j];
        }
    }
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            outChunksMatrix[i][j] = outChunksMatrix[i][j] + (outChunksMatrix[i][j]/(double)coancestrySoFarSum);
        }
    }
}

// For missing recipient assume that its similarity to possible donors is given by the co-ancestry matrix seen so far
void incrementMissingRecipientOnCoancestryLocal(std::vector<std::vector<double> >& outChunksMatrix, std::vector<std::vector<double> >& localChunksMatrix, int i) {
    int Nindividuals = (int)outChunksMatrix[i].size();
    double coancestrySoFarSum = 0;
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            coancestrySoFarSum = coancestrySoFarSum + outChunksMatrix[i][j];
        }
    }
    assert(coancestrySoFarSum == vector_sum(outChunksMatrix[i]));
    for (int j = 0; j < Nindividuals; j++) {
        if (j != i) {
            localChunksMatrix[i][j] = localChunksMatrix[i][j] + (outChunksMatrix[i][j]/(double)coancestrySoFarSum);
        }
    }
}




void checkInputType(const std::vector<std::string>& fields, std::string& inputType) {
    if (fields[2] == "Chr") {
        inputType = "Stacks";
    } else if (fields[0] == "Chr") {
        inputType = "Matrix";
    } else {
        std::cerr << "No location (\"Chr\") info found - assuming the input is a simple data matrix" << std::endl;
        inputType = "SimpleMatrix";
    }
}

std::vector<double> calulateConcastryTetraploid(const std::vector<std::string>& recipientHaps,const std::vector<std::string>& fields, int i, int numIndividuals, std::vector<std::vector<double> >& outChunksMatrix, int tagsRead) {
    assert(recipientHaps.size() == 4); // Handle make sure this is a correctly formatted tatraploid site
    std::vector<double> recipientSimVector;
    string h1 = recipientHaps[0];
    std::vector<double> recipientSimVectorH1 = calculateSimilarity(fields, h1, i, numIndividuals,outChunksMatrix,tagsRead,4);
    string h2 = recipientHaps[1];
    std::vector<double> recipientSimVectorH2 = calculateSimilarity(fields, h2, i, numIndividuals,outChunksMatrix,tagsRead,4);
    string h3 = recipientHaps[2];
    std::vector<double> recipientSimVectorH3 = calculateSimilarity(fields, h3, i, numIndividuals,outChunksMatrix,tagsRead,4);
    string h4 = recipientHaps[3];
    std::vector<double> recipientSimVectorH4 = calculateSimilarity(fields, h4, i, numIndividuals,outChunksMatrix,tagsRead,4);
    for (int j = 0; j < recipientSimVectorH1.size(); j++) {
        recipientSimVector.push_back((recipientSimVectorH1[j] + recipientSimVectorH2[j] + recipientSimVectorH3[j] + recipientSimVectorH3[j])/4);
    }
    return recipientSimVector;
}

struct AllPerChrData {
    std::map<std::string, int> numTagsPerChr;
    std::map<std::string,std::vector<std::vector<double> > > chunksMatrixPerChr;
};

void checkChr(const std::string& thisChr, AllPerChrData& chrTagsAndChunks, int numIndividuals) {
    if (thisChr != "") {
        assert(chrTagsAndChunks.numTagsPerChr.count(thisChr) == chrTagsAndChunks.chunksMatrixPerChr.count(thisChr));
        if (chrTagsAndChunks.numTagsPerChr.count(thisChr) == 1) {
            chrTagsAndChunks.numTagsPerChr[thisChr]++;
        } else {
            chrTagsAndChunks.numTagsPerChr[thisChr] = 1;
            std::vector<std::vector<double> > newThisChrChunksMatrix;
            initialize_matrix_double(newThisChrChunksMatrix, numIndividuals);
            chrTagsAndChunks.chunksMatrixPerChr[thisChr] = newThisChrChunksMatrix;
        }
    } else {
        
    }
}

void addChunksToChr(const std::string& thisChr, AllPerChrData& chrTagsAndChunks, const std::vector<double>& recipientSimVector, int numIndividuals, int i) {
    if (thisChr != "") {
        assert(chrTagsAndChunks.chunksMatrixPerChr.count(thisChr) == 1);
        for (int j = 0; j < recipientSimVector.size(); j++) {
            chrTagsAndChunks.chunksMatrixPerChr[thisChr][i][j] = chrTagsAndChunks.chunksMatrixPerChr[thisChr][i][j] + recipientSimVector[j];
        }
    } else {
        
    }
}


void addChunksToLocalMatrix(const std::string& thisChr, AllPerChrData& chrTagsAndChunks, const std::vector<double>& recipientSimVector, int numIndividuals, int i) {
    if (thisChr != "") {
        assert(chrTagsAndChunks.chunksMatrixPerChr.count(thisChr) == 1);
        for (int j = 0; j < recipientSimVector.size(); j++) {
            chrTagsAndChunks.chunksMatrixPerChr[thisChr][i][j] = chrTagsAndChunks.chunksMatrixPerChr[thisChr][i][j] + recipientSimVector[j];
        }
    } else {
        
    }
}



int paintSqlMain(int argc, char** argv) {
    parsePaintSqlOptions(argc, argv);
    string inputType;
    
    std::cerr << "Painting RAD tags from: " << opt::sqlFileName << std::endl;
    std::ifstream* sqlFile = new std::ifstream(opt::sqlFileName.c_str());
    string fileRoot = stripExtension(opt::sqlFileName);
    
    // Create output files:
    string outChunksMatrixFileName = fileRoot + "_chunks.out";
    std::ofstream* outChunksMatrixFile = new std::ofstream(outChunksMatrixFileName.c_str());
    string outMissingnessFileName = fileRoot + "_missingness.out";
    std::ofstream* outMissingnessFile = new std::ofstream(outMissingnessFileName.c_str());
    string outMissingnessMatrixFileName = fileRoot + "_missingnessMatrix.out";
    std::ofstream* outMissingnessMatrixFile = new std::ofstream(outMissingnessMatrixFileName.c_str());

    std::vector<std::vector<double> > outChunksMatrix;
    double runtime; double runStart = clock();
    
    // Some variables for c calculations
    string thisChr; AllPerChrData chrTagsAndChunks;
    std::vector<std::vector<double>> localChunkMatrix;
    std::vector<std::vector<double>> empiricalVariancesMatrix; std::vector<std::vector<double>> empiricalVariancesPaper;
    std::vector<std::vector<double>> theoreticalVariancesMatrix;
    std::vector<std::vector<double>> c_ij;
    std::vector<std::vector<double>> p_ij_full;
    std::vector<std::vector<double>> s2_ij;
    std::vector<std::vector<double>> my_empiricalVar;
    int blockSize = 100;
    int ploidy = 2;
    int notInformative = 0;
    
    int tagsRead = 0;
    string line;
    getline(*sqlFile, line);
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
   
    //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    
    std::vector<std::string> individuals = split(line, '\t'); checkInputType(individuals, inputType);
    if (inputType == "Stacks") { individuals.erase(individuals.begin(), individuals.begin()+12); }
    else if (inputType == "Matrix") { individuals.erase(individuals.begin());}
    std::cerr << "The file seems to be in a " << inputType << " format" << std::endl;
    int numIndividuals = (int)individuals.size();
    initialize_matrix_double(outChunksMatrix, numIndividuals);
    initialize_matrix_double(outChunksNoMissing, numIndividuals);
    initialize_matrix_double(missingnessMatrix, numIndividuals);
    BlockCoancestries bC(numIndividuals,numIndividuals);
    initialize_matrix_double(localChunkMatrix, numIndividuals);
    initialize_matrix_double(localMissingnessMatrix, numIndividuals);
    initialize_matrix_double(localChunksNoMissing, numIndividuals);
    initialize_matrix_double(empiricalVariancesMatrix, numIndividuals);
    initialize_matrix_double(empiricalVariancesPaper, numIndividuals);
    initialize_matrix_double(theoreticalVariancesMatrix, numIndividuals);
    initialize_matrix_double(c_ij, numIndividuals);
    initialize_matrix_double(p_ij_full, numIndividuals);
    initialize_matrix_double(s2_ij, numIndividuals);
    initialize_matrix_double(my_empiricalVar, numIndividuals);
    std::vector<double> missingness(numIndividuals,0.0);
    //std::cerr << "Number of columns: " << fields[11] << std::endl;
    // Lets start going through the tags:
    while (getline(*sqlFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows

        bool bIndTriallelic = false;
        
        if (tagsRead % 100 == 0 && tagsRead > 0) {
            std::cerr << "Processed: " << tagsRead << " tag loci" << std::endl;
        }
        
        std::vector<std::string> fields;
        
        if (inputType == "Stacks") {
            if (line.length() > 0) {
                fields = split(line, '\t');
            } else {
                if (tagsRead > 1) break;
                else { std::cerr << "The input file seems malformed (empty first line)" << std::endl; exit(1); }
            }
            thisChr = fields[2];
            if (fields.size() < 12) {
                if (tagsRead == 0){
                    std::cerr << "The input file seems malformed (less than 12 columns)" << std::endl; exit(1);
                } else {
                    break; // Assume all tags have been read and these are just the additional lines at the bottom of the file
                }
            }
            fields.erase(fields.begin(), fields.begin()+12);
        } else if (inputType == "Matrix") {
            fields = split(line, '\t');
            thisChr = fields[0];
            fields.erase(fields.begin());
        } else if (inputType == "SimpleMatrix") {
            fields = split(line, '\t');
            thisChr = "";
        }
        
        tagsRead++;
        if (numIndividuals == fields.size() + 1) {
            fields.push_back("");
        }
        
        if (fields.size() != numIndividuals) {
            std::cerr << "fields.size(): " << fields.size() << " numIndividuals:" << numIndividuals << std::endl;
            std::cerr << "tagsRead: " << tagsRead << std::endl;
           // std::cerr << line << std::endl;
        }
        assert(fields.size() == numIndividuals);
        for (int i = 0; i < fields.size(); i++) {
            std::vector<std::string> recipientHaps = split(fields[i], '/');
            if (recipientHaps.size() > 2) {
                bIndTriallelic = true;
            }
        }
        std::regex Ns("N+");
        std::regex NsHet("N+/N+");
        if (!bIndTriallelic) {
            if (!checkIfTagInformative(fields)) {
                notInformative++;
            }
            checkChr(thisChr, chrTagsAndChunks, numIndividuals);
            for (int i = 0; i < fields.size(); i++) {
                if (fields[i] == "" || fields[i] == " " || std::regex_match(fields[i], Ns) || std::regex_match(fields[i], NsHet)) {
                    missingness[i]++; incrementMissingessMatrixBothHaplotypes(missingnessMatrix, i);
                    incrementMissingessMatrixBothHaplotypes(localMissingnessMatrix, i);
                    if (tagsRead < initMissing) {
                        incrementMissingRecipient(outChunksMatrix, i);
                    } else {
                        incrementMissingRecipientOnCoancestry(outChunksMatrix, i);
                    }
                    if (tagsRead % blockSize < 2) {
                        incrementMissingRecipient(localChunkMatrix, i);
                    } else {
                        incrementMissingRecipientOnCoancestry(localChunkMatrix, i);
                    }
                } else {
                    std::vector<std::string> recipientHaps = split(fields[i], '/');
                    std::vector<double> recipientSimVector;
                    if (opt::bTetraploid) {
                        recipientSimVector = calulateConcastryTetraploid(recipientHaps, fields, i, numIndividuals,outChunksMatrix,tagsRead);
                    } else {
                        assert(recipientHaps.size() < 4);
                        // std::cerr << "recipientHaps.size(): " << recipientHaps.size() << std::endl;
                         // Calculate the coancestry values for this tag
                        if (recipientHaps.size() == 2) {
                            string h1 = recipientHaps[0];
                            std::vector<double> recipientSimVectorH1 = calculateSimilarity(fields, h1, i, numIndividuals,outChunksMatrix, tagsRead,2);
                            string h2 = recipientHaps[1];
                            std::vector<double> recipientSimVectorH2 = calculateSimilarity(fields, h2, i, numIndividuals,outChunksMatrix,tagsRead,2);
                            for (int j = 0; j < recipientSimVectorH1.size(); j++) {
                                recipientSimVector.push_back((recipientSimVectorH1[j] + recipientSimVectorH2[j])/2);
                            }
                        } else {
                            string homHap = recipientHaps[0]; // This individual recipient is homozygous
                            recipientSimVector = calculateSimilarity(fields, homHap, i, numIndividuals,outChunksMatrix,tagsRead,1);
                            //std::cerr << "recipientSimVector problem here? " << recipientSimVector[0] << std::endl;
                        }
                    }
                    // Add the coancestry values for this tag to the overall matrix
                    if(!checkSimilaritiesSumToOne(recipientSimVector)) {
                        std::cerr << "recipientSimVector problem: " << recipientSimVector[0] << std::endl;
                        print_vector_stream(recipientSimVector, std::cerr);
                    } assert(checkSimilaritiesSumToOne(recipientSimVector));
                    
                    for (int j = 0; j < recipientSimVector.size(); j++) {
                        outChunksMatrix[i][j] = outChunksMatrix[i][j] + recipientSimVector[j];
                        localChunkMatrix[i][j] = localChunkMatrix[i][j] + recipientSimVector[j];
                    }
                    addChunksToChr(thisChr, chrTagsAndChunks, recipientSimVector, numIndividuals, i);
                }
            }
        }
        if (tagsRead % blockSize == 0) {
            for (int i = 0; i < numIndividuals; i++) {
                for (int j = 0; j < numIndividuals; j++) {
                    double local_p_ij = localChunkMatrix[i][j]/vector_sum(localChunkMatrix[i]);
                    if (tagsRead == blockSize) {
                        bC.matrix[i][j][0] = localChunkMatrix[i][j];
                        bC.p_ij_matrix[i][j][0] = local_p_ij;
                    } else {
                        bC.matrix[i][j].push_back(localChunkMatrix[i][j]);
                        bC.p_ij_matrix[i][j].push_back(local_p_ij);
                    }
                    if (i != j) {
                        if (localChunksNoMissing[i][j] != 0) {
                            // std::cerr << "localChunksNoMissing[i][j] == 0; "  << localChunksNoMissing[i][j] * (1.0/(1-(localMissingnessMatrix[i][j]/blockSize))) << std::endl;
                            if (localMissingnessMatrix[i][j] == blockSize) {
                                std::cerr << "localChunksNoMissing[i][j] = "  << localChunksNoMissing[i][j] << std::endl;
                                std::cerr << "localMissingnessMatrix[i][j] = "  << localMissingnessMatrix[i][j] << std::endl;
                                std::cerr << "tagsRead; [i]; [j] = "  << tagsRead << "; " << i << "; " << j << std::endl;
                                std::cerr << "BUG? Need to sort this out!!! Contact: " << BUGREPORT << std::endl;
                            } else {
                            // std::cerr << "localMissingnessMatrix[i][j] = "  << localMissingnessMatrix[i][j] << std::endl;
                                localChunksNoMissing[i][j] = localChunksNoMissing[i][j] * (1.0/(1-(localMissingnessMatrix[i][j]/blockSize)));
                                localChunksNoMissing[i][j] = (localChunksNoMissing[i][j]/vector_sum(localChunksNoMissing[i])) * blockSize;
                            }
                        }
                        if (tagsRead == blockSize) {
                            bC.matrixNewMissing[i][j][0] = localChunksNoMissing[i][j];
                        } else {
                            bC.matrixNewMissing[i][j].push_back(localChunksNoMissing[i][j]);
                        }
                    } else {
                        if (tagsRead == blockSize) {
                            bC.matrixNewMissing[i][j][0] = 0.0;
                        } else {
                            bC.matrixNewMissing[i][j].push_back(0.0);
                        }
                    }
                }
            }
            // Reset the local matrices to all zeros
            for (int i = 0; i < localChunkMatrix.size(); i++) {
                std::fill(localChunkMatrix[i].begin(), localChunkMatrix[i].end(), 0.0);
                std::fill(localMissingnessMatrix[i].begin(), localMissingnessMatrix[i].end(), 0.0);
                std::fill(localChunksNoMissing[i].begin(), localChunksNoMissing[i].end(), 0.0);
            }
        }
    }
    // Add all the missing data:
    std::vector<std::vector<double> > chunksNoMissingRescaled; initialize_matrix_double(chunksNoMissingRescaled, numIndividuals);
    for (int i = 0; i < missingnessMatrix.size(); i++) {
        // Rescale the outChunksNoMissing matrix so that it reflects relative co-ancestry inferred from observed data
        // independent of missingess
        for (int j = 0; j < missingnessMatrix.size(); j++) {
            if (i != j) {
                chunksNoMissingRescaled[i][j] = outChunksNoMissing[i][j] * (1.0/(1-(missingnessMatrix[i][j]/tagsRead)));
                // old version - BUG!!
                ///chunksNoMissingRescaled[i][j] = outChunksNoMissing[i][j] * (1 + (missingnessMatrix[i][j])/tagsRead);
            }
        }
        for (int j = 0; j < missingnessMatrix.size(); j++) {
            outChunksNoMissing[i][j] = (chunksNoMissingRescaled[i][j]/vector_sum(chunksNoMissingRescaled[i])) * tagsRead;
        }
    }
    std::cerr << individuals[0] << "\t" << vector_sum(chunksNoMissingRescaled[0]) << std::endl;
    std::cerr << individuals[0] << "\t" << vector_sum(outChunksNoMissing[0]) << std::endl;
    
    *outMissingnessMatrixFile << "Recipient" << "\t"; print_vector(individuals, *outMissingnessMatrixFile);
    print_matrix_wNames(missingnessMatrix, *outMissingnessMatrixFile,individuals);
    
    //print_matrix_wNames(chunksNoMissingRescaled, *outChunksMatrixFile,individuals);
    
    // Print missingess:
    for (int i = 0; i < missingness.size(); i++) {
        missingness[i] = missingness[i]/tagsRead;
    }
    std::cerr << "Printing missingness per individual to: " << outMissingnessFileName << std::endl;
    print_vector(individuals, *outMissingnessFile);
    print_vector(missingness, *outMissingnessFile);
    
    // Estimate theoretical variances:
    double R_i = tagsRead/blockSize;
    for (int i = 0; i < numIndividuals; i++) {
        for (int j = 0; j < numIndividuals; j++) {
            if (i != j) {
                double p_ij = outChunksMatrix[i][j]/vector_sum(outChunksMatrix[i]);
                // std::cerr << individuals[i] << "\t" << vector_sum(outChunksMatrix[i]) << std::endl;
                p_ij_full[i][j] = p_ij;
                theoreticalVariancesMatrix[i][j] = (vector_sum(outChunksMatrix[i])*p_ij*(1-p_ij))/R_i;
            }
        }
    }
    
    
    for (int i = 0; i < numIndividuals; i++) {
        for (int j = 0; j < numIndividuals; j++) {
            if (i != j) {
                bC.c_prime_matrix[i][j].resize(bC.matrix[0][0].size(),0.0);
                for (int k = 0; k < bC.matrix[0][0].size(); k++) {
                    bC.c_prime_matrix[i][j][k] = outChunksMatrix[i][j] - bC.matrix[i][j][k];
                    s2_ij[i][j] = s2_ij[i][j] + pow(bC.matrixNewMissing[i][j][k], 2);
                    my_empiricalVar[i][j] = my_empiricalVar[i][j] + pow(bC.matrix[i][j][k] - (vector_sum(bC.matrix[i][j])/R_i), 2);
                }
            }
        }
    }
    
    
    // Estimate empirical variances by jackknife:
    double sumC_ij = 0; double sumC_ij_my = 0;
    double sumTV_ij = 0; double sumEVpaper = 0;
    for (int i = 0; i < numIndividuals; i++) {
        for (int j = 0; j < numIndividuals; j++) {
            if (i != j) {
                
                empiricalVariancesPaper[i][j] = (s2_ij[i][j] - (pow(vector_sum(bC.matrixNewMissing[i][j]), 2)/R_i))/(R_i - 1);
                my_empiricalVar[i][j] = my_empiricalVar[i][j]/(R_i-1);
                //empiricalVariancesMatrix[i][j] = pow(jackknive_std_err(bC.matrix[i][j]),2);
                // empiricalVariancesMatrix[i][j] = pow(jackknive_std_err(bC.c_prime_matrix[i][j]),2);
                empiricalVariancesMatrix[i][j] = pow(jackknive_std_err_sum(bC.matrixNewMissing[i][j]),2)/(R_i - 1);
                c_ij[i][j] = empiricalVariancesMatrix[i][j]/theoreticalVariancesMatrix[i][j];
                sumC_ij = sumC_ij + c_ij[i][j];
                sumC_ij_my = sumC_ij_my + empiricalVariancesPaper[i][j]/theoreticalVariancesMatrix[i][j];
                sumTV_ij = sumTV_ij + theoreticalVariancesMatrix[i][j];
                sumEVpaper = sumEVpaper + empiricalVariancesPaper[i][j];
            }
        }
    }
    // print_vector_stream(bC.matrix[0][1], std::cerr,',');
    // std::cerr << "s2_ij[0][1] = " << s2_ij[0][1] << std::endl;
    // std::cerr << "outChunksMatrix[0][1] = " << outChunksMatrix[0][1] << std::endl;
    // std::cerr << "my_empiricalVar[0][1] = " << my_empiricalVar[0][1] << std::endl;
    // std::cerr << "empiricalVariancesPaper[0][1] = " << empiricalVariancesPaper[0][1] << std::endl;
    // std::cerr << "jackknive_std_err(bC.matrixNewMissing[0][1]) = " << jackknive_std_err(bC.matrix[0][1]) << std::endl;
    // string outTVarName = fileRoot + "_theoreticalVariances.out";
    // string outEmpVarName = fileRoot + "_empiricalVariances.out";
    // std::ofstream* outEmpVarFile = new std::ofstream(outEmpVarName.c_str());
    // std::ofstream* outTVarFile = new std::ofstream(outTVarName.c_str());
    // print_vector(individuals, *outEmpVarFile);
    // print_matrix_wNames(empiricalVariancesMatrix, *outEmpVarFile,individuals);
    // print_vector(individuals, *outTVarFile);
    // print_matrix_wNames(theoreticalVariancesMatrix, *outTVarFile,individuals);
    // std::cerr << "Mean theoretical T_V per block = " << sumTV_ij/(numIndividuals*(numIndividuals-1)) << std::endl;
    
    double jackknifeC = sumC_ij/(numIndividuals*(numIndividuals-1));
    std::cerr << "Theoretical c = " << ploidy * (1.0/(numIndividuals-1)) << std::endl;
    std::cerr << "Jackknife c = " << jackknifeC << std::endl;
    std::cerr << "2012 Manuscript c = " << sumC_ij_my/(numIndividuals*(numIndividuals-1)) << std::endl;
    
    // Print results:
    std::cerr << "notInformative = " << notInformative << std::endl;
    std::cerr << "Printing the final coancestry matrix to " << outChunksMatrixFileName << std::endl;
    *outChunksMatrixFile << "#Cfactor " << jackknifeC << std::endl;
    *outChunksMatrixFile << "Recipient" << "\t"; print_vector(individuals, *outChunksMatrixFile);
    print_matrix_wNames(outChunksNoMissing, *outChunksMatrixFile,individuals);
    
    
    // If requested, print the alternative results with missigng data treated differently
    if (opt::bMissing2) {
        string outChunksMatrixMissingess2FileName = fileRoot + "_missingness2_chunks.out";
        std::ofstream* outChunksMatrixMissingess2File = new std::ofstream(outChunksMatrixMissingess2FileName.c_str());
        *outChunksMatrixMissingess2File << "#Cfactor " << jackknifeC << std::endl;
        *outChunksMatrixMissingess2File << "Recipient" << "\t"; print_vector(individuals, *outChunksMatrixFile);
        print_matrix_wNames(outChunksMatrix, *outChunksMatrixMissingess2File,individuals);
    }

    
    
    // If requested, print the per-chromosome results
    if (opt::bOutputChr) {
        string outChunksMatrixPerChrFileName = fileRoot + "_perChr_chunks.out";
        std::ofstream* outChunksMatrixPerChrFile = new std::ofstream(outChunksMatrixPerChrFileName.c_str());
        std::cerr << "Printing the per chromosome coancestry matrices to " << outChunksMatrixPerChrFileName << std::endl;
        for(std::map<std::string, int>::iterator it = chrTagsAndChunks.numTagsPerChr.begin(); it != chrTagsAndChunks.numTagsPerChr.end(); it++) {
            if (it->second > 50) { // if there are more than 50 tags per this chromosome/scaffold
                *outChunksMatrixPerChrFile << "# " << it->first << std::endl;
                *outChunksMatrixPerChrFile << "# " << it->second << " tags" << std::endl;
                *outChunksMatrixPerChrFile << "Recipient" << "\t";
                print_vector(individuals, *outChunksMatrixPerChrFile);
                print_matrix_wNames(chrTagsAndChunks.chunksMatrixPerChr[it->first], *outChunksMatrixPerChrFile,individuals);
            }
        }
    }
    runtime = ( std::clock() - runStart ) / (double) CLOCKS_PER_SEC;
    std::cout << "Analysis completed in: " << runtime << " seconds (" << runtime/tagsRead << " seconds per RAD locus)" << std::endl;
    return 0;
}



void parsePaintSqlOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
                case 'n': arg >> opt::runName; break;
                case 'c': opt::bOutputChr = true; break;
                case 't': opt::bTetraploid = true; break;
                case 'm': opt::bMissing2 = true; break;
                case '?': die = true; break;
                case 'h': std::cout << PAINTSQL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind == 0)
    {
        std::cerr << "you need to specify an input file\n";
        die = true;
    }
    if (argc - optind > 1)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << PAINTSQL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::sqlFileName = argv[optind++];
}
