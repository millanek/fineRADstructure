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

static const char *PAINTSQL_USAGE_MESSAGE =
"Usage: " BIN " " SUBPROGRAM " [OPTIONS] INPUT.txt\n"
"Generate a co-ancestry matrix from RAD data\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -t, --tetraploid                        the tetraploid mode (see manual for input format)\n"
"       -c, --chr                               output per-chromosome coancestry matrices\n"
"       -n, --run-name                          run-name will be included in the output file name(s)\n\n"
"\n\n"
"\nReport bugs to " BUGREPORT "\n\n";

static const char* shortopts = "hn:";

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "chr",   no_argument, NULL, 'c' },
    { "tetraploid",   no_argument, NULL, 't' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string runName = "";
    static string sqlFileName = "";
    static bool bTetraploid = false;
    static bool bOutputChr = false;
}

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

// Divide similarity equally between samples that are the closest - to equal 1 in total
// For a missing donor, assign similarity 1/(N-1) (N the number of samples), and take a corresponding amount away from all the closest ones
std::vector<double> calculateSimilarity(const std::vector<std::string>& allHaps, const std::string& recipientHap, int thisIndI, int nSamples) {
    if (std::regex_match(recipientHap, std::regex("N+"))) {
	//std::cerr << "recipientHap is N: " << recipientHap << std::endl;
        return getMissingRecipientSimVector(nSamples, thisIndI);
    }
    std::vector<double> diffVector(nSamples*2, 0.0); // Proportion of difference between all haplotypes and the recipient
    std::vector<double> simVectorPerInd(nSamples, 0.0);
    double numMissing = 0;
    for (int i = 0; i < allHaps.size(); i++) {
        if (i != thisIndI) {
            if (allHaps[i] == "" || std::regex_match(allHaps[i], std::regex("N+"))) {
                numMissing++;
                diffVector[2*i] = 1; diffVector[(2*i)+1] = 1;
                simVectorPerInd[i] = 1/(double)(nSamples-1);   // Similarity between the recipent and this donor is 1/(N-1)
            } else {
                std::vector<std::string> donorHaps = split(allHaps[i], '/');
                assert(donorHaps.size() < 3);
                if (donorHaps.size() == 2) {
                    if (std::regex_match(donorHaps[0], std::regex("N+"))) {
                        numMissing = numMissing + 0.5;
                        diffVector[2*i] = 0.25;
                        simVectorPerInd[i] =+ (1/(double)(nSamples-1))/2;
                    } else {
                        diffVector[2*i] = compareSeqs(recipientHap, donorHaps[0]);
                    }
                    if (std::regex_match(donorHaps[1], std::regex("N+"))) {
                        numMissing = numMissing + 0.5;
                        // std::cerr << "numMissing: " << numMissing << std::endl;
                        diffVector[(2*i)+1] = 0.25;
                        simVectorPerInd[i] =+ (1/(double)(nSamples-1))/2;
                    } else {
                        diffVector[(2*i)+1] = compareSeqs(recipientHap, donorHaps[1]);
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
    //std::cerr << "thisIndI: " << thisIndI << std::endl;
    //std::cerr << "numClosest: " << numClosest << std::endl;
    std::vector<int> indicesMin = findMinDiffIndices(diffVector, minDiff);
    assert(numClosest == indicesMin.size());
    //std::cerr << "numMissing: " << numMissing << std::endl;
    // Fill in the per-individual co-ancestry vector
    double totalToMissing = numMissing/(double)(nSamples-1);
    //std::cerr << "totalToMissing: " << totalToMissing << std::endl;
    double subtractFromClosest = totalToMissing/(double)numClosest; // Take this away from the closest (assigned to missing)
    double add = (1/(double)numClosest) - subtractFromClosest;
    //std::cerr << "add: " << add << std::endl;
    for (int i = 0; i < simVectorPerInd.size(); i++) {
        if (std::find(indicesMin.begin(), indicesMin.end(), 2*i) != indicesMin.end()) {
            simVectorPerInd[i] = simVectorPerInd[i] + add;
        }
        if (std::find(indicesMin.begin(), indicesMin.end(), (2*i)+1) != indicesMin.end()) {
            simVectorPerInd[i] = simVectorPerInd[i] + add;
        }
    }
    if(checkSimilaritiesSumToOne(simVectorPerInd) == false) {
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

void checkFieldsStacks(const std::vector<std::string>& fields, int tagsRead) {
    if (fields.size() < 12) {
        if (tagsRead == 0){
            std::cerr << "The input file seems malformed (less than 12 columns)" << std::endl; exit(1);
        } else {
            exit(0); // Assume all tags have been read and these are just the additional lines at the bottom of the file
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

std::vector<double> calulateConcastryTetraploid(const std::vector<std::string>& recipientHaps,const std::vector<std::string>& fields, int i, int numIndividuals) {
    assert(recipientHaps.size() == 4); // Handle make sure this is a correctly formatted tatraploid site
    std::vector<double> recipientSimVector;
    string h1 = recipientHaps[0];
    std::vector<double> recipientSimVectorH1 = calculateSimilarity(fields, h1, i, numIndividuals);
    string h2 = recipientHaps[1];
    std::vector<double> recipientSimVectorH2 = calculateSimilarity(fields, h2, i, numIndividuals);
    string h3 = recipientHaps[2];
    std::vector<double> recipientSimVectorH3 = calculateSimilarity(fields, h3, i, numIndividuals);
    string h4 = recipientHaps[3];
    std::vector<double> recipientSimVectorH4 = calculateSimilarity(fields, h4, i, numIndividuals);
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


int paintSqlMain(int argc, char** argv) {
    parsePaintSqlOptions(argc, argv);
    string inputType;
    
    std::cerr << "Painting RAD tags from: " << opt::sqlFileName << std::endl;
    std::ifstream* sqlFile = new std::ifstream(opt::sqlFileName.c_str());
    string fileRoot = stripExtension(opt::sqlFileName);
    
    // Create output files:
    string outChunksMatrixFileName = fileRoot + "_chunks.out";
    std::ofstream* outChunksMatrixFile = new std::ofstream(outChunksMatrixFileName.c_str());
    string outChunksMatrixPerChrFileName = fileRoot + "_perChr_chunks.out";
    std::ofstream* outChunksMatrixPerChrFile = new std::ofstream(outChunksMatrixPerChrFileName.c_str());

    std::vector<std::vector<double> > outChunksMatrix;
    double runtime; double runStart = clock();
    
    // Some variables for c calculations
    string thisChr; AllPerChrData chrTagsAndChunks;
    
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
    
    //std::cerr << "Number of columns: " << fields[11] << std::endl;
    // Lets start going through the tags:
    while (getline(*sqlFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows

        bool bIndTriallelic = false;
        
        if (tagsRead % 100 == 0 && tagsRead > 0) {
            std::cerr << "Processed: " << tagsRead << " tag loci" << std::endl;
        }
        std::vector<std::string> fields = split(line, '\t');
        if (inputType == "Stacks") {
            thisChr = fields[2];
            checkFieldsStacks(fields,tagsRead);
            fields.erase(fields.begin(), fields.begin()+12);
        } else if (inputType == "Matrix") {
            thisChr = fields[0];
            fields.erase(fields.begin());
        } else if (inputType == "SimpleMatrix") {
            thisChr = "";
        }
        
        tagsRead++;
        if (numIndividuals == fields.size() + 1) {
            fields.push_back("");
        }
        // std::cerr << "fields.size(): " << fields.size() << " numIndividuals:" << numIndividuals << std::endl;
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
            checkChr(thisChr, chrTagsAndChunks, numIndividuals);
            for (int i = 0; i < fields.size(); i++) {
                if (fields[i] == "" || std::regex_match(fields[i], Ns) || std::regex_match(fields[i], NsHet)) {
                    incrementMissingRecipient(outChunksMatrix, i);
                } else {
                    std::vector<std::string> recipientHaps = split(fields[i], '/');
                    std::vector<double> recipientSimVector;
                    if (opt::bTetraploid) {
                        recipientSimVector = calulateConcastryTetraploid(recipientHaps, fields, i, numIndividuals);
                    } else {
                        assert(recipientHaps.size() < 4);
                        // std::cerr << "recipientHaps.size(): " << recipientHaps.size() << std::endl;
                         // Calculate the coancestry values for this tag
                        if (recipientHaps.size() == 2) {
                            string h1 = recipientHaps[0];
                            std::vector<double> recipientSimVectorH1 = calculateSimilarity(fields, h1, i, numIndividuals);
                            string h2 = recipientHaps[1];
                            std::vector<double> recipientSimVectorH2 = calculateSimilarity(fields, h2, i, numIndividuals);
                            for (int j = 0; j < recipientSimVectorH1.size(); j++) {
                                recipientSimVector.push_back((recipientSimVectorH1[j] + recipientSimVectorH2[j])/2);
                            }
                        } else {
                            string homHap = recipientHaps[0]; // This individual recipient is homozygous
                            recipientSimVector = calculateSimilarity(fields, homHap, i, numIndividuals);
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
                    }
                    addChunksToChr(thisChr, chrTagsAndChunks, recipientSimVector, numIndividuals, i);
                }
            }
        }
    }
    std::cerr << "Printing the final coancestry matrix to " << outChunksMatrixFileName << std::endl;
    *outChunksMatrixFile << "#Cfactor " << (1.0/(numIndividuals-1))*((tagsRead/5000)+1) << std::endl;
    *outChunksMatrixFile << "Recipient" << "\t"; print_vector(individuals, *outChunksMatrixFile);
    print_matrix_wNames(outChunksMatrix, *outChunksMatrixFile,individuals);
    
    if (opt::bOutputChr) {
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
