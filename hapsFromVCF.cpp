//
//  hapsFromVCF.cpp
//  fineRADstructure
//
//  Created by Milan Malinsky on 05/03/2018.
//  Copyright © 2018 Milan Malinsky. All rights reserved.
//

#include "hapsFromVCF.h"
#include "utils.h"

#define SUBPROGRAM "hapsFromVCF"

static const char *HAPS_USAGE_MESSAGE =
"Usage: " BIN " " SUBPROGRAM " [OPTIONS] INPUT.vcf\n"
"Generate a co-ancestry matrix from RAD data\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -H,   --het-treatment <r|p>             r: assign het bases randomly (default); p: use the phase information in the VCF\n"
"\n\n"
"\nReport bugs to " BUGREPORT "\n\n";


// Options
static const char* shortopts = "hn:H:";
static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "het-treatment",   required_argument, NULL, 'H' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string runName = "";
    static string VCFfileName = "";
    static char hetTreatment = 'r';
    static bool bOutputChr = false;
}

int VCFhapsMain(int argc, char** argv) {
    string line;
    parseVCFoptions(argc, argv);
    string vcfFileRoot = stripExtension(opt::VCFfileName);
    
    std::cerr << "Generating a haplotype file for RADpainter with variants from: " << opt::VCFfileName << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::istream* vcfFile = createReader(opt::VCFfileName.c_str());
    
    string currentScaffoldNum = "";
    size_t numSamples;
    std::vector<string> sampleNames;
    std::vector<string> scaffoldStrings; std::vector<string> scaffoldStringsH2;
    unsigned int processedVariantCounter = 0; unsigned int usedVariantCounter = 0;
    std::vector<int> hetCounters; bool skipRestDueToUnphasedHets = false;
    int numLoci = 0;
    int missingLociNum = 0; int missingDueToUnphasedHets = 0;
    std::regex NsHet("N+/N+");
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            numSamples = fields.size()-NUM_NON_GENOTYPE_COLUMNS;
            std::cerr << "numSamples: " << numSamples << std::endl;
            // Initialize vectors
            scaffoldStrings.resize(numSamples); scaffoldStringsH2.resize(numSamples); hetCounters.resize(numSamples);
            for (std::vector<string>::size_type i = 0; i != scaffoldStrings.size(); i++) {
                scaffoldStrings[i] = ""; scaffoldStringsH2[i] = ""; hetCounters[i] = 0;
            }
            // Open get sample names
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                sampleNames.push_back(fields[i]);
            }
            print_vector_stream(sampleNames, std::cout);
        } else {
            processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (fields[0] != currentScaffoldNum) {
                if (currentScaffoldNum != "") {
                   // std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << std::endl;
                    for (std::vector<string>::size_type i = 0; i != scaffoldStrings.size(); i++) {
                        if (scaffoldStrings[i] != "") {
                            scaffoldStrings[i] = scaffoldStrings[i] + "/" + scaffoldStringsH2[i];
                        }
                        if (std::regex_match(scaffoldStrings[i], NsHet)) {
                            scaffoldStrings[i] = ""; missingLociNum++;
                        }
                        if(hetCounters[i] > 1 && opt::hetTreatment == 'r') {
                            scaffoldStrings[i] = ""; missingDueToUnphasedHets++;
                        }
                    }
                    print_vector_stream(scaffoldStrings,std::cout); numLoci = numLoci + (int)numSamples;
                    for (std::vector<string>::size_type i = 0; i != scaffoldStrings.size(); i++) {
                        scaffoldStrings[i] = ""; scaffoldStringsH2[i] = ""; hetCounters[i] = 0;
                    }
                }
                // get to the next "chromosome"
                processedVariantCounter = 1; usedVariantCounter = 0;
                currentScaffoldNum = fields[0];
                skipRestDueToUnphasedHets = false;
            }
            if (info[0] != "INDEL") {
                //int lengthToAppend = (atoi(fields[1].c_str()) - 1) - (int)inStrPos;
                // make sure the length is non-negative (can happen
                // if two consecutive variants have the same coordinate)
                // for now we just ignore the additional variant
               // if (lengthToAppend >= 0) {
                std::vector<int> appendVectorInt(numSamples,0);
                std::vector<std::string> appendVector(numSamples,"0");
                std::vector<std::string> appendVectorH2(numSamples,"0");
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    //std::cerr << "Going through genotypes1:" << i << std::endl;
                    //std::cerr << scaffoldStrings.size() << " " << inStrPos << " " << fields[1] << " " << currentScaffoldReference.size() << std::endl;
                    std::vector<string> genotypeFields = split(fields[i], ':');
                    std::vector<char> genotype; genotype.push_back(genotypeFields[0][0]); genotype.push_back(genotypeFields[0][2]);
                    std::vector<std::string> genotypeAndZeroOne = returnGenotypeBasesAndZeroOneTwo(fields[3], fields[4], genotype, opt::hetTreatment);
                    appendVector[i- NUM_NON_GENOTYPE_COLUMNS] = genotypeAndZeroOne[0];
                    appendVectorH2[i- NUM_NON_GENOTYPE_COLUMNS] = genotypeAndZeroOne[1];
                    appendVectorInt[i- NUM_NON_GENOTYPE_COLUMNS] = (int)stringToDouble(genotypeAndZeroOne[2].c_str());
                   /* if (fields[0] == "MC00000568" && i- NUM_NON_GENOTYPE_COLUMNS == 0) {
                        std::cerr << "genotypeAndZeroOne[0] " << genotypeAndZeroOne[0] << "\tgenotypeAndZeroOne[1] " << genotypeAndZeroOne[1] << "\tgenotypeFields[0][0] = " << genotypeFields[0][0] << std::endl;
                    } */
                }
                
                if(vector_sum(appendVectorInt) > 0) {
                    double F = calculateInbreedingCoefficient(appendVectorInt);
                    if (F < -0.3) {
                        std::cerr << "Strongly negative inbreeding coefficient; Variant: " << fields[0] << "\t" << fields[1] << "\tF = " << F << std::endl;
                    } else {
                        if (opt::hetTreatment == 'r') {
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                if (appendVectorInt[i] == 1) {
                                    hetCounters[i]++;
                                }
                                if (hetCounters[i] > 1 && appendVectorInt[i] == 1) {
                                    //std::cerr << "More than one unphased het in variant: " << fields[0] << "\t" << fields[1] << "\tindividual = " << sampleNames[i] << std::endl;
                                    //std::cerr << "hetCounters[i]: " << hetCounters[i] << "\tusedVariantCounter " << usedVariantCounter << std::endl;
                                    appendVector[i] = "N";
                                    appendVectorH2[i] = "N";
                                    //skipRestDueToUnphasedHets = true;
                                }
                            }
                        }

                        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                            scaffoldStrings[i].append(appendVector[i]); scaffoldStringsH2[i].append(appendVectorH2[i]);
                        }
                        usedVariantCounter++;
                    }
                }
               // }
                
                // Double check that all looks OK
                if((int)scaffoldStrings[0].length() != usedVariantCounter) {
                    std::cerr << "usedVariantCounter: " << usedVariantCounter << std::endl;
                    std::cerr << "scaffoldStrings[0].length(): " << scaffoldStrings[0].length() << std::endl;
                    std::cerr << "scaffoldStrings[1].length(): " << scaffoldStrings[1].length() << std::endl;
                    std::cerr << fields[3] << " " << fields[4] << std::endl;
                    std::cerr << "scaffoldStrings[0]: " << scaffoldStrings[0] << std::endl;
                    std::cerr << "scaffoldStrings[1]: " << scaffoldStrings[1] << std::endl;
                }
                assert((int)scaffoldStrings[0].length() == usedVariantCounter);
                
            }
        }
    }
    // Also the final chr
  //  std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << std::endl;
    for (std::vector<string>::size_type i = 0; i != scaffoldStrings.size(); i++) {
        scaffoldStrings[i] = scaffoldStrings[i] + "/" + scaffoldStringsH2[i];
    }
    print_vector_stream(scaffoldStrings,std::cout); numLoci = numLoci + (int)numSamples;
    std::cerr << "Missing loci total: " << missingLociNum+missingDueToUnphasedHets << "; " << (double)(missingLociNum+missingDueToUnphasedHets)/numLoci << "% of total loci"<< std::endl;
    std::cerr << "Due to missing genotypes: " << missingLociNum << "; " << (double)missingLociNum/(missingLociNum+missingDueToUnphasedHets) << "% of missing" << std::endl;
    std::cerr << "Due to multiple unphased hets: " << missingDueToUnphasedHets << "; " << (double)missingDueToUnphasedHets/(missingLociNum+missingDueToUnphasedHets) << "% of missing" << std::endl;
    std::cerr << "DONE!" << std::endl;
    return 0;
}


void parseVCFoptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'n': arg >> opt::runName; break;
            case '?': die = true; break;
            case 'h': std::cout << HAPS_USAGE_MESSAGE;
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
        std::cout << "\n" << HAPS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::VCFfileName = argv[optind++];
}
