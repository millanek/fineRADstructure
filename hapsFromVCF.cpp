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

static int maxLocusGapBp = 1000;

static const char *HAPS_USAGE_MESSAGE =
"Usage: " BIN " " SUBPROGRAM " [OPTIONS] INPUT.vcf\n"
"Generate a co-ancestry matrix from RAD data\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -H,   --het-treatment <r|p>             r: assign the first het base randomly, subsequent as Ns (default); p: use the phase information in the VCF\n"
"       -m,   --maxNproportion <0.5>            maximum proportion of N bases at a locus before it is set to missing entirely (default = 0.5)\n"
"       -p,   --printNproportion FILE           print the proportion of Ns for each individual at each locus\n"
"       -F MIN_F                                minimum acceptable inbreeding coefficient (default: F >= -0.3)\n"
"\n\n"
"\nReport bugs to " BUGREPORT "\n\n";


// Options
static const char* shortopts = "hn:H:F:m:p:";
static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "het-treatment",   required_argument, NULL, 'H' },
    { "maxNproportion",   required_argument, NULL, 'm' },
    { "printNproportion",   required_argument, NULL, 'p' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string runName = "";
    static string VCFfileName = "";
    static string NproportionFileName = "";
    static char hetTreatment = 'r';
    static bool bOutputChr = false;
    static double minF = -0.3;
    static double maxNs = 0.5;
}

class VCFprocessCounts {
public:
    VCFprocessCounts(size_t ns) : numLoci(0), missingDueToUnphasedHets(0), missingLociNum(0),
                                processedVariantCounter(0), usedVariantCounter(0), missingDueToTooManyNs(0) {
        numSamples = ns;
        hetCounters.resize(numSamples);
        Ncounters.resize(numSamples);
    };
    
    int numLoci; int missingDueToUnphasedHets; int missingLociNum;
    int missingDueToTooManyNs;
    unsigned int processedVariantCounter; unsigned int usedVariantCounter;
    size_t numSamples;
    std::vector<int> hetCounters;
    std::vector<int> Ncounters;
};

class Alleles {
public:
    Alleles(size_t numSamples) {
        H1.resize(numSamples); H2.resize(numSamples);
    };
    
    std::vector<string> H1; std::vector<string> H2;
};



void printAllelesIfNotMissing(Alleles* alleles, VCFprocessCounts* counts, std::regex NsHet, std::ofstream* nProportionOutFile) {
    int numAlleles = (int)alleles->H1.size(); std::vector<string> allelesH1H2;
    int numMissingAlleles = 0;
    std::vector<double> Nproportions;
    for (std::vector<string>::size_type i = 0; i != numAlleles; i++) {
        if (alleles->H1[i] != "") {
            double Nproportion = (double)counts->Ncounters[i]/(alleles->H1[i].length()*2);
            Nproportions.push_back(Nproportion);
            allelesH1H2.push_back(alleles->H1[i] + "/" + alleles->H2[i]);
            if (std::regex_match(allelesH1H2[i], NsHet)) {
                allelesH1H2[i] = ""; counts->missingLociNum++;
                numMissingAlleles++;
            } else if(Nproportion > opt::maxNs) {
                //std::cerr << "allelesH1H2[i]: " << allelesH1H2[i] << std::endl;
                //std::cerr << "counts->Ncounters[i]: " << counts->Ncounters[i] << std::endl;
                //std::cerr << "(0.5*alleles->H1[i].length()*2): " << (0.5*alleles->H1[i].length()*2) << std::endl;
                //exit(1);
                allelesH1H2[i] = ""; counts->missingDueToTooManyNs++;
                numMissingAlleles++;
            }
            //} else if(counts->hetCounters[i] > 1 && opt::hetTreatment == 'r') {
            //    allelesH1H2[i] = ""; counts->missingDueToUnphasedHets++;
            //    numMissingAlleles++;
            // }
        } else {
            allelesH1H2.push_back("");
            numMissingAlleles++;
        }
    }
    if (numAlleles > numMissingAlleles) {
        print_vector_stream(allelesH1H2,std::cout);
        if (opt::NproportionFileName != "") {
            print_vector_stream(Nproportions,*nProportionOutFile);
        }
    }
    counts->numLoci = counts->numLoci + (int)counts->numSamples;
}


int VCFhapsMain(int argc, char** argv) {
    string line;
    parseVCFoptions(argc, argv);
    string vcfFileRoot = stripExtension(opt::VCFfileName);
    
    std::cerr << "Generating a haplotype file for RADpainter with variants from: " << opt::VCFfileName << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::istream* vcfFile = createReader(opt::VCFfileName.c_str());
    
    std::ofstream* nProportionOutFile;
    if (opt::NproportionFileName != "") {
        nProportionOutFile = new std::ofstream(opt::NproportionFileName);
    }
    
    string currentScaffoldNum = "";
    int currentCoord = 0;
    std::vector<string> sampleNames;
    Alleles* alleles = nullptr; VCFprocessCounts* counts = nullptr;
    std::regex NsHet("N+/N+");
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            size_t ns = fields.size()-NUM_NON_GENOTYPE_COLUMNS; std::cerr << "numSamples: " << ns << std::endl;
            
            alleles = new Alleles(ns); counts = new VCFprocessCounts(ns);
            // Initialize vectors
            for (std::vector<string>::size_type i = 0; i != counts->numSamples; i++) {
                alleles->H1[i] = ""; alleles->H2[i] = ""; counts->hetCounters[i] = 0; counts->Ncounters[i] = 0;
            }
            // Get sample names
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                sampleNames.push_back(fields[i]);
            }
            print_vector_stream(sampleNames, std::cout);
            if (opt::NproportionFileName != "") {
                print_vector_stream(sampleNames,*nProportionOutFile);
            }
            
        } else {
            counts->processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (fields[0] != currentScaffoldNum || atoi(fields[1].c_str()) - maxLocusGapBp > currentCoord) {
                if (currentScaffoldNum != "") {
                   // std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << std::endl;
                    printAllelesIfNotMissing(alleles,counts,NsHet,nProportionOutFile);
                    for (std::vector<string>::size_type i = 0; i != counts->numSamples; i++) {
                        alleles->H1[i] = ""; alleles->H2[i] = ""; counts->hetCounters[i] = 0; counts->Ncounters[i] = 0;
                    }
                }
                // get to the next "chromosome"
                counts->processedVariantCounter = 1; counts->usedVariantCounter = 0;
                currentScaffoldNum = fields[0]; currentCoord = atoi(fields[1].c_str());
            }
            
            if (info[0] != "INDEL" && fields[3].length() == 1 && fields[4].length() == 1
                && isDNAonly(fields[3][0]) && isDNAonly(fields[4][0])) {
                //int lengthToAppend = (atoi(fields[1].c_str()) - 1) - (int)inStrPos;
                // make sure the length is non-negative (can happen
                // if two consecutive variants have the same coordinate)
                // for now we just ignore the additional variant
               // if (lengthToAppend >= 0) {
                std::vector<int> appendVectorInt(counts->numSamples,0);
                std::vector<std::string> appendVector(counts->numSamples,"0");
                std::vector<std::string> appendVectorH2(counts->numSamples,"0");
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
                    if (F < opt::minF) {
                        std::cerr << "Strongly negative inbreeding coefficient; Variant: " << fields[0] << "\t" << fields[1] << "\tF = " << F << std::endl;
                    } else {
                        if (opt::hetTreatment == 'r') {
                            for (std::vector<std::string>::size_type i = 0; i != counts->numSamples; i++) {
                                if (appendVectorInt[i] == 1) {
                                    counts->hetCounters[i]++;
                                }
                                if (counts->hetCounters[i] > 1 && appendVectorInt[i] == 1) {
                                    //std::cerr << "More than one unphased het in variant: " << fields[0] << "\t" << fields[1] << "\tindividual = " << sampleNames[i] << std::endl;
                                    //std::cerr << "hetCounters[i]: " << hetCounters[i] << "\tusedVariantCounter " << usedVariantCounter << std::endl;
                                    appendVector[i] = "N";
                                    appendVectorH2[i] = "N";
                                }
                            }
                        }

                        for (std::vector<std::string>::size_type i = 0; i != counts->numSamples; i++) {
                            alleles->H1[i].append(appendVector[i]); alleles->H2[i].append(appendVectorH2[i]);
                            if (appendVector[i] == "N") {
                                counts->Ncounters[i]++;
                            }
                            if (appendVectorH2[i] == "N") {
                                counts->Ncounters[i]++;
                            }
                            
                        }
                        counts->usedVariantCounter++;
                        currentCoord = atoi(fields[1].c_str());
                    }
                }
               // }
                
                // Double check that all looks OK
                if((int)alleles->H1[0].length() != counts->usedVariantCounter) {
                    std::cerr << "usedVariantCounter: " << counts->usedVariantCounter << std::endl;
                    std::cerr << "scaffoldStrings[0].length(): " << alleles->H1[0].length() << std::endl;
                    std::cerr << "scaffoldStrings[1].length(): " << alleles->H1[1].length() << std::endl;
                    std::cerr << fields[3] << " " << fields[4] << std::endl;
                    std::cerr << "scaffoldStrings[0]: " << alleles->H1[0] << std::endl;
                    std::cerr << "scaffoldStrings[1]: " << alleles->H1[1] << std::endl;
                }
                assert((int)alleles->H1[0].length() == counts->usedVariantCounter);
                
            }
        }
    }
    // Also the final chr
  //  std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << std::endl;
    printAllelesIfNotMissing(alleles,counts,NsHet,nProportionOutFile);
    
    int missingTotal = counts->missingLociNum + counts->missingDueToUnphasedHets + counts->missingDueToTooManyNs;
    std::cerr << "Missing alleles total: " << missingTotal << "; " << (double)missingTotal/counts->numLoci << "% of total alleles"<< std::endl;
    std::cerr << "Due to missing genotypes: " << counts->missingLociNum << "; " << (double)counts->missingLociNum/missingTotal << "% of missing" << std::endl;
    std::cerr << "Due to too many Ns in the individual: " << counts->missingDueToTooManyNs << "; " << (double)counts->missingDueToTooManyNs/missingTotal << "% of missing" << std::endl;
    std::cerr << "Due to multiple unphased hets: " << counts->missingDueToUnphasedHets << "; " << (double)counts->missingDueToUnphasedHets/missingTotal << "% of missing" << std::endl;
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
            case 'F': arg >> opt::minF; break;
            case 'm': arg >> opt::maxNs; break;
            case 'p': arg >> opt::NproportionFileName; break;
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
