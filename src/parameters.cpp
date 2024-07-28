/**
 * Parameters.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *
 */

// Includes from STL:
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;
#include <algorithm> // std::sort
using std::sort;
#include <assert.h>
#include <memory>
using std::shared_ptr;
#include <math.h>

// YAML Parser
#include "yaml-cpp/yaml.h"

// Includes from our files:
#include "typedefs.h"
#include "parameters.h"
#include "chromosome.h"
#include "ran_mk.h"

//======================================================================

// Constants
#define NUM_RUNS "nRuns"
#define KINGMAN "kingman"
#define DRIFT "drift"
#define MS_OUTPUT "msOutput"
#define POP_SIZES "popSizes"
#define INVERSION_FREQS "inversionFreqs"
#define SPECIATION "speciation"
#define SIMULATE "simulate"
#define TIME "time"
#define INVERSION_FREQ "inversionFreq"
#define DEMOGRAPHY_CHANGE "demographyChange"
#define COEFFICIENT "coefficient"
#define INVERSION_AGE "inversionAge"
#define MIG_RATE "migRate"
#define BP "bp"
#define IS_PHI_RANDOM "isPhiRand"
#define PHI "phi"
#define INVERSION_BREAKPOINTS "inversionBreakpoints"
#define LEFT "L"
#define RIGHT "R"
#define SNP "SNP"
#define IS_FIXED "isFixed"
#define NUM_BASES "nBases"
#define MUTATION_RATE "mutationRate"
#define RANDOM_POS "randPos"
#define RANGE "range"
#define BEGIN "begin"
#define END "end"
#define CARRIERS "carriers"
#define RANDOM "random"
#define DATA "data"

Parameters::Parameters(const char *insstring)
{ // Constructor.  Reads parameter values from the file parameters.txt

    DBG("Parameters:: Constructing...")

    paramData.reset(new ParameterData);

    ifstream fin(insstring);
    string line;
    int temp;
    double doubtemp;

    // Read in different params
    //  number of runs
    getline(fin, line);
    if (!fin)
    {
        std::cerr << "No pars file. \n";
        exit(1);
    }

    istringstream iss(line);
    while (iss >> temp)
    {
        paramData->nRuns = temp;
    }
    std::cerr << "Number of runs: " << paramData->nRuns << '\n';

    // Should we simulate kingman coalescent? (i.e. exponential waiting times). If false, does generation by generation
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->kingman = (bool)temp;
    }
    std::cerr << "Assuming Kingman? " << paramData->kingman << '\n';

    // Simulate drift trajectory?
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->drift = (bool)temp;
    }
    std::cerr << "Simulating drift? " << paramData->drift << '\n';

    if (paramData->drift)
    {
        std::cerr << "Sorry, drift not implemented in this version yet\n";
        exit(1);
    }

    // Output ms-formatted haplotypes?
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->msOutput = (bool)temp;
    }
    std::cerr << "Output each set in ms format? " << paramData->msOutput << '\n';

    // Population Sizes
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->popSizeVec.push_back(temp);
    }
    std::cerr << "Population Sizes: ";
    for (int n = 0; n < (int)paramData->popSizeVec.size(); n++)
    {
        std::cerr << (n) << " = " << paramData->popSizeVec.at(n) << " | ";
    }
    std::cerr << '\n';
    paramData->totalPopSize = std::accumulate(paramData->popSizeVec.begin(), paramData->popSizeVec.end(), 0);

    unsigned int pops = paramData->popSizeVec.size();

    // Inversion initial frequencies
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        paramData->initialFreqs.push_back(doubtemp);
    }

    if (paramData->initialFreqs.size() != pops)
    {
        std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
        exit(1);
    }

    std::cerr << "Frequencies of inverted chromosomes: ";
    for (int n = 0; n < (int)paramData->initialFreqs.size(); n++)
    {
        std::cerr << (n) << " = " << paramData->initialFreqs.at(n) << " | ";
    }
    std::cerr << '\n';

    // Simulate speciation?
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        paramData->speciation.push_back(doubtemp);
    }
    if (paramData->speciation.at(0) == 1)
    {
        if (paramData->speciation.size() < 3)
        {
            std::cerr << "Error; Simulating speciation but no time or invFreq parameters found\n";
            exit(1);
        }
        std::cerr << "Speciation time: " << paramData->speciation.at(1) << "\n";
    }
    else
    {
        std::cerr << "No speciation.\n";
    }

    // Simulate demography?
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        paramData->demography.push_back(doubtemp);
    }
    if (paramData->demography.at(0) == 1)
    {
        if (paramData->speciation.size() < 3)
        {
            std::cerr << "Error; Simulating demography but no time parameter found\n";
            exit(1);
        }
        std::cerr << "Demography change time: " << paramData->demography.at(1) << "\n";
        std::cerr << "Reminder: this coeff multiplies all pop sizes in the preceding epoch\n";
    }
    else
    {
        std::cerr << "No demographic event.\n";
    }

    // Age of inversion in number of generations (EX: 1500)
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->inv_age = temp;
    }
    std::cerr << "Age of inversion = " << paramData->inv_age << '\n';

    if (paramData->inv_age > 0 && paramData->inv_age < std::max(paramData->demography.at(1), paramData->speciation.at(1)))
    {
        std::cerr << "Error: Current implementations only allows for inversions that predate speciation and demographic changes\n";
        exit(1);
    }

    // Migration rate (4Nm)
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        paramData->migRate.push_back(doubtemp);
    }
    std::cerr << "Migration rate (pop0<-->1) = " << paramData->migRate[0] << '\n';

    // bases per Morgan in the hh genotype.
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        paramData->BasesPerMorgan = doubtemp;
    }
    std::cerr << "Bases per Morgan in homokaryotypic recombination " << paramData->BasesPerMorgan << '\n';

    // Random phi?
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->randPhi = (bool)temp;
    }
    std::cerr << "Random phi values? (if 1, range below. if 0, single value read below) " << paramData->randPhi << '\n';

    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        paramData->phi_range.push_back(doubtemp);
    }
    if (paramData->randPhi)
    {
        std::cerr << "Random gene flux (phi) range = " << paramData->phi_range[0] << " - " << paramData->phi_range[1] << '\n';
    }
    else
    {
        std::cerr << "Gene flux (phi) = " << paramData->phi_range[0] << '\n';
    }

    // the range of the inversion
    vector<double> invtemp;
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        invtemp.push_back(doubtemp);
    }
    paramData->invRange.L = invtemp[0] / paramData->BasesPerMorgan;
    paramData->invRange.R = invtemp[1] / paramData->BasesPerMorgan;
    std::cerr << "Inversion from: " << invtemp[0] << " to " << invtemp[1] << " (" << paramData->invRange.L << " - " << paramData->invRange.R << " recUnits)\n";

    // Fixed S? Value of S or sequence length (bases), Theta
    vector<double> stemp;
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> doubtemp)
    {
        stemp.push_back(doubtemp);
    }
    paramData->fixedS = (bool)stemp[0];
    paramData->n_SNPs = (int)stemp[1];
    if (paramData->fixedS)
    {
        paramData->theta = 0;
        std::cerr << "Number of markers to simulate: " << paramData->n_SNPs << '\n';
    }
    else
    {
        paramData->theta = stemp[2];
        std::cerr << "Number of bases (non-recombining) to simulate: " << paramData->n_SNPs << ", with mutation rate: " << paramData->theta << '\n';
    }

    // random positions of SNPs?
    getline(fin, line);
    iss.clear();
    iss.str(line);
    while (iss >> temp)
    {
        paramData->randSNP = (bool)temp;
    }
    std::cerr << "Markers in random locations? " << paramData->randSNP << '\n';

    // // range (positions of begin-end) where the SNPs are -- in bases
    // getline(fin, line);
    // iss.clear();
    // iss.str(line);
    // while (iss >> doubtemp)
    // {
    //     paramData->snpPositions.push_back(doubtemp / paramData->BasesPerMorgan);
    // }
    // std::cerr << paramData->snpPositions.size() << " Site positions read.\n";

    // if (paramData->randSNP)
    // {
    //     paramData->snpRange.push_back(paramData->snpPositions[0]);
    //     paramData->snpRange.push_back(paramData->snpPositions[paramData->snpPositions.size() - 1]);

    //     if (paramData->snpPositions.size() != 2)
    //     {
    //         std::cerr << "Warning: Asking for random SNPs, but number of SNPs read !=2. Will only take first and last SNP positions as the range for random SNP locations\n";
    //     }
    //     if (!paramData->fixedS)
    //     {
    //         std::cerr << "Warning: Number of segregating sites (S) is not fixed. SNPRange will be ignored. Random SNPs will be simulated in non-recombining region, using theta parameter.\n";
    //     }
    //     std::cerr << "Window of SNP locations (in recUnits): " << paramData->snpRange[0] << " - " << paramData->snpRange[paramData->snpRange.size() - 1] << '\n';
    // }

    // if (!paramData->fixedS)
    // {
    //     paramData->snpPositions.resize(1);
    //     paramData->randSNP = false;
    // }

    // // should the sample be random?
    // bool randomSample = false;
    // getline(fin, line);
    // iss.clear();
    // iss.str(line);
    // while (iss >> temp)
    // {
    //     randomSample = (bool)temp;
    // }
    // std::cerr << "Random sample of carriers? " << randomSample << '\n';

    // paramData->nCarriers.resize(pops);

    // if (randomSample)
    // {
    //     vector<int> tempRead;
    //     getline(fin, line);
    //     iss.clear();
    //     iss.str(line);
    //     while (iss >> temp)
    //     {
    //         tempRead.push_back(temp);
    //     }
    //     if (tempRead.size() != pops)
    //     {
    //         std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
    //         exit(1);
    //     }

    //     for (int p = 0; p < pops; ++p)
    //     {
    //         int invcount = randbinom(tempRead[p], paramData->initialFreqs.at(p));
    //         paramData->nCarriers.at(p).push_back(tempRead[p] - invcount);
    //         paramData->nCarriers.at(p).push_back(invcount);
    //     }
    // }
    // else
    // { // number of S, I carriers in each pop
    //     int totalSample = 0;
    //     for (int p = 0; p < pops; ++p)
    //     {
    //         std::cerr << "Sample in Pop " << p << ": ";
    //         getline(fin, line);
    //         iss.clear();
    //         iss.str(line);
    //         while (iss >> temp)
    //         {
    //             paramData->nCarriers.at(p).push_back(temp);
    //             totalSample += temp;
    //             std::cerr << temp << " ";
    //         }
    //         std::cerr << '\n';
    //     }
    //     if (totalSample == 2)
    //         std::cerr << "Sample size = 2 || Warning: Informative sites length function won't work.\n";
    // }

    fin.close();
    
    YAML::Node input = YAML::LoadFile("inLABP.yaml");

    // // Number of runs
    // auto nRuns = input[NUM_RUNS].as<unsigned int>();
    // paramData->nRuns = nRuns;
    // std::cerr << "Number of runs: " << paramData->nRuns << '\n';

    // // Should we simulate kingman coalescent? (i.e. exponential waiting times). If false, does generation by generation
    // bool kingman = input[KINGMAN].as<bool>();
    // paramData->kingman = kingman;
    // std::cerr << "Assuming Kingman? " << paramData->kingman << '\n';

    // // Simulate drift trajectory?
    // bool drift = input[DRIFT].as<bool>();
    // paramData->drift = drift;
    // std::cerr << "Simulating drift? " << paramData->drift << '\n';

    // if (drift)
    // {
    //     std::cerr << "Sorry, drift not implemented in this version yet\n";
    //     exit(1);
    // }

    // // Output ms-formatted haplotypes?
    // bool msOutput = input[MS_OUTPUT].as<bool>();
    // paramData->msOutput = msOutput;
    // std::cerr << "Output each set in ms format? " << paramData->msOutput << '\n';

    // // Population sizes
    // vector<unsigned int> popSizes = input[POP_SIZES].as<vector<unsigned int>>();
    // paramData->popSizeVec = popSizes;
    // std::cerr << "Population Sizes: ";

    // for (int n = 0; n < (int)popSizes.size(); n++)
    // {
    //     std::cerr << (n) << " = " << paramData->popSizeVec.at(n) << " | ";
    // }
    // std::cerr << '\n';

    // paramData->totalPopSize = std::accumulate(popSizes.begin(), popSizes.end(), 0);

    // // Inversion initial frequencies
    // vector<double> initialFreqs = input[INVERSION_FREQS].as<vector<double>>();
    // paramData->initialFreqs = initialFreqs;

    // if (paramData->initialFreqs.size() != paramData->popSizeVec.size())
    // {
    //     std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
    //     exit(1);
    // }

    // std::cerr << "Frequencies of inverted chromosomes: ";
    // for (int n = 0; n < (int) initialFreqs.size(); n++) {
    //     std::cerr << (n) << " = " << paramData->initialFreqs.at(n) << " | ";
    // }
    // std::cerr << '\n';

    // // Simulate speciation?
    // auto speciation = input[SPECIATION];
    // bool simulateSpeciation = speciation[SIMULATE].as<bool>();
    // paramData->speciation.push_back(simulateSpeciation ? 1 : 0);
    // if (simulateSpeciation) {
    //     if (!speciation[TIME] || !speciation[INVERSION_FREQ]) {
    //         std::cerr << "Error; Simulating speciation but no time or invFreq parameters found\n";
    //         exit(1);
    //     }
    //     paramData->speciation.push_back(speciation[TIME].as<double>());
    //     paramData->speciation.push_back(speciation[INVERSION_FREQ].as<double>());

    //     std::cerr << "Speciation time: " << speciation[TIME] << "\n";
    // } else {
    //     std::cerr << "No speciation.\n";
    // }

    // // Simulate demography?
    // auto demography = input[DEMOGRAPHY_CHANGE];
    // bool simulateDemography = demography[SIMULATE].as<bool>();
    // paramData->demography.push_back(simulateDemography ? 1 : 0);
    // if (simulateDemography) {
    //     if (!demography[TIME] || !demography[COEFFICIENT]) {
    //         std::cerr << "Error; Simulating demography but no time parameter found\n";
    //         exit(1);
    //     }
    //     paramData->demography.push_back(demography[TIME].as<double>());
    //     paramData->demography.push_back(demography[COEFFICIENT].as<double>());
    //     std::cerr << "Demography change time: " << demography[TIME] << "\n";
    //     std::cerr << "Reminder: this coeff multiplies all pop sizes in the preceding epoch\n";
    // } else {
    //     std::cerr << "No demographic event.\n";
    // }

    // // Age of inversion in number of generations (EX: 1500)
    // unsigned int inversionAge = input[INVERSION_AGE].as<unsigned int>();
    // if (inversionAge > 0 && inversionAge < std::max(demography[TIME].as<double>(), speciation[TIME].as<double>())) {
    //     std::cerr << "Error: Current implementations only allows for inversions that predate speciation and demographic changes\n";
    //     exit(1);
    // }

    // paramData->inv_age = inversionAge;
    // std::cerr << "Age of inversion = " << paramData->inv_age << '\n';

    // // Migration rate (4Nm)
    // double migrationRate = input[MIG_RATE].as<double>();
    // paramData->migRate.push_back(migrationRate);
    // std::cerr << "Migration rate (pop0<-->1) = " << paramData->migRate[0] << '\n';

    // // bases per Morgan in the hh genotype.
    // double bpm = input[BP].as<double>();
    // paramData->BasesPerMorgan = bpm;
    // std::cerr << "Bases per Morgan in homokaryotypic recombination " << paramData->BasesPerMorgan << '\n';

    // // Random phi?
    // bool randPhi = input[IS_PHI_RANDOM].as<bool>();
    // paramData->randPhi = randPhi;
    // std::cerr << "Random phi values? (if 1, range below. if 0, single value read below) " << paramData->randPhi << '\n';

    // vector<double> phiRange = input[PHI].as<vector<double>>();
    // paramData->phi_range = phiRange;
    // if (randPhi) {
    //     std::cerr << "Random gene flux (phi) range = " << paramData->phi_range[0] << " - " << phiRange[1] << '\n';
    // } else {
    //     std::cerr << "Gene flux (phi) = " << paramData->phi_range[0] << '\n';
    // }

    // // the range of the inversion
    // double left = input[INVERSION_BREAKPOINTS][LEFT].as<double>();
    // double right = input[INVERSION_BREAKPOINTS][RIGHT].as<double>();
    // paramData->invRange.L = left;
    // paramData->invRange.R = right;
    // std::cerr << "Inversion from: " << paramData->invRange.L << " to " << paramData->invRange.R << " (" << paramData->invRange.L << " - " << paramData->invRange.R << " recUnits)\n";

    // // Fixed S? Value of S or sequence length (bases), Theta
    // auto snp = input[SNP];
    // bool SNPFixed = snp[IS_FIXED].as<bool>();
    // paramData->fixedS = SNPFixed;

    // unsigned int numBases = snp[NUM_BASES].as<unsigned int>();
    // paramData->n_SNPs = numBases;

    // double mutationRate = snp[MUTATION_RATE].as<double>();
    // paramData->theta = SNPFixed ? 0 : mutationRate;
    // std::cerr << (SNPFixed ? "Number of markers to simulate: " : "Number of bases (non-recombining) to simulate: ") << paramData->n_SNPs << '\n';

    // // random positions of SNPs?
    // bool randSNP = snp[RANDOM_POS].as<bool>();
    // paramData->randSNP = randSNP;
    // std::cerr << "Markers in random locations? " << paramData->randSNP << '\n';

    // range (positions of begin-end) where the SNPs are -- in bases
    if (!input[SNP][RANGE]) {
        exit(1);
    }

    const vector<double>& snpPositions = input[SNP][RANGE].as<vector<double>>();

    // Reserve space in paramData->snpPositions to avoid repeated allocations
    paramData->snpPositions.reserve(snpPositions.size());

    for (size_t i = 0; i < snpPositions.size(); ++i) {
        paramData->snpPositions.push_back(snpPositions[i] / paramData->BasesPerMorgan);
    }
    std::cerr << paramData->snpPositions.size() << " Site positions read.\n";

    if (paramData->randSNP) {
        paramData->snpRange.push_back(snpPositions[0]);
        paramData->snpRange.push_back(snpPositions.back());

        if (snpPositions.size() != 2) {
            std::cerr << "Warning: Asking for random SNPs, but number of SNPs read != 2. Will only take first and last SNP positions as the range for random SNP locations\n";
        }
        if (!paramData->fixedS) {
            std::cerr << "Warning: Number of segregating sites (S) is not fixed. SNPRange will be ignored. Random SNPs will be simulated in non-recombining region, using theta parameter.\n";
        }
        std::cerr << "Window of SNP locations (in recUnits): " << paramData->snpRange[0] << " - " << paramData->snpRange.back() << '\n';
    } else {
        paramData->snpPositions.resize(1);
    }

    // should the sample be random?
    bool randomSample = input[CARRIERS][RANDOM].as<bool>();
    std::cerr << "Random sample of carriers? " << randomSample << '\n';

    paramData->nCarriers.resize(paramData->popSizeVec.size());

    if (randomSample) {
        vector<int> carriers = input[CARRIERS][DATA][0].as<vector<int>>();
        if (carriers.size() != paramData->popSizeVec.size()) {
            std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
            exit(1);
        }

        for (int p = 0; p < paramData->popSizeVec.size(); ++p) {
            int invcount = randbinom(carriers[p], paramData->initialFreqs.at(p));
            paramData->nCarriers.at(p).push_back(carriers[p] - invcount);
            paramData->nCarriers.at(p).push_back(invcount);
        }
    } else {
        YAML::Node carriersDataRaw = input[CARRIERS][DATA];
        int totalSample = 0;
        for (int p = 0; p < paramData->popSizeVec.size(); ++p) {
            std::cerr << "Sample in Pop " << p << ": ";
            vector<int> row = carriersDataRaw[p].as<vector<int>>();
            for (int j = 0; j < row.size(); j++){
                paramData->nCarriers.at(p).push_back(row[j]);
                totalSample += row[j];
                std::cerr << row[j] << " ";
            }
            std::cerr << '\n';
        }
        if (totalSample == 2) {
            std::cerr << "Sample size = 2 || Warning: Informative sites length function won't work.\n";
        }
    }

} // end constructor

Parameters::~Parameters(){// Destructor
                          DBG("Parameters:: Destructing...")}

vector<unsigned int> Parameters::getPopulationSizes()
{ // Returns vector of population sizes (over populations)
    return paramData->popSizeVec;
}

vector<unsigned> Parameters::getSamplePerPop()
{
    vector<unsigned> popVal;
    for (unsigned i = 0; i < paramData->nCarriers.size(); ++i)
    {
        unsigned int ncarr = accumulate(paramData->nCarriers.at(i).begin(), paramData->nCarriers.at(i).end(), 0);
        for (unsigned j = 0; j < ncarr; ++j)
            popVal.push_back(i);
    }
    return popVal;
}

vector<shared_ptr<Chromosome>> Parameters::getChromVec()
{ // Returns the prototype chromosome for the initial carriers
    return paramData->initChr;
}

shared_ptr<Parameters::ParameterData> Parameters::getpData()
{
    return paramData;
}

void Parameters::setCarriers()
{

    vector<Segment> initSegVec;
    for (int j = 0; j < paramData->neut_site.size(); ++j)
    {
        Segment initSeg;
        initSeg.L = paramData->neut_site[j];
        initSeg.R = paramData->neut_site[j];
        initSegVec.push_back(initSeg);
    }
    paramData->initChr.clear();

    for (int p = 0; p < paramData->popSizeVec.size(); ++p)
    {
        for (int inv = 0; inv < 2; ++inv)
        {

            Context initCtxt(p, inv);
            shared_ptr<ARGNode> nulle;

            for (int n = 0; n < paramData->nCarriers.at(p).at(inv); ++n)
            {
                // Make a prototype chromosome:
                shared_ptr<Chromosome> chr(new Chromosome(initCtxt, initSegVec, paramData->invRange, nulle));
                paramData->initChr.push_back(chr);
            }
        }
    }
}

void Parameters::setPhi()
{
    if (paramData->randPhi)
    {
        double phiExponent = randreal(paramData->phi_range[0], paramData->phi_range[1]);
        paramData->phi = pow(10, phiExponent);
    }
    else
    {
        paramData->phi = paramData->phi_range[0];
    }
}

void Parameters::setSNPs()
{

    paramData->neut_site.clear();
    if (paramData->randSNP)
    {
        for (int p = 0; p < paramData->n_SNPs; ++p)
        {
            paramData->neut_site.push_back(randreal(paramData->snpRange[0], paramData->snpRange[1]));
        }
        sort(paramData->neut_site.begin(), paramData->neut_site.end());
    }
    else
    {
        paramData->neut_site = paramData->snpPositions;
    }
    // std::cerr.precision(12);
    // std::cerr<<"\nSetting SNPs at: ";for(int p=0;p< paramData->n_SNPs;++p){std::cerr<<paramData->neut_site[p]<<" ";}
    // std::cerr<<'\n';
}

void Parameters::setTimeStamp(double t)
{
    paramData->timestamp = t;
}
