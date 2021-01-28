// functions and methods for analysing lattice QCD simulation results via the jackknife method and 2D correlated function fits
// hadron resonance gas model related functions were also implemented

// used headers/libraries
#include <Eigen/Dense>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include <complex>
#include <tuple>
#include <random>

// include own hadron header
#include "Hadron.hh"

// lambda to calculate squares
auto sq = [](auto const &x) {
    return x * x;
};

// some small epsilon for comparisons
int const eps = 1e-6;

// complex unit
std::complex<double> const ci(0., 1.);

//
//
// READING GIVEN DATASET FOR FURTHER ANALYSIS
//
//

// read file with dataset into a raw matrix
auto ReadFile = [](std::string const &fileName) {
    // start reading
    std::ifstream fileToRead;
    fileToRead.open(fileName);

    // determine number of columns
    std::string firstLine;
    std::getline(fileToRead, firstLine);
    std::stringstream firstLineStream(firstLine);

    // number of columns in given file
    int numOfCols = 0;
    std::string tmpString;
    // count number of writes to a temporary string container
    while (firstLineStream >> tmpString)
    {
        numOfCols++;
    }
    fileToRead.close();

    // string for all the lines
    std::string line;

    // data structure (raw matrix) to store the data
    Eigen::MatrixXd rawDataMat(0, numOfCols);

    // reopen file
    fileToRead.open(fileName);
    // check if open
    if (fileToRead.is_open())
    {
        // read line by line
        int i = 0;
        while (std::getline(fileToRead, line))
        {
            // using stringstream to write matrix
            std::stringstream dataStream(line);
            rawDataMat.conservativeResize(i + 1, numOfCols);
            for (int j = 0; j < numOfCols; j++)
            {
                dataStream >> rawDataMat(i, j);
            }
            i++;
        }
        // close file
        fileToRead.close();
    }
    // error check
    else
    {
        std::cout << "ERROR\nProblem occured while reading given file." << std::endl;
        std::exit(-1);
    }

    // return raw data matrix
    return (Eigen::MatrixXd)rawDataMat;
};

//
//
// CREATING LIST OF HADRONS (from PDG file)
//
//

// container for all the hadrons in the PDG list
auto HadronList = [](std::string const &PDGList) {
    // start reading
    std::ifstream fileToRead;

    // string for all the lines
    std::string line;

    // data structure to store hadronic data
    std::vector<Hadron> hadronDataContainer;

    // reopen file
    fileToRead.open(PDGList);
    // check if open
    if (fileToRead.is_open())
    {
        // read line by line
        int i = 0;
        while (std::getline(fileToRead, line))
        {
            /*
                IN INPUT FILE: 
                particle ID
                name
                mass (in GeV)
                width
                spin degeneracy
                Baryon number
                Strangeness
                Charmness
                Bottomness
                a column of zeros for everybody
                electric charge
                number of decay channels 
            */
            // create particle data variables
            // trash (not used right now)
            std::string tmp;
            // particle name;
            std::string particleName;
            // particle mass (in GeV)
            double particleMass;
            // spin degeneracy
            int spinDeg;
            // baryon number
            int B;
            // strangness
            int S;
            // electric charge
            int Q;
            // process lines using stringstream
            std::stringstream dataStream(line);
            dataStream >> tmp;
            dataStream >> particleName;
            dataStream >> particleMass;
            dataStream >> tmp;
            dataStream >> spinDeg;
            dataStream >> B;
            dataStream >> S;
            dataStream >> tmp;
            dataStream >> tmp;
            dataStream >> tmp;
            dataStream >> Q;
            // determine particle type (boson / fermion)
            std::string particleType = "none";
            if ((spinDeg % 2) > eps)
                particleType = "fermion";
            else
                particleType = "boson";
            hadronDataContainer.push_back(Hadron(particleName, particleMass, particleType, B, Q, S, spinDeg));
            i++;
        }
        // close file
        fileToRead.close();
    }
    // error check
    else
    {
        std::cout << "ERROR\nProblem occured while reading given file." << std::endl;
        std::exit(-1);
    }

    // return raw data matrix
    return hadronDataContainer;
};

//
//
// CALCULATING SUSCEPTIBILITIES (with jackknife samples using a "vector" form)
// labeling
// susceptibility   ZContainer[index]
// imZu ----------- 0
// imZs ----------- 1
// Zuu ------------ 2
// Zud ------------ 3
// Zus ------------ 4
// Zss ------------ 5
// Zuuuu ---------- 6
// Zuuud ---------- 7
// Zuuus ---------- 8
// Zuudd ---------- 9
// Zuuds ---------- 10
// Zuuss ---------- 11
// Zudss ---------- 12
// Zusss ---------- 13
// Zssss ---------- 14
// Zuc ------------ 15
// Zsc ------------ 16
// Zcc ------------ 17
// Zuuuc ---------- 18
// Zuudc ---------- 19
// Zuusc ---------- 20
// Zudsc ---------- 21
// Zussc ---------- 22
// Zsssc ---------- 23
// Zuucc ---------- 24
// Zudcc ---------- 25
// Zuscc ---------- 26
// Zsscc ---------- 27
// Zuccc ---------- 28
// Zsccc ---------- 29
// Zcccc ---------- 30
//
//

// imZB
auto imZBCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (2 * Z[0] + Z[1]) / 3;
};

// ------------------------------------------------------------------------------------------------------------

// imZQ
auto imZQCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (Z[0] - Z[1]) / 3;
};

// ------------------------------------------------------------------------------------------------------------

// imZS
auto imZSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return -Z[1];
};

// ------------------------------------------------------------------------------------------------------------

// ZBB
auto ZBBCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (2 * Z[2] + Z[5] + 4 * Z[4] + 2 * Z[3]) / 9;
};

// ------------------------------------------------------------------------------------------------------------

// ZQQ
auto ZQQCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (5 * Z[2] + Z[5] - 2 * Z[4] - 4 * Z[3]) / 9;
};

// ------------------------------------------------------------------------------------------------------------

// ZSS
auto ZSSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return Z[5];
};

// ------------------------------------------------------------------------------------------------------------

// ZII
auto ZIICalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (Z[2] - Z[3]) / 2;
};

// ------------------------------------------------------------------------------------------------------------

// ZBQ
auto ZBQCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (Z[2] - Z[5] - Z[4] + Z[3]) / 9;
};

// ------------------------------------------------------------------------------------------------------------

// ZBS
auto ZBSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return -(Z[5] + 2 * Z[4]) / 3;
};

// ------------------------------------------------------------------------------------------------------------

// ZQS
auto ZQSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (Z[5] - Z[4]) / 3;
};

// ------------------------------------------------------------------------------------------------------------

// ZBBBB
auto ZBBBBCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (2. * Z[6] + Z[14] + 8. * Z[7] + 8. * Z[8] + 8. * Z[13] + 6. * Z[9] + 12. * Z[11] + 24. * Z[10] + 12. * Z[12]) / 81.;
};

// ------------------------------------------------------------------------------------------------------------

// ZSSSS
auto ZSSSSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return Z[14];
};

// ------------------------------------------------------------------------------------------------------------

// ZBSSS
auto ZBSSSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return -(2. * Z[13] + Z[14]) / 3.;
};

// ------------------------------------------------------------------------------------------------------------

// ZBBSS
auto ZBBSSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return (2. * Z[11] + 2. * Z[12] + 4. * Z[13] + Z[14]) / 9.;
};

// ------------------------------------------------------------------------------------------------------------

// ZBBBS
auto ZBBBSCalc = [](std::vector<Eigen::VectorXd> const &Z) {
    return -(2. * Z[8] + Z[14] + 6. * Z[10] + 6. * Z[11] + 6. * Z[13] + 6. * Z[12]) / 27.;
};

//
//
// STATISTICAL FUNCTIONS (ERROR, VARIANCE, JACKKNIFE, BOOTSTRAP, ETC...)
// including sample number reduction methods
//
//

// calculate variance (for jackknife samples)
auto JCKVariance = [](Eigen::VectorXd const &JCKSamples) {
    // size of vector
    int N = JCKSamples.size();
    // pre-factor
    double preFactor = (double)(N - 1) / N;
    // estimator / mean
    double estimator = JCKSamples.mean();
    // calculate variance
    double var = 0.;
    for (int i = 0; i < N; i++)
    {
        double val = JCKSamples(i) - estimator;
        var += sq(val);
    }
    // return variance
    return preFactor * var;
};

// ------------------------------------------------------------------------------------------------------------

// general jackknife error calculator for susceptibilities
auto ZError = [](Eigen::VectorXd const &Z) {
    return std::sqrt(JCKVariance(Z.segment(2, Z.size() - 2)));
};

// ------------------------------------------------------------------------------------------------------------

// calculate original block means (and reducing their number by averaging) from jackknife samples
auto ReducedBlocks = [](Eigen::VectorXd const &JCKSamplesOld, int const &divisor) {
    // number of samples
    int NOld = JCKSamplesOld.size();
    // test if divisor is correct for the original sample number
    if (std::abs(NOld % divisor) > eps)
    {
        std::cout << "ERROR\nIncorrect divisor during Jackknife sample reduction." << std::endl;
        std::exit(-1);
    }
    // empty vector for block values
    Eigen::VectorXd blockVals(NOld);
    // sum of (original) samples
    double sum = JCKSamplesOld.sum();
    // calculate block values and add to vector
    for (int i = 0; i < NOld; i++)
    {
        blockVals(i) = sum - (NOld - 1) * JCKSamplesOld(i);
    }

    // create new blocks
    // old blocks to add up for new blocks
    int reduced = NOld / divisor;
    // vector for new blocks (reduced)
    Eigen::VectorXd newBlocks = Eigen::VectorXd::Zero(reduced);
    // calculate new blocks
    for (int i = 0; i < reduced; i++)
    {
        newBlocks(i) = blockVals.segment(i * divisor, divisor).sum() / divisor;
    }
    // return new blocks
    return newBlocks;
};

// ------------------------------------------------------------------------------------------------------------

// calculate jackknife samples from block means
auto JCKSamplesCalculation = [](Eigen::VectorXd const &blocks) {
    // number of blocks
    int lengthBlocks = blocks.size();
    // vector for jackknife samples
    Eigen::VectorXd Samples(lengthBlocks);
    // copy data to std::vector
    std::vector<double> tmpVec(blocks.data(), blocks.data() + lengthBlocks);
    // create jackknife samples
    std::vector<double> tmpJCKVec{};
    for (int i = 0; i < lengthBlocks; i++)
    {
        // copy data
        tmpJCKVec = tmpVec;
        // delete ith element
        tmpJCKVec.erase(tmpJCKVec.begin() + i);
        // calculate mean
        Samples[i] = std::accumulate(tmpJCKVec.begin(), tmpJCKVec.end(), 0.) / (lengthBlocks - 1);
    }
    // return new jackknife samples
    return Samples;
};

// ------------------------------------------------------------------------------------------------------------

// general jackknife error calculator for susceptibilities with sample number reductions (according to divisors)
auto ZErrorJCKReduced = [](Eigen::VectorXd const &Z, int const &ZDivisor) {
    // number of jackknife samples
    int NOld = Z.size() - 2;
    // get new jackknife samples via calculating old blocks and reducing their number by averaging
    // return jackknfife error
    return std::sqrt(JCKVariance(JCKSamplesCalculation(ReducedBlocks(Z.segment(2, NOld), ZDivisor))));
};

// ------------------------------------------------------------------------------------------------------------

// estimation of function fit error from jackknife sample fits
auto JCKFitErrorEstimation = [](Eigen::VectorXd const &coeffVector, std::vector<Eigen::VectorXd> const &JCK_coeffVector) {
    // number of coefficients
    int numOfCoeffs = coeffVector.size();
    // number of jackknife samples
    int N = static_cast<int>(JCK_coeffVector.size());
    // pre-factor
    double preFactor = (double)(N - 1) / N;

    // vector for estimated errors
    Eigen::VectorXd errorVec(numOfCoeffs);

    // calculate errors via jackknife variance
    for (int coeffIndex = 0; coeffIndex < numOfCoeffs; coeffIndex++)
    {
        double tmp = 0.;
        for (int jckIndex = 0; jckIndex < N; jckIndex++)
        {
            tmp += sq(JCK_coeffVector[jckIndex](coeffIndex) - coeffVector(coeffIndex));
        }
        errorVec(coeffIndex) = std::sqrt(preFactor * tmp);
    }

    // return coefficient errors
    return errorVec;
};

//
//
// FUNCTION FITTING METHODS (2D and/or correlated)
// (not yet generalized entirely)
//
//

// calculate correlation coefficients of two datasets with given means (better this way) (jackknife)
auto CorrCoeffJCK = [](Eigen::VectorXd const &vec1, Eigen::VectorXd const &vec2, double const &mean1, double const &mean2) {
    // number of jackknife samples
    double NJck = (double)vec1.size();

    // calculate correlation (not normed)
    double corr = 0;
    for (int i = 0; i < NJck; i++)
    {
        corr += (vec1(i) - mean1) * (vec2(i) - mean2);
    }

    // return normed correlation coefficient
    return corr * (NJck - 1) / NJck;
};

// ------------------------------------------------------------------------------------------------------------

// block from the blockdiagonal covariance matrix (jackknife)
auto BlockCInverseJCK = [](Eigen::MatrixXd const &JCKs, int const &numOfQs, int const &qIndex, int const &jckNum) {
    // choose appropriate jackknife samples from given JCK matrix
    Eigen::MatrixXd JCKsQ(numOfQs, jckNum);
    for (int i = 0; i < numOfQs; i++)
    {
        JCKsQ.row(i) = JCKs.row(qIndex * numOfQs + i);
    }

    // means to calculate correlations
    Eigen::VectorXd means = Eigen::VectorXd::Zero(numOfQs);
    for (int i = 0; i < numOfQs; i++)
    {
        means(i) = JCKsQ.row(i).mean();
    }

    // covariance matrix block
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(numOfQs, numOfQs);

    for (int i = 0; i < numOfQs; i++)
    {
        for (int j = i; j < numOfQs; j++)
        {
            // check if mean and jackknife samples are zero --> not measured and set to zero artificially
            if (means(i) == 0 && means(j) == 0 && JCKsQ.row(i).isZero() && JCKsQ.row(j).isZero())
            {
                // set to identity matrix
                if (i == j)
                    C(j, i) = 1.;
                else
                    C(j, i) = 0.;
            }
            // triangular part
            else
            {
                C(j, i) = CorrCoeffJCK(JCKsQ.row(i), JCKsQ.row(j), means(i), means(j));
                // using symmetries
                if (i != j)
                    C(i, j) = C(j, i);
            }
        }
    }

    // return inverse covariance matrix block
    return (Eigen::MatrixXd)C.inverse();
};

// ------------------------------------------------------------------------------------------------------------

// general basis function with given ansatz --> ** determines the fit **
auto BasisFunc = [](int const &B,
                    int const &S,
                    int const &BOrder,
                    int const &SOrder,
                    Eigen::VectorXd const &muB,
                    Eigen::VectorXd const &muS,
                    int const &index) {
    // total number of partial derivations
    int FullOrder = BOrder + SOrder;
    // first derivative
    if (FullOrder % 2 == 1)
    {
        return std::pow(B, BOrder) * std::pow(-S, SOrder) * std::sin(B * muB(index) - S * muS(index));
    }
    // second derivative
    else if (FullOrder % 2 == 0)
    {
        return std::pow(B, BOrder) * std::pow(-S, SOrder) * std::cos(B * muB(index) - S * muS(index));
    }
    else
    {
        std::cout << "ERROR\nInvalid derivative order." << std::endl;
        std::exit(-1);
    }
};

// ------------------------------------------------------------------------------------------------------------

// LHS matrix element for given fit
auto MatElement = [](int const &i,
                     int const &j,
                     std::vector<std::pair<int, int>> const &DOrders,
                     std::vector<std::pair<int, int>> const &DOrdersMuZero,
                     std::vector<std::pair<int, int>> const &BSNumbers,
                     Eigen::VectorXd const &muB,
                     Eigen::VectorXd const &muS,
                     std::vector<Eigen::MatrixXd> const &CInvContainer,
                     Eigen::MatrixXd const &CInvMuZero) {
    // number of quantites to fit
    int numOfQs = static_cast<int>(DOrders.size());
    int numOfQsMuZero = static_cast<int>(DOrdersMuZero.size());
    // vectors to store base function data
    Eigen::VectorXd baseFunc_i(numOfQs), baseFunc_j(numOfQs);
    Eigen::VectorXd baseFuncMuZero_i(numOfQsMuZero), baseFuncMuZero_j(numOfQsMuZero);

    // helper variables
    int B_i = BSNumbers[i].first, S_i = BSNumbers[i].second;
    int B_j = BSNumbers[j].first, S_j = BSNumbers[j].second;

    // calculate matrix element
    double sum = 0.;
    // vector element for mu = 0
    for (int qIndex = 0; qIndex < numOfQsMuZero; qIndex++)
    {
        // derivation orders
        int BOrder = DOrdersMuZero[qIndex].first;
        int SOrder = DOrdersMuZero[qIndex].second;
        // fill basis function vectors
        baseFuncMuZero_i(qIndex) = BasisFunc(B_i, S_i, BOrder, SOrder, muB, muS, 0);
        baseFuncMuZero_j(qIndex) = BasisFunc(B_j, S_j, BOrder, SOrder, muB, muS, 0);
    }

    sum += baseFuncMuZero_i.transpose() * CInvMuZero * baseFuncMuZero_j;

    for (int m = 1; m < muB.size(); m++)
    {
        // create vector elements
        for (int qIndex = 0; qIndex < numOfQs; qIndex++)
        {
            // derivation orders
            int BOrder = DOrders[qIndex].first;
            int SOrder = DOrders[qIndex].second;
            // fill basis function vectors
            baseFunc_i(qIndex) = BasisFunc(B_i, S_i, BOrder, SOrder, muB, muS, m);
            baseFunc_j(qIndex) = BasisFunc(B_j, S_j, BOrder, SOrder, muB, muS, m);
        }

        // add to sum the proper covariance matrix contribution
        sum += baseFunc_i.transpose() * CInvContainer[m - 1] * baseFunc_j;
    }

    // return calculated matrix element
    return sum;
};

// ------------------------------------------------------------------------------------------------------------

// LHS matrix for the linear equation system
auto MatLHS = [](std::vector<std::pair<int, int>> const &BSNumbers,
                 std::vector<std::pair<int, int>> const &DOrders,
                 std::vector<std::pair<int, int>> const &DOrdersMuZero,
                 Eigen::VectorXd const &muB,
                 Eigen::VectorXd const &muS,
                 std::vector<Eigen::MatrixXd> const &CInvContainer,
                 Eigen::MatrixXd const &CInvMuZero) {
    // size of matrix
    int size = static_cast<int>(BSNumbers.size());

    // square matrix with the above size
    Eigen::MatrixXd LHS(size, size);

    // fill matrix
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            LHS(i, j) = MatElement(i, j, DOrders, DOrdersMuZero, BSNumbers, muB, muS, CInvContainer, CInvMuZero);
        }
    }

    // return LHS matrix
    return (Eigen::MatrixXd)LHS;
};

// ------------------------------------------------------------------------------------------------------------

// RHS vector element for given fit
auto VecElement = [](int const &i,
                     std::vector<std::pair<int, int>> const &BSNumbers,
                     std::vector<std::pair<int, int>> const &DOrders,
                     std::vector<std::pair<int, int>> const &DOrdersMuZero,
                     Eigen::MatrixXd const &yMat,
                     Eigen::MatrixXd const &yMatMuZero,
                     Eigen::VectorXd const &muB,
                     Eigen::VectorXd const &muS,
                     std::vector<Eigen::MatrixXd> const &CInvContainer,
                     Eigen::MatrixXd const &CInvMuZero) {
    // number of quantites to fit
    int numOfQs = static_cast<int>(DOrders.size());
    int numOfQsMuZero = static_cast<int>(DOrdersMuZero.size());
    // vectors to store base function data
    Eigen::VectorXd baseFunc_i(numOfQs);
    Eigen::VectorXd baseFuncMuZero_i(numOfQsMuZero);

    // vector to store given y values
    Eigen::VectorXd yVec(numOfQs);
    Eigen::VectorXd yVecMuZero(numOfQsMuZero);

    // helper variables
    int B_i = BSNumbers[i].first, S_i = BSNumbers[i].second;

    // calculate vector element
    double sum = 0.;
    for (int qIndex = 0; qIndex < numOfQsMuZero; qIndex++)
    {
        // derivation orders
        int BOrder = DOrdersMuZero[qIndex].first;
        int SOrder = DOrdersMuZero[qIndex].second;
        // fill basis function vectors
        baseFuncMuZero_i(qIndex) = BasisFunc(B_i, S_i, BOrder, SOrder, muB, muS, 0);
        // fill y vectors
        yVecMuZero(qIndex) = yMatMuZero.row(qIndex)(0);
    }

    sum += yVecMuZero.transpose() * CInvMuZero * baseFuncMuZero_i;

    for (int m = 1; m < muB.size(); m++)
    {
        // create vector elements
        for (int qIndex = 0; qIndex < numOfQs; qIndex++)
        {
            // derivation orders
            int BOrder = DOrders[qIndex].first;
            int SOrder = DOrders[qIndex].second;
            // fill basis function vectors
            baseFunc_i(qIndex) = BasisFunc(B_i, S_i, BOrder, SOrder, muB, muS, m);
            // fill y vectors
            yVec(qIndex) = yMat.row(qIndex)(m - 1);
        }

        // add to sum the covariance matrix contribution
        sum += yVec.transpose() * CInvContainer[m - 1] * baseFunc_i;
    }

    // return calculated matrix element
    return sum;
};

// ------------------------------------------------------------------------------------------------------------

// RHS vector for the linear equation system
auto VecRHS = [](std::vector<std::pair<int, int>> const &BSNumbers,
                 std::vector<std::pair<int, int>> const &DOrders,
                 std::vector<std::pair<int, int>> const &DOrdersMuZero,
                 Eigen::MatrixXd const &yMat,
                 Eigen::MatrixXd const &yMatMuZero,
                 Eigen::VectorXd const &muB,
                 Eigen::VectorXd const &muS,
                 std::vector<Eigen::MatrixXd> const &CInvContainer,
                 Eigen::MatrixXd const &CInvMuZero) {
    // size of vector
    int size = static_cast<int>(BSNumbers.size());

    // empty vector with given size
    Eigen::VectorXd RHS(size);

    // fill vector
    for (int i = 0; i < size; i++)
    {
        RHS(i) = VecElement(i, BSNumbers, DOrders, DOrdersMuZero, yMat, yMatMuZero, muB, muS, CInvContainer, CInvMuZero);
    }

    // return RHS vector
    return (Eigen::VectorXd)RHS;
};

//
//
// HADRON RESONANCE GAS (HRG) FUNCTIONS
//
//

// determine eta function for given particle type (boson / fermion)
auto EtaDetermination = [](Hadron const &H) {
    // get particle type
    std::string particleType = H.getType();

    int eta = 0;
    if (particleType == "boson")
        eta = -1;
    else if (particleType == "fermion")
        eta = 1;
    else
    {
        std::cout << "ERROR\nGiven particle type is not appropriate." << std::endl;
        H.hadronData();
        std::exit(-1);
    }

    // return appropriate eta-value
    return eta;
};

// ------------------------------------------------------------------------------------------------------------

// partial pressure calculator (for dimension = 3) at mu = 0
// (kCut + 1) index is not included in the final summation
auto iPartialPressure = [](double const &temperature, Hadron const &H, int const &kCut) {
    // determine hadron type (boson / fermion)
    int eta = EtaDetermination(H);
    // determine spin degeneracy
    int iSpinDeg = H.getSpinDegeneracy();
    // determine mass
    double iHadronMass = H.getMass();

    // pre-factor
    double preFactor = iSpinDeg * sq(temperature * iHadronMass / M_PI) / 2;

    // summation of Macdonald function
    double sumBessel = 0.;
    for (int k = 1; k <= kCut; k++)
    {
        double argumentBessel = k * iHadronMass / temperature;
        sumBessel += std::pow(-eta, k + 1) / sq(k) * gsl_sf_bessel_Kn(2, argumentBessel);
    }

    // return partial pressure
    return preFactor * sumBessel;
};

// ------------------------------------------------------------------------------------------------------------

// partial energy density calculator (for dimension = 3) at mu = 0
// (kCut + 1) index is not included in the final summation
auto iPartialEnergyDensity = [](double const &temperature, Hadron const &H, int const &kCut) {
    // determine hadron type (boson / fermion)
    int eta = EtaDetermination(H);
    // determine spin degeneracy
    int iSpinDeg = H.getSpinDegeneracy();
    // determine mass
    double iHadronMass = H.getMass();

    // pre-factor
    double preFactor = iSpinDeg * sq(temperature * iHadronMass / M_PI) / 2;

    // summation of Macdonald function
    double sumBessel = 0.;
    for (int k = 1; k <= kCut; k++)
    {
        double argumentBessel = k * iHadronMass / temperature;
        sumBessel += std::pow(-eta, k + 1) / sq(k) * (3 * gsl_sf_bessel_Kn(2, argumentBessel) + argumentBessel * gsl_sf_bessel_Kn(1, argumentBessel));
    }

    // return partial energy density
    return preFactor * sumBessel;
};

// ------------------------------------------------------------------------------------------------------------

// partial trace anomaly (interaction measure; for dimension = 3) at mu = 0
auto iPartialTraceAnomaly = [](double const &partialPressure, double const &partialEnergyDensity) {
    // return trace anomaly
    return partialEnergyDensity - 3 * partialPressure;
};

// ------------------------------------------------------------------------------------------------------------

// partial (even) suscebtibility calculator (for dimension = 3) at mu = 0 (pressure and chemical potentials are reduced)
// (kCut + 1) index is not included in the final summation
auto iPartialSusceptibility = [](int const &orderB, int const &orderS, int const &orderQ, double const &temperature, Hadron const &H, int const &kCut) {
    // check if orders are even
    if ((orderB + orderS + orderQ) % 2 > eps)
    {
        std::cout << "ERROR\nGiven susceptibility orders are not appropriate." << std::endl;
        std::exit(-1);
    }

    // determine hadron type (boson / fermion)
    int eta = EtaDetermination(H);
    // determine spin degeneracy
    int iSpinDeg = H.getSpinDegeneracy();
    // determine mass
    double iHadronMass = H.getMass();
    // determine baryon number
    int iBaryonNumber = H.getB();
    // determine electric charge
    int iElectricCharge = H.getQ();
    // determine strangeness
    int iStrangeness = H.getS();

    // pre-factor
    double preFactor = iSpinDeg * sq(iHadronMass / M_PI / temperature) / 2;

    // summation of Macdonald function
    double sumBessel = 0.;
    for (int k = 1; k <= kCut; k++)
    {
        double argumentBessel = k * iHadronMass / temperature;
        sumBessel += std::pow(-eta, k + 1) / sq(k) * std::pow(k * iBaryonNumber, orderB) * std::pow(k * iStrangeness, orderS) * std::pow(k * iElectricCharge, orderQ) * gsl_sf_bessel_Kn(2, argumentBessel);
    }

    // return partial susceptibility
    return preFactor * sumBessel;
};

// ------------------------------------------------------------------------------------------------------------

// pressure in HRG
auto PressureHRG = [](std::vector<Hadron> const &hadronList, double const &temperature, int const &kCut) {
    // number of hadrons to consider
    int numOfHadrons = static_cast<int>(hadronList.size());
    // calculate pressure
    double pressure = 0.;
    for (int i = 0; i < numOfHadrons; i++)
    {
        pressure += iPartialPressure(temperature, hadronList[i], kCut);
    }

    // return pressure value
    return pressure;
};

// ------------------------------------------------------------------------------------------------------------

// energy density in HRG
auto EnergyDensityHRG = [](std::vector<Hadron> const &hadronList, double const &temperature, int const &kCut) {
    // number of hadrons to consider
    int numOfHadrons = static_cast<int>(hadronList.size());
    // calculate energy density
    double energyDensity = 0.;
    for (int i = 0; i < numOfHadrons; i++)
    {
        energyDensity += iPartialEnergyDensity(temperature, hadronList[i], kCut);
    }

    // return energy density value
    return energyDensity;
};

// ------------------------------------------------------------------------------------------------------------

// susceptibility in HRG
auto SusceptibilityHRG = [](std::vector<Hadron> const &hadronList,
                            int const &orderB,
                            int const &orderS,
                            int const &orderQ,
                            double const &temperature,
                            int const &kCut) {
    // number of hadrons to consider
    int numOfHadrons = static_cast<int>(hadronList.size());
    // calculate given susceptibility
    double susceptibility = 0.;
    for (int i = 0; i < numOfHadrons; i++)
    {
        susceptibility += iPartialSusceptibility(orderB, orderS, orderQ, temperature, hadronList[i], kCut);
    }

    // return susceptibility
    return susceptibility;
};

//
//
// GOODNESS OF FIT TESTS
//
//

// chiSquared fit quality
auto ChiSq = [](std::vector<std::pair<int, int>> const &BSNumbers,
                std::vector<std::pair<int, int>> const &DOrders,
                Eigen::MatrixXd const &yMat,
                Eigen::VectorXd const &muB,
                Eigen::VectorXd const &muS,
                std::vector<Eigen::MatrixXd> const &CInvContainer,
                Eigen::VectorXd const &coeffVector) {
    // number of quantites fitted on
    int numOfQs = static_cast<int>(DOrders.size());
    // number of fitted parameters
    int nParams = coeffVector.size();
    // sum over blocks
    double sum = 0.;
    for (int i = 0; i < muB.size(); i++)
    {
        // delta vector for given block ~ size is equal to number of measured quantities
        Eigen::VectorXd deltaVec(numOfQs);
        // filling delta and y data vectors
        for (int j = 0; j < numOfQs; j++)
        {
            // calculate fitted function for data points
            double deltaSum = 0.;
            for (int k = 0; k < nParams; k++)
            {
                // helper variables
                int B_k = BSNumbers[k].first, S_k = BSNumbers[k].second;
                // choose y data
                deltaSum += coeffVector(k) * BasisFunc(B_k, S_k, DOrders[j].first, DOrders[j].second, muB, muS, i);
            }

            // fill delta vector
            deltaVec(j) = yMat(j, i) - deltaSum;
        }

        // add contribution to chiSquared (matrix multiplication block by block)
        sum += deltaVec.transpose() * CInvContainer[i] * deltaVec;
    }

    // return chiSquared
    return sum;
};

// ------------------------------------------------------------------------------------------------------------

// number of degrees of freedom
auto NDoF = [](std::vector<Eigen::MatrixXd> const &CInvContainer, Eigen::VectorXd const &coeffVector) {
    // determine number of datapoints
    int dataSize = 0;
    // number of quantitites in the fit
    int numOfQs = CInvContainer[0].rows();
    for (int iContainer = 0; iContainer < static_cast<int>(CInvContainer.size()); iContainer++)
    {
        for (int q = 0; q < numOfQs; q++)
        {
            // must be exactly 1 (possible mistakes in the future)
            if (CInvContainer[iContainer](q, q) != 1)
            {
                dataSize += 1;
            }
        }
    }

    // return ndof
    return dataSize - coeffVector.size();
};

// ------------------------------------------------------------------------------------------------------------

// AIC weight
auto AIC_weight = [](double const &chiSq, int const &ndof) {
    return std::exp(-0.5 * (chiSq - 2.0 * ((double)ndof)));
};

// ------------------------------------------------------------------------------------------------------------

// Q weight
auto Q_weight = [](double const &chiSq, int const &ndof) {
    return gsl_cdf_chisq_Q(chiSq, (double)ndof);
};
