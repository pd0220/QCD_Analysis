// extrapolation to real chemical potentials for lattice QCD simulations
// using previously determined sector coefficients

// used headers and libraries
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include "AnalysisTools.hh"

// ---------------------------------------------------------------------------------------

// constants
int const maxIteration = 10000;
int const muBNum = 200;
double const epsExtrapolation = 1e-10;
double const muBMin = 0.01;
double const muBMax = 2.00;

// main function
// argv[1] --> file name to use
// argv[2] --> where to cut sectors (B = 2 or 3)
int main(int argc, char **argv)
{
    // check length of argument list
    int const argcExpected = 3;
    if (argc > argcExpected)
    {
        std::cout << "ERROR\nNot enough arguments: " << argc << " is given instead of " << argcExpected << "." << std::endl;
        std::exit(-1);
    }

    // reading argument list
    // sector coefficients at fixed temperature
    std::string const fileName = argv[1];
    // where to cut sector (B = 2 or 3)
    int const SectorCut = std::atoi(argv[2]);
    if (SectorCut != 2 && SectorCut != 3)
    {
        std::cout << "ERROR\nSector cut parameter is not appropriate: instead of " << SectorCut << " it should be 2 or 3" << std::endl;
        std::exit(-2);
    }

    // raw data matrix: sector coefficients and its error
    Eigen::MatrixXd PMatrix = ReadFile(fileName);
    // number of data points
    Eigen::VectorXd muB(muBNum);
    for (int i = 0; i < muBNum; i++)
    {
        muB(i) = muBMin + i * muBMax / muBNum;
    }
    Eigen::MatrixXd muS = Eigen::MatrixXd::Zero(muBNum, PMatrix.cols());

    // used sectors
    std::vector<std::pair<int, int>> SectorNumbers;
    if (SectorCut == 2)
        SectorNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}};
    else
        SectorNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}, {3, 0}, {3, 1}, {3, 2}, {3, 3}};

    // determine muS(muB) dependence via making the expactation value of strangeness to be zero
    // using Newton iteration
    for (int iMuB = 0; iMuB < muBNum; iMuB++)
    {
        for (int iSample = 0; iSample < PMatrix.cols(); iSample++)
        {
            // Newton iteration to find roots
            // 0th step
            double guess = 1.;
            // 1st step
            double res = guess - ZFuncReal(muB(iMuB), guess, 0, 1, SectorNumbers, PMatrix.col(iSample)) / ZFuncReal(muB(iMuB), guess, 0, 2, SectorNumbers, PMatrix.col(iSample));
            // start iteration
            int nIter = 1;
            while (std::abs(guess - res) > epsExtrapolation && nIter < maxIteration)
            {
                guess = res;
                res = guess - ZFuncReal(muB(iMuB), guess, 0, 1, SectorNumbers, PMatrix.col(iSample)) / ZFuncReal(muB(iMuB), guess, 0, 2, SectorNumbers, PMatrix.col(iSample));
                nIter++;
            }
            muS(iMuB, iSample) = res;
        }
    }

    // susceptibility ratios and their estimated errors
    Eigen::VectorXd Z12 = Eigen::VectorXd::Zero(muBNum);
    Eigen::VectorXd Z31 = Eigen::VectorXd::Zero(muBNum);
    Eigen::VectorXd Z42 = Eigen::VectorXd::Zero(muBNum);
    Eigen::MatrixXd Z12Err = Eigen::MatrixXd::Zero(muBNum, PMatrix.cols() - 1);
    Eigen::MatrixXd Z31Err = Eigen::MatrixXd::Zero(muBNum, PMatrix.cols() - 1);
    Eigen::MatrixXd Z42Err = Eigen::MatrixXd::Zero(muBNum, PMatrix.cols() - 1);

    // extra data 1
    Eigen::VectorXd Diff31 = Eigen::VectorXd::Zero(muBNum);
    Eigen::VectorXd Diff42 = Eigen::VectorXd::Zero(muBNum);
    Eigen::MatrixXd Diff31Err = Eigen::MatrixXd::Zero(muBNum, PMatrix.cols() - 1);
    Eigen::MatrixXd Diff42Err = Eigen::MatrixXd::Zero(muBNum, PMatrix.cols() - 1);
    for (int i = 0; i < muBNum; i++)
    {
        // results
        Z12(i) = ZFuncReal(muB(i), muS(i, 0), 1, 0, SectorNumbers, PMatrix.col(0)) / ZFuncReal(muB(i), muS(i, 0), 2, 0, SectorNumbers, PMatrix.col(0));
        Z31(i) = ZFuncReal(muB(i), muS(i, 0), 3, 0, SectorNumbers, PMatrix.col(0)) / ZFuncReal(muB(i), muS(i, 0), 1, 0, SectorNumbers, PMatrix.col(0));
        Z42(i) = ZFuncReal(muB(i), muS(i, 0), 4, 0, SectorNumbers, PMatrix.col(0)) / ZFuncReal(muB(i), muS(i, 0), 2, 0, SectorNumbers, PMatrix.col(0));

        // extra 1 results
        Diff31(i) = Z31(i) - Z31(0);
        Diff42(i) = (Z42(i) - Z42(0)) / 3.;

        // errors from jackknife samples
        for (int iSample = 0; iSample < PMatrix.cols() - 1; iSample++)
        {
            Z12Err(i, iSample) = ZFuncReal(muB(i), muS(i, iSample + 1), 1, 0, SectorNumbers, PMatrix.col(iSample + 1)) / ZFuncReal(muB(i), muS(i, iSample + 1), 2, 0, SectorNumbers, PMatrix.col(iSample + 1));
            Z31Err(i, iSample) = ZFuncReal(muB(i), muS(i, iSample + 1), 3, 0, SectorNumbers, PMatrix.col(iSample + 1)) / ZFuncReal(muB(i), muS(i, iSample + 1), 1, 0, SectorNumbers, PMatrix.col(iSample + 1));
            Z42Err(i, iSample) = ZFuncReal(muB(i), muS(i, iSample + 1), 4, 0, SectorNumbers, PMatrix.col(iSample + 1)) / ZFuncReal(muB(i), muS(i, iSample + 1), 2, 0, SectorNumbers, PMatrix.col(iSample + 1));

            // extra 1 results
            Diff31Err(i, iSample) = Z31Err(i, iSample) - Z31Err(0, iSample);
            Diff42Err(i, iSample) = (Z42Err(i, iSample) - Z42Err(0, iSample)) / 3.;
        }
    }

    for (int i = 0; i < muBNum; i++)
    {
        // write to screen
        std::cout << muB(i) << " "
                  << muS.col(0)(i) << " " << std::sqrt(JCKVariance(muS.row(i).segment(1, PMatrix.cols() - 1))) << " "
                  << Z12(i) << " " << std::sqrt(JCKVariance(Z12Err.row(i))) << " "
                  << Z31(i) << " " << std::sqrt(JCKVariance(Z31Err.row(i))) << " "
                  << Z42(i) << " " << std::sqrt(JCKVariance(Z42Err.row(i))) << std::endl;
        //<< Diff31(i) << " " << std::sqrt(JCKVariance(Diff31Err.row(i))) << " "
        //<< Diff42(i) << " " << std::sqrt(JCKVariance(Diff42Err.row(i))) << std::endl;
    }

    /*
    for (int i = 0; i < muBNum; i++)
    {
        std::cout << muB(i) << " ";
        std::cout << Z12(i) << " ";
        for (int j = 0; j < Z12Err.cols(); j++)
        {
            std::cout << Z12Err(i, j) << " ";
        }
        std::cout << Z31(i) << " ";
        for (int j = 0; j < Z31Err.cols(); j++)
        {
            std::cout << Z31Err(i, j) << " ";
        }
        std::cout << Z42(i) << " ";
        for (int j = 0; j < Z42Err.cols(); j++)
        {
            std::cout << Z42Err(i, j) << " ";
        }
        std::cout << std::endl;
    }
    */
}