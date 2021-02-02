// continuum limes calculation for sector coefficients

// including used header
#include "AnalysisTools.hh"

// ------------------------------------------------------------------------------------------------------------

// helper lambda for reading files
auto FileName = [](std::string const &data) {
    return "../FitResult/Sectors." + data + ".txt";
};

// ------------------------------------------------------------------------------------------------------------

// main function
// argv[1] --> temperature (145, 150, 155 or 160 [MeV])
// argv[2] --> where to cut sectors: B = 2 or 3
// argv[3] --> number of jackknife samples
int main(int argc, char **argv)
{
    // check length of argument list
    int const argcExpected = 4;
    if (argc > argcExpected)
    {
        std::cout << "ERROR\nNumber of arguments is not appropriate: " << argc << " is given instead of " << argcExpected << "." << std::endl;
        std::exit(-1);
    }

    // reading argument list
    //
    // fixing temperature value
    std::string const TString = argv[1];

    // where to cut baryon numbers (B = 2 or 3)
    int const baryonCut = std::atoi(argv[2]);
    std::string const baryonCutString = argv[2];
    if (baryonCut != 2 && baryonCut != 3)
    {
        std::cout << "ERROR\nSector cut parameter is not appropriate: instead of " << baryonCut << " it should be 2 or 3" << std::endl;
        std::exit(-2);
    }
    int const jckNum = std::atoi(argv[3]);

    // number of coefficients
    int coeffNum = 0;
    if (baryonCut == 2)
        coeffNum = 12;
    else
        coeffNum = 16;

    // file specifiers
    // rule: smallest lattice goes FIRST
    std::vector<std::string> dataNames{
        "T" + TString + ".36x12.B" + baryonCutString,
        "T" + TString + ".32x10.B" + baryonCutString,
        "T" + TString + ".24x8.B" + baryonCutString};
    // number of different lattice sizes
    int const latticeNum = static_cast<int>(dataNames.size());

    // container for estimated sector coefficients
    std::vector<Eigen::VectorXd> SectorCoeffsVector(coeffNum);
    // container for sector coefficient square errors (jackknife)
    std::vector<Eigen::VectorXd> SectorErrsSqVector(coeffNum);
    // container for jackknife samples for each coefficients
    std::vector<Eigen::MatrixXd> SectorJCKVector(coeffNum);
    // set vector sizes and values (to zero)
    for (int iVec = 0; iVec < coeffNum; iVec++)
    {
        SectorCoeffsVector[iVec] = Eigen::VectorXd::Zero(latticeNum);
        SectorErrsSqVector[iVec] = Eigen::VectorXd::Zero(latticeNum);
        SectorJCKVector[iVec] = Eigen::MatrixXd::Zero(latticeNum, jckNum);
    }

    // reading files
    for (int iFile = 0; iFile < latticeNum; iFile++)
    {
        // raw data matrix
        Eigen::MatrixXd RawDataMat = ReadFile(FileName(dataNames[iFile]));

        // saving coefficients at different temperatures
        for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
        {
            SectorCoeffsVector[iCoeff](iFile) = RawDataMat(iCoeff, 0);
            SectorJCKVector[iCoeff].row(iFile) = RawDataMat.row(iCoeff).segment(1, RawDataMat.cols() - 1);
            SectorErrsSqVector[iCoeff](iFile) = JCKVariance(RawDataMat.row(iCoeff).segment(1, RawDataMat.cols() - 1));
        }
    }

    // data points for the linear fit
    Eigen::VectorXd NtInvSq(latticeNum);
    NtInvSq << 1. / 144., 1. / 100., 1. / 64.;

    // basis functions for linear fits
    Eigen::VectorXd basisConstant = Eigen::VectorXd::Constant(latticeNum, 1);
    Eigen::VectorXd basisLinear = NtInvSq;
    std::vector<Eigen::VectorXd> basisFunctions{basisConstant, basisLinear};

    // LHS matrices
    std::vector<Eigen::MatrixXd> LHSMatConatiner(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        LHSMatConatiner[iCoeff] = Eigen::MatrixXd::Zero(2, 2);
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
                LHSMatConatiner[iCoeff](i, j) = basisFunctions[i].transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * basisFunctions[j];
        }
    }

    // RHS vectors
    std::vector<Eigen::VectorXd> RHSVecContainer(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        RHSVecContainer[iCoeff] = Eigen::VectorXd::Zero(2);
        for (int i = 0; i < 2; i++)
            RHSVecContainer[iCoeff](i) = SectorCoeffsVector[iCoeff].transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * basisFunctions[i];
    }

    // RHS vectors for jackknife error estimation
    std::vector<std::vector<Eigen::VectorXd>> RHSVecJCKContainer(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        RHSVecJCKContainer[iCoeff] = std::vector<Eigen::VectorXd>(jckNum, Eigen::VectorXd::Zero(2));
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
        {
            for (int i = 0; i < 2; i++)
                RHSVecJCKContainer[iCoeff][iJCK](i) = SectorJCKVector[iCoeff].col(iJCK).transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * basisFunctions[i];
        }
    }

    // calculate continnum limes results and estimated errors
    std::vector<Eigen::VectorXd> continuumLimesRes(coeffNum);
    std::vector<Eigen::VectorXd> continuumLimesErr(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        continuumLimesRes[iCoeff] = Eigen::VectorXd::Zero(2);
        continuumLimesErr[iCoeff] = Eigen::VectorXd::Zero(2);
        continuumLimesRes[iCoeff] = (LHSMatConatiner[iCoeff]).fullPivLu().solve(RHSVecContainer[iCoeff]);

        Eigen::MatrixXd jckLimes = Eigen::MatrixXd::Zero(2, jckNum);
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
            jckLimes.col(iJCK) = (LHSMatConatiner[iCoeff]).fullPivLu().solve(RHSVecJCKContainer[iCoeff][iJCK]);
        for (int i = 0; i < 2; i++)
            continuumLimesErr[iCoeff](i) = std::sqrt(JCKVariance(jckLimes.row(i)));
    }

    // sectors
    std::vector<std::pair<int, int>> BSNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}, {3, 0}, {3, 1}, {3, 2}, {3, 3}};

    // write results to screen
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        std::cout << "{" << BSNumbers[iCoeff].first << " , " << BSNumbers[iCoeff].second << "} " << continuumLimesRes[iCoeff](0) << " +/- " << continuumLimesErr[iCoeff](0) << std::endl;
    }
}