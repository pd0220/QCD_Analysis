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
// argv[1] --> where to cut sectors: B = 2 or 3
// argv[2] --> number of jackknife samples
// argv[3] --> number of parameters in the continuum fit
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
    // where to cut baryon numbers (B = 2 or 3)
    int const baryonCut = std::atoi(argv[1]);
    std::string const baryonCutString = argv[1];
    if (baryonCut != 2 && baryonCut != 3)
    {
        std::cout << "ERROR\nSector cut parameter is not appropriate: instead of " << baryonCut << " it should be 2 or 3" << std::endl;
        std::exit(-2);
    }

    // number of jackknife samples
    int const jckNum = std::atoi(argv[2]);
    // number of paramteres in continuum fit
    int const paramNum = std::atoi(argv[3]);
    if (paramNum != 5 && paramNum != 6)
    {
        std::cout << "ERROR\nContinuum fit parameter number is not appropriate: must be 5 or 6." << std::endl;
        std::exit(-3);
    }

    // number of coefficients
    int coeffNum = 0;
    if (baryonCut == 2)
        coeffNum = 12;
    else
        coeffNum = 16;

    // file specifiers
    // rule: smallest lattice goes FIRST
    std::vector<std::string> dataNames{"T145.36x12.B" + baryonCutString, "T145.32x10.B" + baryonCutString, "T145.24x8.B" + baryonCutString,
                                       "T150.36x12.B" + baryonCutString, "T150.32x10.B" + baryonCutString, "T150.24x8.B" + baryonCutString,
                                       "T155.36x12.B" + baryonCutString, "T155.32x10.B" + baryonCutString, "T155.24x8.B" + baryonCutString,
                                       "T160.36x12.B" + baryonCutString, "T160.32x10.B" + baryonCutString, "T160.24x8.B" + baryonCutString};

    // number of datapoints for one coefficient ~ T datapoints * 1/Nt^2 datapoints
    int const dataSize = static_cast<int>(dataNames.size());

    // container for estimated sector coefficients
    std::vector<Eigen::VectorXd> SectorCoeffsVector(coeffNum);
    // container for sector coefficient square errors (jackknife)
    std::vector<Eigen::VectorXd> SectorErrsSqVector(coeffNum);
    // container for jackknife samples for each coefficients
    std::vector<Eigen::MatrixXd> SectorJCKVector(coeffNum);
    // set vector sizes and values (to zero)
    for (int iVec = 0; iVec < coeffNum; iVec++)
    {
        SectorCoeffsVector[iVec] = Eigen::VectorXd::Zero(dataSize);
        SectorErrsSqVector[iVec] = Eigen::VectorXd::Zero(dataSize);
        SectorJCKVector[iVec] = Eigen::MatrixXd::Zero(dataSize, jckNum);
    }

    // reading files
    for (int iFile = 0; iFile < dataSize; iFile++)
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
    Eigen::VectorXd NtInvSq(dataSize);
    Eigen::VectorXd T(dataSize);
    // setting up data point vectors to be compatible with data
    NtInvSq << 1. / 144., 1. / 100., 1. / 64.,
        1. / 144., 1. / 100., 1. / 64.,
        1. / 144., 1. / 100., 1. / 64.,
        1. / 144., 1. / 100., 1. / 64.;

    T << 145., 145., 145.,
        150., 150., 150.,
        155., 155., 155.,
        160., 160., 160.;

    // basis functions for linear fits
    // constant ~ 1
    Eigen::VectorXd basisConstant = Eigen::VectorXd::Constant(dataSize, 1);
    // linear in temperature ~ T
    Eigen::VectorXd basisLinearT = T;
    // quadratic in temperature ~ T^2
    Eigen::VectorXd basisQuadraT(dataSize);
    for (int iT = 0; iT < dataSize; iT++)
        basisQuadraT(iT) = sq(T(iT));
    // linear in inverse time lattice squared ~ 1 / Nt^2
    Eigen::VectorXd basisLinearNt = NtInvSq;
    // ~ T / Nt^2
    Eigen::VectorXd basisLinearNtLinearT(dataSize);
    for (int iData = 0; iData < dataSize; iData++)
        basisLinearNtLinearT(iData) = T(iData) * NtInvSq(iData);
    // ~ T^2 / Nt^2
    Eigen::VectorXd basisLinearNtQuadraT(dataSize);
    for (int iData = 0; iData < dataSize; iData++)
        basisLinearNtQuadraT(iData) = sq(T(iData)) * NtInvSq(iData);

    // container for basis functions
    std::vector<Eigen::VectorXd> basisFunctions;
    // set number of parameters / basis functions in continuum fit
    if (paramNum == 5)
        basisFunctions = std::vector<Eigen::VectorXd>{basisConstant, basisLinearT, basisQuadraT, basisLinearNt, basisLinearNtLinearT};
    else
        basisFunctions = std::vector<Eigen::VectorXd>{basisConstant, basisLinearT, basisQuadraT, basisLinearNt, basisLinearNtLinearT, basisLinearNtQuadraT};

    // LHS matrices
    std::vector<Eigen::MatrixXd> LHSMatContainer(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        LHSMatContainer[iCoeff] = Eigen::MatrixXd::Zero(paramNum, paramNum);
        for (int i = 0; i < paramNum; i++)
        {
            for (int j = 0; j < paramNum; j++)
                LHSMatContainer[iCoeff](i, j) = basisFunctions[i].transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * basisFunctions[j];
        }

        // checking something else...
        Eigen::MatrixXd FMat = Eigen::MatrixXd::Zero(dataSize, paramNum);
        for (int iParam = 0; iParam < paramNum; iParam++)
            FMat.col(iParam) = basisFunctions[iParam];

        Eigen::MatrixXd tempLHS = FMat.transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * FMat;

        LHSMatContainer[iCoeff] = tempLHS;
    }

    // RHS vectors
    std::vector<Eigen::VectorXd> RHSVecContainer(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        RHSVecContainer[iCoeff] = Eigen::VectorXd::Zero(paramNum);
        for (int i = 0; i < paramNum; i++)
            RHSVecContainer[iCoeff](i) = SectorCoeffsVector[iCoeff].transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * basisFunctions[i];

        // checking something else...
        Eigen::MatrixXd FMat = Eigen::MatrixXd::Zero(dataSize, paramNum);
        for (int iParam = 0; iParam < paramNum; iParam++)
            FMat.col(iParam) = basisFunctions[iParam];

        //Eigen::VectorXd tempRHS = FMat.transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * SectorCoeffsVector[iCoeff];
        Eigen::VectorXd tempRHS = SectorCoeffsVector[iCoeff].transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * FMat;

        RHSVecContainer[iCoeff] = tempRHS;
    }

    // RHS vectors for jackknife error estimation
    std::vector<std::vector<Eigen::VectorXd>> RHSVecJCKContainer(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        RHSVecJCKContainer[iCoeff] = std::vector<Eigen::VectorXd>(jckNum, Eigen::VectorXd::Zero(paramNum));
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
        {
            for (int i = 0; i < paramNum; i++)
                RHSVecJCKContainer[iCoeff][iJCK](i) = SectorJCKVector[iCoeff].col(iJCK).transpose() * SectorErrsSqVector[iCoeff].asDiagonal().inverse() * basisFunctions[i];
        }
    }

    // calculate continnum limes results and estimated errors
    std::vector<Eigen::VectorXd> continuumLimesRes(coeffNum);
    std::vector<Eigen::VectorXd> continuumLimesErr(coeffNum);
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        continuumLimesRes[iCoeff] = Eigen::VectorXd::Zero(paramNum);
        continuumLimesErr[iCoeff] = Eigen::VectorXd::Zero(paramNum);
        continuumLimesRes[iCoeff] = (LHSMatContainer[iCoeff]).fullPivLu().solve(RHSVecContainer[iCoeff]);

        Eigen::MatrixXd jckLimes = Eigen::MatrixXd::Zero(paramNum, jckNum);
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
            jckLimes.col(iJCK) = (LHSMatContainer[iCoeff]).fullPivLu().solve(RHSVecJCKContainer[iCoeff][iJCK]);
        for (int i = 0; i < paramNum; i++)
            continuumLimesErr[iCoeff](i) = std::sqrt(JCKVariance(jckLimes.row(i)));
    }

    // sectors
    std::vector<std::pair<int, int>> BSNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}, {3, 0}, {3, 1}, {3, 2}, {3, 3}};

    // write results to screen
    for (int iCoeff = 0; iCoeff < coeffNum; iCoeff++)
    {
        std::cout << "{" << BSNumbers[iCoeff].first << " , " << BSNumbers[iCoeff].second << "} " << std::endl;

        for (int iParam = 0; iParam < paramNum; iParam++)
        {
            std::cout << "         " << continuumLimesRes[iCoeff](iParam) << " +/- " << continuumLimesErr[iCoeff](iParam) << std::endl;
        }
    }

    /*
    for (int i = 0; i < coeffNum; i++)
    {
        std::cout << "{" << BSNumbers[i].first << " , " << BSNumbers[i].second << "}" << std::endl;
        std::cout << (LHSMatContainer[i]).fullPivLu().solve(RHSVecContainer[i]) << "\n" << std::endl;
    }
    */
}