// continuum limes calculation for sector coefficients

// including used headers
#include "AnalysisTools.hh"
//#include "MPRealSupport.hh"
#include <mpfr.h>

// ------------------------------------------------------------------------------------------------------------

// helper lambda for reading files
auto FileName = [](std::string const &ZRatio, std::string const &data) {
    return "../RESULTS/Extrapolation/" + ZRatio + "/" + ZRatio + "." + data + ".txt";
    //return "../ExtrapolationResult/" + ZRatio + "/" + ZRatio + "." + data + ".txt";
};

// ------------------------------------------------------------------------------------------------------------

// main function
// argv[1] --> type of ratio
// argv[2] --> where to cut sectors: B = 2 or 3
// argv[3] --> number of jackknife samples
// argv[4] --> number of parameters in the continuum fit
// argv[5] --> number of muB values where the extrapolation is evaluated
int main(int argc, char **argv)
{
    // declare matrix and vector types with multi-precision scalar type
    typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXld;
    typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> VectorXld;

    // check length of argument list
    int const argcExpected = 6;
    if (argc > argcExpected)
    {
        std::cout << "ERROR\nNumber of arguments is not appropriate: " << argc << " is given instead of " << argcExpected << "." << std::endl;
        std::exit(-1);
    }

    // reading argument list
    //
    // ratio type
    std::string const ZRatio = argv[1];
    // where to cut baryon numbers (B = 2 or 3)
    int const baryonCut = std::atoi(argv[2]);
    std::string const baryonCutString = argv[2];
    if (baryonCut != 2 && baryonCut != 3)
    {
        std::cout << "ERROR\nSector cut parameter is not appropriate: instead of " << baryonCut << " it should be 2 or 3" << std::endl;
        std::exit(-2);
    }

    // number of jackknife samples
    int const jckNum = std::atoi(argv[3]);
    // number of paramteres in continuum fit
    int const paramNum = std::atoi(argv[4]);
    if (paramNum != 5 && paramNum != 6)
    {
        std::cout << "ERROR\nContinuum fit parameter number is not appropriate: must be 5 or 6." << std::endl;
        std::exit(-3);
    }
    // number muB datapoints
    int const muBNum = std::atoi(argv[5]);

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

    // container for the given fluctuation ratio at different muB-s
    std::vector<Eigen::VectorXd> RatioVector(muBNum);
    // container for the given flucttation ratio errors (squared) at different muB-s
    std::vector<Eigen::VectorXd> RatioErrSqVector(muBNum);
    // container for jackknife samples for given fluctutaion ratio at different muB-s
    std::vector<Eigen::MatrixXd> RatioJCKVector(muBNum);
    // container for muB values
    Eigen::VectorXd muB(muBNum);
    // set vector sizes and values (to zero)
    for (int iVec = 0; iVec < muBNum; iVec++)
    {
        RatioVector[iVec] = Eigen::VectorXd::Zero(dataSize);
        RatioErrSqVector[iVec] = Eigen::VectorXd::Zero(dataSize);
        RatioJCKVector[iVec] = Eigen::MatrixXd::Zero(dataSize, jckNum);
    }

    // reading files
    for (int iFile = 0; iFile < dataSize; iFile++)
    {
        // raw data matrix
        Eigen::MatrixXd RawDataMat = ReadFile(FileName(ZRatio, dataNames[iFile]));

        // muB values (simple)
        muB = RawDataMat.col(0);

        // saving ratios
        for (int iMuB = 0; iMuB < muBNum; iMuB++)
        {
            RatioVector[iMuB](iFile) = RawDataMat(iMuB, 1);
            RatioErrSqVector[iMuB](iFile) = JCKVariance(RawDataMat.row(iMuB).segment(2, RawDataMat.cols() - 2));
            RatioJCKVector[iMuB].row(iFile) = RawDataMat.row(iMuB).segment(2, RawDataMat.cols() - 2);
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
    //
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
    std::vector<Eigen::MatrixXd> LHSMatContainer(muBNum);
    for (int iMuB = 0; iMuB < muBNum; iMuB++)
    {
        LHSMatContainer[iMuB] = Eigen::MatrixXd::Zero(paramNum, paramNum);
        Eigen::MatrixXd FMat = Eigen::MatrixXd::Zero(dataSize, paramNum);
        for (int iParam = 0; iParam < paramNum; iParam++)
            FMat.col(iParam) = basisFunctions[iParam];

        Eigen::MatrixXd tempLHS = FMat.transpose() * RatioErrSqVector[iMuB].asDiagonal().inverse() * FMat;

        LHSMatContainer[iMuB] = tempLHS;
    }

    // RHS vectors
    std::vector<Eigen::VectorXd> RHSVecContainer(muBNum);
    for (int iMuB = 0; iMuB < muBNum; iMuB++)
    {
        RHSVecContainer[iMuB] = Eigen::VectorXd::Zero(paramNum);
        Eigen::MatrixXd FMat = Eigen::MatrixXd::Zero(dataSize, paramNum);
        for (int iParam = 0; iParam < paramNum; iParam++)
            FMat.col(iParam) = basisFunctions[iParam];

        Eigen::VectorXd tempRHS = RatioVector[iMuB].transpose() * RatioErrSqVector[iMuB].asDiagonal().inverse() * FMat;

        RHSVecContainer[iMuB] = tempRHS;
    }

    // RHS vectors for jackknife error estimation
    std::vector<std::vector<Eigen::VectorXd>> RHSVecJCKContainer(muBNum);
    for (int iMuB = 0; iMuB < muBNum; iMuB++)
    {
        RHSVecJCKContainer[iMuB] = std::vector<Eigen::VectorXd>(jckNum, Eigen::VectorXd::Zero(paramNum));
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
        {
            for (int i = 0; i < paramNum; i++)
                RHSVecJCKContainer[iMuB][iJCK](i) = RatioJCKVector[iMuB].col(iJCK).transpose() * RatioErrSqVector[iMuB].asDiagonal().inverse() * basisFunctions[i];
        }
    }

    // calculate continnum limes results and estimated errors
    std::vector<VectorXld> continuumLimesRes(muBNum);
    std::vector<VectorXld> continuumLimesErr(muBNum);
    // jackknife samples
    std::vector<MatrixXld> continuumLimesJCK(muBNum);
    for (int iMuB = 0; iMuB < muBNum; iMuB++)
    {
        // initialise to zero
        continuumLimesRes[iMuB] = VectorXld::Zero(paramNum);
        continuumLimesErr[iMuB] = VectorXld::Zero(paramNum);
        continuumLimesJCK[iMuB] = MatrixXld::Zero(paramNum, jckNum);
        // temporary LHS matrix and RHS vector
        MatrixXld tempLHS = LHSMatContainer[iMuB].cast<long double>();
        VectorXld tempRHS = RHSVecContainer[iMuB].cast<long double>();
        // calculating coefficients
        continuumLimesRes[iMuB] = (tempLHS).partialPivLu().solve(tempRHS);

        // repeat calculation with jackknife samples and error estimation
        MatrixXld jckLimes = MatrixXld::Zero(paramNum, jckNum);
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
        {
            VectorXld tempRHSJCK = RHSVecJCKContainer[iMuB][iJCK].cast<long double>();
            jckLimes.col(iJCK) = (tempLHS).partialPivLu().solve(tempRHSJCK);
        }
        // save new jackknife samples
        continuumLimesJCK[iMuB] = jckLimes;
        // estimate jackknife errors
        for (int i = 0; i < paramNum; i++)
            continuumLimesErr[iMuB](i) = std::sqrt(JCKVariance(jckLimes.row(i).cast<double>()));
    }

    // sectors
    std::vector<std::pair<int, int>> BSNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}, {3, 0}, {3, 1}, {3, 2}, {3, 3}};

    // results to screen
    for (int iMuB = 0; iMuB < muBNum; iMuB++)
    {
        std::cout << muB(iMuB) << " ";
        for (int iParam = 0; iParam < 3; iParam++)
            std::cout << continuumLimesRes[iMuB](iParam) << " " << continuumLimesJCK[iMuB].row(iParam) << " ";
        std::cout << std::endl;
    }
}