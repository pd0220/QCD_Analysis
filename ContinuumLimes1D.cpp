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
int main(int argc, char **argv)
{
    // check length of argument list
    int const argcExpected = 3;
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

    // data points for the linear fit
    Eigen::VectorXd Nt(3);
    Nt << 1. / 144., 1. / 100., 1. / 64.;

    // container for estimated sector coefficients
    std::vector<Eigen::VectorXd> SectorCoeffsVector(coeffNum);
    // container for sector coefficient square errors (jackknife)
    std::vector<Eigen::VectorXd> SectorErrsSqVector(coeffNum);
    // set vector sizes and values (to zero)
    for (int iVec = 0; iVec < coeffNum; iVec++)
    {
        SectorCoeffsVector[iVec] = Eigen::VectorXd::Zero(latticeNum);
        SectorErrsSqVector[iVec] = Eigen::VectorXd::Zero(latticeNum);
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
            SectorErrsSqVector[iCoeff](iFile) = JCKVariance(RawDataMat.row(iCoeff).segment(1, RawDataMat.cols() - 1));
        }
    }

    std::cout << SectorCoeffsVector[0] << std::endl;
}