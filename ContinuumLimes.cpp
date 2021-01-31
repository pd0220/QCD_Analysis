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

    // file specifiers
    // rule: smallest lattice goes FIRST
    std::vector<std::string> dataNames{
        "T" + TString + ".24x8.B" + baryonCutString,
        "T" + TString + ".32x10.B" + baryonCutString,
        "T" + TString + ".36x12.B" + baryonCutString,
        "T" + TString + ".48x16.B" + baryonCutString,
    };
    // number of different lattice sizes
    int const numLattice = static_cast<int>(dataNames.size());

    // reading files
    // container for estimated sector coefficients
    std::vector<Eigen::VectorXd> SectorCoeffsVector(numLattice);
    // container for sector coefficient square errors (jackknife)
    std::vector<Eigen::VectorXd> SectorErrsSqVector(numLattice);
    for (int iFile = 0; iFile < numLattice; iFile++)
    {
        // raw data matrix
        Eigen::MatrixXd RawDataMat = ReadFile(FileName(dataNames[iFile]));
        // sector coefficient
        SectorCoeffsVector[iFile] = RawDataMat.col(0);
        // sector coefficient error
        SectorErrsSqVector[iFile] = Eigen::VectorXd::Zero(RawDataMat.rows());
        for (int iRow = 0; iRow < RawDataMat.rows(); iRow++)
        {
            SectorErrsSqVector[iFile](iRow) = JCKVariance(RawDataMat.row(iRow).segment(1, RawDataMat.cols() - 1));
        }
    }

    
}