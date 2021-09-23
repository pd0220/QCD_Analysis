// including used headers
#include "AnalysisTools.hh"

// PDG file name (hadron list)
std::string const PDG = "../PDG.txt";

// main function
int main()
{
    // hadrons from PDG
    std::vector<Hadron> hadronList = HadronList(PDG);

    // where to cut the summation in the partial pressure
    int kCut = 1000;

    // temperature values (MeV)
    std::vector<double> Ts{110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0,
                           190.0, 195.0, 200.0};

    // container for ZBQ values
    std::vector<double> ZBQ_HRG(Ts.size(), 0.);

    // loop for temperatures and write to file
    std::ofstream file("ZBQ_HRG.txt");
    for (int Ti = 0; Ti < static_cast<int>(Ts.size()); Ti++)
    {
        // calculate susceptibilities
        ZBQ_HRG[Ti] = SusceptibilityHRG(hadronList, 1, 0, 1, Ts[Ti] / 1000, kCut);
        file << Ts[Ti] / 1000 << " " << ZBQ_HRG[Ti] << std::endl;
    }
    file.close();
}