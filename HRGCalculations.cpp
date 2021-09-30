// including used headers
#include "AnalysisTools.hh"

// PDG file name (hadron list)
const std::string PDG = "../PDG.txt";

// phenomenological constant in GeV^-2
//const double KPhenom = 58.566954;
const double KPhenom = 58.86;

// main function
int main()
{
    // hadrons from PDG
    std::vector<Hadron> hadronList = HadronList(PDG);

    // where to cut the summation in the partial pressure
    int kCut = 10;

    // temperature values (MeV)
    std::vector<double> Ts{145., 146., 147., 148., 149., 150., 151., 152., 153., 154., 155., 156., 157., 158., 159., 160.};

    // determining contributions in ideal HRG and in repulsive mean field HRG
    std::vector<double> P21_T(Ts.size(), 0.);
    std::vector<double> P20_T(Ts.size(), 0.);
    std::vector<double> P22_T(Ts.size(), 0.);
    std::vector<double> repP21_T(Ts.size(), 0.);
    std::vector<double> repP20_T(Ts.size(), 0.);
    std::vector<double> repP22_T(Ts.size(), 0.);

    // loop for temperatures
    for (int iT = 0; iT < static_cast<int>(Ts.size()); iT++)
    {
        // temperature in GeV
        double T = Ts[iT] / 1000.;

        // contribution containers
        double Ps0 = 0.;
        double Ps1 = 0.;
        double Ps2 = 0.;

        // loop for hadrons
        for (int iH = 0; iH < static_cast<int>(hadronList.size()); iH++)
        {
            // hadron
            Hadron hadron = hadronList[iH];

            // get contributions from different hadrons
            P20_T[iT] += SectorContribution(2, 0, T, hadron, kCut);
            P21_T[iT] += SectorContribution(2, 1, T, hadron, kCut);
            P22_T[iT] += SectorContribution(2, 2, T, hadron, kCut);

            // contribution for the repulsive mean field HRG
            Ps0 += SectorContribution(1, 0, T, hadron, kCut);
            Ps1 += SectorContribution(1, 1, T, hadron, kCut);
            Ps2 += SectorContribution(1, 2, T, hadron, kCut);
        }
        repP20_T[iT] += -KPhenom / sq(T) * sq(Ps0) / sq(sq(T)) / 2.;
        repP21_T[iT] += -KPhenom / sq(T) * Ps0 * Ps1 / sq(sq(T));
        repP22_T[iT] += -KPhenom / sq(T) * (sq(Ps1) + 2 * Ps0 * Ps2) / sq(sq(T)) / 2.;
    }

    // write results to screen
    for (int iT = 0; iT < static_cast<int>(Ts.size()); iT++)
    {
        std::cout << Ts[iT] << " " << P20_T[iT] << " " << repP20_T[iT] << std::endl;
        //std::cout << Ts[iT] << " " << P22_T[iT] << " " << repP22_T[iT] << std::endl;
        //std::cout << Ts[iT] << " " << P21_T[iT] << " " << repP21_T[iT] << std::endl;
    }

    /*
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
*/
}