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

    // determining contributions in ideal HRG to the P11 sector
    // container for P11 at different temperatures
    std::vector<double> P11_T(Ts.size(), 0.);

    // loop for temperatures
    for (int iT = 0; iT < static_cast<int>(Ts.size()); iT++)
    {
        // temperature in GeV
        double T = Ts[iT] / 1000.;

        // loop for hadrons
        for (int iH = 0; iH < static_cast<int>(hadronList.size()); iH++)
        {
            // hadron
            Hadron hadron = hadronList[iH];
            // baryon number and strangeness
            int B = hadron.getB();
            int S = -hadron.getS();

            // particle data: mass, eta (boson / fermion), spin degeneracy
            double mass = hadron.getMass();
            int eta = EtaDetermination(hadron);
            int spinDeg = hadron.getSpinDegeneracy();

            // prefactor of Macdonald function
            double prefactor = spinDeg * sq(T * mass / M_PI) / 2;

            // loop for k index in partial pressure (bit too general...)
            for (int k = 1; k < kCut; k++)
            {
                // determine what sector to update
                int sectorB = k * B;
                int sectorS = k * S;

                // only consider B = 1 and S = 1
                if (sectorB != 1 || sectorS != 1)
                    break;

                // argument of Macdonald function
                double argument = k * mass / T;

                P11_T[iT] += prefactor * std::pow(-eta, k + 1) / sq(k) * gsl_sf_bessel_Kn(2, argument);
            }
        }
    }

    // write results to screen
    for (int iT = 0; iT < static_cast<int>(Ts.size()); iT++)
    {
        std::cout << Ts[iT] << " " << P11_T[iT] << std::endl;
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