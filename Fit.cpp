// fitting sector coefficients using lattice QCD simulation results

// including used header
#include "AnalysisTools.hh"

// ------------------------------------------------------------------------------------------------------------

// PDG file name (hadron list)
std::string const PDG = "../PDG.txt";

// ------------------------------------------------------------------------------------------------------------

// main function
// argv[1] --> name of given file with dataset
// argv[2] --> name of given file with dataset (zero chemical potential)
// argv[3] --> number of jackknife samples
// argv[4] --> number of used susceptibilities (Zu, Zs, etc.)
// argv[5] --> divisor for jackknife sample number reduction
// argv[6] --> where to cut sectors --> B = 2 or 3
int main(int, char **argv)
{
    // prepare for reading given file
    // string for file name (im mu)
    std::string const fileName = argv[1];
    // string for file name (mu = 0)
    std::string const fileNameMuZero = argv[2];
    // matrix container for raw data for analysis
    Eigen::MatrixXd const rawDataMat = ReadFile(fileName);
    // data at mu = 0
    Eigen::MatrixXd const rawDataMatMuZero = ReadFile(fileNameMuZero);

    // number of jackknife samples
    int jckNum = std::atoi(argv[3]);
    // number of used susceptibilities
    int const ZNum = std::atoi(argv[4]);
    // divisor for jackknife sample number reduction
    int const divisor = std::atoi(argv[5]);
    // check if the number of jackknife samples can be divided by the divisor
    if (jckNum % divisor > eps)
    {
        std::cout << "ERROR\nThe ,,jckNum'' and ,,divisor'' pair is not appropriate." << std::endl;
        std::exit(-1);
    }
    // where to cut sectors
    int const BCut = std::atoi(argv[6]);
    // check input: B = 2 or 3
    if (BCut < 2 || BCut > 3)
    {
        std::cout << "ERROR\nSector cut must be 2 or 3." << std::endl;
        std::exit(-1);
    }

    // new matrix for the fit data
    Eigen::MatrixXd FitMat = Eigen::MatrixXd::Zero(rawDataMat.rows() + 1, rawDataMat.cols());
    FitMat.bottomRows(rawDataMat.rows()) = rawDataMat;
    // beta and Nt (not relevant now)
    FitMat(0, 0) = rawDataMatMuZero(0, 0);
    FitMat(0, 1) = rawDataMatMuZero(0, 1);
    // muB = 0
    FitMat(0, 2) = 0;
    // muS = 0
    FitMat(0, 3) = 0;
    // non-zero data (imZu = 0 and imZs = 0 at mu = 0)
    FitMat.row(0).segment(4 + 2 * (jckNum + 2), (ZNum - 2) * (jckNum + 2)) = rawDataMatMuZero.row(0).segment(2, (ZNum - 2) * (jckNum + 2));

    // number of rows of fit matrix
    int const rows = FitMat.rows();

    // chemical potentials for baryon numbers and strangeness
    Eigen::VectorXd const muB = FitMat.col(2);
    Eigen::VectorXd const muS = FitMat.col(3);

    // susceptibilities (regarding u, d, s flavours) with error and jackknife samples
    // size of vectors (val + err + jck samples)
    int const ZSize = 2 + jckNum;
    // overwrite jackknife sample number for sample number reduction
    jckNum = jckNum / divisor;
    // vectors to store imZB values and their estimated errors (results) + JCK
    Eigen::VectorXd imZBVals(rows);
    Eigen::VectorXd imZBErrs(rows);
    Eigen::MatrixXd imZBJCKs(rows, jckNum);
    // vectors to store imZQ values and their estimated errors (results) + JCK
    Eigen::VectorXd imZQVals(rows);
    Eigen::VectorXd imZQErrs(rows);
    Eigen::MatrixXd imZQJCKs(rows, jckNum);
    // vectors to store imZS values and their estimated errors (results) + JCK
    Eigen::VectorXd imZSVals(rows);
    Eigen::VectorXd imZSErrs(rows);
    Eigen::MatrixXd imZSJCKs(rows, jckNum);
    // vectors to store ZBB values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBBVals(rows);
    Eigen::VectorXd ZBBErrs(rows);
    Eigen::MatrixXd ZBBJCKs(rows, jckNum);
    // vectors to store ZQQ values and their estimated errors (results) + JCK
    Eigen::VectorXd ZQQVals(rows);
    Eigen::VectorXd ZQQErrs(rows);
    Eigen::MatrixXd ZQQJCKs(rows, jckNum);
    // vectors to store ZSS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZSSVals(rows);
    Eigen::VectorXd ZSSErrs(rows);
    Eigen::MatrixXd ZSSJCKs(rows, jckNum);
    // vectors to store ZBQ values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBQVals(rows);
    Eigen::VectorXd ZBQErrs(rows);
    Eigen::MatrixXd ZBQJCKs(rows, jckNum);
    // vectors to store ZBS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBSVals(rows);
    Eigen::VectorXd ZBSErrs(rows);
    Eigen::MatrixXd ZBSJCKs(rows, jckNum);
    // vectors to store ZQS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZQSVals(rows);
    Eigen::VectorXd ZQSErrs(rows);
    Eigen::MatrixXd ZQSJCKs(rows, jckNum);
    // vectors to store ZII values and their estimated errors (results) + JCK
    Eigen::VectorXd ZIIVals(rows);
    Eigen::VectorXd ZIIErrs(rows);
    Eigen::MatrixXd ZIIJCKs(rows, jckNum);
    // vectors to store ZBBBB values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBBBBVals(rows);
    Eigen::VectorXd ZBBBBErrs(rows);
    Eigen::MatrixXd ZBBBBJCKs(rows, jckNum);
    // vectors to store ZSSSS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZSSSSVals(rows);
    Eigen::VectorXd ZSSSSErrs(rows);
    Eigen::MatrixXd ZSSSSJCKs(rows, jckNum);
    // vectors to store ZBSSS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBSSSVals(rows);
    Eigen::VectorXd ZBSSSErrs(rows);
    Eigen::MatrixXd ZBSSSJCKs(rows, jckNum);
    // vectors to store ZBBSS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBBSSVals(rows);
    Eigen::VectorXd ZBBSSErrs(rows);
    Eigen::MatrixXd ZBBSSJCKs(rows, jckNum);
    // vectors to store ZBBBS values and their estimated errors (results) + JCK
    Eigen::VectorXd ZBBBSVals(rows);
    Eigen::VectorXd ZBBBSErrs(rows);
    Eigen::MatrixXd ZBBBSJCKs(rows, jckNum);
    // container for "flavour vectors"
    std::vector<Eigen::VectorXd> ZContainer(ZNum);
    // loop for every row
    for (int i = 0; i < rows; i++)
    {
        // filling up vectors
        for (int j = 0; j < ZNum; j++)
        {
            // if sample number reduction is ON
            if (std::abs(1 - divisor) > eps)
            {
                // create temporary vector for sample number reduction
                Eigen::VectorXd tmpVec = FitMat.row(i).segment(4 + j * ZSize, ZSize);

                // temporary JCK vectors
                Eigen::VectorXd tmpJCKVec_OLD = tmpVec.segment(2, ZSize - 2);
                // calculate original blocks and perform the sample number reduction
                Eigen::VectorXd tmpJCKVec_NEW = JCKSamplesCalculation(ReducedBlocks(tmpJCKVec_OLD, divisor));

                // JCK result (+ val + err)
                // new size for Z
                int ZSize_NEW = 2 + jckNum;
                // results
                Eigen::VectorXd tmpResult(ZSize_NEW);
                // value
                tmpResult(0) = tmpVec(0);
                // error
                tmpResult(1) = tmpVec(1);
                // jackknife samples
                tmpResult.segment(2, jckNum) = tmpJCKVec_NEW;

                ZContainer[j] = tmpResult;
            }
            else
                ZContainer[j] = FitMat.row(i).segment(4 + j * ZSize, ZSize);
        }

        // calculate (in vector form with jackknife samples)
        Eigen::VectorXd imZB = imZBCalc(ZContainer);
        Eigen::VectorXd imZQ = imZQCalc(ZContainer);
        Eigen::VectorXd imZS = imZSCalc(ZContainer);
        Eigen::VectorXd ZBB = ZBBCalc(ZContainer);
        Eigen::VectorXd ZQQ = ZQQCalc(ZContainer);
        Eigen::VectorXd ZSS = ZSSCalc(ZContainer);
        Eigen::VectorXd ZBQ = ZBQCalc(ZContainer);
        Eigen::VectorXd ZBS = ZBSCalc(ZContainer);
        Eigen::VectorXd ZQS = ZQSCalc(ZContainer);
        Eigen::VectorXd ZII = ZIICalc(ZContainer);
        Eigen::VectorXd ZBBBB = ZBBBBCalc(ZContainer);
        Eigen::VectorXd ZSSSS = ZSSSSCalc(ZContainer);
        Eigen::VectorXd ZBSSS = ZBSSSCalc(ZContainer);
        Eigen::VectorXd ZBBSS = ZBBSSCalc(ZContainer);
        Eigen::VectorXd ZBBBS = ZBBBSCalc(ZContainer);

        // save results
        // imZB
        imZBVals(i) = imZB(0);
        // imZQ
        imZQVals(i) = imZQ(0);
        // imZS
        imZSVals(i) = imZS(0);
        // ZBB
        ZBBVals(i) = ZBB(0);
        // ZQQ
        ZQQVals(i) = ZQQ(0);
        // ZSS
        ZSSVals(i) = ZSS(0);
        // ZBQ
        ZBQVals(i) = ZBQ(0);
        // ZBS
        ZBSVals(i) = ZBS(0);
        // ZQS
        ZQSVals(i) = ZQS(0);
        // ZII
        ZIIVals(i) = ZII(0);
        // ZBBBB
        ZBBBBVals(i) = ZBBBB(0);
        // ZSSSS
        ZSSSSVals(i) = ZSSSS(0);
        // ZBSSS
        ZBSSSVals(i) = ZBSSS(0);
        // ZBBSS
        ZBBSSVals(i) = ZBBSS(0);
        // ZBBBS
        ZBBBSVals(i) = ZBBBS(0);

        // save jackknife samples
        // imZB
        imZBJCKs.row(i) = imZB.segment(2, jckNum);
        // imZQ
        imZQJCKs.row(i) = imZQ.segment(2, jckNum);
        // imZS
        imZSJCKs.row(i) = imZS.segment(2, jckNum);
        // ZBB
        ZBBJCKs.row(i) = ZBB.segment(2, jckNum);
        // ZQQ
        ZQQJCKs.row(i) = ZQQ.segment(2, jckNum);
        // ZSS
        ZSSJCKs.row(i) = ZSS.segment(2, jckNum);
        // ZBQ
        ZBQJCKs.row(i) = ZBQ.segment(2, jckNum);
        // ZBS
        ZBSJCKs.row(i) = ZBS.segment(2, jckNum);
        // ZQS
        ZQSJCKs.row(i) = ZQS.segment(2, jckNum);
        // ZII
        ZIIJCKs.row(i) = ZII.segment(2, jckNum);
        // ZBBBB
        ZBBBBJCKs.row(i) = ZBBBB.segment(2, jckNum);
        // ZSSSS
        ZSSSSJCKs.row(i) = ZSSSS.segment(2, jckNum);
        // ZBSSS
        ZBSSSJCKs.row(i) = ZBSSS.segment(2, jckNum);
        // ZBBSS
        ZBBSSJCKs.row(i) = ZBBSS.segment(2, jckNum);
        // ZBBBS
        ZBBBSJCKs.row(i) = ZBBBS.segment(2, jckNum);

        // save errors
        // imZB
        imZBErrs(i) = ZError(imZB);
        // imZQ
        imZQErrs(i) = ZError(imZQ);
        // imZS
        imZSErrs(i) = ZError(imZS);
        // ZBB
        ZBBErrs(i) = ZError(ZBB);
        // ZQQ
        ZQQErrs(i) = ZError(ZQQ);
        // ZSS
        ZSSErrs(i) = ZError(ZSS);
        // ZBQ
        ZBQErrs(i) = ZError(ZBQ);
        // ZBS
        ZBSErrs(i) = ZError(ZBS);
        // ZQS
        ZQSErrs(i) = ZError(ZQS);
        // ZII
        ZIIErrs(i) = ZError(ZII);
        // ZBBBB
        ZBBBBErrs(i) = ZError(ZBBBB);
        // ZSSSS
        ZSSSSErrs(i) = ZError(ZSSSS);
        // ZBSSS
        ZBSSSErrs(i) = ZError(ZBSSS);
        // ZBBSS
        ZBBSSErrs(i) = ZError(ZBBSS);
        // ZBBBS
        ZBBBSErrs(i) = ZError(ZBBBS);
    }

    //
    // START FIT
    // --> imZB and imZS (correlated)
    // --> ZBB, ZBS, ZSS, ZBBBB, ZSSSS, ZBSSS, ZBBSS and ZBBBS at mu = 0
    //

    // number of x-values (muB and muS)
    int const N = muB.size();

    // what quantities we are fitting on (imZB and imZS)
    std::vector<std::pair<int, int>> DOrders{{1, 0}, {0, 1}};
    // what quantities we are fittin at mu = 0 (ZBB, ZBS, ZSS, ZBBBB, ZSSSS, ZBSSS, ZBBSS and ZBBBS)
    std::vector<std::pair<int, int>> DOrdersMuZero{{2, 0}, {1, 1}, {0, 2}, {4, 0}, {0, 4}, {1, 3}, {2, 2}, {3, 1}};
    // number of quantitites
    int const numOfQs = static_cast<int>(DOrders.size());
    int const numOfQsMuZero = static_cast<int>(DOrdersMuZero.size());

    // y data matrix to calculate RHS vector
    Eigen::MatrixXd yMat(numOfQs, N - 1);
    yMat.row(0) = imZBVals.segment(1, N - 1);
    yMat.row(1) = imZSVals.segment(1, N - 1);

    // y data matrix to calculate RHS vector at mu = 0
    Eigen::MatrixXd yMatMuZero(numOfQsMuZero, 1);
    yMatMuZero(0, 0) = ZBBVals(0);
    yMatMuZero(1, 0) = ZBSVals(0);
    yMatMuZero(2, 0) = ZSSVals(0);
    yMatMuZero(3, 0) = ZBBBBVals(0);
    yMatMuZero(4, 0) = ZSSSSVals(0);
    yMatMuZero(5, 0) = ZBSSSVals(0);
    yMatMuZero(6, 0) = ZBBSSVals(0);
    yMatMuZero(7, 0) = ZBBBSVals(0);

    // JCK samples with ordered structure (required for covariance matrix estimation)
    Eigen::MatrixXd JCKSamplesForFit(numOfQs * (N - 1), jckNum);
    for (int i = 0; i < N - 1; i++)
    {
        for (int q = 0; q < numOfQs; q++)
        {
            if (q == 0)
                JCKSamplesForFit.row(numOfQs * i + q) = imZBJCKs.row(i + 1);
            else if (q == 1)
                JCKSamplesForFit.row(numOfQs * i + q) = imZSJCKs.row(i + 1);
        }
    }
    // JCK samples with ordered structure (required for covariance matrix estimation) at mu = 0
    Eigen::MatrixXd JCKSamplesForFitMuZero(numOfQsMuZero, jckNum);
    JCKSamplesForFitMuZero.row(0) = ZBBJCKs.row(0);
    JCKSamplesForFitMuZero.row(1) = ZBSJCKs.row(0);
    JCKSamplesForFitMuZero.row(2) = ZSSJCKs.row(0);
    JCKSamplesForFitMuZero.row(3) = ZBBBBJCKs.row(0);
    JCKSamplesForFitMuZero.row(4) = ZSSSSJCKs.row(0);
    JCKSamplesForFitMuZero.row(5) = ZBSSSJCKs.row(0);
    JCKSamplesForFitMuZero.row(6) = ZBBSSJCKs.row(0);
    JCKSamplesForFitMuZero.row(7) = ZBBBSJCKs.row(0);

    // inverse covariance matrix blocks
    std::vector<Eigen::MatrixXd> CInvContainer(N - 1, Eigen::MatrixXd(numOfQs, numOfQs));
    for (int i = 0; i < N - 1; i++)
    {
        CInvContainer[i] = BlockCInverseJCK(JCKSamplesForFit, numOfQs, i, jckNum);
    }
    // inverse covariance matrix blocks at mu = 0
    Eigen::MatrixXd CInvMuZero = BlockCInverseJCK(JCKSamplesForFitMuZero, numOfQsMuZero, 0, jckNum);

    // what basis functions shall be included in the fit {B, S} ~ sectors
    std::vector<std::pair<int, int>> BSNumbers;
    if (BCut == 2)
        BSNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}};
    else
        BSNumbers = {{1, 0}, {0, 1}, {1, -1}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {0, 2}, {0, 3}, {3, 0}, {3, 1}, {3, 2}, {3, 3}};
    // number of sectors
    //int sectorNumber = static_cast<int>(BSNumbers.size());

    // LHS matrix for the linear equation system
    Eigen::MatrixXd LHS = MatLHS(BSNumbers, DOrders, DOrdersMuZero, muB, muS, CInvContainer, CInvMuZero);

    // RHS vector for the linear equation system
    Eigen::VectorXd RHS = VecRHS(BSNumbers, DOrders, DOrdersMuZero, yMat, yMatMuZero, muB, muS, CInvContainer, CInvMuZero);

    // solving the linear equqation system for fitted coefficients
    Eigen::VectorXd coeffVector = (LHS).fullPivLu().solve(RHS);

    // chi squared value
    //double chiSq = ChiSq(BSNumbers, DOrders, yMat, muB, muS, CInvContainer, coeffVector);
    // number of degrees of freedom
    //int ndof = NDoF(CInvContainer, coeffVector);

    // error estimation via jackknife method
    std::vector<Eigen::VectorXd> JCK_RHS(jckNum);
    for (int i = 0; i < jckNum; i++)
    {
        // y data matrix for jackknife fits
        Eigen::MatrixXd yMatJCK(numOfQs, N - 1);
        yMatJCK.row(0) = imZBJCKs.col(i).segment(1, N - 1);
        yMatJCK.row(1) = imZSJCKs.col(i).segment(1, N - 1);

        // y data matrix for jackknife fits at mu = 0
        Eigen::MatrixXd yMatJCKMuZero(numOfQsMuZero, 1);
        yMatJCKMuZero(0, 0) = ZBBJCKs.col(i)(0);
        yMatJCKMuZero(1, 0) = ZBSJCKs.col(i)(0);
        yMatJCKMuZero(2, 0) = ZSSJCKs.col(i)(0);
        yMatJCKMuZero(3, 0) = ZBBBBJCKs.col(i)(0);
        yMatJCKMuZero(4, 0) = ZSSSSJCKs.col(i)(0);
        yMatJCKMuZero(5, 0) = ZBSSSJCKs.col(i)(0);
        yMatJCKMuZero(6, 0) = ZBBSSJCKs.col(i)(0);
        yMatJCKMuZero(7, 0) = ZBBBSJCKs.col(i)(0);

        // RHS vectors from jackknife samples
        JCK_RHS[i] = VecRHS(BSNumbers, DOrders, DOrdersMuZero, yMatJCK, yMatJCKMuZero, muB, muS, CInvContainer, CInvMuZero);
    }
    // fit with jackknife samples
    std::vector<Eigen::VectorXd> JCK_coeffVector(jckNum);
    for (int i = 0; i < jckNum; i++)
    {
        JCK_coeffVector[i] = (LHS).fullPivLu().solve(JCK_RHS[i]);
    }
    // estimate error from jackknife fits
    Eigen::VectorXd errorVec = JCKFitErrorEstimation(coeffVector, JCK_coeffVector);

    // results and tests
    // fit quality tests
    //std::cout << "\nchiSq = " << chiSq << std::endl;
    //std::cout << "ndof = " << ndof << std::endl;
    //std::cout << "AIC = " << AIC_weight(chiSq, ndof) << std::endl;
    //std::cout << "Q = " << Q_weight(chiSq, ndof) << std::endl;

    // write result coefficients to screen
    /*
    std::cout << "\nFitted parameters:" << std::endl;
    for (int coeffIndex = 0; coeffIndex < sectorNumber; coeffIndex++)
    {
        std::cout << "{" << BSNumbers[coeffIndex].first << " , " << BSNumbers[coeffIndex].second << "}: " << coeffVector(coeffIndex) << " +/- " << errorVec(coeffIndex) << std::endl;
    }
    */

    for (int i = 0; i < coeffVector.size(); i++)
    {
        std::cout << coeffVector[i] << " ";
        for (int iJCK = 0; iJCK < jckNum; iJCK++)
        {
            std::cout << JCK_coeffVector[iJCK](i) << " ";
        }
        std::cout << std::endl;
    }
}