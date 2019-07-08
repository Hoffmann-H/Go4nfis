#include "PuFC.h"
using namespace std;

PuFC::PuFC(Bool_t draw)
{
    Name = "PuFC";
    cout << endl << "Creating fission chamber " << Name << endl;
    CommentFlag = kTRUE;
    InitVar(draw);
    InitPuVar();

/////////////////////////////////////////////////////////////////////////////////////////87976.7);//
//    pR[0] = new Run("PuFC", "Open", 0);
//    FgMon = pR[0]->SetNeutronField(27492078.81+33478369.95+30916955.27+54792678.94, 0.0014, 16794.37+20286.24+18736.55+34490.15);
//    UseHist("NIF", "Open", 0);
//    FgRuns = 1;
//    pR[1] = new Run("PuFC", "SB", 0);
//    BgMon = pR[1]->SetNeutronField(3623068.621401+28646613.7419+4757385.125897, 0.0015, 2287.18+18331.64+3098.89);
//    UseHist("SB", "SB", 1);
//    BgRuns = 1;
/////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
    Int_t k = 0;
    FgMon = 0;
    pR[k] = new Run("PuFC", "Open", k);
    FgMon += pR[k]->SetNeutronField(27492078.81, 0.0014, 16794.37);
//    UseHists(159, 164, "Open", k);
    UseHist("PuFC_FG_MS4", "Open", k);
    k++;
    pR[k] = new Run("PuFC", "Open", k);
    FgMon += pR[k]->SetNeutronField(33478369.95, 0.0014, 20286.24);
//    UseHists(165, 170, "Open", k);
    UseHist("PuFC_FG_MS5", "Open", k);
    k++;
    pR[k] = new Run("PuFC", "Open", k);
    FgMon += pR[k]->SetNeutronField(30916955.27, 0.0014, 18736.55);
//    UseHists(172, 177, "Open", k);
    UseHist("PuFC_FG_MS6", "Open", k);
    k++;
    pR[k] = new Run("PuFC", "Open", k);
    FgMon += pR[k]->SetNeutronField(54792678.94, 0.0014, 34490.15);
//    UseHists(179, 190, "Open", k);
    UseHist("PuFC_FG_MS7", "Open", k);
    k++;
    FgRuns = k;

    BgMon = 0;
//    pR[k] = new Run("PuFC", "SB", k - FgRuns);
//    FgMon += pR[k]->SetNeutronField(219560.406192, 0.0025, 140.45);
////    UseHists(191, 191, "SB", k);
//    UseHist("PuFC_BG_MS8", "SB", k);
//    k++;
    pR[k] = new Run("PuFC", "SB", k - FgRuns);
    BgMon += pR[k]->SetNeutronField(3623068.621401, 0.0015, 2287.18);
//    UseHists(192, 192, "SB", k);
    UseHist("PuFC_BG_MS9", "SB", k);
    k++;
    pR[k] = new Run("PuFC", "SB", k - FgRuns);
    FgMon += pR[k]->SetNeutronField(28646613.7419, 0.0014, 18331.64);
//    UseHists(194, 201, "SB", k);
    UseHist("PuFC_BG_MS10", "SB", k);
    k++;
    pR[k] = new Run("PuFC", "SB", k - FgRuns);
    FgMon += pR[k]->SetNeutronField(4757385.125897, 0.0014, 3098.89);
//    UseHists(202, 203, "SB", k);
    UseHist("PuFC_BG_MS11", "SB", k);
    k++;
    BgRuns = k - FgRuns;
/////////////////////////////////////////////////////////////////////////////////////////

    pRSF = new Run("PuFC","SF", 0);
    pHSF = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root", "SF", Name, "SF", FgRuns + BgRuns);

    RegisterHists();
    cout << " Run " << pRSF->Name << endl;
    cout << "  Hist " << pHSF->Name << endl;

    cout << endl << "Created " << Name << endl;
}


//void PuFC::plt()
//{// So funktioniert die Übergabe an Plot!
//    Double_t x1[NumCh];
//    Double_t y1[NumCh];
//    Double_t x2[NumCh];
//    Double_t y2[NumCh];
//    for (int i = 0; i < NumCh; i++)
//    {
//        x1[i] = 4*i;
//        y1[i] = i;
//        x2[i] = 4-i;
//        y2[i] = 0.5;
//    }
//    plot->ExpT(x1, y1, x2, y2);
//}


PuFC::~PuFC()
{

}


void PuFC::InitPuVar()
{
    cout << endl << "Initializing Pu variables.." << endl;
    PuSFT2 = 6.76E10 * 365.24*24*60*60; // Pu-242 spontaneaus fission half-life period in s
    DPuSFT2 = 7E8 * 365.24*24*60*60;
//    MonitorFG = 146680082.972163; // MS#4 - MS#7
//    DMonitorFG = 0.0014;
//    MonitorBG = 3623069+28646614+4757385;
//    DMonitorBG = 0.0015;

    Double_t L = 1500, DL = 1;
    for (int i = 0; i < NumHist; i++)
    {// Distances source-deposit
        sd[i] = L + 191 - 20.8 * i;
        Dsd[i] = DL;
    }
    cout << "Done: Pu variables" << endl;
}


void PuFC::HardCodedThresholds()
{
    Double_t avCut[] = {907, 853, 896, 849, 979, 892, 904, 837};

    cout << endl << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " <<
                       avCut[3] << " " << avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
}


//void PuFC::AnalyzeQDC()
//{// Calculate the average QDC minimum position
//    if (CommentFlag)
//        cout << endl << "Calculate QDC thesholds" << endl;
//    Double_t avCut[NumHist];
//    Double_t DavCut[NumHist];

//    for (Int_t j = 0; j < FgRuns; j++)
//        pHBG[j]->DoAnalyzeQDC();
//    for (Int_t j = 0; j < BgRuns; j++)
//        pHFG[j]->DoAnalyzeQDC();
//    pHSF->DoAnalyzeQDC();

//    for (int i = 0; i < NumHist; i++)
//    {
//        Double_t w[FgRuns + BgRuns + 1]; // weights
//        Double_t wSum = 0;
//        avCut[i] = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//        {
//            w[j] = pHFG[j]->GetNevents(i);
//            avCut[i] += w[j] * pHFG[j]->CutQDC[i];
//            wSum += w[j];
//        }
//        for (Int_t j = 0; j < BgRuns; j++)
//        {
//            w[FgRuns + j] = pHBG[j]->GetNevents(i);
//            avCut[i] += w[FgRuns + j] * pHBG[j]->CutQDC[i];
//            wSum += w[FgRuns + j];
//        }
//        w[FgRuns + BgRuns] = pHSF->GetNevents(i);
//        avCut[i] += w[FgRuns + BgRuns] * pHSF->CutQDC[i];
//        wSum += w[FgRuns + BgRuns];
//        avCut[i] /= wSum;
//        DavCut[i] = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//            DavCut[i] += w[j] * pow(pHFG[j]->CutQDC[i] - avCut[i], 2);
//        for (Int_t j = 0; j < BgRuns; j++)
//            DavCut[i] += w[FgRuns + j] * pow(pHBG[j]->CutQDC[i] - avCut[i], 2);
//        DavCut[i] += w[FgRuns + BgRuns] * pow(pHSF->CutQDC[i] - avCut[i], 2);
//        DavCut[i] = sqrt(DavCut[i] / wSum);
//        if (CommentFlag)
//            cout << " ch " << i+1 << ": " << avCut[i] << "+-" << DavCut[i] << endl;
//    }
//    cout << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " << avCut[3] << " " <<
//                                  avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
//    cout << "Done: QDC thresholds" << endl;
//    DoneQDC = kTRUE;
//}


//void PuFC::AnalyzeDtBG()
//{   // Calculate average fission background per livetime and Dt bin
//    if (!DoneLimits)
//        GetLimits();
//    cout << endl << "Analyzing ToF background..." << endl;
//    if (CommentFlag)
//        cout << " Ch   FG   BG   SF   avBg" << endl;
//    for (int i = 0; i < NumHist; i++)
//    {
//        // Integrate background
//        Double_t counts[4];
//        counts[0] = 0; counts[1] = 0; counts[2] = 0;
//        Double_t tBin[4];
//        tBin[0] = 0; tBin[1] = 0; tBin[2] = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//        {
//            counts[0] += pHFG[j]->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHFG[j]->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
//            tBin[0] += pHFG[j]->t_live * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]);
//        }
//        for (Int_t j = 0; j < BgRuns; j++)
//        {
//            counts[1] += pHBG[j]->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHBG[j]->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
//            tBin[1] += pHBG[j]->t_live * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]);
//        }
//        counts[2] += pHSF->pHDtG[i]->Integral(lim[0][i], lim[3][i]);
//        tBin[2] += pHSF->t_live * (lim[3][i] - lim[0][i] + 1); // Counted bins * live time
//        counts[3] = counts[0] + counts[1] + counts[2];
//        tBin[3] = tBin[0] + tBin[1] + tBin[2];
//        avBg[i] = counts[3] / tBin[3];
//        DavBg[i] = sqrt(counts[3]) / tBin[3];
//        if (CommentFlag)
//            cout << " " << i+1 << "  " << counts[0]/tBin[0] << "+-" << sqrt(counts[0])/tBin[0] << "  " << counts[1]/tBin[1] << "+-" << sqrt(counts[1])/tBin[1] << "  " << counts[2]/tBin[2] << "  " << avBg[i] << "+-" << DavBg[i] << endl;
//    }
//    DoneDtBG = kTRUE;
//    cout << "Done: ToF background" << endl;
//}


void PuFC::GetNatoms()
{
//    if (!DoneDt)
//        AnalyzeDt();
    cout << endl << "Calculating effective number of Pu atoms..." << endl;
    cout << " T2SF(Pu-242) = " << PuSFT2 << "+-" << DPuSFT2 << endl;
    Double_t eMinimum = 0.986;
    Double_t DeMinimum = 0.010;
    Double_t cAtoms[NumCh];
    Double_t DcAtoms[NumCh];
    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
    cout << "Ch   SF-rate   eff.Atoms   Atoms" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        Double_t nSF = 0;
        Double_t t_live = 0;
        for (Int_t j = 0; j < FgRuns + BgRuns; j++)
        {
            nSF += pH[j]->nSF[i];
            t_live += pH[j]->t_live;
        }
        nSF += pHSF->nSF[i];
        t_live += pHSF->t_live;
        nAtoms[i] = rSF_Toni[i] * PuSFT2 / log(2.0);
//        nAtoms[i] = nSF / t_live * PuSFT2 / log(2.0);
        DnAtoms[i] = sqrt(nSF) / t_live * PuSFT2 / log(2.0);
        cAtoms[i] = nAtoms[i] / eMinimum;
        DcAtoms[i] = cAtoms[i] * sqrt( pow(DnAtoms[i] / nAtoms[i], 2) + pow(DeMinimum / eMinimum, 2) );
        cout << " " << i+1 << "   " << nSF / t_live << "+-" << sqrt(nSF) / t_live << "   " << nAtoms[i] << "+-" << DnAtoms[i] << "   " << cAtoms[i] << "+-" << DcAtoms[i] << endl;
    }
    for (Int_t k = 0; k < FgRuns + BgRuns; k++)
        pR[k]->SetNatoms(nAtoms, DnAtoms);
    cout << "Done: effective number of Pu atoms" << endl;
    DoneNatoms = kTRUE;
}


void PuFC::IsoVec()
{
    cout << endl << "PuFC::IsoVec() defining isotope vector correction factor..." << endl;
    // PuFC: Do nothing.
    fIsoVec = 1;
    DfIsoVec = 0;
    sIsoVec = 0;
    DsIsoVec = 0;
    cout << "Done: IsoVec" << endl;
    DoneIso = kTRUE;
}


//void PuFC::GetExpT()
//{
//    cout << endl << "PuFC::GetExpT() experimental transmission" << endl;
//    for (int i = 0; i < NumCh; i++)
//    {
//        Double_t nSF = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//            nSF += pHFG[j]->nSF[i];
//        uT[i] = nFG[i] / nSF; // ExpT := N(n,f) / N(sf)
//        DuT[i] = sqrt( pow(DnFG[i] / nFG[i], 2) +
//                   1.0 / nSF ) * uT[i];
//    }
//}
