#include "PuFC.h"
using namespace std;

PuFC::PuFC(Plot *p)
{
    Name = "PuFC";
    cout << endl << "Creating fission chamber " << Name << endl;
    CommentFlag = kTRUE;
    InitVar();
    InitPuVar();
    SimPath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_5E6.root";
    SetDraw(p);

    pHFG = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root", "NIF");
    pHBG = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root", "SB");
    pHSF = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root", "SF");
    pHFG->SetNeutronField(Yield, DYield, MonitorFG, DMonitorFG, 1500, 1);
    pHBG->SetNeutronField(Yield, DYield, MonitorBG, DMonitorBG, 1500, 1);

    cout << "Created " << Name << endl;
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
    PuSFT2 = 6.77E10 * 365.24*24*60*60; // Pu-242 spontaneaus fission half-life period in s
    DPuSFT2 = 7E8 * 365.24*24*60*60;
    MonitorFG = 27492079+33478370+30916955+54792679;
    DMonitorFG = MonitorFG * 0.0014;
    MonitorBG = 3623069+28646614+4757385;
    DMonitorBG = MonitorBG * 0.0015;

    Double_t L = 1500, DL = 1;
//    Double_t T[] = {0.936084, 0.942565, 0.949147, 0.955535, 0.962355, 0.968932, 0.974894, 0.981514};
//    Double_t DT[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//    Double_t S[] = {0.0999258, 0.0918293, 0.0852001, 0.0811008, 0.0807139, 0.0770468, 0.0685954, 0.0604568};
//    Double_t DS[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < NumHist; i++)
    {// Distances source-deposit
        sd[i] = L + 191 - 20.8 * i;
        Dsd[i] = DL;
    }
    cout << "Done: Pu variables" << endl;
}


void PuFC::SetDraw(Plot *p)
{
    if (p == 0)
        return;
    plot = p;
    pHFG->SetDraw(p);
    pHBG->SetDraw(p);
    pHSF->SetDraw(p);
    Draw = kTRUE;
}


void PuFC::HardCodedThresholds()
{
    Double_t avCut[] = {907, 853, 896, 849, 979, 892, 904, 837};

    cout << endl << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " <<
                       avCut[3] << " " << avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
}


void PuFC::AnalyzeQDC()
{// Calculate the average QDC minimum position
    if (CommentFlag)
        cout << endl << "Calculate QDC thesholds" << endl;
    Double_t avCut[NumHist];
    Double_t DavCut[NumHist];

    pHFG->DoAnalyzeQDC();
    pHBG->DoAnalyzeQDC();
    pHSF->DoAnalyzeQDC();

    for (int i = 0; i < NumHist; i++)
    {
        Double_t w[3]; // weights
        w[0] = pHFG->GetNevents(i);
        w[1] = pHBG->GetNevents(i);
        w[2] = pHSF->GetNevents(i);
        Double_t wSum = w[0] + w[1] + w[2];
        avCut[i] = (w[0] * pHFG->CutQDC[i] +
                    w[1] * pHBG->CutQDC[i] +
                    w[2] * pHSF->CutQDC[i]) / wSum;
        DavCut[i] = (w[0] * pow(pHFG->CutQDC[i] - avCut[i], 2) +
                     w[1] * pow(pHBG->CutQDC[i] - avCut[i], 2) +
                     w[2] * pow(pHSF->CutQDC[i] - avCut[i], 2) )/ wSum;
        if (CommentFlag)
            cout << " ch " << i+1 << ": " << avCut[i] << "+-" << DavCut[i] << " (" <<
                    pHFG->CutQDC[i] << ", " << pHBG->CutQDC[i] << ", " << pHSF->CutQDC[i] << ")" << endl;
    }
    cout << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " << avCut[3] << " " <<
                                  avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
    cout << "Done: QDC thresholds" << endl;
    DoneQDC = kTRUE;
}


void PuFC::AnalyzeDtBG()
{   // Calculate average fission background per livetime and Dt bin
    if (!DoneLimits)
        GetLimits();
    cout << endl << "Analyzing TimeDiff background..." << endl;

    Double_t tFG = pHFG->t_live;
    Double_t tBG = pHBG->t_live;
    Double_t tSF = pHSF->t_live;

    for (int i = 0; i < NumHist; i++)
    {
        // Integrate background
        Double_t counts = pHFG->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHFG->pHDtG[i]->Integral(lim[1][i], lim[2][i]) +
                          pHBG->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHBG->pHDtG[i]->Integral(lim[1][i], lim[2][i]) +
                          pHSF->pHDtG[i]->Integral(lim[0][i], lim[3][i]);
        Double_t tBin = (tFG + tBG) * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]) +
                        tSF * (lim[3][i] - lim[0][i] + 1); // Counted bins * live time
        avBg[i] = counts / tBin;
        DavBg[i] = sqrt(counts) / tBin;
        if (CommentFlag)
            cout << " ch " << i+1 << "  t " << tFG << " " << tBG << " " << tSF << "  counts " << counts << "  bins " << lim[1][i]- lim[0][i] + lim[3][i] - lim[2][i] << " " << lim[3][i] - lim[0][i] + 1 << endl <<
                    "  avBg " << avBg[i] << "+-" << DavBg[i] << endl;
    }
    DoneDtBG = kTRUE;
    cout << "Done: TimeDiff background" << endl;
}


void PuFC::GetNatoms()
{
    if (!DoneDt)
        AnalyzeDt();
    cout << endl << "Calculating effective number of Pu atoms..." << endl;
    cout << " T2SF(Pu-242) = " << PuSFT2 << "+-" << DPuSFT2 << endl;
    for (int i = 0; i < NumHist; i++)
    {
        Double_t t_live = pHFG->t_live + pHBG->t_live + pHSF->t_live;
        Double_t nSF = pHFG->pHDtG[i]->Integral() - nFG[i] +
                       pHBG->pHDtG[i]->Integral() - nBG[i] +
                       pHSF->pHDtG[i]->Integral();
        nAtoms[i] = nSF / t_live * PuSFT2 / log(2.0);
        DnAtoms[i] = sqrt(nSF) / t_live * PuSFT2 / log(2.0);
        cout << " Ch " << i+1 << ", eff. atoms " << nAtoms[i] << "+-" << DnAtoms[i] << endl;
    }
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


void PuFC::ExpTrans()
{
    if (!DoneDt)
        AnalyzeDt();
    if (!DoneDtBG)
        AnalyzeDtBG();
    cout << endl << "Experimental transmission" << endl;

    Double_t uT[NumCh];
    Double_t DuT[NumCh];
    Double_t T[NumCh];
    Double_t DT[NumCh];
    for (int i = 0; i < NumCh; i++)
    {
        Double_t nSF = pHFG->pHDtG[i]->Integral() - nFG[i];
        uT[i] = nFG[i] / nSF; // ExpT := N(n,f) / N(sf)
        DuT[i] = sqrt( pow(DnFG[i] / nFG[i], 2) +
                       1.0 / nSF ) * uT[i];
        T[i] = uT[i] * pow(sd[i] / sd[7], 2);
        DT[i] = DuT[i] * pow(sd[i] / sd[7], 2);
        if (CommentFlag)
            cout << " Ch " << i+1 << ", T " << uT[i] << "+-" << DuT[i] << " / " << T[i] << "+-" << DT[i] << endl;
    }
    if (!Draw)
        return;
    plot->ExpT(uT, DuT, T, DT);
}
