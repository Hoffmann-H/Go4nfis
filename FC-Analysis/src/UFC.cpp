#include "UFC.h"
using namespace std;

UFC::UFC(Plot *p)
{
    Name = "UFC";
    cout << endl << "Creating fission chamber " << Name << endl;
    CommentFlag = kTRUE;
    InitVar();
    InitUVar();
    SimPath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/UFC_Open_06052019_1730.root";
    SetDraw(p);
    
    pHFG = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root", "NIF");
    pHBG = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root", "SB");
    pHFG->SetNeutronField(Yield, DYield, MonitorFG, DMonitorFG, 1500, 1);
    pHBG->SetNeutronField(Yield, DYield, MonitorBG, DMonitorBG, 1500, 1);
    cout << "Created " << Name << endl;
}


UFC::~UFC()
{
    
}


void UFC::InitUVar()
{
    cout << endl << "Initializing U variables.." << endl;
    frac235 = 0.904;
    Dfrac235 = 0.005;
    frac238 = 0.0912;
    Dfrac238 = 0.0006;
    sigma238 = 1.25E-22; // in mm^2
    Dsigma238 = 0.0307 * sigma238;
    MonitorFG = 9509138+12371470+12876188+27941726+9962706;
    DMonitorFG = MonitorFG * 0.0014;
    MonitorBG = 13815794;
    DMonitorBG = MonitorBG * 0.0014;
    
    Double_t L = 1500, DL = 1;
//    Double_t T[] = {0.941864, 0.947976, 0.954109, 0.960458, 0.96626, 0.97273, 0.979405, 0.985676};
//    Double_t DT[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//    Double_t S[] = {0.085898, 0.0790897, 0.0745957, 0.0703662, 0.0661441, 0.0613817, 0.05615, 0.0506034};
//    Double_t DS[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < NumHist; i++)
    {// Distances source-deposit
        sd[i] = L + 119 - 10.8 * i;
        Dsd[i] = DL;
    }
    cout << "Done: U variables" << endl;
}


void UFC::SetDraw(Plot *p)
{
    if (p == 0)
        return;
    plot = p;
    pHFG->SetDraw(p);
    pHBG->SetDraw(p);
    Draw = kTRUE;
}


void UFC::HardCodedThresholds(){
    Double_t avCut[] = {212, 470, 207, 206, 199, 152, 184, 124};

    cout << endl << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " <<
                       avCut[3] << " " << avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
}


void UFC::AnalyzeQDC()
{// Calculate the average QDC minimum position
    if (CommentFlag)
        cout << endl << "Calculate QDC thesholds" << endl;
    Double_t avCut[NumHist];
    Double_t DavCut[NumHist];

    pHFG->DoAnalyzeQDC();
    pHBG->DoAnalyzeQDC();

    for (int i = 0; i < NumHist; i++)
    {
        Double_t w[2]; // weights
        w[0] = pHFG->GetNevents(i);
        w[1] = pHBG->GetNevents(i);
        Double_t wSum = w[0] + w[1];
        avCut[i] = (w[0] * pHFG->CutQDC[i] +
                    w[1] * pHBG->CutQDC[i]) / wSum;
        DavCut[i] = (w[0] * pow(pHFG->CutQDC[i] - avCut[i], 2) +
                     w[1] * pow(pHBG->CutQDC[i] - avCut[i], 2) )/ wSum;
        if (CommentFlag)
            cout << " ch " << i+1 << ": " << avCut[i] << "+-" << DavCut[i] << " (" <<
                    pHFG->CutQDC[i] << ", " << pHBG->CutQDC[i] << ")" << endl;
    }
    cout << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " << avCut[3] << " " <<
                                  avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
    cout << "Done: QDC thresholds" << endl;
    DoneQDC = kTRUE;
}


void UFC::AnalyzeDtBG()
{   // Calculate average fission background per livetime and Dt bin
    if (!DoneLimits)
        GetLimits();
    cout << endl << "Analyzing TimeDiff background..." << endl;

    Double_t tFG = pHFG->t_live;
    Double_t tBG = pHBG->t_live;

    for (int i = 0; i < NumHist; i++)
    {
        // Integrate background
        Double_t counts = pHFG->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHFG->pHDtG[i]->Integral(lim[1][i], lim[2][i]) +
                          pHBG->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHBG->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
        Double_t tBin = (tFG + tBG) * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]); // Counted bins * live time
        avBg[i] = counts / tBin;
        DavBg[i] = sqrt(counts) / tBin;
        if (CommentFlag)
            cout << " ch " << i+1 << "  t " << tFG << " " << tBG << "  counts " << counts << "  bins " << lim[1][i]- lim[0][i] + lim[3][i] - lim[2][i] << " " << endl <<
                    "  avBg " << avBg[i] << "+-" << DavBg[i] << endl;
    }
    DoneDtBG = kTRUE;
    cout << "Done: TimeDiff background" << endl;
}


void UFC::GetNatoms()
{
    if (!DoneDt)
        AnalyzeDt();
    cout << endl << "Calculating effective number of U atoms..." << endl;

    Double_t frac234 = 0.00459; // isotope fractions
    Double_t Dfrac234 = 0.00005;
    Double_t frac236 = 0.00401;
    Double_t Dfrac236 = 0.00005;
    Double_t emA[] = {365.2, 396.2, 400.4, 393, 403.9, 403.7, 406.6, 397.1}; // deposits' efficient Areal mass density. Unit: 10^-6 g / cm^2
    Double_t DemA[] = {1.7, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};
    Double_t M = 235 - 1*frac234 + 1*frac236 + 3*frac238;
    Double_t DM = sqrt( pow(Dfrac234 * 1, 2) +
                        pow(Dfrac235 * 0, 2) +
                        pow(Dfrac236 * 1, 2) +
                        pow(Dfrac238 * 3, 2) );

    for (int i = 0; i < NumHist; i++)
    {
        nAtoms[i] = Area * emA[i] / (M * u) * 1.E-8;
        DnAtoms[i] = nAtoms[i] * sqrt( pow(DArea / Area, 2) +
                                       pow(DemA[i] / emA[i], 2) +
                                       pow(DM / M, 2) );
    }
    cout << "Done: effective number of U atoms" << endl;
    DoneNatoms = kTRUE;
}


void UFC::IsoVec()
{
    if (!DoneNatoms)
        GetNatoms();
    if (!DoneDt)
        AnalyzeDt();
    cout << endl << "UFC::IsoVec() calculating isotope vector correction..." << endl;

    fIsoVec = 1.0 / frac235;
    DfIsoVec = Dfrac235 / frac235 * fIsoVec;
    cout << "factor " << fIsoVec << "+-" << DfIsoVec << endl;
    sIsoVec = sigma238 * frac238 / frac235;
    DsIsoVec = sIsoVec * sqrt( pow(Dsigma238 / sigma238, 2) +
                               pow(Dfrac238 / frac238, 2) +
                               pow(Dfrac235 / frac235, 2) );
    cout << "subtrahend " << sIsoVec * 1.E22 << "+-" << DsIsoVec * 1.E22 << " barn" << endl;

    cout << "Done: IsoVec" << endl;
    DoneIso = kTRUE;
}


void UFC::ExpTrans()
{
    cout << "UFC::ExpTrans() experimental transmission" << endl;
}


//void UFC::CalcCS(Int_t i)
//{
//    // nInc[]: transmission corrected incident fluence
//    cout << " UFC::CalcCS(" << i << ")" << endl;
//    CS[i] = nfDirect[i] / (frac235 * nAtoms[i] * nInc[i]) - sigma238 * frac238 / frac235;
//    DCS[i] = sqrt( pow(DnfDirect[i] / (frac235 * nAtoms[i] * nInc[i]), 2) +
//                   pow(CS[i] * Dfrac235 / frac235, 2) +
//                   pow(nfDirect[i] * DnAtoms[i] / (frac235 * nAtoms[i] * nAtoms[i] * nInc[i]), 2) +
//                   pow(nfDirect[i] * DnInc[i] / (nAtoms[i] * nInc[i] * nInc[i]), 2) +
//                   pow(Dsigma238 * frac238 / frac235, 2) +
//                   pow(sigma238 * Dfrac238 / frac235, 2) );
//}

























