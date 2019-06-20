#include "UFC.h"
using namespace std;

UFC::UFC(Bool_t draw)
{
    Name = "UFC";
    cout << endl << "Creating fission chamber " << Name << endl;
    CommentFlag = kFALSE;
    InitVar(draw);
    InitUVar();
/*
    12371470
    12876188
    14226964
    27941726
    9962706

    13815794
*/

//    FgRuns = 1;
//    pHFG[0] = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root", "Open", "UFC", 0);
//    pHFG[0]->SetNeutronField(72661227.5811, 0.0014, 61012.07);

    FgRuns = 3;
    pHFG[0] = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_FG_MS20_2.root", "Open", "UFC", 0);
    pHFG[0]->SetNeutronField(9509137.9391, 0.0014, 7880.71);
    pHFG[1] = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_FG_MS20_3.root", "Open", "UFC", 1);
    pHFG[1]->SetNeutronField(12371469.5909, 0.0014, 11193.57);
    pHFG[2] = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_FG_MS20_4.root", "Open", "UFC", 2);
    pHFG[2]->SetNeutronField(12876188.3535, 0.0014, 10940.15);
    pHFG[3] = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_FG_MS21_2.root", "Open", "UFC", 3);
    pHFG[3]->SetNeutronField(27941726.1118, 0.0014, 22636.3);

    BgRuns = 1;
    pHBG[0] = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root", "SB", "UFC", 0);
    pHBG[0]->SetNeutronField(14226964.0879, 0.0014, 11456.44);
    cout << "Created " << Name << endl;
}


UFC::~UFC()
{
    
}


void UFC::InitUVar()
{
    cout << endl << "Initializing U variables.." << endl;
    // isotope fractions
    Double_t frac233 = 0.0000011;
    Double_t Dfrac233 = 0.0000005;
    Double_t frac234 = 0.00459;
    Double_t Dfrac234 = 0.00005;
    frac235 = 0.904;
    Dfrac235 = 0.005;
    Double_t frac236 = 0.00401;
    Double_t Dfrac236 = 0.00005;
    frac238 = 0.0912;
    Dfrac238 = 0.0006;
    Double_t norm = frac233 + frac234 + frac235 + frac236 + frac238;
    // average atomic mass
    M = (frac233 * 233.039634367 + frac234 * 234.040950370 + frac235 * 235.043928190 + frac236 * 236.045566201 + frac238 * 238.050786996) / norm;
    DM = sqrt( pow(Dfrac233 * -2, 2) +
               pow(Dfrac234 * -1, 2) +
               pow(Dfrac235 * 0.043928190, 2) +
               pow(Dfrac236 * 1, 2) +
               pow(Dfrac238 * 3, 2) );
    cout << "Average atomic mass: " << M << "+-" << DM << endl << "norm: " << norm << endl;

    // cross section standard
    sigma238 = 1.235; // in barns
    Dsigma238 = 0.01 * sigma238;

//    MonitorFG = 9509138+12371470+12876188+27941726+9962706;
//    DMonitorFG = MonitorFG * 0.0014;
//    MonitorBG = 14226964.0879;
//    DMonitorBG = MonitorBG * 0.0014;
    
    Double_t L = 1500, DL = 1;
    Double_t ema[] = {365.2, 396.2, 400.4, 393, 403.9, 403.7, 406.6, 397.1}; // deposits' efficient Areal mass density. Unit: 10^-6 g / cm^2
    Double_t Dema[] = {1.7, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};
    for (int i = 0; i < NumHist; i++)
    {// Distances source-deposit
        sd[i] = L + 119 - 10.8 * i;
        Dsd[i] = DL;
        emA[i] = ema[i] * 1.E-8; // Unit: g/mm^2
        DemA[i] = Dema[i] * 1.E-8;
    }
    cout << "Done: U variables" << endl;
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

    for (Int_t j = 0; j < FgRuns; j++)
        pHBG[j]->DoAnalyzeQDC();
    for (Int_t j = 0; j < BgRuns; j++)
        pHFG[j]->DoAnalyzeQDC();

    for (int i = 0; i < NumHist; i++)
    {
        Double_t w[FgRuns + BgRuns]; // weights
        Double_t wSum = 0;
        avCut[i] = 0;
        for (Int_t j = 0; j < FgRuns; j++)
        {
            w[j] = pHFG[j]->GetNevents(i);
            avCut[i] += w[j] * pHFG[j]->CutQDC[i];
            wSum += w[j];
        }
        for (Int_t j = 0; j < BgRuns; j++)
        {
            w[FgRuns + j] = pHBG[j]->GetNevents(i);
            avCut[i] += w[FgRuns + j] * pHBG[j]->CutQDC[i];
            wSum += w[FgRuns + j];
        }
        avCut[i] /= wSum;
        DavCut[i] = 0;
        for (Int_t j = 0; j < FgRuns; j++)
            DavCut[i] += w[j] * pow(pHFG[j]->CutQDC[i] - avCut[i], 2);
        for (Int_t j = 0; j < BgRuns; j++)
            DavCut[i] += w[FgRuns + j] * pow(pHBG[j]->CutQDC[i] - avCut[i], 2);
        DavCut[i] = sqrt(DavCut[i] / wSum);
        if (CommentFlag)
            cout << " ch " << i+1 << ": " << avCut[i] << "+-" << DavCut[i] << endl;
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
    cout << endl << "Analyzing ToF background..." << endl;

    if (CommentFlag)
        cout << "Ch   counts   tBin   avBg" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        // Integrate background
        Double_t counts = 0;
        Double_t tBin = 0;
        for (Int_t j = 0; j < FgRuns; j++)
        {
            counts += pHFG[j]->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHFG[j]->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
            tBin += pHFG[j]->t_live * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]);
        }
        for (Int_t j = 0; j < BgRuns; j++)
        {
            counts += pHBG[j]->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHBG[j]->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
            tBin += pHBG[j]->t_live * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]);
        }
        avBg[i] = counts / tBin;
        DavBg[i] = sqrt(counts) / tBin;

        if (CommentFlag)
            cout << " " << i+1 << "  " << counts << "  " << tBin << "  " << avBg[i] << "+-" << DavBg[i] << endl;
    }
    DoneDtBG = kTRUE;
    cout << "Done: ToF background" << endl;
}


void UFC::GetNatoms()
{
    if (!DoneDt)
        AnalyzeDt();
    cout << endl << "Calculating effective number of U atoms..." << endl;

    for (int i = 0; i < NumHist; i++)
    {
        nAtoms[i] = Area * emA[i] / (M * u); // units: mm^2 * g/mm^2 / g == 1
        DnAtoms[i] = nAtoms[i] * sqrt( pow(DArea / Area, 2) +
                                       pow(DemA[i] / emA[i], 2) +
                                       pow(DM / M, 2) );
        n235[i] = frac235 * nAtoms[i];
        Dn235[i] = n235[i] * sqrt( pow(Dfrac235 / frac235, 2) +
                                   pow(DemA[i] / emA[i], 2) );
        n238[i] = frac238 * nAtoms[i];
        Dn238[i] = n238[i] * sqrt( pow(Dfrac238 / frac238, 2) +
                                   pow(DemA[i] / emA[i], 2) ); // Statistical uncertainties
    }
    for (Int_t j = 0; j < FgRuns; j++)
        pHFG[j]->SetNatoms(nAtoms, DnAtoms);
    for (Int_t j = 0; j < BgRuns; j++)
        pHBG[j]->SetNatoms(nAtoms, DnAtoms);
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
    sIsoVec = sigma238 * frac238 / frac235;
    DsIsoVec = sIsoVec * sqrt( pow(Dsigma238 / sigma238, 2) +
                               pow(Dfrac238 / frac238, 2) +
                               pow(Dfrac235 / frac235, 2) );
//    fIsoVec = 1.0;
//    DfIsoVec = 0;
//    sIsoVec = 0;
//    DsIsoVec = 0;
    cout << "factor " << fIsoVec << "+-" << DfIsoVec << endl;
    cout << "subtrahend " << sIsoVec << "+-" << DsIsoVec << " barn" << endl;

    cout << "Done: IsoVec" << endl;
    DoneIso = kTRUE;
}


void UFC::GetExpT()
{
    cout << endl << "UFC::GetExpT() experimental transmission" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        uT[i] = nFG[i] / emA[i];
        DuT[i] = sqrt( pow(DnFG[i] / nFG[i], 2) +
                       pow(DemA[i] / emA[i], 2) ) * uT[i];
    }
}


void UFC::Calibrate()
{
    if (!DoneDt)
        AnalyzeDt();
    if (!DoneSimFg)
        GetSimFg();
    cout << endl << "U atom number Calibration..." << endl;

    // IAEA cross section standards
    Double_t sCS235 = 2.12623431;
    Double_t DsCS235 = 0.015082;
    Double_t sCS238 = 1.23462382;
    Double_t DsCS238 = 0.0098885;

    // cross section of isotope mixture
    Double_t effCS = frac235 * sCS235 + frac238 * sCS238;
    Double_t DeffCS = sqrt( pow(Dfrac235 * sCS235, 2) +
                            pow(frac235 * DsCS235, 2) +
                            pow(Dfrac238 * sCS238, 2) +
                            pow(frac238 * DsCS238, 2) );

    Double_t area[NumCh];
    Double_t Darea[NumCh];

    cout << " Ch   raw atoms   corrected atoms   eff mass density [10^-6 g/cm^2]" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        Double_t rawAtoms = nFG[i] / (tFG * effCS * nFlux[i]) * 1.E22;
        nAtoms[i] = fTS[i] * rawAtoms;
        DnAtoms[i] = nAtoms[i] * sqrt( pow(DnFG[i] / nFG[i], 2) +
                              pow(DeffCS / effCS, 2) +
                              pow(DnFlux[i] / nFlux[i], 2) +
                              pow(DfTS[i] / fTS[i], 2) );
        area[i] = nAtoms[i] * M * u / emA[i];
        Darea[i] = DnAtoms[i] / nAtoms[i] * area[i];
        emA[i] = nAtoms[i] * M * u / Area * 1.E8;
        DemA[i] = emA[i] * DnAtoms[i] / nAtoms[i];
        Double_t Dsyst = emA[i] * sqrt( pow(DM / M, 2) +
                                       pow(DArea / Area, 2) );

        cout << " " << i+1 << "   " << rawAtoms << "   " << nAtoms[i] << "+-" << DnAtoms[i] << "   " << emA[i] << "+-" << DemA[i] << "(stat)+-" << Dsyst << "(syst)"<< endl;
    }
    cout << "Done: U atom number calibration" << endl;
    if (!DrawSingle)
        return;
    plot->CalibrateUFC(emA, DemA, area, Darea);
}
