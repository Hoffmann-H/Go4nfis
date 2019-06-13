#include "Xsection.h"

using namespace std;

Xsection::Xsection(Bool_t draw)
{
    CommentFlag = kFALSE;
    Draw = draw;
    Pu = new PuFC(Draw);
    U = new UFC(Draw);
    simPu = new AnaSim("PuFC", 1, 0);
    simU = new AnaSim("UFC", 1, 0);
}


void Xsection::RelativeCS()
{
    cout << endl << "Relative Cross section 242Pu/235U..." << endl;
    Double_t cPu[NumCh];
    Double_t DcPu[NumCh];
    Double_t cU[NumCh];
    Double_t DcU[NumCh];
    Double_t nPu[NumCh];
    Double_t DnPu[NumCh];
    Double_t nU[NumCh];
    Double_t DnU[NumCh];
    Double_t flPu[NumCh];
    Double_t DflPu[NumCh];
    Double_t flU[NumCh];
    Double_t DflU[NumCh];
    Pu->AnalyzeDt();
    Pu->GetNatoms();
    simPu->Corrections();
    Double_t sum_cPu = 0;
    Double_t D2sum_cPu = 0;
    Double_t sum_nflPu = 0;
    Double_t D2sum_nflPu = 0;
    Double_t sum_sigmaPu = 0;
    Double_t D2sum_sigmaPu = 0;
    Double_t sigma, D2sigma, f, Df;
    for (Int_t i = 0; i < NumCh; i++)
    {
        cout << "Channel " << i+1 << endl;
        cPu[i] = Pu->nFG[i]; DcPu[i] = Pu->DnFG[i];
        sum_cPu += cPu[i]; D2sum_cPu += pow(DcPu[i], 2);
        cout << " Pu(n,f) counts: " << cPu[i] << "+-" << DcPu[i] << endl;
        nPu[i] = Pu->nAtoms[i]; DnPu[i] = Pu->DnAtoms[i];
        cout << " 242Pu atoms: " << nPu[i] << "+-" << DnPu[i] << endl;
        flPu[i] = Pu->pHFG->NeutronFluence[i]; DflPu[i] = Pu->pHFG->DstatNeutronFluence[i];
        sum_nflPu += nPu[i] * flPu[i]; D2sum_nflPu += pow(flPu[i] * DnPu[i], 2) + pow(nPu[i] * DflPu[i], 2);
        cout << " Pu n fluence: " << flPu[i] << "+-" << DflPu[i] << " mm^-2" << endl;
        f = simPu->F[i]; Df = simPu->DF[i];
        cout << " Pu T&S correction factor: " << f << "+-" << Df << endl;
        sigma = f * cPu[i] / (nPu[i] * flPu[i]) * 1.E22;
        D2sigma = pow(sigma, 2) * (pow(Df / f, 2) +
                                   pow(DcPu[i] / cPu[i], 2) +
                                   pow(DnPu[i] / nPu[i], 2) +
                                   pow(DflPu[i] / flPu[i], 2));
        sum_sigmaPu += sigma;
        D2sum_sigmaPu += D2sigma;
        cout << " Cross section: " << sigma << "+-" << sqrt(D2sigma) << endl;
    }
    U->AnalyzeDt();
    U->GetNatoms();
    simU->Corrections();
    Double_t sum_cU = 0;
    Double_t D2sum_cU = 0;
    Double_t sum_nflU = 0;
    Double_t D2sum_nflU = 0;
    Double_t sum_sigmaU = 0;
    Double_t D2sum_sigmaU = 0;
    for (Int_t i = 0; i < NumCh; i++)
    {
        cout << "Channel " << i+1 << endl;
        cU[i] = U->nFG[i]; DcU[i] = U->DnFG[i];
        sum_cU += cU[i]; D2sum_cU += pow(DcU[i], 2);
        cout << " U(n,f) counts: " << cU[i] << "+-" << DcU[i] << endl;
        nU[i] = U->n235[i]; DnU[i] = U->Dn235[i];
        cout << " 235U atoms: " << nU[i] << "+-" << DnU[i] << endl;
        flU[i] = U->pHFG->NeutronFluence[i]; DflU[i] = U->pHFG->DstatNeutronFluence[i];
        sum_nflU += nU[i] * flU[i]; D2sum_nflU += pow(flU[i] * DnU[i], 2) + pow(nU[i] * DflU[i], 2);
        cout << " U n fluence: " << flU[i] << "+-" << DflU[i] << " mm^-2" << endl;
        f = simU->F[i]; Df = simU->DF[i];
        cout << " U T&S correction factor: " << f << "+-" << Df << endl;
        Double_t sigma_raw = f * cU[i] / (nU[i] * flU[i]) * 1.E22 / U->frac235;
        sigma = sigma_raw - U->frac238 * U->sigma238 / U->frac235;
        D2sigma = pow(sigma, 2) * (pow(Df / f, 2) +
                                   pow(DcU[i] / cU[i], 2) +
                                   pow(DnU[i] / nU[i], 2) +
                                   pow(DflU[i] / flU[i], 2)) +
                  pow(U->Dfrac238 * U->sigma238 / U->frac235, 2) +
                  pow(U->Dfrac235 * U->frac238 * U->sigma238 / (U->frac235*U->frac235), 2);
        sum_sigmaU += sigma;
        D2sum_sigmaU += D2sigma;
        cout << " Cross section: " << sigma << "+-" << sqrt(D2sigma) << endl;
    }
    Double_t relSigmaRaw = sum_sigmaPu / sum_sigmaU;
    Double_t DrelSigmaRaw = relSigmaRaw * sqrt( D2sum_sigmaPu / pow(sum_sigmaPu, 2) + D2sum_sigmaU / pow(sum_sigmaU, 2));
    cout << "Quotient of averages: " << relSigmaRaw << "+-" << DrelSigmaRaw << endl;
}
/*
void Xsection::SetParam()
{
    // physics parameters
    u = 1.660539E-24; // atomic mass unit in g
    PuLit = 2120E-25; // 242Pu neutron-induced fission cross section in mm^2, original in mb
    DPuLit = 35E-25;
    ULit = 2.1; // 235U (n,f) cs in b
    DULit = 0.0315;
    PuSFT2 = 6.77E10 * 365.24*24*60*60; // Pu-242 spontaneaus fission half-life period in s
    DPuSFT2 = 7E8 * 365.24*24*60*60;
    Yield = 2.15882E4;
    DYield = 543.3;
    MonitorNIF = 27492079+33478370+30916955+54792679;
    DMonitorNIF = MonitorNIF * 0.0014;
    MonitorSB = 3623069+28646614+4757385;
    DMonitorSB = MonitorSB * 0.0015;
    MonitorUNIF = 9509138+12371470+12876188+27941726+9962706;
    DMonitorUNIF = MonitorUNIF * 0.0014;
    MonitorUSB = 13815794;
    DMonitorUSB = MonitorUSB * 0.0014;
//    AreaPuFC = 4300; // in mm^2
//    DAreaPuFC = 120;
    AreaUFC = 4300; // in mm^2
    DAreaUFC = 120;
    mPu = 0.03724; // in g
    DmPu = 0.00001;
    mU = 0.15801;
    DmU = 0.00001;
//    eSimGayther = 0.988;
//    DeSimGayther = 0.005;
    eSimMinimum = 0.986;
    DeSimMinimum = 0.010;
}

void Xsection::CalculateThresholds()
{ // Merge QDC spectra analyses. Calculate average extrema positions
//    cout << "Common cut positions" << endl;
    Double_t PuFC_Cut[NumHist];
    Double_t PuFC_DCut[NumHist];
    Double_t UFC_Cut[NumHist];
    Double_t UFC_DCut[NumHist];
    Double_t NIF_Ped, NIF_Cut, NIF_Max,
             SB_Ped, SB_Cut, SB_Max,
             SF_Ped, SF_Cut, SF_Max;
    /// method 1: Allow different relative minimum positions in deposits
    for (int i = 0; i < NumHist; i++)
    {
        // PuFC
        NIF_Cut = pHNIF->CutQDC[i];
        SB_Cut = pHSB->CutQDC[i];
        SF_Cut = pHUG->CutQDC[i];
        Double_t counts[] = {pHNIF->GetNevents(i), pHSB->GetNevents(i), pHUG->GetNevents(i)};
//        {pHNIF->t_live*(pHNIF->SFRate[i]+pHNIF->NIFRate[i]), // number of FF events for relative weighting
//                             pHSB->t_live*(pHSB->NIFRate[i]+pHSB->NIFRate[i]),
//                             pHUG->t_live*(pHUG->SFRate[i]+pHUG->NIFRate[i])};//
        Double_t sum = counts[0] + counts[1] + counts[2];
        Double_t w[] = {counts[0]/sum, counts[1]/sum, counts[2]/sum}; // weights
//        cout << "ch " << i+1 << "  " << w[0] << "  " << w[1] << "  " << w[2] << endl;
        PuFC_Cut[i] = w[0] * NIF_Cut + w[1] * SB_Cut + w[2] * SF_Cut;
        PuFC_DCut[i] = sqrt( w[0] * pow(NIF_Cut - PuFC_Cut[i], 2) +
                             w[1] * pow(SB_Cut - PuFC_Cut[i], 2) +
                             w[2] * pow(SF_Cut - PuFC_Cut[i], 2) );
//        cout << "      " << PuFC_Cut[i] << "+-" << PuFC_DCut[i] << endl;
        // UFC...
    }
    // Draw
    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, PuFC_Cut, xerr, PuFC_DCut);
    g1->SetNameTitle("gCut", "PuFC combined QDC minimum position vs channel");
    SaveToFile("PuFC/QDC", g1);
}


void Xsection::DoAnalyzeDt(string FC, Bool_t method)
{
    cout << endl << "=== " << FC << " in-scattering correction ===" << endl;
    cout << "SB peak integration method: ";
    switch (method) {
//    case 2:
//        cout << "Simulation" << endl;
//        break;
    case 1:
        cout << "Fit" << endl;
        break;
    case 0:
    default:
        cout << "Subtraction" << endl;
        break;
    }

    Bool_t PuFC = strcmp(FC.c_str(), "UFC");
    Hist *hNIF, *hSB;
    if (PuFC)
    {
        hNIF = pHNIF;
        hSB = pHSB;
    } else {
        hNIF = pHUNIF;
        hSB = pHUSB;
    }


    Double_t Dt_min = 63600;
    Double_t Dt_max = 78500;
    Double_t ToF_low, ToF_up, ChPerBin, BinOffset;
//    Double_t TimePerCh = 25.E-12; // 25 ps in seconds
    Double_t tNIF = hNIF->t_live;
    Double_t tSB = hSB->t_live;
    Double_t tUG; // only necessary for PuFC
    Double_t MonNIF = PuFC ? MonitorNIF : MonitorUNIF;
    Double_t DMonNIF = PuFC ? DMonitorNIF : DMonitorUNIF;
    Double_t MonSB = PuFC ? MonitorSB : MonitorUSB;
    Double_t DMonSB = PuFC ? DMonitorSB : DMonitorUSB;
    cout << "Times: " << tNIF << ", " << tSB << endl;
    cout << "Neutron monitor: " << MonNIF << ", " << MonSB << endl;
    Int_t lim[4];
    Double_t nNIF[NumHist];
    Double_t DnNIF[NumHist];
    Double_t D2nNIF[NumHist];
    Double_t nNIFfit[NumHist];
    Double_t DnNIFfit[NumHist];
    Double_t D2nNIFfit[NumHist];
    Double_t nSBsub[NumHist];
    Double_t DnSBsub[NumHist];
    Double_t D2nSBsub[NumHist];
    Double_t nSBfit[NumHist];
    Double_t DnSBfit[NumHist];
    Double_t D2nSBfit[NumHist];
    Double_t nSB[NumHist];
    Double_t DnSB[NumHist];
    Double_t D2nSB[NumHist];
    Double_t nSignal[NumHist];
    Double_t DnSignal[NumHist];
    Double_t nNIFtoMon[NumHist];
    Double_t DnNIFtoMon[NumHist];
    Double_t nSBtoMon[NumHist];
    Double_t DnSBtoMon[NumHist];
    Double_t InScatPart[NumHist];
    Double_t DInScatPart[NumHist];
    Double_t sum = 0;
    Double_t D2sum = 0;
    TH1I *pH1NIF, *pH1SB, *pH1UG;
    char name[64] = "";
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

//    cout << "ch   average Ug   Ug-rate   nNIF   nSB   corrected NIF   in-scat ratio" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        // get hand-made peak limits
        ToF_low = hNIF->GetPeakLow(i);
        ToF_up = hNIF->GetPeakUp(i);

        // set limits from 3 sigma to n sigma
        Double_t n = 3;
        Double_t center = 0.5 * (ToF_low + ToF_up);
        Double_t width = (ToF_up - ToF_low) / 6.0;
        ToF_low = center - n * width;
        ToF_up = center + n * width;

        // Get TimeDiff histograms
        pH1NIF = hNIF->pHDtG[i];
        pH1SB = hSB->pHDtG[i];
        if (PuFC)
            pH1UG = pHUG->pHDtG[i];
        // Calculate integration limit bins
        ChPerBin = pH1NIF->GetBinWidth(0);
        BinOffset = pH1NIF->GetBinCenter(0);
        lim[0] = (Dt_min - BinOffset) / ChPerBin;
        lim[1] = (ToF_low - BinOffset) / ChPerBin;
        lim[2] = (ToF_up - BinOffset) / ChPerBin;
        lim[3] = (Dt_max - BinOffset) / ChPerBin;
        if(CommentFlag)
        {
            cout << "Channel " << i+1 << endl << " Limits: ";
            for(int j = 0; j<4;j++)
                cout << pH1NIF->GetBinCenter(lim[j]) << " ";
            cout << endl;
        }

        //// Peak integration
        Int_t NIFPeakCount = pH1NIF->Integral(lim[1], lim[2]); // D2PeakCount = PeakCount
        Int_t SBPeakCount = pH1SB->Integral(lim[1], lim[2]);

        //// Fit constant background
        Double_t UgCount, DUgCount, avUg, DavUg;
        if (PuFC)
        {
            tUG = pHUG->t_live;
            UgCount = pH1NIF->Integral(lim[0], lim[3]) - NIFPeakCount +
                    pH1SB->Integral(lim[0], lim[3]) - SBPeakCount +
                    pH1UG->Integral(lim[0], lim[3]);
            DUgCount = sqrt(UgCount);
            // avUg: average underground count per bin and t_live
            avUg = UgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) + tUG * (lim[3] - lim[0] + 1) );
            DavUg = DUgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) + tUG * (lim[3] - lim[0] + 1) );
        } else {
            UgCount = pH1NIF->Integral(lim[0], lim[3]) - NIFPeakCount +
                               pH1SB->Integral(lim[0], lim[3]) - SBPeakCount;
            DUgCount = sqrt(UgCount);
            // avUg: average underground count per bin and t_live
            avUg = UgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) );
            DavUg = DUgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) );
        }
        if(CommentFlag)
            cout << " Underground per bin: " << avUg << "+-" << DavUg << endl;

        //// save underground to file
        // NIF
        Double_t x[] = {Dt_min, Dt_max};
        Double_t y[] = {avUg * tNIF, avUg * tNIF};
        TGraph* g0 = new TGraph(2, x, y);
        sprintf(name, "f%sNIFUg_%i", FC.c_str(), i+1);
        g0->SetNameTitle(name, "Pulse-heigt gated TimeDiff Underground Fit");
        sprintf(name, "%s/TimeDiff/Underground/NIF", FC.c_str());
        SaveToFile(name, g0);
        x[0] = ToF_low; x[1] = ToF_up;
        TGraph* g1 = new TGraph(2, x, y);
        sprintf(name, "f%sNIFUgPeak_%i", FC.c_str(), i+1);
        g1->SetNameTitle(name, "TimeDiff Peak area Underground");
        sprintf(name, "%s/TimeDiff/Underground/NIF", FC.c_str());
        SaveToFile(name, g1);
        // SB
        y[0] = avUg * tSB; y[1] = avUg * tSB;
        TGraph* g2 = new TGraph(2, x, y);
        sprintf(name, "f%sSBUgPeak_%i", FC.c_str(), i+1);
        g2->SetNameTitle(name, "TimeDiff Peak area Underground");
        sprintf(name, "%s/TimeDiff/Underground/SB", FC.c_str());
        SaveToFile(name, g2);
        x[0] = Dt_min; x[1] = Dt_max;
        TGraph* g3 = new TGraph(2, x, y);
        sprintf(name, "fPuFCSBUg_%i", i+1);
        g3->SetNameTitle(name, "Pulse-heigt gated TimeDiff Underground Fit");
        sprintf(name, "%s/TimeDiff/Underground/SB", FC.c_str());
        SaveToFile(name, g3);
        if (PuFC)
        { // UG
            y[0] = avUg * tUG; y[1] = avUg * tUG;
            TGraph* g4 = new TGraph(2, x, y);
            sprintf(name, "fPuFCUGUg_%i", i+1);
            g4->SetNameTitle(name, "Pulse-heigt gated TimeDiff Underground Fit");
            SaveToFile("PuFC/TimeDiff/Underground/UG", g4);
        }

        // number of time-correlated FF events in NIF measurement
        nNIF[i] = NIFPeakCount - (lim[2] - lim[1] + 1) * tNIF * avUg;
        D2nNIF[i] = NIFPeakCount + pow((lim[2] - lim[1] + 1) * tNIF * DavUg, 2);
        DnNIF[i] = sqrt(D2nNIF[i]);

        //// number of time-correlated FF events in SB measurement
        //// 1st method: subtraction of time-independent underground
        // number of time-correlated FF events in SB measurement
        nSBsub[i] = SBPeakCount - (lim[2] - lim[1] + 1) * tSB * avUg;
        D2nSBsub[i] = SBPeakCount + pow((lim[2] - lim[1] + 1) * tSB * DavUg, 2);
        DnSBsub[i] = sqrt( D2nSBsub[i] );
        if(CommentFlag)
        {
            cout << " NIF peak integral: " << NIFPeakCount << "+-" << sqrt(NIFPeakCount) << endl <<
                    " NIF peak underground: " << (lim[2] - lim[1] + 1) * tNIF * avUg << "+-" << (lim[2] - lim[1] + 1) * tNIF * DavUg << endl <<
                    " NIF peak content: " << nNIF[i] << "+-" << sqrt(D2nNIF[i]) << endl;
            cout << " SB peak integral: " << SBPeakCount << "+-" << sqrt(SBPeakCount) << endl <<
                    " SB peak underground: " << (lim[2] - lim[1] + 1) * tSB * avUg << "+-" << (lim[2] - lim[1] + 1) * tSB * DavUg << endl <<
                    " SB peak content (1st method): " << nSBsub[i] << "+-" << sqrt(D2nSBsub[i]) << endl;
        }

        //// 2nd method: Fit the peaks
        Double_t Dwidth, ampl, Dampl;
        Double_t pi = 3.1415926535898;
        sprintf(name, "f%sNIFDtPeak_%i", FC.c_str(), i+1);
        TF1* fNIF = new TF1(name, func_peak, Dt_min, Dt_max, 4);
        fNIF->SetNpx(1000);
//        TF1* fNIF = new TF1(name, "gaus", Dt_min, Dt_max, 4);
//        fNIF->SetRange(Dt_min, Dt_max);
        fNIF->SetParameters(100, center, width, tNIF * avUg);
        pH1NIF->Fit(name, "LQR");
        sprintf(name, "%s/TimeDiff/Signal/NIF/Fit", FC.c_str());
        SaveToFile(name, fNIF);
        ampl = (Double_t)fNIF->GetParameter(0);
        Dampl = (Double_t)fNIF->GetParError(0);
        center = (Double_t)fNIF->GetParameter(1);
        width = (Double_t)fNIF->GetParameter(2);
        Dwidth = (Double_t)fNIF->GetParError(2);
        nNIFfit[i] = sqrt(pi) / ChPerBin * ampl * width;
        D2nNIFfit[i] = pi / pow(ChPerBin, 2) * (pow(Dampl * width, 2) + pow(ampl * Dwidth, 2));
        DnNIFfit[i] = sqrt(D2nNIFfit[i]);

        sprintf(name, "f%sSBDtPeak_%i", FC.c_str(), i+1);
        TF1* fSB = new TF1(name, func_peak, Dt_min, Dt_max, 4);
        fSB->SetNpx(1000);
        fSB->SetRange(Dt_min, Dt_max);
        fSB->SetParameters(100, center, width, tSB * avUg);
//        fSB->FixParameter(1, center);
        fSB->FixParameter(2, width);
        pH1SB->Fit(name, "BLQR");
        sprintf(name, "%s/TimeDiff/Signal/SB/Fit", FC.c_str());
        SaveToFile(name, fSB);
        ampl = (Double_t)fSB->GetParameter(0);
        Dampl = (Double_t)fSB->GetParError(0);
        nSBfit[i] = sqrt(pi) / ChPerBin * ampl * width; // using width from NIF fit
        D2nSBfit[i] = pi / pow(ChPerBin, 2) * (pow(Dampl * width, 2) + pow(ampl * Dwidth, 2));
        DnSBfit[i] = sqrt(D2nSBfit[i]);
        if(CommentFlag)
            cout << " SB peak content (2nd method): " << nSBfit[i] << "+-" << sqrt(D2nSBfit[i]) << endl;

        //// In-Scattering correction
        // from data
        // Choose a method from subtraction, fit
        if (method == 0) {
            nSB[i] = nSBsub[i];
            D2nSB[i] = D2nSBsub[i];
        } else {
            nSB[i] = nSBfit[i];
            D2nSB[i] = D2nSBfit[i];
        }
        DnSB[i] = sqrt(D2nSB[i]);

        // norm signals to monitor counts
        nNIFtoMon[i] = nNIF[i] / MonNIF;
        DnNIFtoMon[i] = sqrt( pow(DnNIF[i] / MonNIF, 2) +
                              pow(nNIF[i] * DMonNIF, 2) / pow(MonNIF, 4) );
        nSBtoMon[i] = nSB[i] / MonSB;
        DnSBtoMon[i] = sqrt( pow(DnSB[i] / MonSB, 2) +
                              pow(nSB[i] * DMonSB, 2) / pow(MonSB, 4) );

        // subtract normed SB signals
        nSignal[i] = nNIF[i] - nSB[i] * (hNIF->NeutronFluence[i] / hSB->NeutronFluence[i]);
        DnSignal[i] = sqrt( D2nNIF[i] +
                            D2nSB[i] * pow(hNIF->NeutronFluence[i] / hSB->NeutronFluence[i], 2) +
                            pow(nSB[i] * (hNIF->DNeutronFluence[i] / hSB->NeutronFluence[i]), 2) +
                            pow(nSB[i] * (hNIF->NeutronFluence[i] * hSB->DNeutronFluence[i] / pow(hSB->NeutronFluence[i], 2)), 2) );

        // in-scattering portion
        InScatPart[i] = 100 * nSB[i] / nNIF[i] * MonNIF / MonSB;
        DInScatPart[i] = InScatPart[i] * sqrt( D2nSB[i] / pow(nSB[i], 2) +
                                                  D2nNIF[i] / pow(nNIF[i], 2) +
                                                  pow(DMonNIF / MonNIF, 2) +
                                                  pow(DMonSB / MonSB, 2) );
        // Average in-scattering portion
        sum += InScatPart[i];
        D2sum += pow(DInScatPart[i], 2);

        // scale SB count
        Double_t scaledSB = nSB[i] * MonNIF / MonSB;
        Double_t DscaledSB = scaledSB * sqrt( D2nSB[i] / pow(nSB[i], 2) +
                                              pow(DMonNIF / MonNIF, 2) +
                                              pow(DMonSB / MonSB, 2) );
        if (PuFC)
        {
            NIFRate[i] = nNIF[i] / tNIF;//nSignal[i] / tNIF; // nNIF: no SB-Correction. nSignal: with SB-Corr. Correct only once!
            DNIFRate[i] = DnNIF[i] / tNIF;//DnSignal[i] / tNIF;

            SFRate[i] = (pH1NIF->Integral() + pH1SB->Integral() + pH1UG->Integral() - nNIF[i] - nSBsub[i]) / (tNIF + tSB + tUG);
            DSFRate[i] = sqrt( pH1NIF->Integral() + pH1SB->Integral() + pH1UG->Integral() + D2nNIF[i] + D2nSBsub[i] ) / (tNIF + tSB + tUG);

            TGraphErrors* g6 = new TGraphErrors(NumHist, x, SFRate, xerr, DSFRate);
            g6->SetNameTitle("SF_Rate", "SF detection rate");
            sprintf(name, "%s/TimeDiff/Underground", FC.c_str());
            SaveToFile(name, g6);

            TGraphErrors* g7 = new TGraphErrors(NumHist, x, NIFRate, xerr, DNIFRate);
            g7->SetNameTitle("NIF_Rate", "In-scattering corrected NIF detection rate");
            sprintf(name, "%s/TimeDiff/Signal", FC.c_str());
            SaveToFile(name, g7);
        } else {
            UNIFRate[i] = nNIF[i] / tNIF;//nSignal[i] / tNIF;
            DUNIFRate[i] = DnNIF[i] / tNIF;//DnSignal[i] / tNIF;

            USFRate[i] = (pH1NIF->Integral() + pH1SB->Integral() - nNIF[i] - nSBsub[i]) / (tNIF + tSB);
            DUSFRate[i] = sqrt( pH1NIF->Integral() + pH1SB->Integral() + D2nNIF[i] + D2nSBsub[i] ) / (tNIF + tSB);

            TGraphErrors* g6 = new TGraphErrors(NumHist, x, USFRate, xerr, DUSFRate);
            g6->SetNameTitle("SF_Rate", "SF detection rate");
            sprintf(name, "%s/TimeDiff/Underground", FC.c_str());
            SaveToFile(name, g6);

            TGraphErrors* g7 = new TGraphErrors(NumHist, x, UNIFRate, xerr, DUNIFRate);
            g7->SetNameTitle("NIF_Rate", "In-scattering corrected NIF detection rate");
            sprintf(name, "%s/TimeDiff/Signal", FC.c_str());
            SaveToFile(name, g7);
        }
        cout << " ch " << i+1 << "   FG " << nNIF[i] << "+-" << sqrt(D2nNIF[i]) <<
                "   SB " << nSB[i] << "+-" << sqrt(D2nSB[i]) <<
                "   scaled SB " << scaledSB << "+-" << DscaledSB <<
                "   Signal " << nSignal[i] << "+-" << DnSignal[i] <<
                "   In-scat " << InScatPart[i] << "+-" << DInScatPart[i] << " %" << endl;
    }
    cout << " Avarage: " << sum / NumHist << " +- " << sqrt(D2sum / NumHist) << " %" << endl;

    TGraphErrors* g8 = new TGraphErrors(NumHist, x, nSignal, xerr, DnSignal);
    g8->SetNameTitle("cNIF", "In-scattering corrected NIF count");
    sprintf(name, "%s/TimeDiff/Signal", FC.c_str());
    SaveToFile(name, g8);

    TGraphErrors* g9 = new TGraphErrors(NumHist, x, nNIF, xerr, DnNIF);
    g9->SetNameTitle("nNIF", "number of detected time-correlated fission events. Method: subtraction");
    sprintf(name, "%s/TimeDiff/Signal/NIF", FC.c_str());
    SaveToFile(name, g9);

    TGraphErrors* g10 = new TGraphErrors(NumHist, x, nSBsub, xerr, DnSBsub);
    g10->SetNameTitle("nNIF", "number of detected time-correlated fission events. Method: subtraction");
    sprintf(name, "%s/TimeDiff/Signal/SB", FC.c_str());
    SaveToFile(name, g10);

    TGraphErrors* g11 = new TGraphErrors(NumHist, x, nNIFfit, xerr, DnNIFfit);
    g11->SetNameTitle("nNIFfit", "number of detected time-correlated fission events. Method: fit");
    sprintf(name, "%s/TimeDiff/Signal/NIF", FC.c_str());
    SaveToFile(name, g11);

    TGraphErrors* g12 = new TGraphErrors(NumHist, x, nSBfit, xerr, DnSBfit);
    g12->SetNameTitle("nNIFfit", "number of detected time-correlated fission events. Method: fit");
    sprintf(name, "%s/TimeDiff/Signal/SB", FC.c_str());
    SaveToFile(name, g12);

    TGraphErrors* g13 = new TGraphErrors(NumHist, x, nNIFtoMon, xerr, DnNIFtoMon);
    g13->SetNameTitle("NIFtoMon", "NIF rate normed to monitor");
    sprintf(name, "%s/TimeDiff/Signal/NIF", FC.c_str());
    SaveToFile(name, g13);

    TGraphErrors* g14 = new TGraphErrors(NumHist, x, nSBtoMon, xerr, DnSBtoMon);
    g14->SetNameTitle("NIFtoMon", "NIF rate normed to monitor");
    sprintf(name, "%s/TimeDiff/Signal/SB", FC.c_str());
    SaveToFile(name, g14);

    TGraphErrors* g15 = new TGraphErrors(NumHist, x, InScatPart, xerr, DInScatPart);
    g15->SetNameTitle("gScat", "In-scattered part in open geometry");
    sprintf(name, "%s/TimeDiff", FC.c_str());
    SaveToFile(name, g15);
}


Xsection::~Xsection()
{

}

/*void Xsection::PrintInScat()
{
    cout << endl << "=== In-scattering corrections ===" << endl;
    cout << "Neutron-induced fissions per monitor count" << endl;
    cout << "Channel nr, Free, Shadow bar, Ratio, scat-corr. NIF rate" << endl;
    for(int i = 0; i < 8; i++)
    {
        cout <<  i+1
             << ",  " << pHNIF->NIFRate[i] / MonitorNIF * pHNIF->t_real
             << " +- " << pHNIF->DNIFRate[i] / MonitorNIF * pHNIF->t_real;
        cout << ",  " << pHSB->NIFRate[i] / MonitorSB * pHSB->t_real
             << " +- " << pHSB->DNIFRate[i] / MonitorSB * pHSB->t_real;
        cout << ",  " << pHSB->NIFRate[i] / MonitorNIF * pHNIF->t_real / pHNIF->NIFRate[i] * MonitorSB / pHSB->t_real
             << " +- " << sqrt( pow(pHSB->DNIFRate[i]/pHNIF->NIFRate[i], 2) +
                                pow(pHNIF->DNIFRate[i]*pHSB->NIFRate[i]/pow(pHNIF->NIFRate[i], 2), 2) )
             << ",  " << NIFRate[i] << " +- " << DNIFRate[i] << endl;
    }
}

void Xsection::ScatCorrNIF()
{ // Calculate In-scattering corrected NIF rate normed to neutron flux
    for (int i = 0; i < NumHist; i++)
    {
        NIFRate[i] = pHNIF->NIFRate[i] - pHSB->NIFRate[i] * pHNIF->NeutronFlux[i] / pHSB->NeutronFlux[i];
        DNIFRate[i] = sqrt( pow(pHNIF->DNIFRate[i], 2) +
                            pow(pHSB->DNIFRate[i] * pHNIF->NeutronFlux[i] / pHSB->NeutronFlux[i], 2) +
                            pow(pHSB->NIFRate[i] * pHNIF->DNeutronFlux[i] / pHSB->NeutronFlux[i], 2) +
                            pow(pHSB->NIFRate[i] * pHNIF->NeutronFlux[i] * pHSB->DNeutronFlux[i] / pow(pHSB->NeutronFlux[i], 2), 2) );
    }
}//

void Xsection::CalculateNPu()
{
    cout << endl << "=== Number of 242Pu atoms ===" << endl;
//    Double_t w[3];
//    Double_t sum;
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {   // Calculate N(242Pu) from SF rate for each plate.
        // Calculate effective N(242Pu)
        nPuSF[i_ch] = SFRate[i_ch] * PuSFT2 / log(2.0); // SFRate, PuSFT2: both in seconds!
        DnPuSF[i_ch] = sqrt( pow(DSFRate[i_ch] * PuSFT2 / (log(2.0) * pHUG->eInt[i_ch]), 2) +
                             pow(SFRate[i_ch] * DPuSFT2 / (log(2.0) * pHUG->eInt[i_ch]), 2) );

        // Calculate N(242Pu) using internal efficiency from SF measurement
        nPu[i_ch] = nPuSF[i_ch] / pHUG->eInt[i_ch];
        DnPu[i_ch] = sqrt( pow(DnPuSF[i_ch] / pHUG->eInt[i_ch], 2) +
                           pow(nPuSF[i_ch] * pHUG->DeInt[i_ch] / pow(pHUG->eInt[i_ch], 2), 2) );

        cout << "    #Pu: " << nPu[i_ch] << " +- " << DnPu[i_ch] << ", eff=" << pHUG->eInt[i_ch] << endl;
    }

    // Draw...
    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

//    TGraphErrors* g1 = new TGraphErrors(NumHist, x, SFRate, xerr, DSFRate);
//    g1->SetNameTitle("SF_Rate", "PuFC combined SF detection rate");
//    SaveToFile("PuFC/N_atoms", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, nPuSF, xerr, DnPuSF);
    g2->SetNameTitle("NPuEff", "Effective number of 242Pu atoms from SF");
    SaveToFile("PuFC/N_atoms", g2);

    TGraphErrors* g3 = new TGraphErrors(NumHist, x, nPu, xerr, DnPu);
    g3->SetNameTitle("NPu", "Number of 242Pu atoms from spontaneaus fission");
    SaveToFile("PuFC/N_atoms", g3);
//
}

//void Xsection::CalculateEfficiency()
//{ // calculate the PuFC neutron-induced detection efficiency via comparing internal efficiencies of NIF and SB setup
//    //// note: method not successful, too large uncertainty.
//    Double_t NIFcSFRate, D2NIFcSFRate; // efficiency-corrected SF rate of NIF setup, squared error
//    Double_t NIFcNIFRate, D2NIFcNIFRate;
//    Double_t eNIF[NumHist], DeNIF[NumHist];
//    for (int i = 0; i < NumHist; i++)
//    {
//        eSF[i] = pHUG->eInt[i];
//        DeSF[i] = pHUG->DeInt[i];

//        NIFcSFRate = pHNIF->SFRate[i] / eSF[i];
//        D2NIFcSFRate = pow(pHNIF->DSFRate[i] / eSF[i], 2) +
//                       pow(pHNIF->SFRate[i] * DeSF[i], 2) / pow(eSF[i], 4);
//        NIFcNIFRate = (pHNIF->NIFRate[i] + pHNIF->SFRate[i]) / pHNIF->eInt[i] - NIFcSFRate;
//        D2NIFcNIFRate = D2NIFcSFRate +
//                        pow(pHNIF->DNIFRate[i] / pHNIF->eInt[i], 2) +
//                        pow(pHNIF->DNIFRate[i] / pHNIF->eInt[i], 2) +
//                        pow((pHNIF->NIFRate[i] + pHNIF->SFRate[i]) * pHNIF->DeInt[i], 2) / pow(pHNIF->eInt[i], 4);
//        eNIF[i] = pHNIF->NIFRate[i] / NIFcNIFRate;
//        DeNIF[i] = sqrt( D2NIFcNIFRate * pHNIF->NIFRate[i] / pow(NIFcNIFRate, 2) +
//                         pow(pHNIF->DNIFRate[i] / NIFcNIFRate, 2) );
//    }
//    // Draw...
//    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
//    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

//    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eNIF, xerr, DeNIF);
//    g1->SetNameTitle("eNIFdiff", "PuFC neutron-induced detection efficiency");
//    SaveToFile("PuFC/Efficiency", g1);
//}


Double_t Xsection::GetCorrectionFactor(string FC, Int_t i)
{// Correction factor for transmission and in-scattering
    if (strcmp(FC.c_str(), "UFC"))
    {
        Double_t Transm[] = {0.936084, 0.942565, 0.949147, 0.955535, 0.962355, 0.968932, 0.974894, 0.981514};
        Double_t InScat[] = {0.0999258, 0.0918293, 0.0852001, 0.0811008, 0.0807139, 0.0770468, 0.0685954, 0.0604568};
        return Transm[i] + InScat[i];
    } else {
        Double_t Transm[] = {0.937442, 0.942856, 0.948595, 0.954762, 0.960783, 0.9665, 0.972764, 0.979301};
        Double_t InScat[] = {0.0893054, 0.0819305, 0.0771868, 0.0732408, 0.0690476, 0.0643461, 0.0581159, 0.0518615};
        return Transm[i] + InScat[i];
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 1: #FF, neutron flux, #SF, T2_SF  -->  #PuEff, sigma_nfis                      ////
////////////////////////////////////////////////////////////////////////////////////////////////////
///   Assume SF and NIF FF detection efficiencies to be equal

void Xsection::Evaluation1()
{
    cout << endl << "=== Evaluation 1 ===" << endl;
//    CalculateNPu();
    CalculateCrossSection();
}


void Xsection::CalculateCrossSection()
{
    cout << "Cross section" << endl
         << " assuming equal efficiencies for NIF and SF" << endl
         << "ch  XS" << endl;
    Double_t CrossSection, DCrossSection;
    for(int i_ch = 0; i_ch < NumHist; i_ch ++)
    {
        Double_t FluxNIF = pHNIF->NeutronFluence[i_ch] / pHNIF->t_live * GetCorrectionFactor("UFC", i_ch);
        Double_t DFluxNIF = pHNIF->DNeutronFluence[i_ch] / pHNIF->t_live * GetCorrectionFactor("UFC", i_ch); // rel. unc. about 2.5 %
//        cout << "  Relative unc. " << DFluxNIF / FluxNIF << endl;

        // cross section = reaction rate / (EFFECTIVE number of target atoms * EFFECTIVE flux)
        Double_t Xs = NIFRate[i_ch] / (nPuSF[i_ch] * FluxNIF); // nPuSF == eff.*nPu
        Double_t DXs = sqrt( pow(DNIFRate[i_ch] / (nPuSF[i_ch] * FluxNIF), 2) +
                               pow(NIFRate[i_ch] * DnPuSF[i_ch] / (pow(nPuSF[i_ch], 2) * FluxNIF), 2) +
                               pow(NIFRate[i_ch] * DFluxNIF / (nPuSF[i_ch] * pow(FluxNIF, 2)), 2) );

        cout << " " << i_ch+1 << "  " << Xs * 1.E22 << " +- " << DXs * 1.E22 << " barn" << endl;
        XSec[i_ch] = Xs;
        DXSec[i_ch] = DXs;
    }
    // Calculate cross section average across channels
    Double_t sum = 0;
    Double_t Dsum = 0;
    for (int i = 0; i < NumHist; i++)
    {
        sum += XSec[i];
        Dsum += pow(DXSec[i], 2);
    }
    CrossSection = sum / NumHist;
    DCrossSection = sqrt(Dsum) / NumHist;
    cout << "Combined:  " << CrossSection * 1.E22 << " +- " << DCrossSection * 1.E22 << " barn" << endl;

    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    TGraphErrors* g3 = new TGraphErrors(NumHist, x, XSec, xerr, DXSec);
    g3->SetNameTitle("XS_SF", "Cross section (PuFC alone)");
    SaveToFile("PuFC/Evaluation", g3);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 2: #FF, neutron flux, sigma_nfis(Lit.), T2_SF  -->  #Pu                        ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation2()
{
    cout << endl << "=== Evaluation 2 ===" << endl;
    Double_t sum = 0;
    Double_t Dsum = 0;
    cout << "Effective number of Pu-242 atoms calculated with literature value" << endl <<
            "sigma_nfis = " << PuLit*1.E22 << " +- " << DPuLit*1.E22 << " barn" << endl <<
            "ch  N" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        nPuNIF[i] = NIFRate[i] / (pHNIF->NeutronFluence[i] / pHNIF->t_live * PuLit);
        DnPuNIF[i] = sqrt( pow(DNIFRate[i] / (pHNIF->NeutronFluence[i] / pHNIF->t_live * PuLit), 2) +
                        pow(NIFRate[i] * pHNIF->DNeutronFluence[i] / pHNIF->t_live / ( pow(pHNIF->NeutronFluence[i] / pHNIF->t_live, 2) * PuLit ), 2) +
                        pow(NIFRate[i] * DPuLit / ( pHNIF->NeutronFluence[i] / pHNIF->t_live * pow(PuLit, 2) ), 2) );
        cout << " " << i+1 << "  " << nPuNIF[i] << " +- " << DnPu[i] << endl;
        sum += nPuNIF[i];
        Dsum += pow(DnPuNIF[i], 2);
    }
    cout << "Sum = " << sum << " +- " << sqrt(Dsum) << endl;

    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, nPuNIF, xerr, DnPuNIF);
    g1->SetNameTitle("NPuLit", "Effective number of 242Pu atoms from NIF");
    SaveToFile("PuFC/N_atoms", g1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 3: #FF, neutron flux, sigma_nfis(Lit.), #SF, T2_SF  -->  efficiency            ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation3()
{
    cout << endl << "=== Evaluation 3 ===" << endl;
    cout << "Detection efficiencies' relation" << endl;

//    cout << "ch  NIF  NIF/SF" << endl;
//    for (int i = 0; i < NumHist; i++)
//    {
//        eRel[i] = nPuNIF[i] / nPuSF[i];
//        DeRel[i] = sqrt( pow(DnPuNIF[i] / nPuSF[i], 2) +
//                         pow(nPuNIF[i] * DnPuSF[i], 2) / pow(nPuSF[i], 4) );

//        eNIF[i] = eRel[i] * eSF[i];
//        DeNIF[i] = sqrt( pow(DeRel[i] * eSF[i], 2) +
//                         pow(eRel[i] * DeSF[i], 2) );

//        cout << " " << i+1 << "  " << eNIF[i] << "+-" << DeNIF[i] << "  " << eRel[i] << "+-" << DeRel[i] << endl;
//    }
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

//    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eNIF, xerr, DeNIF);
//    g1->SetNameTitle("eNIF", "Neutron-induced fission detection efficiency, determined from efficiency ratio");
//    SaveToFile("PuFC/Efficiency", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, eSF, xerr, DeSF);
    g2->SetNameTitle("eSF", "Spontaneaus fission detection efficiency");
    SaveToFile("PuFC/Efficiency", g2);

//    TGraphErrors* g3 = new TGraphErrors(NumHist, x, eRel, xerr, DeRel);
//    g3->SetNameTitle("eRel", "NIF to SF detection efficiency ratio, determined with given cross section");
//    SaveToFile("PuFC/Efficiency", g3);

//    Double_t eSimG[] = {eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther};
//    Double_t DeSimG[] = {DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther};
//    TGraphErrors* g4 = new TGraphErrors(NumHist, x, eSimG, xerr, DeSimG);
//    g4->SetNameTitle("eSimG", "Simulated efficiency (Gayther)");
//    SaveToFile("PuFC/Efficiency", g4);

    Double_t eSimM[] = {eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum};
    Double_t DeSimM[] = {DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum};
    TGraphErrors* g5 = new TGraphErrors(NumHist, x, eSimM, xerr, DeSimM);
    g5->SetNameTitle("eSimM", "Simulated efficiency (Minimum)");
    SaveToFile("PuFC/Efficiency", g5);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 4: #FF(Pu), #FF(U), neutron flux ratio, #SF, T2_SF, #U  -->  sigma_nfis        ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation4()
{
//    cout << endl << "=== Evaluation 4 ===" << endl << "UFC vs PuFC" << endl;
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g0 = new TGraphErrors(NumHist, x, UNIFRate, xerr, DUNIFRate);
    g0->SetNameTitle("cNIF", "UFC in-scattering corrected NIF detection rate");
    SaveToFile("UFC/TimeDiff/Signal", g0);

    //// UFC internal efficiency
    cout << "=== Internal efficiencies ===" << endl
         << "ch  eff" << endl;
    for ( int i = 0; i < NumHist; i++)
    { // i: UFC channel
        // use efficiencies from NIF only, SB too inaccurate
        eU[i] = pHUNIF->eInt[i];
        DeU[i] = pHUNIF->DeInt[i];
        cout << " " << i+1 << "  " << eU[i] << "+-" << DeU[i] << endl;
    }
    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eU, xerr, DeU);
    g1->SetNameTitle("eNIF", "UFC efficiency for NIF setup");
    SaveToFile("UFC/Efficiency", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, pHUSB->eInt, xerr, pHUSB->DeInt);
    g2->SetNameTitle("eSB", "UFC efficiency for SB setup");
    SaveToFile("UFC/Efficiency", g2);

    //// Number of U atoms
    cout << "=== Number of 235U atoms ===" << endl;
    Double_t N235 = 3.64247102724592E+020;  // rel. unc. 0.0074856945
    Double_t DN235 = 2.72664253120661E+018;
    Double_t frac235 = 0.904; // isotope fractions
    Double_t Dfrac235 = 0.005;
    Double_t frac238 = 0.0912;
    Double_t Dfrac238 = 0.0006;
    Double_t emA[] = {365.2, 396.2, 400.4, 393, 403.9, 403.7, 406.6, 397.1}; // deposits' average efficient Areal mass density. Unit: 10^-6 g / cm^2
    Double_t DemA[] = {1.7, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};
    Double_t nUeff[NumHist], DnUeff[NumHist];
    cout << "ch   n(U)   n(235U)eff" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        nUeff[i] = AreaUFC * emA[i] / (235 * u) * 1.E-8;
        DnUeff[i] = nUeff[i] * sqrt( pow(DAreaUFC / AreaUFC, 2) +
                                     pow(DemA[i] / emA[i], 2) );
        nU235eff[i] = frac235 * nUeff[i];
        DnU235eff[i] = sqrt( pow(Dfrac235 * nUeff[i], 2) +
                             pow(frac235 * DnUeff[i], 2) );
        nU238eff[i] = frac238 * nUeff[i];
        DnU238eff[i] = sqrt( pow(Dfrac238 * nUeff[i], 2) +
                             pow(frac238 * DnUeff[i], 2) );
        nU[i] = nUeff[i] / eU[i];
        DnU[i] = sqrt( pow(DnUeff[i] / eU[i], 2) +
                       pow(nUeff[i] * DeU[i], 2) / pow(eU[i], 2) );
        cout << " " << i+1 << "  " << nU[i] << "+-" << DnU[i] << "  " << nU235eff[i] << "+-" << DnU235eff[i] << endl;
    }
    if(CommentFlag)
        cout << "Check: Sum = " << nU[0] + nU[1] + nU[2] + nU[3] + nU[4] + nU[5] + nU[6] + nU[7] << endl;

//    TGraphErrors* g4 = new TGraphErrors(NumHist, x, nUSF, xerr, DnUSF);
//    g4->SetNameTitle("NUEff", "Effective number of 235U atoms from SF");
//    SaveToFile("UFC/N_atoms", g4);

    TGraphErrors* g5 = new TGraphErrors(NumHist, x, nU, xerr, DnU);
    g5->SetNameTitle("NU", "Number of Uranium atoms from H19 with internal efficiencies");
    SaveToFile("UFC/N_atoms", g5);

    //// Cross section
    cout << "=== U235 Cross section ===" << endl
         << "ch  nFlux  U238NIF  U235NIF  sigma  emA  deviation" << endl;
    Double_t CrossSection238 = 1.25; // neutron-induced fission cross section of 238U in barns
    Double_t DCrossSection238 = 0.0307 * CrossSection238;
    Double_t nFlux, DnFlux;
    Double_t U235NIFRate[NumHist];
    Double_t DU235NIFRate[NumHist];
    Double_t U238NIFRate[NumHist];
    Double_t DU238NIFRate[NumHist];
    Double_t Sum = 0, D2Sum = 0;
    Double_t CrossSection, DCrossSection;
    Double_t CSeff = frac235 * ULit + frac238 * CrossSection238; // Uranium average cross section
    Double_t DCSeff = sqrt( pow(Dfrac235 * ULit, 2) +
                            pow(frac235 * DULit, 2) +
                            pow(Dfrac238 * CrossSection238, 2) +
                            pow(frac238 * DCrossSection238,2) );
    Double_t emAcs[NumHist]; // efficient mass density determined from cross section
    Double_t DemAcs[NumHist];
    Double_t m = 235*u + frac238 * 3*u; // Average atom mass (des U-Isotopengemischs) in g
    Double_t Dm = 3*u * Dfrac238;
    for (int i = 0; i < NumHist; i++)
    {
        // determine cross section from mass density
        nFlux = pHUNIF->NeutronFluence[i] / pHNIF->t_live * GetCorrectionFactor("PuFC", i);
        DnFlux = pHUNIF->DNeutronFluence[i] / pHNIF->t_live * GetCorrectionFactor("PuFC", i);
        U238NIFRate[i] = 1.E-22 * nU238eff[i] * nFlux * CrossSection238; // flux: 1/(s*mm^2), CS: barn, factor 1.E-22 -> unit 1/s
        DU238NIFRate[i] = 1.E-22 * sqrt( pow(DnU238eff[i] * nFlux * CrossSection238, 2) +
                                         pow(nU238eff[i] * DnFlux * CrossSection238, 2) +
                                         pow(nU238eff[i] * nFlux * DCrossSection238, 2) );
        U235NIFRate[i] = UNIFRate[i] - U238NIFRate[i];
        DU235NIFRate[i] = sqrt( pow(DUNIFRate[i], 2) +
                                pow(DU238NIFRate[i], 2) );
        UXSec[i] = 1.E22 * U235NIFRate[i] / (nU235eff[i] * nFlux);
        DUXSec[i] = UXSec[i] * sqrt( pow(DU235NIFRate[i] / U235NIFRate[i], 2) +
                                     pow(DnU235eff[i] / nU235eff[i], 2) +
                                     pow(DnFlux / nFlux, 2) );
        Sum += UXSec[i];
        D2Sum += pow(DUXSec[i], 2);

        // Determine mass density from cross section
        emAcs[i] = m * UNIFRate[i] / (CSeff * nFlux * AreaUFC);
        DemAcs[i] = emAcs[i] * sqrt( pow(Dm / m, 2) +
                                     pow(DUNIFRate[i] / UNIFRate[i], 2) +
                                     pow(DCSeff / CSeff, 2) +
                                     pow(DnFlux / nFlux, 2) +
                                     pow(DAreaUFC / AreaUFC, 2) );


        cout << " " << i+1 << "  " << nFlux << "+-" << DnFlux << "/(s*mm^2)  "
             << U238NIFRate[i] << "+-" << DU238NIFRate[i] << "/s  "
             << U235NIFRate[i] << "+-" << DU235NIFRate[i] << "/s  "
             << UXSec[i] << "+-" << DUXSec[i] << " barn  "
             << emA[i] * ULit / UXSec[i] << "mug/cm^2  "
             << ULit / UXSec[i]
//             << emAcs[i] << "+-" << DemAcs[i] << "g/mm^2  "
//             << emAcs[i] / emA[i] << "+-" << sqrt( pow(DemAcs[i] / emA[i], 2) + pow(emAcs[i] * DemA[i] / emA[i] / emA[i], 2) )
             << endl;

    }
    CrossSection = Sum / NumHist;
    DCrossSection = sqrt(D2Sum) / NumHist;
    cout << "Combined: " << CrossSection << " barn +- " << 100 * DCrossSection / CrossSection << " %" << endl;

    TGraphErrors* g6 = new TGraphErrors(NumHist, x, UXSec, xerr, DUXSec);
    g6->SetNameTitle("XS_SF", "Cross section (UFC alone)");
    SaveToFile("UFC/Evaluation", g6);


}


void Xsection::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    TFile* file = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/Evaluation.root", "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder path already exists, otherwise create it
    if (file->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) file->Get(path.c_str());
    else
    {   file->mkdir(path.c_str(), "Folder containing offline Analysis objects");
        EvalDir = file->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    file->Save(); file->Close();
}


Double_t Xsection::func_peak(Double_t *x, Double_t *par)
{
    return par[0] * exp( - pow((x[0] - par[1]) / par[2], 2)) + par[3];
}
*/

/*void Xsection::CompareFiles(string path, Int_t start, Int_t stop)
{
    TCanvas* c1 = new TCanvas("SFRates", "SFRates", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    TMultiGraph* mg1 = new TMultiGraph();
    mg1->SetTitle("Spontaneaus fission detection rates; Deposit; Rate / Hz");
    TLegend* l1 = new TLegend(0.6, 0.2, 0.85, 0.40, "File nr");

    char fname[64] = "";
    Double_t x[] = {1,2,3,4,5,6,7,8};
    Double_t xerr[] = {0,0,0,0,0,0,0,0};
    Color_t color[] = {kRed, kGreen, kBlue, kYellow, kBlue, kYellow, kOrange, kRed, kGreen, kBlue, kGray, kYellow, kMagenta, kOrange};

    for ( int i = start; i < stop; i++ )
    {
        sprintf(fname, "%s/SB%i.root", path.c_str(), i);
        Hist* pH = new Hist(fname, "SB");
        pH->DoAnalyzeDt();
        cout << "nr " << i << ",  t_live/t_real=" << pH->t_live / pH->t_real << endl;
        for (int j = 0; j < 8; j++)
        {
            cout << "  ch " << j+1 << ",  SF rate=" << pH->SFRate[j] << endl;
        }
        TGraphErrors* ge = new TGraphErrors(NumHist, x, pH->SFRate, xerr, pH->DSFRate);
//        sprintf(fname, "SF_Rate_%i", i);
//        ge->SetNameTitle(fname, fname);
        ge->SetLineWidth(2);
        ge->SetLineColor(color[i-start]);
        mg1->Add(ge);
        sprintf(fname, "%i", i);
        l1->AddEntry(ge, fname, "lp");
    }
    mg1->Draw("AP");
    l1->Draw();
//    mg1->GetYaxis()->SetRangeUser(0, 2);
    c1->Modified();
    c1->Update();
}*/
