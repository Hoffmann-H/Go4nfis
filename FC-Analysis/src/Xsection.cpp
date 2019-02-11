#include "Xsection.h"

using namespace std;

Xsection::Xsection()
{
    CommentFlag = kFALSE;
    SetParam();

    pHNIF = new Hist("/home/hoffma93/Go4nfis/offline/results/NIF.root", "NIF");
    pHNIF->SetNeutronField(Yield, DYield, MonitorNIF, DMonitorNIF, 1500, 1);
    pHSB  = new Hist("/home/hoffma93/Go4nfis/offline/results/SB.root", "SB");
    pHSB->SetNeutronField(Yield, DYield, MonitorSB, DMonitorSB, 1500, 1);
    pHUG  = new Hist("/home/hoffma93/Go4nfis/offline/results/SF.root", "SF");
    pHNIF->DoAnalyzeQDC();
    pHSB ->DoAnalyzeQDC();
    pHUG ->DoAnalyzeQDC();
    CalculateThresholds();
    DoAnalyzeDt();
    CalculateNPu();
    CalculateEfficiency();
    Evaluation1();
    Evaluation2();
    Evaluation3();

    pHUNIF = new Hist("/home/hoffma93/Go4nfis/offline/results/UFC_NIF.root", "NIF", "UFC");
    pHUNIF->SetNeutronField(Yield, DYield, MonitorUNIF, DMonitorUNIF, 1500, 1);
    pHUSB = new Hist("/home/hoffma93/Go4nfis/offline/results/UFC_SB.root", "SB", "UFC");
    pHUSB->SetNeutronField(Yield, DYield, MonitorUSB, DMonitorUSB, 1500, 1);

    pHUNIF->DoAnalyzeQDC();
    pHUSB->DoAnalyzeQDC();
    Evaluation4();

//    CompareFiles("/home/hoffma93/Go4nfis/offline/results/test", 191, 204);
}

void Xsection::SetParam()
{
    // physics parameters
    u = 1.660539E-24; // atomic mass unit in g
    PuLit = 2120E-25; // 242Pu neutron-induced fission cross section in mm^2, original in mb
    DPuLit = 35E-25;
    ULit = 2.1E-22; // 235U nif cs in mm^2, original in b
    DULit = 0.0315E-22;
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
    AreaPuFC = 4300; // in mm^2
    DAreaPuFC = 120;
    AreaUFC = 4300; // in mm^2
    DAreaUFC = 120;
    mPu = 0.03724; // in g
    DmPu = 0.00001;
    mU = 0.15801;
    DmU = 0.00001;
    eSimGayther = 0.988;
    DeSimGayther = 0.005;
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
        Double_t counts[] = {pHNIF->GetNumberEvents(i), pHSB->GetNumberEvents(i), pHUG->GetNumberEvents(i)};
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
    SaveToFile("Analysis/PuFC/QDC", g1);
}


void Xsection::DoAnalyzeDt()
{
    cout << "Analyzing TimeDiff spectra. " << endl;
    Double_t Dt_min = 63600;
    Double_t Dt_max = 78500;
    Double_t ToF_low, ToF_up, ChPerBin, BinOffset;
    Double_t TimePerCh = 25.E-12; // 25 ps in seconds
    Double_t tNIF = pHNIF->t_live;
    Double_t tSB = pHSB->t_live;
    Double_t tUG = pHUG->t_live;
    cout << "Times: " << tNIF << ", " << tSB << ", " << tUG << endl;
    Int_t lim[4];
    Double_t nNIF[NumHist];
    Double_t D2nNIF[NumHist];
    Double_t nSBsub[NumHist];
    Double_t D2nSBsub[NumHist];
    Double_t nSBfit[NumHist];
    Double_t D2nSBfit[NumHist];
    Double_t nSB[NumHist];
    Double_t D2nSB[NumHist];
    Double_t nSignal[NumHist];
    Double_t DnSignal[NumHist];
    TH1I *pH1NIF, *pH1SB, *pH1UG;
    char name[64] = "";

//    cout << "ch   average Ug   Ug-rate   nNIF   nSB   corrected NIF   in-scat ratio" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        // get hand-made peak limits
        ToF_low = pHNIF->GetPeakLow(i);
        ToF_up = pHNIF->GetPeakUp(i);

        // set limits from 3 sigma to n sigma
        Double_t n = 1;
        Double_t center = 0.5 * (ToF_low + ToF_up);
        Double_t width = (ToF_up - ToF_low) / 6.0;
        ToF_low = center - n * width;
        ToF_up = center + n * width;

        pH1NIF = pHNIF->pHDtG[i];
        pH1SB = pHSB->pHDtG[i];
        pH1UG = pHUG->pHDtG[i];
        ChPerBin = pH1NIF->GetBinWidth(0);
        BinOffset = pH1NIF->GetBinCenter(0);
        // Calculate integration limit bins
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

        //// Fit constant underground
        Double_t UgCount = pH1NIF->Integral(lim[0], lim[3]) - NIFPeakCount +
                           pH1SB->Integral(lim[0], lim[3]) - SBPeakCount +
                           pH1UG->Integral(lim[0], lim[3]);
        Double_t DUgCount = sqrt(UgCount);
        // avUg: average underground count per bin and t_live
        Double_t avUg = UgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) + tUG * (lim[3] - lim[0] + 1) );
        Double_t DavUg = DUgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) + tUG * (lim[3] - lim[0] + 1) );
        if(CommentFlag)
            cout << " Underground per bin: " << avUg << "+-" << DavUg << endl;

        // number of time-correlated FF events in NIF measurement
        nNIF[i] = NIFPeakCount - (lim[2] - lim[1] + 1) * tNIF * avUg;
        D2nNIF[i] = NIFPeakCount + pow((lim[2] - lim[1] + 1) * tNIF * DavUg, 2);

        //// number of time-correlated FF events in SB measurement
        //// 1st method: subtraction of time-independent underground
        // number of time-correlated FF events in SB measurement
        nSBsub[i] = SBPeakCount - (lim[2] - lim[1] + 1) * tSB * avUg;
        D2nSBsub[i] = SBPeakCount + pow((lim[2] - lim[1] + 1) * tSB * DavUg, 2);
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
        sprintf(name, "fNIFDtPeak_%i", i+1);
        TF1* fNIF = new TF1(name, func_peak, Dt_min, Dt_max, 4);
        fNIF->SetRange(Dt_min, Dt_max);
        fNIF->SetParameters(100, center, width, tNIF * avUg);
        pH1NIF->Fit(name, "RQ");
        center = (Double_t)fNIF->GetParameter(1);
        width = (Double_t)fNIF->GetParameter(2);
        Dwidth = (Double_t)fNIF->GetParError(2);

        sprintf(name, "fSBDtPeak_%i", i+1);
        TF1* fSB = new TF1(name, func_peak, Dt_min, Dt_max, 4);
        fSB->SetRange(Dt_min, Dt_max);
        fSB->SetParameters(100, center, width, tSB * avUg);
        fSB->FixParameter(1, center);
        fSB->FixParameter(2, 3*width);
        pH1SB->Fit(name, "RBQ");
        ampl = (Double_t)fSB->GetParameter(0);
        Dampl = (Double_t)fSB->GetParError(0);
        Double_t pi = 3.1415926535898;
        nSBfit[i] = sqrt(pi) / ChPerBin * ampl * width;
        D2nSBfit[i] = pi / pow(ChPerBin, 2) * (pow(Dampl * width, 2) + pow(ampl * Dwidth, 2));
        if(CommentFlag)
            cout << " SB peak content (2nd method): " << nSBfit[i] << "+-" << sqrt(D2nSBfit[i]) << endl;
        nSB[i] = nSBsub[i]; // choose method
        D2nSB[i] = D2nSBsub[i];

        //// In-Scattering correction
        nSignal[i] = nNIF[i] - nSB[i] * (tNIF / tSB) * (pHNIF->NeutronFlux[i] / pHSB->NeutronFlux[i]);
        DnSignal[i] = sqrt( D2nNIF[i] +
                            D2nSB[i] * pow((tNIF / tSB) * (pHNIF->NeutronFlux[i] / pHSB->NeutronFlux[i]), 2) +
                            pow(nSB[i] * (tNIF / tSB) * (pHNIF->DNeutronFlux[i] / pHSB->NeutronFlux[i]), 2) +
                            pow(nSB[i] * (tNIF / tSB) * (pHNIF->NeutronFlux[i] * pHSB->DNeutronFlux[i] / pow(pHSB->NeutronFlux[i], 2)), 2) );
        NIFRate[i] = nSignal[i] / tNIF;
        DNIFRate[i] = DnSignal[i] / tNIF;

        SFRate[i] = (pH1NIF->Integral() + pH1SB->Integral() + pH1UG->Integral() - nNIF[i] - nSB[i]) / (tNIF + tSB + tUG);
        DSFRate[i] = sqrt( pH1NIF->Integral() + pH1SB->Integral() + pH1UG->Integral() + D2nNIF[i] + D2nSB[i] ) / (tNIF + tSB + tUG);

        cout << "PuFC ch " << i+1 << ", NIF " << nSignal[i] << "+-" << DnSignal[i] <<
                ", rate " << NIFRate[i] << "+-" << DNIFRate[i] << " / s, Underground " << SFRate[i] << "+-" << DSFRate[i] << endl;
    }
}

void Xsection::UFC_AnalyzeDt()
{
    cout << "Analyzing TimeDiff spectra. " << endl;
    Double_t Dt_min = 63600;
    Double_t Dt_max = 78500;
    Double_t ToF_low, ToF_up, ChPerBin, BinOffset;
    Double_t TimePerCh = 25.E-12; // 25 ps in seconds
    Double_t tNIF = pHUNIF->t_live;
    Double_t tSB = pHUSB->t_live;
    cout << "Times: " << tNIF << ", " << tSB << endl;
    Int_t lim[4];
    Double_t nNIF[NumHist];
    Double_t D2nNIF[NumHist];
    Double_t nSB[NumHist];
    Double_t D2nSB[NumHist];
    Double_t nSignal[NumHist];
    Double_t DnSignal[NumHist];
    TH1I *pH1NIF, *pH1SB;
    char name[64] = "";

    for (int i = 0; i < NumHist; i++)
    {
        // get hand-made peak limits
        ToF_low = pHUNIF->GetPeakLow(i);
        ToF_up = pHUNIF->GetPeakUp(i);

        // set limits from 3 sigma to n sigma
        Double_t n = 3;
        Double_t center = 0.5 * (ToF_low + ToF_up);
        Double_t width = (ToF_up - ToF_low) / 6.0;
        ToF_low = center - n * width;
        ToF_up = center + n * width;

        pH1NIF = pHUNIF->pHDtG[i];
        pH1SB = pHUSB->pHDtG[i];
        ChPerBin = pH1NIF->GetBinWidth(0);
        BinOffset = pH1NIF->GetBinCenter(0);
        // Calculate integration limit bins
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

        //// Fit constant underground
        Double_t UgCount = pH1NIF->Integral(lim[0], lim[3]) - NIFPeakCount +
                           pH1SB->Integral(lim[0], lim[3]) - SBPeakCount;
        Double_t DUgCount = sqrt(UgCount);
        // avUg: average underground count per bin and t_live
        Double_t avUg = UgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) );
        Double_t DavUg = DUgCount / ( (tNIF + tSB) * (lim[1] - lim[0] + lim[3] - lim[2]) );
        if(CommentFlag)
            cout << " Underground per bin: " << avUg << "+-" << DavUg << endl;

        // number of time-correlated FF events in NIF measurement
        nNIF[i] = NIFPeakCount - (lim[2] - lim[1] + 1) * tNIF * avUg;
        D2nNIF[i] = NIFPeakCount + pow((lim[2] - lim[1] + 1) * tNIF * DavUg, 2);

        // number of time-correlated FF events in SB measurement
        nSB[i] = SBPeakCount - (lim[2] - lim[1] + 1) * tSB * avUg;
        D2nSB[i] = SBPeakCount + pow((lim[2] - lim[1] + 1) * tSB * DavUg, 2);
        if(CommentFlag)
        {
            cout << " NIF peak integral: " << NIFPeakCount << "+-" << sqrt(NIFPeakCount) << endl <<
                    " NIF peak underground: " << (lim[2] - lim[1] + 1) * tNIF * avUg << "+-" << (lim[2] - lim[1] + 1) * tNIF * DavUg << endl <<
                    " NIF peak content: " << nNIF << "+-" << sqrt(D2nNIF[i]) << endl;
            cout << " SB peak integral: " << SBPeakCount << "+-" << sqrt(SBPeakCount) << endl <<
                    " SB peak underground: " << (lim[2] - lim[1] + 1) * tSB * avUg << "+-" << (lim[2] - lim[1] + 1) * tSB * DavUg << endl <<
                    " SB peak content (1st method): " << nSB << "+-" << sqrt(D2nSB[i]) << endl;
        }

        //// In-Scattering correction
        nSignal[i] = nNIF[i] - nSB[i] * (tNIF / tSB) * (pHUNIF->NeutronFlux[i] / pHUSB->NeutronFlux[i]);
        DnSignal[i] = sqrt( D2nNIF[i] +
                            D2nSB[i] * pow((tNIF / tSB) * (pHUNIF->NeutronFlux[i] / pHUSB->NeutronFlux[i]), 2) +
                            pow(nSB[i] * (tNIF / tSB) * (pHUNIF->DNeutronFlux[i] / pHUSB->NeutronFlux[i]), 2) +
                            pow(nSB[i] * (tNIF / tSB) * (pHUNIF->NeutronFlux[i] * pHUSB->DNeutronFlux[i] / pow(pHUSB->NeutronFlux[i], 2)), 2) );
        UNIFRate[i] = nSignal[i] / tNIF;
        DUNIFRate[i] = DnSignal[i] / tNIF;

        USFRate[i] = TimePerCh * ChPerBin * avUg;
        DUSFRate[i] = TimePerCh * ChPerBin * DavUg;

        cout << "UFC ch " << i+1 << ", NIF " << nSignal[i] << "+-" << DnSignal[i] <<
                ", rate " << UNIFRate[i] << "+-" << DUNIFRate[i] << " / s, Underground " << USFRate[i] << "+-" << DUSFRate[i] << " / s" << endl;
    }
    return;
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
}//*/

void Xsection::CalculateNPu()
{
    cout << endl << "=== Number of 242Pu atoms ===" << endl;
//    Double_t w[3];
//    Double_t sum;
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {   // Calculate N(242Pu) from SF rate for each plate.
//        cout << "Spontaneaus fission rates in s^-1 for channel " << i_ch + 1 << endl;
//        // NIF measurement:
//        cout << "    NIF: " << pHNIF->SFRate[i_ch] << " +- " << pHNIF->DSFRate[i_ch] << endl;
//        cout << "    SB:  " << pHSB->SFRate[i_ch] << " +- " << pHSB->DSFRate[i_ch] << endl;
//        cout << "    SF:  " << pHUG->SFRate[i_ch] << " +- " << pHUG->DSFRate[i_ch] << endl;
//        // Combine 3 measurements
//        sum = 1/pHNIF->DSFRate[i_ch] + 1/pHSB->DSFRate[i_ch] + 1/pHUG->DSFRate[i_ch];
//        w[0] = 1/pHNIF->DSFRate[i_ch]/sum; // weighting rates according to inverse uncertainty
//        w[1] = 1/pHSB->DSFRate[i_ch]/sum;
//        w[2] = 1/pHUG->DSFRate[i_ch]/sum;
//        if(CommentFlag)
//            cout << "    Weighting  " << w[0] << "  " << w[1] << "  " << w[2] << endl;
//        SFRate[i_ch] = w[0] * pHNIF->SFRate[i_ch] +
//                       w[1] * pHSB->SFRate[i_ch] +
//                       w[2] * pHUG->SFRate[i_ch];
//        DSFRate[i_ch] = sqrt( pow(w[0] * pHNIF->DSFRate[i_ch], 2) +
//                              pow(w[1] * pHSB->DSFRate[i_ch], 2) +
//                              pow(w[2] * pHUG->DSFRate[i_ch], 2) );

        // Calculate effective N(242Pu)
        nPuSF[i_ch] = SFRate[i_ch] * PuSFT2 / log(2.0); // SFRate, PuSFT2: both in seconds!
        DnPuSF[i_ch] = sqrt( pow(DSFRate[i_ch] * PuSFT2 / (log(2.0)/* * pHUG->eInt[i_ch]*/), 2) +
                             pow(SFRate[i_ch] * DPuSFT2 / (log(2.0)/* * pHUG->eInt[i_ch]*/), 2) );

//        cout << "Test: pHNIF->DSFRate[3-7]" << endl;
//        for (int i = 0; i < NumHist; i++)
//            cout << i+1 << "  " << pHNIF->SFRate[i] << "+-" << pHNIF->DSFRate[i] << endl;

        // Calculate N(242Pu) using internal efficiency from SF measurement
        nPu[i_ch] = nPuSF[i_ch] / pHUG->eInt[i_ch];

//        cout << "Test: pHNIF->DSFRate[3-7]" << endl;
//        for (int i = 0; i < NumHist; i++)
//            cout << i+1 << "  " << pHNIF->SFRate[i] << "+-" << pHNIF->DSFRate[i] << endl;

        DnPu[i_ch] = sqrt( pow(DnPuSF[i_ch] / pHUG->eInt[i_ch], 2) +
                           pow(nPuSF[i_ch] * pHUG->DeInt[i_ch] / pow(pHUG->eInt[i_ch], 2), 2) );

        cout << "    #Pu: " << nPu[i_ch] << " +- " << DnPu[i_ch] << ", eff=" << pHUG->eInt[i_ch] << endl;
    }

    // Draw...
    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, SFRate, xerr, DSFRate);
    g1->SetNameTitle("SF_Rate", "PuFC combined SF detection rate");
    SaveToFile("Analysis/PuFC/SF", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, nPuSF, xerr, DnPuSF);
    g2->SetNameTitle("NPuEff", "Effective number of 242Pu atoms from SF");
    SaveToFile("Analysis/PuFC/SF", g2);

    TGraphErrors* g3 = new TGraphErrors(NumHist, x, nPu, xerr, DnPu);
    g3->SetNameTitle("NPu", "Number of 242Pu atoms from spontaneaus fission");
    SaveToFile("Analysis/PuFC/SF", g3);
//*/
}

void Xsection::CalculateEfficiency()
{ // calculate the PuFC neutron-induced detection efficiency via comparing internal efficiencies of NIF and SB setup
    Double_t NIFcSFRate, D2NIFcSFRate; // efficiency-corrected SF rate of NIF setup, squared error
    Double_t NIFcNIFRate, D2NIFcNIFRate;
    Double_t eNIF[NumHist], DeNIF[NumHist];
    for (int i = 0; i < NumHist; i++)
    {
        eSF[i] = pHUG->eInt[i];
        DeSF[i] = pHUG->DeInt[i];

        NIFcSFRate = pHNIF->SFRate[i] / eSF[i];
        D2NIFcSFRate = pow(pHNIF->DSFRate[i] / eSF[i], 2) +
                       pow(pHNIF->SFRate[i] * DeSF[i], 2) / pow(eSF[i], 4);
        NIFcNIFRate = (pHNIF->NIFRate[i] + pHNIF->SFRate[i]) / pHNIF->eInt[i] - NIFcSFRate;
        D2NIFcNIFRate = D2NIFcSFRate +
                        pow(pHNIF->DNIFRate[i] / pHNIF->eInt[i], 2) +
                        pow(pHNIF->DNIFRate[i] / pHNIF->eInt[i], 2) +
                        pow((pHNIF->NIFRate[i] + pHNIF->SFRate[i]) * pHNIF->DeInt[i], 2) / pow(pHNIF->eInt[i], 4);
        eNIF[i] = pHNIF->NIFRate[i] / NIFcNIFRate;
        DeNIF[i] = sqrt( D2NIFcNIFRate * pHNIF->NIFRate[i] / pow(NIFcNIFRate, 2) +
                         pow(pHNIF->DNIFRate[i] / NIFcNIFRate, 2) );
    }
    // Draw...
    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eNIF, xerr, DeNIF);
    g1->SetNameTitle("eNIFdiff", "PuFC neutron-induced detection efficiency");
    SaveToFile("Analysis/PuFC/Efficiency", g1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 1: #FF, neutron flux, #SF, T2_SF  -->  #PuEff, sigma_nfis                      ////
////////////////////////////////////////////////////////////////////////////////////////////////////
///   Assume SF and NIF efficiencies to be equal

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
        Double_t FluxNIF = pHNIF->NeutronFlux[i_ch];
        Double_t DFluxNIF = pHNIF->DNeutronFlux[i_ch];

        // cross section = reaction rate / (effective number of target atoms * flux)
        Double_t Xs = NIFRate[i_ch] / (nPu[i_ch] * FluxNIF);
        Double_t DXs = sqrt( pow(DNIFRate[i_ch] / (nPu[i_ch] * FluxNIF), 2) +
                               pow(NIFRate[i_ch] * DnPu[i_ch] / (pow(nPu[i_ch], 2) * FluxNIF), 2) +
                               pow(NIFRate[i_ch] * DFluxNIF / (nPu[i_ch] * pow(FluxNIF, 2)), 2) );

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
    cout << "Combined:  " << CrossSection * 1.E22 << " barn +- " << 100 * DCrossSection / CrossSection << " %" << endl;

    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    TGraphErrors* g3 = new TGraphErrors(NumHist, x, XSec, xerr, DXSec);
    g3->SetNameTitle("XS_SF", "Cross section (PuFC alone)");
    SaveToFile("Analysis/Evaluation", g3);
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
        nPuNIF[i] = NIFRate[i] / (pHNIF->NeutronFlux[i] * PuLit);
        DnPuNIF[i] = sqrt( pow(DNIFRate[i] / (pHNIF->NeutronFlux[i] * PuLit), 2) +
                        pow(NIFRate[i] * pHNIF->DNeutronFlux[i] / ( pow(pHNIF->NeutronFlux[i], 2) * PuLit ), 2) +
                        pow(NIFRate[i] * DPuLit / ( pHNIF->NeutronFlux[i] * pow(PuLit, 2) ), 2) );
        cout << " " << i+1 << "  " << nPuNIF[i] << /*" +- " << DnPu[i] << */endl;
        sum += nPuNIF[i];
        Dsum += pow(DnPuNIF[i], 2);
    }
    cout << "Sum = " << sum << " +- " << sqrt(Dsum) << endl;

    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, nPuNIF, xerr, DnPuNIF);
    g1->SetNameTitle("NPuLit", "Effective number of 242Pu atoms from NIF");
    SaveToFile("Analysis/Evaluation", g1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 3: #FF, neutron flux, sigma_nfis(Lit.), #SF, T2_SF  -->  efficiency            ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation3()
{
    cout << endl << "=== Evaluation 3 ===" << endl;
    cout << "Detection efficiencies' relation" << endl;

//    Double_t sumSF = 0;
//    Double_t D2sumSF = 0;
//    for (int i = 0; i < NumHist; i++)
//    {
//        sumSF += SFRate[i];
//        D2sumSF += pow(DSFRate[i], 2);
//    }

//    eSF = sumSF * PuSFT2 / (N * log(2.0));
//    DeSF = sqrt( D2sumSF*pow(PuSFT2 / (N * log(2.0)), 2) +
//                         pow(sumSF * DPuSFT2 / (N * log(2.0)), 2) +
//                         pow(sumSF * PuSFT2 * DN / (N * N * log(2.0)), 2) );
//    cout << "average SF eff.: " << eSF << " +- " << DeSF << endl;

    cout << "ch  NIF  NIF/SF" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        eRel[i] = nPuNIF[i] / nPuSF[i];
        DeRel[i] = sqrt( pow(DnPuNIF[i] / nPuSF[i], 2) +
                         pow(nPuNIF[i] * DnPuSF[i], 2) / pow(nPuSF[i], 4) );

        eNIF[i] = eRel[i] * eSF[i];
        DeNIF[i] = sqrt( pow(DeRel[i] * eSF[i], 2) +
                         pow(eRel[i] * DeSF[i], 2) );

        cout << " " << i+1 << "  " << eNIF[i] << "+-" << DeNIF[i] << "  " << eRel[i] << "+-" << DeRel[i] << endl;
    }
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eNIF, xerr, DeNIF);
    g1->SetNameTitle("eNIF", "Neutron-induced fission detection efficiency");
    SaveToFile("Analysis/PuFC/Efficiency", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, eSF, xerr, DeSF);
    g2->SetNameTitle("eSF", "Spontaneaus fission detection efficiency");
    SaveToFile("Analysis/PuFC/Efficiency", g2);

    TGraphErrors* g3 = new TGraphErrors(NumHist, x, eRel, xerr, DeRel);
    g3->SetNameTitle("eRel", "NIF to SF detection efficiency ratio");
    SaveToFile("Analysis/PuFC/Efficiency", g3);

    Double_t eSimG[] = {eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther, eSimGayther};
    Double_t DeSimG[] = {DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther, DeSimGayther};
    TGraphErrors* g4 = new TGraphErrors(NumHist, x, eSimG, xerr, DeSimG);
    g4->SetNameTitle("eSimG", "Simulated efficiency (Gayther)");
    SaveToFile("Analysis/PuFC/Efficiency", g4);

    Double_t eSimM[] = {eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum, eSimMinimum};
    Double_t DeSimM[] = {DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum, DeSimMinimum};
    TGraphErrors* g5 = new TGraphErrors(NumHist, x, eSimM, xerr, DeSimM);
    g5->SetNameTitle("eSimM", "Simulated efficiency (Minimum)");
    SaveToFile("Analysis/PuFC/Efficiency", g5);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////  Evaluation 4: #FF(Pu), #FF(U), neutron flux ratio, #SF, T2_SF, #U  -->  sigma_nfis        ////
////////////////////////////////////////////////////////////////////////////////////////////////////
void Xsection::Evaluation4()
{
    cout << endl << "=== Evaluation 4 ===" << endl << "UFC vs PuFC" << endl;
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    //// UFC In-scattering correction
    UFC_AnalyzeDt();

    TGraphErrors* g0 = new TGraphErrors(NumHist, x, UNIFRate, xerr, DUNIFRate);
    g0->SetNameTitle("cNIF", "UFC in-scattering corrected NIF detection rate");
    SaveToFile("Analysis/UFC", g0);

    //// UFC internal efficiency
    cout << "Internal efficiencies" << endl
         << "ch  eff" << endl;
    for ( int i = 0; i < NumHist; i++)
    { // i: UFC channel
        eU[i] = pHUNIF->eInt[i];
        DeU[i] = pHUNIF->DeInt[i];
        cout << " " << i+1 << "  " << eU[i] << "+-" << DeU[i] << endl;
    }
    TGraphErrors* g1 = new TGraphErrors(NumHist, x, eU, xerr, DeU);
    g1->SetNameTitle("eNIF", "UFC efficiency for NIF setup");
    SaveToFile("Analysis/UFC/Efficiency", g1);

    TGraphErrors* g2 = new TGraphErrors(NumHist, x, pHUSB->eInt, xerr, pHUSB->DeInt);
    g2->SetNameTitle("eSB", "UFC efficiency for SB setup");
    SaveToFile("Analysis/UFC/Efficiency", g2);

    //// Number of U atoms
    cout << "Number of 235U atoms" << endl;
    Double_t N235 = 3.64247102724592E+020;  // rel. unc. 0.0074856945
    Double_t DN235 = 2.72664253120661E+018;
    cout << "Total number fixed to N=" << N235 << " +- " << DN235 << endl;
    Double_t emA[] = {365.2, 396.2, 400.4, 393, 403.9, 403.7, 406.6, 397.1}; // deposits' average efficient Areal mass density. Unit: 10^-6 g / cm^2
    Double_t DemA[] = {1.7, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};
    cout << "ch  <N>  N_eff" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        nU[i] = N235 / NumHist;
        DnU[i] = DN235 / NumHist;
        nUNIF[i] = AreaUFC * emA[i] / (235 * u) * 1.E-8;
        DnUNIF[i] = nUNIF[i] * sqrt( pow(DAreaUFC / AreaUFC, 2) +
                                     pow(DemA[i] / emA[i], 2) );
        cout << " " << i+1 << "  " << nU[i] << "+-" << DnU[i] << "  " << nUNIF[i] << "+-" << DnUNIF[i] << endl;
    }
    if(CommentFlag)
        cout << "Check: Sum = " << nUNIF[0] + nUNIF[1] + nUNIF[2] + nUNIF[3] + nUNIF[4] + nUNIF[5] + nUNIF[6] + nUNIF[7] << endl;

//    TGraphErrors* g4 = new TGraphErrors(NumHist, x, nUSF, xerr, DnUSF);
//    g4->SetNameTitle("NUEff", "Effective number of 235U atoms from SF");
//    SaveToFile("Analysis/UFC/SF", g4);

    TGraphErrors* g5 = new TGraphErrors(NumHist, x, nU, xerr, DnU);
    g5->SetNameTitle("NU", "Number of 235U atoms from spontaneaus fission");
    SaveToFile("Analysis/UFC/SF", g5);

    //// Cross section
    cout << "Cross section" << endl
         << "ch  sigma" << endl;
    Double_t nFlux, DnFlux;
    Double_t Sum = 0, D2Sum = 0;
    Double_t CrossSection, DCrossSection;
    for (int i = 0; i < NumHist; i++)
    {
        nFlux = pHUNIF->NeutronFlux[i];
        DnFlux = pHUNIF->DNeutronFlux[i];
        UXSec[i] = UNIFRate[i] / (nUNIF[i] * nFlux);
        DUXSec[i] = UXSec[i] * sqrt( pow(DUNIFRate[i] / UNIFRate[i], 2) +
                                     pow(DnUNIF[i] / nUNIF[i], 2) +
                                     pow(DnFlux / nFlux, 2) );
        Sum += UXSec[i];
        D2Sum += pow(DUXSec[i], 2);
        cout << " " << i+1 << "  " << UXSec[i] * 1.E22 << "+-" << DUXSec[i] * 1.E22 << " barn" << endl;
    }
    CrossSection = Sum / NumHist;
    DCrossSection = sqrt(D2Sum) / NumHist;
    cout << "Combined: " << CrossSection * 1.E22 << " barn +- " << 100 * DCrossSection / CrossSection << " %" << endl;

    TGraphErrors* g6 = new TGraphErrors(NumHist, x, UXSec, xerr, DUXSec);
    g6->SetNameTitle("XS_SF", "Cross section (UFC alone)");
    SaveToFile("Analysis/UFC", g6);
}


void Xsection::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    TFile* file = TFile::Open("/home/hoffma93/Go4nfis/offline/results/Evaluation.root", "UPDATE");
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
