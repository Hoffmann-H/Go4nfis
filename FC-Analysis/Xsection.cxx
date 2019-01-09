#include "Xsection.h"
#include "Hist.h"
using namespace std;

Xsection::Xsection()
{
    CommentFlag = kFALSE;

    // physics
    PuSFT2 = 6.77E10; // Pu-242 spontaneaus fission half-life period
    DPuSFT2 = 7E8;
    NeutronFlux = 1.0;
    //define surrounding peak limits
    ToF_low = 66500;
    ToF_up  = 68500;
    //define time window for underground integration
    Dt_min = 63300;
    Dt_max = 78500;
    // Load histograms
    char hist_name[64] = "";

    pHNIF = new Hist("/home/hoffma93/Go4nfis/offline/results/NIF.root");
    pHSB  = new Hist("/home/hoffma93/Go4nfis/offline/results/SB.root");

    pNIF_t_live = pHNIF->GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pNIF_t_real = pHNIF->GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_48");
    pSB_t_live = pHSB->GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pSB_t_real = pHSB->GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_48");

    for(int i_ch = 0; i_ch < NumHist; i_ch++) {
        sprintf(hist_name, "/Histograms/Analysis/PuFC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i_ch + 1);
        pNIF_DtG[i_ch] = pHNIF->GetTH1I(hist_name);
        pSB_DtG[i_ch] = pHSB->GetTH1I(hist_name);
        sprintf(hist_name, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i_ch + 1);
        pHRawQDCl[i_ch] = pHNIF->GetTH1I(hist_name);
        sprintf(hist_name, "/Histograms/Analysis/PuFC/QDC/low/trig/H1AnaQDCl_trig_%i", i_ch + 1);
        pHAnaQDCl[i_ch] = pHNIF->GetTH1I(hist_name);
    }
}

void Xsection::DoAnalyze()
{
    /// 1. Divide QDC-gated TimeDiff into peak and underground
    char sname[32] = "";
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {
        sprintf(sname, "DtG_%i", i_ch+1);
        AnalyzePeak(pNIF_DtG[i_ch], &N1NIF_raw[i_ch], &DN1NIF_raw[i_ch], &N1SF_raw[i_ch], &DN1SF_raw[i_ch]);
        pHNIF->Save(pNIF_DtG[i_ch], "Analysis/NIF", sname);
        AnalyzePeak(pSB_DtG[i_ch], &N2NIF_raw[i_ch], &DN2NIF_raw[i_ch], &N2SF_raw[i_ch], &DN2SF_raw[i_ch]);
        pHSB->Save(pSB_DtG[i_ch], "Analysis/NIF", sname);
    }
    /// 2. Get times
    GetTimes();
    /// 3. Dead-time correction
    DeadTimeCorrection();
    /// 4. Calculate 242Pu amount
    CalculateNPu();
    /// 5. Calculate Cross Section
    CalculateCrossSection();

    /// Test
//    cout << pHRawQDCl[0]->GetBinWidth(0) << endl;
//    Int_t i_min;
//    i_min = GetMinBin(pHRawQDCl[0], 700, 800);
//    cout << i_min << endl;
//    i_min = GetMinBin(pHRawQDCl[0], 700, 1300);
//    cout << i_min << endl;
}

void Xsection::AnalyzePeak(TH1I *pH, Double_t *pNIF, Double_t *pDNIF, Double_t *pSF, Double_t *pDSF)
{
    Double_t ChPerBin = pH->GetBinWidth(0);
    Double_t BinOffset = pH->GetBinCenter(0);
    // Calculate integration limit bins
    Double_t lim_0 = (Dt_min - BinOffset) / ChPerBin;
    Double_t lim_1 = (ToF_low - BinOffset) / ChPerBin;
    Double_t lim_2 = (ToF_up - BinOffset) / ChPerBin;
    Double_t lim_3 = (Dt_max - BinOffset) / ChPerBin;
    if(CommentFlag)
        cout << "Analyzing peak. Ch per bin: " << ChPerBin << ", Offset: " << BinOffset << endl <<
                "Integration limits: Bin " << lim_0 << " " << lim_1 << " " << lim_2 << " " << lim_3 << endl;
    // Integrations
    Int_t PeakCount = pH->Integral(lim_1, lim_2); // Sum within peak region
    Double_t UgCount = pH->Integral(lim_0, lim_3) - PeakCount; // Sum within underground region; follows Poisson-stat
    Double_t DUgCount = sqrt(UgCount); // Underground deviation
    Double_t UgPerBin = 1.0 * UgCount / (lim_1 - lim_0 + lim_3 - lim_2);
    Double_t DUgPerBin = 1.0 * DUgCount / (lim_1 - lim_0 + lim_3 - lim_2);
    *pNIF = PeakCount - (lim_2 - lim_1) * UgPerBin;
    *pDNIF = (lim_2 - lim_1) * DUgPerBin;
    *pSF = pH->Integral() - *pNIF; // Calculate total underground
    *pDSF = sqrt((*pDNIF)*(*pDNIF) + *pSF); // Error according to Gauss
    if(CommentFlag)
        cout << "NIF: " << *pNIF << " +- " << *pDNIF << endl <<
                "SF:  " << *pSF  << " +- " << *pDSF  << endl;
    // Draw
    pH->GetXaxis()->SetRangeUser(62000, 80000);
    pH->Draw();
    TLine *lH = new TLine(Dt_min, UgPerBin, Dt_max, UgPerBin);
    lH->SetLineColor(kRed);
    lH -> Draw();
    Double_t max = pH->GetMaximum();
    TLine *lV0 = new TLine(ToF_low, 0, ToF_low, max);
    lV0->Draw();
    TLine *lV1 = new TLine(ToF_up , 0, ToF_up , max);
    lV1->Draw();
    char message[32] = "";
    sprintf(message, "NIF counts: %i", (int)*pNIF);
    TText *t = new TText();
    t->SetNDC();
    t->DrawText(0.4, 0.8, message);
    pH->GetListOfFunctions()->Add(lH);
    pH->GetListOfFunctions()->Add(lV0);
    pH->GetListOfFunctions()->Add(lV1);
}

void Xsection::AnalyzeUnderground(TH1I *pH, Double_t *pSF, Double_t *pDSF)
{
    // Analyze underground analogue to AnalyzePeak(), but without peak area exclusion.
}

void Xsection::GetTimes()
{
    t_live_NIF = pNIF_t_live->Integral();
    t_real_NIF = pNIF_t_real->Integral();
    t_live_SB = pSB_t_live->Integral();
    t_real_SB = pSB_t_real->Integral();
//    t_live_SF = pSF_t_live->Integral();
//    t_real_SF = pSF_t_real->Integral();
}

void Xsection::DeadTimeCorrection()
{
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {
        N1SF_dtc[i_ch] = N1SF_raw[i_ch] * t_real_NIF / t_live_NIF;
        DN1SF_dtc[i_ch] = DN1SF_raw[i_ch] * t_real_NIF / t_live_NIF;
        N1NIF_dtc[i_ch] = N1NIF_raw[i_ch] * t_real_NIF / t_live_NIF;
        DN1NIF_dtc[i_ch] = DN1NIF_raw[i_ch] * t_real_NIF / t_live_NIF;
        N2SF_dtc[i_ch] = N2SF_raw[i_ch] * t_real_SB / t_live_SB;
        DN2SF_dtc[i_ch] = DN2SF_raw[i_ch] * t_real_SB / t_live_SB;
        N2NIF_dtc[i_ch] = N2NIF_raw[i_ch] * t_real_SB / t_live_SB;
        DN2NIF_dtc[i_ch] = DN2NIF_raw[i_ch] * t_real_SB / t_live_SB;
//        N3SF_dtc[i_ch] = N3SF_raw[i_ch] * t_real_SF / t_live_SF;
//        DN3SF_dtc[i_ch] = DN3SF_raw[i_ch] * t_real_SF / t_live_SF;
    }
}

void Xsection::CalculateNPu()
{
    Double_t SFRate, DSFRate;
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {   // Calculate N(242Pu) for each plate.
        // NIF measurement:
        SFRate = N1SF_dtc[i_ch] / t_real_NIF;
        NPuNIF[i_ch] = SFRate * PuSFT2 / log(2.0);
        // Uncertainties
        DSFRate = DN1SF_dtc[i_ch] / t_real_NIF; // neglect time error
        DNPuNIF[i_ch] = sqrt( pow(DPuSFT2 * SFRate, 2) + pow(DSFRate * PuSFT2, 2) );
        // SB measurement:
        SFRate = N2SF_dtc[i_ch] / t_real_SB;
        NPuSB[i_ch] = SFRate * PuSFT2 / log(2.0);
        // Uncertainties
        DSFRate = DN2SF_dtc[i_ch] / t_real_SB; // neglect time error
        DNPuSB[i_ch] = sqrt( pow(DPuSFT2 * SFRate, 2) + pow(DSFRate * PuSFT2, 2) );
        NPu[i_ch] = NPuNIF[i_ch];   // Here one could make an average over all undergrounds
        DNPu[i_ch] = DNPuNIF[i_ch];
//        cout << NPu[i_ch] << "+-" << DNPu[i_ch] << endl;
//        cout << NPu[i_ch] << "+-" << DNPu[i_ch] << endl;
    }

    double x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};

    pNIFNPu = new TCanvas("NPu","Number of Pu-242 atoms vs Channel", 0, 0, 1, 1);
    TGraphErrors* g1 = new TGraphErrors(NumHist, x, N1NIF_dtc, xerr, DN1NIF_dtc);
    g1->GetYaxis()->SetTitle("N");
    g1->GetXaxis()->SetTitle("Channel");
    g1->Draw("ap");
    pHNIF->SaveCanvas(pNIFNPu, "Analysis/SF", "NPu");

    pSBNPu = new TCanvas("NPu","Number of Pu-242 atoms vs Channel", 600, 400);
    TGraphErrors* g2 = new TGraphErrors(NumHist, x, N2NIF_dtc, xerr, DN2NIF_dtc);
    g2->GetYaxis()->SetTitle("N");
    g2->GetXaxis()->SetTitle("Channel");
    g2->Draw("ap");
    pHSB->SaveCanvas(pSBNPu, "Analysis/SF", "NPu");

}

void Xsection::CalculateCrossSection()
{
    cout << "Cross section" << endl;
    Double_t CrossSection[NumHist];
    Double_t DCrossSection[NumHist];
    for(int i_ch = 0;i_ch < NumHist; i_ch ++)
    {
        Double_t NIFRate = N1NIF_dtc[i_ch] / t_real_NIF - N2NIF_dtc[i_ch] / t_real_SB; // Subtract scattered NIF underground
        Double_t DNIFRate = 0;
        CrossSection[i_ch] = NIFRate / ( NPu[i_ch] * Neutron_flux );
        DCrossSection[i_ch] = 0;
        cout << NIFRate << "  " << NPu[i_ch] << "  " << CrossSection[i_ch] << endl;
    }
}

Int_t Xsection::GetMinBin(TH1I *pH, Int_t low, Int_t up)
{
    Int_t Min_bin = pH->GetBinContent(pH->GetBin(low));
    for (int i_bin = pH->GetBin(low); i_bin < pH->GetBin(up); i_bin++) {
        if(pH->GetBinContent(i_bin) < pH->GetBinContent(Min_bin))
            Min_bin = i_bin;
    }
    return Min_bin;
}
