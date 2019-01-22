#include "Xsection.h"

using namespace std;

Xsection::Xsection()
{
    CommentFlag = kFALSE;

    // physics parameters
    PuSFT2 = 6.77E10 * 365.24*24*60*60; // Pu-242 spontaneaus fission half-life period in s
    DPuSFT2 = 7E8 * 365.24*24*60*60;
    Yield = 2.15882E4;
    DYield = 543.3;
    MonitorNIF = 27492079+33478370+30916955+54792679;
    DMonitorNIF = 108342;
    MonitorSB = 3623069+28646614+4757385;
    DMonitorSB = 41697;
    Area = 4300; // in mm^2
    DArea = 120;


    pHNIF = new Hist("/home/hoffma93/Go4nfis/offline/results/NIF.root");
    pHSB  = new Hist("/home/hoffma93/Go4nfis/offline/results/SB.root");
    pHSF  = new Hist("/home/hoffma93/Go4nfis/offline/results/SF.root");
    pHNIF->DoAnalyzeDt(kTRUE);
//    pHNIF->DoAnalyzeQDC();
    pHSB ->DoAnalyzeDt(kTRUE);
//    pHSB ->DoAnalyzeQDC();
    pHSF ->DoAnalyzeDt(kTRUE);
//    pHSF ->DoAnalyzeQDC();

    PrintInScat();
    CalculateNPu();
    CalculateCrossSection();
}

Xsection::~Xsection()
{

}

void Xsection::CalculateNPu()
{
//    Double_t SFRate, DSFRate;
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {   // Calculate SF rate in 1/s and N(242Pu) for each plate.
        cout << "Spontaneaus fission rates in s^-1 for channel " << i_ch + 1 << endl;
        // NIF measurement:
        cout << "    NIF: " << pHNIF->SFRate[i_ch] << " +- " << pHNIF->DSFRate[i_ch] << endl;
        cout << "    SB:  " << pHSB->SFRate[i_ch] << " +- " << pHSB->DSFRate[i_ch] << endl;
        cout << "    SF:  " << pHSF->SFRate[i_ch] << " +- " << pHSF->DSFRate[i_ch] << endl;
        // Combine 3 measurements
        Double_t sum = 1/pHNIF->DSFRate[i_ch] + 1/pHSB->DSFRate[i_ch] + 1/pHSF->DSFRate[i_ch];
        Double_t w[] = {1/pHNIF->DSFRate[i_ch]/sum, // weighting rates acc. to inverse uncertainty
                        1/pHSB->DSFRate[i_ch]/sum,
                        1/pHSF->DSFRate[i_ch]/sum};
        SFRate[i_ch] = w[0] * pHNIF->SFRate[i_ch] +
                 w[1] * pHSB->SFRate[i_ch] +
                 w[2] * pHSF->SFRate[i_ch];
        DSFRate[i_ch] = sqrt( pow(w[0] * pHNIF->DSFRate[i_ch], 2) +
                              pow(w[1] * pHSB->DSFRate[i_ch], 2) +
                              pow(w[2] * pHSF->DSFRate[i_ch], 2) );
        cout << "    All: " << SFRate[i_ch] << " +- " << DSFRate[i_ch] << endl;

        // Calculate N(242Pu)
        nPu[i_ch] = SFRate[i_ch] * PuSFT2 / log(2.0); // SFRate, PuSFT2: both in seconds!
        DnPu[i_ch] = sqrt( pow(DSFRate[i_ch] * PuSFT2 / log(2.0), 2) +
                           pow(SFRate[i_ch] * DPuSFT2 / log(2.0), 2) );
        cout << "    #Pu: " << nPu[i_ch] << " +- " << DnPu[i_ch] << " (efficiency=100%)" << endl;
    }
/*
    // Alter Code zum Zeichnen...
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
//*/
}

void Xsection::CalculateCrossSection()
{
    for(int i_ch = 0;i_ch < NumHist; i_ch ++)
    {
        cout << "Cross section channel " << i_ch+1 << endl;

        Double_t distance = 1500 + 80 - 10 * i_ch; // distance from target to deposit i_ch in mm
        Double_t Ddistance = 1; // mm
        // neutron flux at deposit i_ch, SB measurement:
        Double_t FluxSB = Yield * MonitorSB / (distance*distance * pHSB->t_real);
//        cout << Yield << " " << MonitorSB << " " << pHSB->t_real << endl;
        // neutron flux uncertainty assuming exact time
        Double_t DFluxSB = sqrt( pow(DYield * MonitorSB / (pow(distance, 2) * pHSB->t_real), 2) +
                                 pow(Yield * DMonitorSB / (pow(distance, 2) * pHSB->t_real), 2) +
                                 pow(Yield * MonitorSB * 2*Ddistance / (pow(distance, 3) * pHSB->t_real), 2) );
        // cross section = reaction rate / (number of target atoms * flux)
        Double_t Scat = pHSB->NIFRate[i_ch] / (nPu[i_ch] * FluxSB);
//        cout << pHSB->NIFRate[i_ch] << " " << nPu[i_ch] << " " << FluxSB << endl;
        Double_t DScat = sqrt( pow(pHSB->DNIFRate[i_ch] / (nPu[i_ch] * FluxSB), 2) +
                               pow(pHSB->NIFRate[i_ch] * DnPu[i_ch] / (pow(nPu[i_ch], 2) * FluxSB), 2) +
                               pow(pHSB->NIFRate[i_ch] * DFluxSB / (nPu[i_ch] * pow(FluxSB, 2)), 2) );
        cout << "In-scattering cross section: " << Scat << " +- " << DScat << " mm^2" << endl;

        // neutron flux at deposit i_ch, NIF measurement:
        Double_t FluxNIF = Yield * MonitorNIF / (distance*distance * pHNIF->t_real);
        // neutron flux uncertainty assuming exact time
        Double_t DFluxNIF = sqrt( pow(DYield * MonitorNIF / (pow(distance, 2) * pHNIF->t_real), 2) +
                                 pow(Yield * DMonitorNIF / (pow(distance, 2) * pHNIF->t_real), 2) +
                                 pow(Yield * MonitorNIF * 2*Ddistance / (pow(distance, 3) * pHNIF->t_real), 2) );
        // cross section = reaction rate / (number of target atoms * flux)
        Double_t XsRaw = pHNIF->NIFRate[i_ch] / (nPu[i_ch] * FluxNIF);
        Double_t DXsRaw = sqrt( pow(pHNIF->DNIFRate[i_ch] / (nPu[i_ch] * FluxNIF), 2) +
                               pow(pHNIF->NIFRate[i_ch] * DnPu[i_ch] / (pow(nPu[i_ch], 2) * FluxNIF), 2) +
                               pow(pHNIF->NIFRate[i_ch] * DFluxNIF / (nPu[i_ch] * pow(FluxNIF, 2)), 2) );
        cout << "uncorrected cross section: " << XsRaw << " +- " << DXsRaw << " mm^2" << endl;

        // In-scattering underground subtraction
        Double_t Xs = XsRaw - Scat;
        Double_t DXs = sqrt( pow(DXsRaw, 2) + pow(DScat, 2) );
        cout << "cross section: " << Xs << " +- " << DXs << " mm^2" << endl;
        cout << "             = " << Xs * 1.E22 << " +- " << DXs * 1.E22 << " barn" << endl;
        XSec[i_ch] = Xs;
        DXSec[i_ch] = DXs;
    }
    Double_t sum = 0;
    Double_t Dsum = 0;
    for (int i = 0; i < NumHist; i++)
    {
        sum += XSec[i];
        Dsum += pow(DXSec[i], 2);
    }
    CrossSection = sum / NumHist;
    DCrossSection = sqrt(Dsum/NumHist);
    cout << "Combined: " << CrossSection * 1.E22 << " +- " << DCrossSection * 1.E22 << " barn" << endl;
}

void Xsection::PrintInScat()
{
    cout << "=== In-scattering corrections ===" << endl;
    cout << "Neutron-induced fissions per monitor count" << endl;
    cout << "Channel nr, Free, Shadow bar, Ratio " << endl;
    for(int i = 0; i < 8; i++)
    {
        cout << /*"  Channel " <<*/ i+1 /*<< endl*/;
        cout << ",  " << pHNIF->NIFRate[i] / MonitorNIF * pHNIF->t_real/*
             << " +- " << pHNIF->DNIFRate[i] / MonitorNIF * pHNIF->t_real << endl*/;
        cout << ",  " << pHSB->NIFRate[i] / MonitorSB * pHSB->t_real/*
             << " +- " << pHSB->DNIFRate[i] / MonitorSB * pHSB->t_real << endl*/;
//        cout << ",  " << pHNIF->NIFRate[i] / MonitorNIF * pHNIF->t_real - pHSB->NIFRate[i] / MonitorSB * pHSB->t_real
//             << " +- " << sqrt( pow(pHNIF->DNIFRate[i] / MonitorNIF * pHNIF->t_real, 2) +
//                                pow(pHSB->DNIFRate[i] / MonitorSB * pHSB->t_real, 2) ) /*<< endl*/;
        cout << ",  " << pHSB->NIFRate[i] / MonitorNIF * pHNIF->t_real / pHNIF->NIFRate[i] * MonitorSB / pHSB->t_real
             << " +- " << sqrt( pow(pHSB->DNIFRate[i]/pHNIF->NIFRate[i], 2) +
                                pow(pHNIF->DNIFRate[i]*pHSB->NIFRate[i]/pow(pHNIF->NIFRate[i], 2), 2) ) << endl;
    }
}

/*void Xsection::DoAnalyze()
{
    /// 1. Analyze TimeDiff spectra, Dead-time correction
//    pHNIF->DoAnalyzeDt(kTRUE);
    /// 3. Analyze pulse height spectra
//    pHNIF->DoAnalyzeQDC();
    /// 4. Calculate 242Pu amount
//    CalculateNPu();
    /// 5. Calculate Cross Section
//    CalculateCrossSection();

}//*/
/*




*/
