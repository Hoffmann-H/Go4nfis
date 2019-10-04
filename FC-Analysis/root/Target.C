#ifndef TARGET_H
#define TARGET_H
#include "SaveToFile.C"
#include "/home/hoffma93/Programme/ROOT/nToF.C"
#include "TRandom.h"

TH1D* GetToF(TH1D *pE, Double_t FlightPath, Double_t N = 100000)
{
    TH1D *pToF = new TH1D("ToF_E", "ToF(E); #font[12]{t} [ns]; #frac{d#font[12]{n}}{d#font[12]{t}} [1/100ps]", 2000, 0, 200);
    Double_t E, t;
    Int_t BinE, BinT;
    TRandom *r = new TRandom();
    for (Int_t k = 0; k < N; k++)
    {
        E = r->Uniform(16.0);
        BinE = pE->FindBin(E);
        t = nToF(FlightPath, E) - 0.4; // Hard-coded adjusting. Direct neutons are not delayed.
        BinT = pToF->FindBin(t);
        pToF->AddBinContent(BinT, pE->GetBinContent(BinE));
    }
    pToF->Scale(pE->GetNbinsX() / N);
    return pToF;
}

void TargetToF(string File, string Subfolder = "")
{
    TGraph *geToF = new TGraph(File.c_str(), "%*lg %lg %*lg %*lg %lg");
    TH1D *hToF = new TH1D("ToF_Target", "Total TARGET ToF; #font[12]{t} [ns]; #frac{d#font[12]{n}}{d#font[12]{t}} [1/100ps]", 2000, 0, 200);
    Int_t N = geToF->GetN();
    for (Int_t k = 0; k < N; k++)
    {
        Double_t x, y;
        geToF->GetPoint(k, x, y);
        Int_t bin = hToF->FindBin(x-0.0001);
        hToF->SetBinContent(bin, y);
    }
//    hToF->Scale(1/hToF->Integral());
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Save(fAna, "Simulation/Target/"+Subfolder, hToF);
    fAna->Save();
    fAna->Close();//*/
}

void TargetE(string File, string Subfolder = "", Double_t Np = 1)
{
    TGraph *geDir = new TGraph(File.c_str(), "%*lg %lg %lg");
    TGraph *geSc = new TGraph(File.c_str(), "%*lg %lg %*lg %lg");
    TGraph *geTot = new TGraph(File.c_str(), "%*lg %lg %*lg %*lg %lg");

    TH1D *hDir = new TH1D("Direct", "Direct TARGET spectrum; #font[12]{E}_{n} [MeV]; #frac{d#font[12]{n}}{d#font[12]{E}} [1/10keV]", 1600, 0, 16);
    TH1D *hSc = new TH1D("Scattered", "Scattered TARGET spectrum; #font[12]{E}_{n} [MeV]; #frac{d#font[12]{n}}{d#font[12]{E}} [1/10keV]", 1600, 0, 16);
    TH1D *hTot = new TH1D("Total", "Total TARGET spectrum; #font[12]{E}_{n} [MeV]; #frac{d#font[12]{n}}{d#font[12]{E}} [1/10keV]", 1600, 0, 16);
    Int_t N = geDir->GetN();
    for (Int_t k = 0; k < N; k++)
    {
        Double_t x, y;
        geDir->GetPoint(k, x, y);
        Int_t bin = hDir->FindBin(x-0.0001);
        hDir->SetBinContent(bin, y);
        geSc->GetPoint(k, x, y);
        hSc->SetBinContent(bin, y);
        geTot->GetPoint(k, x, y);
        hTot->SetBinContent(bin, y);
    }
    TH1D *hT = GetToF(hTot, 1.65, 1000000);

    cout << "Simulated ions: " << Np << endl;
    cout << "Integral E-Spektrum: " << hTot->Integral() << endl;
//    Double_t Scale = 1.0 / Np;
//    hDir->Scale(Scale);
//    hSc->Scale(Scale);
//    hTot->Scale(Scale);
//    cout << "Direct: " << 100 * hDir->Integral() << " % " << endl;
//    cout << "Scattered in target: " << 100 * hSc->Integral() << " % " << endl;

//    pT->Draw();

    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Save(fAna, "Simulation/Target/"+Subfolder, hDir);
    Save(fAna, "Simulation/Target/"+Subfolder, hSc);
    Save(fAna, "Simulation/Target/"+Subfolder, hTot);
    Save(fAna, "Simulation/Target/"+Subfolder, hT);
    fAna->Save();
    fAna->Close();//*/
}

void TargetAng(string File, string Subfolder = "")
{
    
}

void Target()
{
    TargetE("/home/hoffma93/Programme/TARGET/Results/STARGET_15MEV_Variation/400keV_1E8_ENE", "", 100000000);
    TargetToF("/home/hoffma93/Programme/TARGET/Results/STARGET_15MEV_Variation/400keV_1E8_TOF", "");
    
//    TargetE("/home/hoffma93/Programme/TARGET/Results/FC_STARGET_15MEENE", "3.3m", 500000);
//    TargetE("/home/hoffma93/Programme/TARGET/Results/STARGET_15MEV_Variation/400keV_1E6", "400keV/", 1000000);

//    Target("/home/hoffma93/Programme/TARGET/Results/STARGET_15MEV_Variation/Konstant_400keV_1E6", "Konstant/");
//    Target("/home/hoffma93/Programme/TARGET/Results/STARGET_15MEV_Variation/1-0_400keV_1E6", "1-0/");
}

#endif
