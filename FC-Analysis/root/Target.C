#ifndef TARGET_H
#define TARGET_H
#include "SaveToFile.C"
#include "/home/hoffma93/Programme/ROOT/nToF.C"
#include "TRandom.h"

TH1D* GetToF(TH1D *pE, Double_t FlightPath, string sample, Double_t N = 100000)
{ // convert a TARGET energy distribution into ToF.
  // N: Monte-Carlo parameter
    char name[64] = "";
    sprintf(name, "ToF_%s", pE->GetName());
    TH1D *pToF = new TH1D(name, name, 2000, 0, 200);
    sprintf(name, "%s. ToF", sample.c_str());
    pToF->SetTitle(name);
    pToF->GetXaxis()->SetTitle("#font[12]{t} [ns]");
    pToF->GetYaxis()->SetName(pE->GetYaxis()->GetName());
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

TH1D* ScaleWithSigma(TH1D *hN, TGraphErrors *gSigma, string sample = "Dir")
{
    char name[64] = "";
    Int_t nBins = hN->GetNbinsX();
    Double_t xmin = hN->GetBinLowEdge(1);
    Double_t xmax = hN->GetBinLowEdge(nBins+1);
    sprintf(name, "Fission_%s", sample.c_str());
    TH1D *hFis = new TH1D(name, name, nBins, xmin, xmax);
    sprintf(name, "%s. Fission", sample.c_str());
    hFis->SetTitle(name);
    hFis->GetXaxis()->SetTitle(hN->GetXaxis()->GetTitle());
    sprintf(name, "#font[12]{N}_{(n,f)%s}", sample.c_str());
    hFis->GetYaxis()->SetTitle(name);

    for (Int_t j = 1; j < nBins; j++)
        hFis->SetBinContent(j, gSigma->Eval(hN->GetBinCenter(j)) * hN->GetBinContent(j));

    hFis->Scale(hFis->Integral());
    return hFis;
}

void TargetToF(string File, string Subfolder = "")
{ // copy TARGET ToF histograms into Analysis.root/Subfolder
    TGraph *geDir = new TGraph(File.c_str(), "%*lg %lg %lg");
    TGraph *geSc  = new TGraph(File.c_str(), "%*lg %lg %*lg %lg");
    TGraph *geTot = new TGraph(File.c_str(), "%*lg %lg %*lg %*lg %lg");

    TH1D *hDir = new TH1D("ToF_Dir", "Direct TARGET ToF; #font[12]{t} [ns]; #frac{d#font[12]{n}}{d#font[12]{t}} [1/100ps]", 2000, 0, 200);
    TH1D *hSc = new TH1D("ToF_Sc", "Scattered TARGET ToF; #font[12]{t} [ns]; #frac{d#font[12]{n}}{d#font[12]{t}} [1/100ps]", 2000, 0, 200);
    TH1D *hTot = new TH1D("ToF_Tot", "Total TARGET ToF; #font[12]{t} [ns]; #frac{d#font[12]{n}}{d#font[12]{t}} [1/100ps]", 2000, 0, 200);
    Int_t N = geDir->GetN();
    for (Int_t k = 0; k < N; k++)
    {
        Double_t x, y;
        geDir->GetPoint(k, x, y);
        Int_t bin = hDir->FindBin(x-0.0001); // upper bin edges
        hDir->SetBinContent(bin, y);
        geSc->GetPoint(k, x, y);
        hSc->SetBinContent(bin, y);
        geTot->GetPoint(k, x, y);
        hTot->SetBinContent(bin, y);
    }
//    hToF->Scale(1/hToF->Integral());
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Save(fAna, "Simulation/Target/"+Subfolder, hDir);
    Save(fAna, "Simulation/Target/"+Subfolder, hSc);
    Save(fAna, "Simulation/Target/"+Subfolder, hTot);
    fAna->Save();
    fAna->Close();//*/
}

void TargetE(string File, string Subfolder = "", Double_t Np = 1)
{ // Copy TARGET energy histograms from File into Analysis.root/Subfolder
    TGraph *geDir = new TGraph(File.c_str(), "%*lg %lg %lg");
    TGraph *geSc = new TGraph(File.c_str(), "%*lg %lg %*lg %lg");
    TGraph *geTot = new TGraph(File.c_str(), "%*lg %lg %*lg %*lg %lg");

    TH1D *hDir = new TH1D("E_Dir", "Direct TARGET spectrum; #font[12]{E}_{n} [MeV]; #frac{d#font[12]{n}}{d#font[12]{E}} [1/10keV]", 1600, 0, 16);
    TH1D *hSc = new TH1D("E_Sc", "Scattered TARGET spectrum; #font[12]{E}_{n} [MeV]; #frac{d#font[12]{n}}{d#font[12]{E}} [1/10keV]", 1600, 0, 16);
    TH1D *hTot = new TH1D("E_Tot", "Total TARGET spectrum; #font[12]{E}_{n} [MeV]; #frac{d#font[12]{n}}{d#font[12]{E}} [1/10keV]", 1600, 0, 16);
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
    TH1D *hT_Tot = GetToF(hTot, 1.65, "Tot", 1000000);
    TH1D *hT_Sc = GetToF(hSc, 1.65, "Sc", 1000000);

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
    Save(fAna, "Simulation/Target/"+Subfolder, hT_Tot);
    fAna->Save();
    fAna->Close();//*/
}

void TargetAng(string File, string Subfolder = "", Double_t Np = 1)
{ // copy TARGET angular distribution from File into Analysis.root/Subfolder
    TGraph *gAng = new TGraph(File.c_str(), "%lg %lg");
    gAng->SetNameTitle("Angular", "Angular neutron distribution; Angel / deg; #sigma / mb/sr");
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Save(fAna, "Simulation/Target/"+Subfolder, gAng);
    fAna->Save();
    fAna->Close();
}

void Target()
{
    TargetE("/home/hoffma93/Programme/TARGET/Results/Variation/400keV_1E8_ENE", "", 100000000);
    TargetToF("/home/hoffma93/Programme/TARGET/Results/Variation/400keV_1E8_TOF", "");

//    TargetE("/home/hoffma93/Programme/TARGET/Results/Variation/Double_Distance/FC_STARGET_15MEENE", "3.3m", 100000000);

//    TargetAng("/home/hoffma93/Programme/TARGET/Results/Variation/FC_STARGET_15MEANG", "", 100000000);
}

#endif
