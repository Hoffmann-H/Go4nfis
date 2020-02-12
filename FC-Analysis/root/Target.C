#ifndef TARGET_H
#define TARGET_H
#include "SaveToFile.C"
#include "FC.C"
#include "DrawPics.C"
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

Double_t TargetFactor(TFile *fAna)
{
    char name[32] = "";
    sprintf(name, "Simulation/Target/E_Dir");
    TH1D *hTargetDir = (TH1D*)fAna->Get(name);
    if (!hTargetDir) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/Target/E_Tot");
    TH1D *hTargetTot = (TH1D*)fAna->Get(name);
    if (!hTargetTot) cout << "Could not get " << name << endl;
    return hTargetDir->Integral() / hTargetTot->Integral();
}

void TargetGateU(Bool_t Draw = 0)
{
    char name[64] = "";
    TFile *fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    if (!fAna) cout << "Could not open " << "fAna" << endl;
    Double_t fTarget = TargetFactor(fAna);
    TFile *fFG = TFile::Open("~/Programme/Geant4-Work/results/UFC_ideal_c_FG.root");
    if (!fFG) cout << "Could not open " << "fFG" << endl;
    TFile *fTot = TFile::Open("~/Programme/Geant4-Work/results/UFC_ideal_c_FG+BG.root");
    if (!fTot) cout << "Could not open " << "fTot" << endl;
    TGraphErrors *geU = new TGraphErrors(8);
    geU->SetName("UFC_Target_Gate");
    geU->SetTitle("UFC Target scattering and Gating correction; Deposit; #font[12]{k}_{T}");
    TGraph *gLim1 = new TGraph(8);
    gLim1->SetNameTitle("UFC_Geant4_lim1", "Left ToF gate bin; Deposit; bin nr.");
    TGraph *gLim2 = new TGraph(8);
    gLim2->SetNameTitle("UFC_Geant4_lim2", "Right ToF gate bin; Deposit; bin nr.");
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "UFC/ToF/UFC_ToF_Ch.%i", i+1);
        TH1F *hFG = (TH1F*) fFG->Get(name); if (!hFG) cout << "Could not get FG " << name << endl;
        hFG->Scale(fTarget);
        TH1F *hTot = (TH1F*) fTot->Get(name); if (!hTot) cout << "Could not get Tot " << name << endl;
        Int_t tMax = (Int_t) hTot->GetBinCenter(hTot->GetMaximumBin());
        Int_t bl = hTot->FindBin(tMax - Left("UFC"));
        Int_t br = hTot->FindBin(tMax + Right("UFC") + 1.0) - 1;
        gLim1->SetPoint(i, i+1, bl / 10);
        gLim2->SetPoint(i, i+1, br / 10);
        Double_t k = hFG->Integral() * fTarget / hTot->Integral(bl, br);
        cout << hTot->GetBinWidth(1) << " " << hTot->GetMaximumBin() << " " << tMax << " " << bl << "-" << br << " " << k << endl;
        geU->SetPoint(i, i+1, k);
        if (Draw)
        {
            sprintf(name, "UFC_Tot_%i", i+1);
            hTot->SetName(name);
            sprintf(name, "UFC_FG_%i", i+1);
            hFG->SetName(name);
            hTot->SetStats(0);
            SetSize(hTot);
            hTot->GetYaxis()->SetTitle("#font[12]{N}_{(n,f)} [a.u.]");
            hTot->SetLineColor(kBlue);
            hFG->SetLineColor(kGreen);
            TLegend *l = new TLegend(0.6, 0.6, 0.85, 0.8);
            l->AddEntry(hTot, "FG+BG");
            l->AddEntry(hFG, "FG");

            sprintf(name, "cT_UFC_%i", i+1);
            TCanvas *cT = new TCanvas(name);
            gPad->SetTicks(1, 1);
//            gPad->SetLogx(1);
            gPad->SetLogy(1);
            hTot->Draw("hist");
            hFG->Draw("same hist");
            l->Draw();
            Double_t y0 = cT->GetUymin(), y1 = hTot->GetMaximum();
            Double_t xval = hTot->GetBinLowEdge(bl);
            TLine *ll = new TLine(xval, y0, xval, y1);
            ll->SetLineWidth(2);
            ll->SetLineStyle(3);
            ll->Draw("same");
            xval = hTot->GetBinLowEdge(br);
            TLine *lr = new TLine(xval, y0, xval, y1);
            lr->SetLineWidth(2);
            lr->SetLineStyle(3);
            lr->Draw("same");
        }
    }
    Save(fAna, "UFC/Correction", geU);
    Save(fAna, "Simulation/Geant4", gLim1);
    Save(fAna, "Simulation/Geant4", gLim2);
    fAna->Save();
    fAna->Close();
}

void TargetGatePu(Bool_t Draw = 0)
{
    char name[64] = "";
    TFile *fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    if (!fAna) cout << "Could not open " << "fAna" << endl;
    Double_t fTarget = TargetFactor(fAna);
    TFile *fFG = TFile::Open("~/Programme/Geant4-Work/results/PuFC_ideal_c_FG.root");
    if (!fFG) cout << "Could not open " << "fFG" << endl;
    TFile *fTot = TFile::Open("~/Programme/Geant4-Work/results/PuFC_ideal_c_FG+BG.root");
    if (!fTot) cout << "Could not open " << "fTot" << endl;
    TGraphErrors *gePu = new TGraphErrors(8);
    gePu->SetName("PuFC_Target_Gate");
    gePu->SetTitle("PuFC Target scattering and Gating correction; Deposit; #font[12]{k}_{T}");
    TGraph *gLim1 = new TGraph(8);
    gLim1->SetNameTitle("PuFC_Geant4_lim1", "Left ToF gate bin; Deposit; bin nr.");
    TGraph *gLim2 = new TGraph(8);
    gLim2->SetNameTitle("PuFC_Geant4_lim2", "Right ToF gate bin; Deposit; bin nr.");
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "PuFC/ToF/PuFC_ToF_Ch.%i", i+1);
        TH1F *hFG = (TH1F*) fFG->Get(name); if (!hFG) cout << "Could not get FG " << name << endl;
        hFG->Scale(fTarget);
        TH1F *hTot = (TH1F*) fTot->Get(name); if (!hTot) cout << "Could not get Tot " << name << endl;
        Int_t tMax = (Int_t) hTot->GetBinCenter(hTot->GetMaximumBin());
        Int_t bl = hTot->FindBin(tMax - Left("PuFC"));
        Int_t br = hTot->FindBin(tMax + Right("PuFC") + 1.0) - 1;
        gLim1->SetPoint(i, i+1, bl / 10);
        gLim2->SetPoint(i, i+1, br / 10);
        Double_t k = hFG->Integral() / hTot->Integral(bl, br);
        cout << hTot->GetBinWidth(1) << " " << tMax << " " << bl << "-" << br << " " << k << endl;
        gePu->SetPoint(i, i+1, k);
        if (Draw)
        {
            sprintf(name, "PuFC_Tot_%i", i+1);
            hTot->SetName(name);
            sprintf(name, "PuFC_FG_%i", i+1);
            hFG->SetName(name);
            hTot->SetStats(0);
            SetSize(hTot);
            hTot->GetYaxis()->SetTitle("Spaltereignisse / ns");
            hTot->SetLineColor(kBlue);
            hFG->SetLineColor(kGreen);
            hTot->GetXaxis()->SetRangeUser(0, 80);
            TLegend *l = new TLegend(0.55, 0.7, 0.75, 0.9);
            l->AddEntry(hTot, "FG+BG");
            l->AddEntry(hFG, "FG");

            sprintf(name, "cT_PuFC_%i", i+1);
            TCanvas *cT = new TCanvas(name);
            gPad->SetTicks(1, 1);
//            gPad->SetLogx(1);
            gPad->SetLogy(1);
            hTot->Draw("hist");
            hFG->Draw("same hist");
            l->Draw();
            Double_t y0 = cT->GetUymin(), y1 = hTot->GetMaximum();
            Double_t xval = hTot->GetBinLowEdge(bl);
            TLine *ll = new TLine(xval, y0, xval, y1);
            ll->SetLineWidth(2);
            ll->SetLineStyle(3);
            ll->Draw("same");
            xval = hTot->GetBinLowEdge(br);
            TLine *lr = new TLine(xval, y0, xval, y1);
            lr->SetLineWidth(2);
            lr->SetLineStyle(3);
            lr->Draw("same");
        }
    }
    Save(fAna, "PuFC/Correction", gePu);
    Save(fAna, "Simulation/Geant4", gLim1);
    Save(fAna, "Simulation/Geant4", gLim2);
    fAna->Save();
    fAna->Close();
}

void Target()
{
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle(kTRUE);
    gStyle->SetLegendFont(132);
//    TargetE("/home/hoffma93/Programme/TARGET/Results/Variation/400keV_1E8_ENE", "", 100000000);
//    TargetToF("/home/hoffma93/Programme/TARGET/Results/Variation/400keV_1E8_TOF", "");
//    TargetGateU(0);
    TargetGatePu(1);

//    TargetE("/home/hoffma93/Programme/TARGET/Results/Variation/Double_Distance/FC_STARGET_15MEENE", "3.3m", 100000000);

//    TargetAng("/home/hoffma93/Programme/TARGET/Results/Variation/FC_STARGET_15MEANG", "", 100000000);
}

#endif
