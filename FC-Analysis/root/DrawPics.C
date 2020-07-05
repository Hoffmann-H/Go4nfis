#ifndef DRAWPICS_H
#define DRAWPICS_H

#include "TMathBase.h"
#include "TLine.h"
#include "TLatex.h"
#include "TColor.h"
#include "TTree.h"
#include "TStyle.h"
#include "FC.C"
#include "Runs.C"
#include "MCNPtoROOT.C"
#include "/home/hoffma93/StyleSheets/StyleSheet.C"
#include "TStyle.h"
#define lw 2
//#define a 0.25 // alpha
using namespace std;

void DrawSigma(Int_t isotope = 242)
{
    TGraphErrors *pSigma;
    if (isotope == 242)
        pSigma = new TGraphErrors("/home/hoffma93/Programme/ROOT/Data/Pu242.dat", "%lg %lg %lg");
    if (isotope == 235)
        pSigma = new TGraphErrors("/home/hoffma93/Programme/ROOT/Data/U235.dat", "%lg %lg %lg");
    if (isotope == 238)
        pSigma = new TGraphErrors("/home/hoffma93/Programme/ROOT/Data/U238.dat", "%lg %lg %lg");
    if (!pSigma)
        cout << "Isotope " << isotope << " not found." << endl;

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("^{242}Pu(n,f) Wirkungsquerschnitt; E_{n}; #sigma [b]");
    pSigma->SetFillColorAlpha(1, 0.5);
    mg->Add(pSigma);
    TCanvas *cS = new TCanvas("cS", "Evaluated Cross section");
    gPad->SetTicks(1, 1);
    mg->GetXaxis()->SetRangeUser(0, 20);
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(-0.0);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(-0.0);
    mg->Draw("a3");
    cS->Update();
}

void BiasX(TGraphErrors* pG, Double_t dx, Double_t scale = 1)
{
    Double_t x = 0, y = 0, Dx = 0, Dy = 0;
    for (int i = 0; i < pG->GetN(); i++)
    {
        pG->GetPoint(i, x, y);
        Dx = pG->GetErrorX(i);
        Dy = pG->GetErrorY(i);
        pG->SetPointError(i, Dx, scale * Dy);
        pG->SetPoint(i, x + dx, scale * y);
    }
}

TH1I* BiasX(TH1I *pH, Double_t dx)
{
    char name[64] = "";
    Int_t N = pH->GetNbinsX();
    Double_t Min = pH->GetXaxis()->GetBinLowEdge(1) + dx;
    Double_t Max = pH->GetXaxis()->GetBinLowEdge(N+1) + dx;
    sprintf(name, "%s_cal", pH->GetName());
    TH1I *pHdx = new TH1I(name, pH->GetTitle(), N, Min, Max);
    for (Int_t k = 0; k < N + 2; k++)
        pHdx->SetBinContent(k, pH->GetBinContent(k));
    return pHdx;
}

void DrawCrossSectionRuns(TFile *f, string FC = "PuFC")
{
    char name[64] = "";
    TGraphErrors *geCS[8];
    TMultiGraph *mg = new TMultiGraph();
    sprintf(name, "%s unkorrigierter WQ; Run nr.; #sigma [b]", FC.c_str());
    mg->SetTitle(name);
    mg->GetXaxis()->SetLabelSize(0.07);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(-0.0);
    mg->GetYaxis()->SetLabelSize(0.07);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(-0.0);
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.9, 0.4, 1.0, 0.8, name);
    l->SetTextFont(132);
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/CrossSection/%s_CS_raw_%i", FC.c_str(), FC.c_str(), i+1);
        geCS[i] = (TGraphErrors*)f->Get(name); if (!geCS[i]) cout << "Could not get " << name << endl;
        BiasX(geCS[i], -0.175 + i * 0.05);
        geCS[i]->SetMarkerStyle(20);
        geCS[i]->SetMarkerSize(1);
        geCS[i]->SetMarkerColor(i+2);
        geCS[i]->SetLineColor(i+2);
        mg->Add(geCS[i]);
        sprintf(name, "%i", i+1);
        l->AddEntry(geCS[i], name);
    }
    TCanvas *c = new TCanvas("c", "Runs Cross Section", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    //mg->GetXaxis()->SetNdivisions(9, kTRUE);
    mg->Draw("AP");
    l->Draw("same");
    c->Modified();
    c->Update();
}

TH1F* TH1ItoTH1F(TH1I* pH)
{
    // copy manually...
    Int_t nbins = pH->GetNbinsX();
    Int_t xmin = pH->GetBinLowEdge(1);
    Int_t xmax = pH->GetBinLowEdge(nbins + 1);
    TH1F* pH1F = new TH1F(pH->GetName(), "ToF; t / ns; counts", nbins, xmin, xmax);
    pH1F->Add(pH, 1);
    return pH1F;
}

TH1F* GetDtRate(TFile* f, string FC, string Setup, Int_t channel)
{
    char name[64] = "";
    char title[64] = "";
    Double_t chMin, chMax, tMin, tMax;
    
    // Get live time
    TH1D *pHtLive = (TH1D*) f->Get("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t t_live = pHtLive->Integral();

    // Open PH-gated TimeDiff spectrum
    sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", channel+1);
    TH1I *pHDtG = (TH1I*) f->Get(name);

    // Create rate over time histogram
    TH1F *pHDtRate = TH1ItoTH1F(pHDtG);
    sprintf(name, "%s_%s_DtG_%i", FC.c_str(), Setup.c_str(), channel+1);
    pHDtRate->SetName(name);

    // Scale
    pHDtRate->Scale(1.0 / t_live);
    
    return pHDtRate;
}

void SetSize(TH1 *h)
{
//    gPad->SetTicks(1, 1);
    h->SetStats(0);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.9);
}

void SetSize(TGraph *g)
{
//    gPad->SetTicks(1, 1);
    g->GetXaxis()->SetLabelSize(0.06);
    g->GetXaxis()->SetTitleSize(0.07);
    g->GetXaxis()->SetTitleOffset(0.9);
//    g->GetXaxis()->SetNdivisions(109);
    g->GetYaxis()->SetLabelSize(0.06);
    g->GetYaxis()->SetTitleSize(0.07);
    g->GetYaxis()->SetTitleOffset(0.9);
//    g->GetYaxis()->SetNdivisions(1005);
}

void DrawDtComp4(Int_t RebinFactor)
{
    char name[64] = "";
    TFile* fPuNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    TFile* fPuSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
    TFile* fUNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");
    TFile* fUSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root");
    TCanvas *pC[8];
//    Int_t PuFg[] = {0.03483, 0.02902, 0.0258, 0.0265, 0.0251, 0.0306, 0.0289, 0.0411};
//    Int_t PuBg[] = {61, -18, -54, 23, 8, 98, 24, 147};
//    Int_t UFg[] = {5099, 5659, 5925, 5793, 5781, 6180, 6241, 6226};
//    Int_t UBg[] = {91, 41, 74, 75, 57, 74, 62, 62};

    Double_t PuFg = 0.0412;
    Double_t PuBg = 0.0064;
    Double_t UFg = 0.1054;
    Double_t UBg = 0.0057;
    TH1F* pHDtPuNIF[8];
    TH1F* pHDtPuSB[8];
    TH1F* pHDtUNIF[8];
    TH1F* pHDtUSB[8];

    for (int i = 0; i < 8; i++)
    {
        // Create
        pHDtPuNIF[i] = GetDtRate(fPuNIF, "PuFC", "NIF", i);
        pHDtPuNIF[i]->Rebin(RebinFactor);
        pHDtPuNIF[i]->Scale(1000);
        sprintf(name, "Kanal %i; #font[12]{ToF} [ns]; #frac{#font[12]{SR}}{ns} [10^{-3} s^{-1}]", i+1);
        pHDtPuNIF[i]->SetTitle(name);
        pHDtPuNIF[i]->SetTitleSize(0.08f, "t");
        pHDtPuNIF[i]->GetXaxis()->SetRangeUser(1550, 2000);
        pHDtPuNIF[i]->GetXaxis()->SetLabelSize(0.06);
        pHDtPuNIF[i]->GetXaxis()->SetTitleSize(0.08);
        pHDtPuNIF[i]->GetXaxis()->SetTitleOffset(0.8);
        pHDtPuNIF[i]->GetYaxis()->SetLabelSize(0.06);
        pHDtPuNIF[i]->GetYaxis()->SetTitleSize(0.08);
        pHDtPuNIF[i]->GetYaxis()->SetTitleOffset(0.8);
        pHDtPuNIF[i]->SetLineColor(kMagenta);

        pHDtPuSB[i] = GetDtRate(fPuSB, "PuFC", "SB", i);
        pHDtPuSB[i]->SetName("Dia2");
        pHDtPuSB[i]->Rebin(RebinFactor);
        pHDtPuSB[i]->Scale(1000);
        pHDtPuSB[i]->SetLineColor(14);

        pHDtUNIF[i] = GetDtRate(fUNIF, "UFC", "NIF", i);
        pHDtUNIF[i]->SetName("Dia3");
        pHDtUNIF[i]->Scale(1000);
        pHDtUNIF[i]->Rebin(RebinFactor);
        pHDtUNIF[i]->SetLineColor(kBlue);

        pHDtUSB[i] = GetDtRate(fUSB, "UFC", "SB", i);
        pHDtUSB[i]->SetName("Dia4");
        pHDtUSB[i]->Scale(1000);
        pHDtUSB[i]->Rebin(RebinFactor);
        pHDtUSB[i]->SetLineColor(14);

        // Draw
        sprintf(name, "DtComp_%i", i+1);
        pC[i] = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetTicks(1, 1);
        TLegend *l = new TLegend(0.7, 0.75, 1.0, 1.0);
        l->AddEntry(pHDtPuNIF[i], "Offen");
        l->AddEntry(pHDtPuSB[i], "^{242}Pu, Schattenkegel");
        l->AddEntry(pHDtUSB[i], "^{235}U, Schattenkegel");
        pHDtPuNIF[i]->Draw();
        pHDtUNIF[i]->Draw("same");
        pHDtPuSB[i]->Draw("same");
        pHDtUSB[i]->Draw("same");
        l->Draw();

        TLatex *t1 = new TLatex();
        t1->SetNDC();
        sprintf(name, "%.4f s^{-1}", PuFg);
        t1->SetTextColor(kMagenta);
        t1->DrawLatex(0.4, 0.8, name);

        TLatex *t2 = new TLatex();
        t2->SetNDC();
        sprintf(name, "%.4f s^{-1}", PuBg);
        t2->SetTextColor(14);
        t2->DrawLatex(0.4, 0.6, name);

        TLatex *t3 = new TLatex();
        t3->SetNDC();
        sprintf(name, "%.4f s^{-1}", UFg);
        t3->SetTextColor(kBlue);
        t3->DrawLatex(0.6, 0.4, name);

        TLatex *t4 = new TLatex();
        t4->SetNDC();
        sprintf(name, "%.4f s^{-1}", UBg);
        t4->SetTextColor(14);
        t4->DrawLatex(0.6, 0.2, name);
    }
}

void DrawPuTransmission(TFile *f)
{
	char name[64] = "";
    TGraphErrors *geExp = (TGraphErrors*)f->Get("PuFC/Correction/ExpT");
    TGraphErrors *geSim = (TGraphErrors*)f->Get("Simulation/MCNP/Correction/T");
	Double_t X[] = {1, 2, 3, 4, 5, 6, 7, 8};
	Double_t Xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
	Double_t anaT[] = {094.59, 095.24, 095.89, 096.55, 097.21, 097.87, 098.54, 099.22};
    Double_t DanaT[] = {0.12, 0.10, 0.09, 0.07, 0.05, 0.04, 0.02, 0.01};
    TGraphErrors *geAna = new TGraphErrors(8, X, anaT, Xerr, DanaT);
    TMultiGraph *mgT = new TMultiGraph();

    TCanvas *c1 = new TCanvas("name", "title", 200, 10, 700, 500);

	
}

void DrawUscattering(TFile *f)
{
	char name[64] = "";
    TGraphErrors *geExp = (TGraphErrors*)f->Get("UFC/Correction/ExpS");
        if (!geExp) cout << "Could not get " << "UFC/Correction/ExpS" << endl;
    TGraphErrors *geSim = (TGraphErrors*)f->Get("Simulation/Geant4/UFC_Open/Correction/S");
        if (!geSim) cout << "Could not get " << "Simulation/Geant4/UFC_Open/Correction/S" << endl;
    TGraphErrors *geSimSB = (TGraphErrors*)f->Get("UFC/Correction/SimS");
        if (!geSimSB) cout << "Could not get " << "UFC/Correction/SimS" << endl;
    BiasX(geExp, -0.05);
    geExp->SetLineColor(1);
    geExp->SetMarkerColor(1);
    geExp->SetMarkerStyle(20);
    geExp->SetMarkerSize(2);
    geSim->SetLineColor(2);
    geSim->SetMarkerColor(2);
    geSim->SetMarkerStyle(21);
    geSim->SetMarkerSize(2);
    BiasX(geSimSB, 0.05);
    geSimSB->SetLineColor(3);
    geSimSB->SetMarkerColor(3);
    geSimSB->SetMarkerStyle(22);
    geSimSB->SetMarkerSize(2);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(geExp);
    mg->Add(geSim);
    mg->Add(geSimSB);
    mg->GetXaxis()->SetRangeUser(0.5, 8.5);
    mg->SetTitle("Direkter Anteil");
    mg->GetXaxis()->SetTitle("Deposit");
    mg->GetYaxis()->SetTitle("S");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    TLegend *l = new TLegend(0.5, 0.15, 0.85, 0.4);
    l->SetTextFont(132);
    l->AddEntry(geExp, "Schattenkegel Messung");
    l->AddEntry(geSim, "Simulation");
    l->AddEntry(geSimSB, "Simulation Schattenkegel");
    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c->Update();
}

void DrawPuCorrection(TFile *f)
{ // PuFC. S, T for MCNP, Geant4
    char name[64] = "";
    TGraphErrors *gTM = (TGraphErrors*)f->Get("Simulation/MCNP/Correction/T");
    TGraphErrors *gTG = (TGraphErrors*)f->Get("Simulation/Geant4/PuFC_Open/Correction/T");
    TGraphErrors *gSM = (TGraphErrors*)f->Get("Simulation/MCNP/Correction/S");
    TGraphErrors *gSG = (TGraphErrors*)f->Get("Simulation/Geant4/PuFC_Open/Correction/S");

    BiasX(gTM, -0.15);
    gTM->SetLineColor(1);
    gTM->SetMarkerColor(1);
    gTM->SetMarkerStyle(20);
    gTM->SetMarkerSize(2);
    BiasX(gTG, -0.05);
    gTG->SetLineColor(2);
    gTG->SetMarkerColor(2);
    gTG->SetMarkerStyle(21);
    gTG->SetMarkerSize(2);
    BiasX(gSM, 0.05);
    gSM->SetLineColor(3);
    gSM->SetMarkerColor(3);
    gSM->SetMarkerStyle(22);
    gSM->SetMarkerSize(2);
    BiasX(gSG, 0.15);
    gSG->SetLineColor(4);
    gSG->SetMarkerColor(4);
    gSG->SetMarkerStyle(23);
    gSG->SetMarkerSize(2);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gTM);
    mg->Add(gTG);
    mg->Add(gSM);
    mg->Add(gSG);
    mg->GetXaxis()->SetRangeUser(0.5, 8.5);
    mg->SetTitle("Korrekturfaktoren Plutonium");
    mg->GetXaxis()->SetTitle("Deposit");
    mg->GetYaxis()->SetTitle("Faktor");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    TLegend *l = new TLegend(0.5, 0.15, 0.85, 0.4);
    l->SetTextFont(132);
    l->AddEntry(gTM, "T MCNP");
    l->AddEntry(gTG, "T Geant4");
    l->AddEntry(gSM, "S MCNP");
    l->AddEntry(gSG, "S Geant4");
    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c->Update();
}

void DrawMonitorRate()
{
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "READ");
    TGraphErrors *gU = (TGraphErrors*)fAna->Get("UFC/NeutronField/UFC_MonitorRate");
    if (!gU) cout << "Could not get " << "UFC/NeutronField/UFC_MonitorRate" << ", run NeutronField.C" << endl;
    TGraphErrors *gPu = (TGraphErrors*)fAna->Get("PuFC/NeutronField/PuFC_MonitorRate");
    if (!gPu) cout << "Could not get " << "PuFC/NeutronField/PuFC_MonitorRate" << ", run NeutronField.C" << endl;
    gU->SetLineWidth(1);
    gU->SetMarkerStyle(5);
    gPu->SetLineWidth(1);
    gPu->SetMarkerStyle(5);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gU);
//    mg->Add(gPu);
    mg->GetXaxis()->SetTimeDisplay(1);
    mg->GetXaxis()->SetTimeFormat("#splitline{%d.%m.}{%H:%M} %F1970-01-01");
    mg->GetXaxis()->SetTitle("#font[12]{t}");
    mg->GetYaxis()->SetTitle("Monitor rate / s^{-1}");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetLabelOffset(0.05);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetXaxis()->SetNdivisions(608);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetRangeUser(800, 1800);
    TLegend *l = new TLegend(0.35, 0.2, 0.7, 0.35);
    l->SetTextFont(132);
    l->AddEntry(gU, "Calendar week 21");
//    l->AddEntry(gPu, "Calendar week 22");
    TCanvas *c1 = new TCanvas("name", "title", 1, 1, 1200, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c1->Update();
    SaveToFile(fAna, "UFC/NeutronField", gU);
    SaveToFile(fAna, "PuFC/NeutronField", gPu);
    fAna->Save();
    fAna->Close();
}

Int_t Color(Int_t i)
{
    Int_t color[] = {600,632,418,867,887,820,616,807,432,1};
//    Int_t color[] = {kBlue, kRed, kGreen, kCyan, 9, kSpring, kMagenta, kOrange};
    return color[i];
}

void DrawPeakWidth(TFile *f, string FC = "PuFC")
{
    char name[64] = "";
    sprintf(name, "%s/ToF/Gate/Peak/Right/gPeak_av", FC.c_str());
    TGraphErrors *geAv = (TGraphErrors*)f->Get(name);
    if (!geAv)
        cout << "Could not open " << name << endl;
    TGraphErrors *gePeak[8];
    TMultiGraph *mg = new TMultiGraph();
    sprintf(name, "%s deposit", FC.c_str());
    TLegend *l = new TLegend(0.3, 0.2, 0.9, 0.35, name);//0.5, 0.15, 0.85, 0.4, name); l->SetTextFont(132);
    l->SetNColumns(9);
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/ToF/Gate/Peak/gPeak_%i", FC.c_str(), i+1);
        gePeak[i] = (TGraphErrors*)f->Get(name);
        BiasX(gePeak[i], -0.7 + 0.2*i);
        gePeak[i]->SetLineColorAlpha(Color(i), 0.8);//1.0);//
        gePeak[i]->SetMarkerColorAlpha(Color(i), 0.8);//1.0);//
        gePeak[i]->SetMarkerStyle(20);
        gePeak[i]->SetMarkerSize(2);
        mg->Add(gePeak[i]);
        sprintf(name, "%i", i+1);
        l->AddEntry(gePeak[i], name, "PE");
    }
    geAv->SetLineColor(1);
    geAv->SetMarkerColor(1);
    geAv->SetMarkerStyle(20);
    geAv->SetMarkerSize(2);
    mg->Add(geAv);
    l->AddEntry(geAv, "mean", "PE");
    l->SetTextFont(132);
    mg->GetXaxis()->SetRangeUser(0, 120);
    mg->SetTitle("Peak intervals");
    mg->GetXaxis()->SetTitle("gate length / ns");
    mg->GetYaxis()->SetTitle("events");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c->Update();
}

void DrawSimPeak(TFile *f, string Run, string key, Int_t ch, string Simulation)
{
    char name[64] = "";
    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *fExp = TFile::Open(name); if (!fExp) cout << "Could not open " << name << endl;
    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1F *pH1Exp = (TH1F*)fExp->Get(name); if (!pH1Exp) cout << "Could not get " << name << endl;

    sprintf(name, "Simulation/%s/%s_%s/FitToF/%s_FitT_%s_%i", Simulation.c_str(), FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), ch+1);
    TH1D *pH1Sim = (TH1D*)f->Get(name);
    if (!pH1Sim) cout << "Could not get " << name << endl;

    sprintf(name, "%s Deposit %i", FC.c_str(), ch+1);
    TLegend *l = new TLegend(0.5, 0.6, 0.85, 0.85, name);
    pH1Exp->SetLineWidth(2);
    SetSize(pH1Exp);
    pH1Exp->GetXaxis()->SetRangeUser(0, 100);
    l->AddEntry(pH1Exp, "Messung");

    Int_t Gate[] = {Gate_1(ch, FC), Gate_2(ch, FC)};
    Double_t IntExp = pH1Exp->Integral(Gate[0], Gate[1]);
    Double_t IntSim = pH1Sim->Integral(pH1Sim->FindBin(pH1Exp->GetBinLowEdge(Gate[0])), pH1Sim->FindBin(pH1Exp->GetBinLowEdge(Gate[1]+1)-1));
    Double_t IntSimTotal = pH1Sim->Integral();
    cout << Gate[0] << " " << Gate[1] << " " << IntExp << " " << IntSim << " " << IntSimTotal << endl;
//    Double_t Scale = IntExp / IntSim * pH1Exp->GetXaxis()->GetBinWidth(1) / pH1Sim->GetXaxis()->GetBinWidth(1);
//    pH1Sim->Scale(Scale);
    pH1Sim->SetLineColor(kRed);
    pH1Sim->SetLineWidth(2);
    sprintf(name, "%s-Simulation", Simulation.c_str());
    l->AddEntry(pH1Sim, name);

//    Double_t min = pH1Exp->GetMinimum();
    Double_t max = pH1Exp->GetMaximum();
    Double_t xl = pH1Exp->GetBinLowEdge(Gate_1(ch, FC));
    Double_t xr = pH1Exp->GetBinLowEdge(Gate_2(ch, FC)+1);
    TLine *ll = new TLine(xl, 0, xl, max);
    ll->SetLineStyle(3);
    ll->SetLineWidth(2);
    TLine *lr = new TLine(xr, 0, xr, max);
    lr->SetLineStyle(3);
    lr->SetLineWidth(2);
    Double_t Efficiency = IntSim / IntSimTotal;

    sprintf(name, "c%s_%i", FC.c_str(), ch+1);
    TCanvas *c = new TCanvas(name, "title", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pH1Exp->Sumw2();
    pH1Exp->Draw("p");
    pH1Sim->Draw("same hist");
    l->Draw();
    ll->Draw();
    lr->Draw();;
//    TLatex *t = new TLatex();
//    sprintf(name, "%.1f%% im Zeitfenster", 100*Efficiency);
//    t->SetTextColor(kRed);
//    t->DrawLatexNDC(0.5, 0.5, name);
    c->Update();//*/
}

Double_t GetStartTime(string Run)
{
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *f = TFile::Open(name);
    TH1D *pH = (TH1D*)f->Get("Histograms/Raw/Scaler/Rates/H1RawRate_48");
    Int_t N = pH->GetNbinsX();
    Int_t bin = 0;
    do {
        bin++;
    } while (bin < N && pH->GetBinContent(bin) == 0);
    return pH->GetBinLowEdge(bin);
}

Double_t GetRealTime(string Run)
{
    char name[  128] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *f = TFile::Open(name);
    if (!f)
        cout << "Could not open " << name << endl;
    TH1D *pH = (TH1D*)f->Get("Histograms/Raw/Scaler/Rates/H1RawRate_48");
    if (!pH)
        cout << "Could not open " << "Histograms/Raw/Scaler/Rates/H1RawRate_48" << endl;
    return pH->Integral();
}

void DrawMonitorRatio(TFile *f, Int_t ch)
{
    char name[64] = "";
    sprintf(name, "UFC/CrossSection/UFC_MonitorRatio_%i", ch+1);
    TGraphErrors *geUFC = (TGraphErrors*)f->Get(name);
    if (!geUFC)
        cout << "Could not open " << name << endl;
    sprintf(name, "PuFC/CrossSection/PuFC_MonitorRatio_%i", ch+1);
    TGraphErrors *gePuFC = (TGraphErrors*)f->Get(name);
    if (!gePuFC)
        cout << "Could not open " << name << endl;

    Double_t x, y, yerr, tStart, tLength;
    for (Int_t j = 0; j < 7; j++)
    {
        geUFC->GetPoint(j, x, y);
        yerr = geUFC->GetErrorY(j);
        string Run = GetRunName("UFC", j);
        tStart = GetStartTime(Run);
        tLength = GetRealTime(Run);
        geUFC->SetPoint(j, tStart + 0.5*tLength, y);
        geUFC->SetPointError(j, 0.5*tLength, yerr);

        gePuFC->GetPoint(j, x, y);
        yerr = gePuFC->GetErrorY(j);
        Run = GetRunName("PuFC", j);
        tStart = GetStartTime(Run);
        tLength = GetRealTime(Run);
        gePuFC->SetPoint(j, tStart + 0.5*tLength, y);
        gePuFC->SetPointError(j, 0.5*tLength, yerr);
    }
    geUFC->SetLineWidth(2);
    gePuFC->SetLineWidth(2);
    geUFC->SetLineColor(kBlue);
    gePuFC->SetLineColor(kMagenta);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(geUFC);
    mg->Add(gePuFC);
    mg->GetXaxis()->SetTimeDisplay(1);
    mg->GetXaxis()->SetTimeFormat("%d.%m. %H:%M %F1970-01-01");
    mg->GetXaxis()->SetTitle("#font[12]{t}");
    mg->GetYaxis()->SetTitle("Signal zu Monitor");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    TLegend *l = new TLegend(0.4, 0.2, 0.6, 0.3);
    l->AddEntry(geUFC, "KW 21");
//    l->AddEntry(gePuFC, "KW 22");
    TCanvas *c1 = new TCanvas("name", "title", 1, 1, 1200, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c1->Update();
}//*/

void DrawMonitorRatio(TFile *f, string FC = "PuFC")
{
    char name[64] = "";
    Int_t ch_start = 0, ch_stop = 8;
    TGraphErrors *g[8];
    for (Int_t i = ch_start; i < ch_stop; i++)
    {
        sprintf(name, "%s/CrossSection/%s_MonitorRatio_%i", FC.c_str(), FC.c_str(), i+1);
        g[i] = (TGraphErrors*)f->Get(name);
        if (!g[i])
            cout << "Could not open " << name << endl;
    }
    TGraphErrors *gAv = new TGraphErrors(7);
//    gAv->SetName("gAv");

    Double_t x, y, yerr, tStart, tLength, sumY, sumYerr2;
    for (Int_t j = 0; j < 7; j++)
    {
        string Run = GetRunName(FC, j);
        tStart = GetStartTime(Run);
        tLength = GetRealTime(Run);
        sumY = 0;
        sumYerr2 = 0;
        for (Int_t i = ch_start; i < ch_stop; i++)
        {
            g[i]->GetPoint(j, x, y);
            yerr = g[i]->GetErrorY(j);
            g[i]->SetPoint(j, tStart + 0.5*tLength, y);
            g[i]->SetPointError(j, 0.5*tLength, yerr);
            sumY += y;
            sumYerr2 += pow(yerr, 2);
        }
        gAv->SetPoint(j, tStart + 0.5*tLength, sumY / (ch_stop - ch_start));
        gAv->SetPointError(j, 0.5*tLength, sqrt(sumYerr2) / (ch_stop - ch_start));
    }
    TMultiGraph *mg = new TMultiGraph();
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.4, 0.2, 0.6, 0.3, name);
    for (Int_t i = ch_start; i < ch_stop; i++)
    {
        g[i]->SetLineWidth(1);
        g[i]->SetLineColor(i+2);
        mg->Add(g[i]);
        sprintf(name, "%i", i+1);
        l->AddEntry(g[i], name);
    }
    gAv->SetLineWidth(2);
    gAv->SetLineColor(1);
    mg->Add(gAv);
    l->AddEntry(gAv, "Mittel");
//    Double_t X[] = {140050000, 140075000};
//    Double_t Y[] = {0, 0};
//    TGraph *gInv = new TGraph(2, X, Y);
//    gInv->SetLineColor(1);
//    mg->Add(gInv);
//    mg->GetXaxis()->SetRangeUser(140000000, 140080000);
    mg->GetXaxis()->SetNdivisions(1005);
    mg->GetXaxis()->SetTimeDisplay(1);
    mg->GetXaxis()->SetTimeFormat("%d.%m. %H:%M %F1970-01-01");
    mg->GetXaxis()->SetTitle("#font[12]{t}");
    mg->GetYaxis()->SetTitle("Signal zu Monitor");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    TCanvas *c1 = new TCanvas("name", "title", 1, 1, 1200, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c1->Update();
}

void DrawPuSponFis()
{
    char name[64] = "";
    TFile *fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    TFile *fSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
    TFile *fSF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root");
//    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");
    TGraphErrors *geNIF = (TGraphErrors*)fNIF->Get("Analysis/results/effSponFis"); if (!geNIF) cout << "Could not get " << "Analysis/results/effSponFis" << endl;
    TGraphErrors *geSB = (TGraphErrors*)fSB->Get("Analysis/results/effSponFis"); if (!geSB) cout << "Could not get " << "Analysis/results/effSponFis" << endl;
    TGraphErrors *geSF = new TGraphErrors(8);
    TGraphErrors *geToni = new TGraphErrors(8);

    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
    Double_t DrSF_Toni[] = {0.0013, 0.0012, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0013};
    TH1D *pHt = (TH1D*)fSF->Get("Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t t_live = pHt->Integral();
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I *pHDt = (TH1I*)fSF->Get(name); if (!pHDt) cout << "Could not get " << name << endl;
        Double_t Integral = pHDt->Integral();
        geSF->SetPoint(i, i+1, Integral / t_live);
        geSF->SetPointError(i, 0, sqrt(Integral) / t_live);
        geToni->SetPoint(i, i+1, rSF_Toni[i]);
        geToni->SetPointError(i, 0, DrSF_Toni[i]);
    }
    BiasX(geNIF, -0.15);
    BiasX(geSB, -0.05);
    BiasX(geSF, +0.05);
    BiasX(geToni, +0.15);
    geNIF->SetLineWidth(5);
    geSB->SetLineWidth(5);
    geSF->SetLineWidth(5);
    geToni->SetLineWidth(5);
    geNIF->SetLineColor(2);
    geSB->SetLineColor(3);
    geSF->SetLineColor(4);
    geToni->SetLineColor(1);
    TLegend *l = new TLegend(0.4, 0.2, 0.6, 0.3);
    l->AddEntry(geNIF, "Vordergrund");
    l->AddEntry(geSB, "Schattenkegel");
    l->AddEntry(geSF, "Beam off");
    l->AddEntry(geToni, "Toni");

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(geNIF);
    mg->Add(geSB);
    mg->Add(geSF);
    mg->Add(geToni);
    mg->GetXaxis()->SetTitle("Deposit");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetTitle("Spaltrate [1/s]");
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.7);
    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c->Update();
}

void DrawTimeDiff(TFile *f, Int_t ch, string FC = "PuFC", string Setup = "")
{ // Geht nicht: NIF.root gibt Fehler beim Oeffnen
    char name[64] = "";
    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    cout << name << endl;
    TH1I *pHraw = (TH1I*) f->Get(name);
//    TH1I *pH = BiasX(pHraw, 100);
    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
//    SetSize(pHraw);
    pHraw->GetXaxis()->SetNdivisions(510);
    pHraw->Draw();
    c->Update();
}

void DrawExfor()
{
    TGraphErrors *geManabe = new TGraphErrors("/home/hoffma93/Experiment/Evaluations/Pu-U-Manabe-Plot.txt", "%lg %lg %lg %lg");

    Double_t x, y;
    geManabe->GetPoint(0, x, y);
    cout << x << ", " << y << endl;
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(geManabe);
    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
    mg->Draw();
    c->Update();
}

void DrawSourceE(Double_t E = 14.97, Double_t DE = 0.3)
{
    char name[64] = "";
    TGraph *g = new TGraph("/home/hoffma93/Programme/ROOT/Data/Source_E.dat");

    TCanvas *c = new TCanvas("Source_E", "Source_E", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
//    gPad->SetLogx(1);
    g->SetNameTitle("Source_E", "TARGET Neutron Energy Spectrum");
    g->GetXaxis()->SetTitle("#font[12]{E} / MeV");
    g->GetXaxis()->SetTitleSize(0.08);
    g->GetXaxis()->SetLabelSize(0.06);
    g->GetXaxis()->SetTitleOffset(0.8);
    g->GetYaxis()->SetTitle("#font[12]{n}(#font[12]{E}) [a.u.]");
    g->GetYaxis()->SetTitleSize(0.08);
    g->GetYaxis()->SetLabelSize(0.06);
    g->GetYaxis()->SetTitleOffset(0.7);
    g->SetLineWidth(2);

    sprintf(name, "E = %.2f #pm %.2f MeV", E, DE);
    TLatex *l = new TLatex(14.5, 0.2, name);

    g->Draw();
    l->Draw();
    c->Modified();
    c->Update();
}

void DrawEnELBE()
{
    TGraph *g = new TGraph("/home/hoffma93/Programme/ROOT/Data/nELBE_E.dat");

    TCanvas *c = new TCanvas("E_nELBE", "E_nELBE", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogx(1);
    g->SetNameTitle("E_nELBE", "nELBE Neutron Energy Spectrum");
    g->GetXaxis()->SetTitle("#font[12]{E} / MeV");
    g->GetXaxis()->SetTitleSize(0.08);
    g->GetXaxis()->SetLabelSize(0.06);
    g->GetXaxis()->SetTitleOffset(0.8);
    g->GetYaxis()->SetTitle("#font[12]{n}(#font[12]{E}) [a.u.]");
    g->GetYaxis()->SetTitleSize(0.08);
    g->GetYaxis()->SetLabelSize(0.06);
    g->GetYaxis()->SetTitleOffset(0.7);
    g->SetLineWidth(2);

    g->Draw();
    c->Modified();
    c->Update();
}

void DrawTargetE(TFile *f, string Subfolder, Int_t Population)
{
    if (Population != 1.E8) cout << "Warning: Hard-coded population 1.E8 !" << endl;
    char name[64] = "";
    sprintf(name, "Simulation/Target/%sE_Dir", Subfolder.c_str());
    TH1D *pDir = (TH1D*)f->Get(name);
    if (!pDir)
        cout << "Could not open " << name << endl;
    sprintf(name, "Simulation/Target/%sE_Sc", Subfolder.c_str());
    TH1D *pSc = (TH1D*)f->Get(name);
    sprintf(name, "Simulation/Target/%sE_Tot", Subfolder.c_str());
    TH1D *pTot = (TH1D*)f->Get(name);
    pDir->Scale(Population/pTot->Integral());
    pSc->Scale(Population/pTot->Integral());
    pTot->Scale(Population/pTot->Integral());
    pDir->SetLineColor(kRed);
    pSc->SetLineColor(kBlue);
//    pTot->SetLineColor(kWhite);
    pDir->SetLineWidth(lw);
    pSc->SetLineWidth(lw);
    pTot->SetLineWidth(0);
    pTot->SetLineColor(kBlack);
    pTot->SetTitle("; #font[12]{E}_{n} / MeV; #font[12]{N} / 10keV");

    sprintf(name, "1 #upoint 10^{8} simulated ions");
    gStyle->SetLegendFont(132);
    TLegend *l = new TLegend(0.2, 0.72, 0.7, 0.9, name);
    sprintf(name, "%.2f%% direct neutrons", 100 * pDir->Integral() / pTot->Integral());
    l->AddEntry(pDir, name);
    sprintf(name, "%.2f%% indirect neutrons", 100 * pSc->Integral() / pTot->Integral());
    l->AddEntry(pSc, name);
    SetSize(pTot);
    TCanvas *c = new TCanvas("Target_E", "Target_E", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);
    pTot->Draw("hist");
    pDir->Draw("same hist");
    pSc->Draw("same hist");
    l->Draw();
    c->Update();
}

void DrawTargetToF(TFile *f, string Subfolder)
{
    char name[64] = "";
    sprintf(name, "Simulation/Target/%sToF_Tot", Subfolder.c_str());
    TH1D *pTarget = (TH1D*)f->Get(name); if (!pTarget) cout << "Could not open " << name << endl;
    sprintf(name, "Simulation/Target/%sToF_E_Tot", Subfolder.c_str());
    TH1D *pE = (TH1D*)f->Get(name); if (!pE) cout << "Could not open " << name << endl;
    cout << pTarget->Integral() << "  " << pE->Integral() << endl;
    pTarget->Scale(1.0 / pTarget->Integral());
    pE->Scale(1.0 / pE->Integral());

    pE->SetTitle("; #font[12]{t} / ns; #font[12]{N} / 100ps");
    pE->SetLineColor(kBlue);
    pE->SetStats(0);
    pTarget->SetLineColor(kRed);
    SetSize(pE);
//    pE->GetXaxis()->SetLabelSize(0.06);
//    pE->GetXaxis()->SetTitleSize(0.07);
//    pE->GetXaxis()->SetTitleOffset(0.8);
//    pE->GetYaxis()->SetLabelSize(0.06);
//    pE->GetYaxis()->SetTitleSize(0.07);
//    pE->GetYaxis()->SetTitleOffset(0.7);

    gStyle->SetLegendFont(132);
    TLegend *l = new TLegend(0.6, 0.7, 0.9, 0.85);
    l->AddEntry(pTarget, "TARGET-Ausgabe");
    l->AddEntry(pE, "ToF(#font[12]{E})");
    TCanvas *cT = new TCanvas("cT", "TARGET simulated ToF", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);
    pE->Draw("hist");
    pTarget->Draw("same hist");
    l->Draw();
    cT->Update();
}

void DrawTargetAng(TFile *f, string Subfolder = "")
{
    char name[64] = "";
    sprintf(name, "Simulation/Target/%sAngular", Subfolder.c_str());
    TGraph *gAng = (TGraph*)f->Get(name); if (!gAng) cout << "Could not get " << name << endl;
    Double_t Deviation = gAng->Eval(4.84) / gAng->Eval(0);
    cout << "Max. Abweichung bis 4.84°: " << Deviation << endl;
    new TCanvas();
    SetSize(gAng);
    gAng->GetXaxis()->SetNdivisions();
    gAng->GetXaxis()->SetRangeUser(0, 180);
    gAng->GetYaxis()->SetRangeUser(0, 80);
    gAng->Draw();
    gAng->GetXaxis()->SetTitle("Winkel / #circ");
    gAng->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} / mb/sr");
    gAng->GetYaxis()->SetTitleOffset(0.8);
//    TLatex *t = new TLatex();
}

void DrawTargetAir(TFile *fAna)
{
    char name[64] = "";
    sprintf(name, "Simulation/Target/E_Tot");
    TH1D *h1 = (TH1D*)fAna->Get(name);
    h1->SetName("E_Tot_1");
    h1->SetStats(0);
    SetSize(h1);
    h1->SetLineColor(kBlue);

    sprintf(name, "Simulation/Target/3.3m/E_Tot");
    TH1D *h2 = (TH1D*)fAna->Get(name);
    h2->SetName("E_Tot_2");
    h2->SetLineColor(kRed);
//    h2->Scale(4);
//    cout << h1->Integral() / h2->Integral() << endl;

    TLegend *l = new TLegend(0.2, 0.7, 0.5, 0.9, "Flugstrecke");
    l->AddEntry(h1, "1,65 m");
    l->AddEntry(h2, "3,30 m");
    new TCanvas("cTd");
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);
    h1->Draw("hist");
    h2->Draw("same hist");
    l->Draw();
}

void DrawGEF(Int_t isotope, string spectrum)
{
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Experiment/Carlson-Korrektur/GEF/GEF_92_%i_nf_%s_1e6/GEF_U%i_nf_%s_1e6.root", isotope, spectrum.c_str(), isotope, spectrum.c_str());
    TFile *fGEF = TFile::Open(name);
    if (!fGEF)
        cout << "Could not open " << name << endl;
    TH1F *hMass = (TH1F*)fGEF->Get("Mass");
    sprintf(name, "^{%i}U, %s, Spaltfragmentmasse; A; Spaltfragmente", isotope, spectrum.c_str());
    hMass->SetTitle(name);
    SetSize(hMass);
    hMass->GetYaxis()->SetTitleOffset(1.1);
    TCanvas *c1 = new TCanvas("Mass", "Mass", 200, 10, 700, 500);
    hMass->Draw();
    c1->Update();
    TH1F *hTKE = (TH1F*)fGEF->Get("TKE");
    sprintf(name, "^{%i}U, %s, Totale kinetische Energie; #font[12]{TKE} [MeV]; Spaltereignisse", isotope, spectrum.c_str());
    hTKE->SetTitle(name);
    SetSize(hTKE);
    hTKE->GetYaxis()->SetTitleOffset(1.1);
    TCanvas *c2 = new TCanvas("TKE", "TKE", 200, 10, 700, 500);
    hTKE->Draw();
    c2->Update();
    TH2F *hNvsZ = (TH2F*)fGEF->Get("NvsZ");
    sprintf(name, "^{%i}U, %s, Tochternuklide; Z; A", isotope, spectrum.c_str());
    hNvsZ->SetTitle(name);
    SetSize(hNvsZ);
    hNvsZ->GetZaxis()->SetLabelSize(0.06);
    TCanvas *c3 = new TCanvas("NvsZ", "NvsZ", 200, 10, 700, 500);
    hNvsZ->Draw("colz");
    c3->Update();
}

void DrawGEF(string hist_name = "Mass", string x_title = "#font[12]{A}")
{
    char name[128] = "";
    Int_t Isotopes[] = {235, 235, 238, 238};
    string Spectra[] = {"nELBE", "15MeV", "nELBE", "15MeV"};
    TH1F *pH[4];
    TCanvas *c4 = new TCanvas(hist_name.c_str(), hist_name.c_str(), 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(0.125);
    TLegend *l = new TLegend(0.75, 0.75, 1, 1);
    l->SetTextFont(132);
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(name, "/home/hoffma93/Programme/GEF/results/GEF_U%i_nf_%s_1e6.root", Isotopes[j], Spectra[j].c_str());
        TFile *f = TFile::Open(name);
        if (!f) cout << "Could not open " << name << endl;
        pH[j] = (TH1F*)f->Get(hist_name.c_str());
        if (!pH[j]) cout << "Could not get " << hist_name << endl;
        pH[j]->SetLineWidth(2);
        pH[j]->SetLineColorAlpha(Isotopes[j] == 235 ? kBlue : kOrange, 0.8);
        pH[j]->SetLineStyle(strcmp(Spectra[j].c_str(), "nELBE") ? 2 : 1);
        sprintf(name, "%s ^{%i}U", Spectra[j].c_str(), Isotopes[j]);
        l->AddEntry(pH[j], name);
    }
    pH[0]->SetTitle(hist_name.c_str());
    pH[0]->GetXaxis()->SetTitle(x_title.c_str());
    pH[0]->GetYaxis()->SetTitle("Spaltfragmente");
    SetSize(pH[0]);
    pH[0]->GetYaxis()->SetTitleOffset(1.1);
    pH[0]->GetXaxis()->SetNdivisions(805);
    pH[0]->Draw("hist");
    pH[1]->Draw("same hist"); pH[2]->Draw("same hist"); pH[3]->Draw("same hist");
    l->Draw("same");
    c4->Update();
}

void DrawFFRange(TFile *f)
{
    TH1F *pH0 = (TH1F*)f->Get("Carlson/U235/nELBE/U235_nELBE_FF_Range");
    if (!pH0) cout << "Could not get " << "Carlson/U235/nELBE/U235_nELBE_FF_Range" << endl;
    TH1F *pH1 = (TH1F*)f->Get("Carlson/U235/15MeV/U235_15MeV_FF_Range");
    if (!pH1) cout << "Could not get " << "Carlson/U235/15MeV/U235_15MeV_FF_Range" << endl;
    TH1F *pH2 = (TH1F*)f->Get("Carlson/U238/nELBE/U238_nELBE_FF_Range");
    if (!pH2) cout << "Could not get " << "Carlson/U238/nELBE/U238_nELBE_FF_Range" << endl;
    TH1F *pH3 = (TH1F*)f->Get("Carlson/U238/15MeV/U238_15MeV_FF_Range");
    if (!pH3) cout << "Could not get " << "Carlson/U238/15MeV/U238_15MeV_FF_Range" << endl;

    pH0->SetLineWidth(2);
    pH0->SetLineColorAlpha(kBlue, 0.8);
    pH0->SetLineStyle(1);
    pH0->GetXaxis()->SetTitle("#font[12]{R'} / mg/cm^{2}");
    pH1->SetLineWidth(2);
    pH1->SetLineColorAlpha(kBlue, 0.8);
    pH1->SetLineStyle(2);
    pH2->SetLineWidth(2);
    pH2->SetLineColorAlpha(kOrange, 0.8);
    pH2->SetLineStyle(1);
    pH3->SetLineWidth(2);
    pH3->SetLineColorAlpha(kOrange, 0.8);
    pH3->SetLineStyle(2);
    TLegend *l = new TLegend(0.75, 0.75, 1, 1);
    l->SetTextFont(132);
    l->AddEntry(pH0, "#font[12]{n}ELBE ^{235}U");
    l->AddEntry(pH1, "15MeV ^{235}U");
    l->AddEntry(pH2, "#font[12]{n}ELBE ^{238}U");
    l->AddEntry(pH3, "15MeV ^{238}U");

    pH0->SetTitle("FF Range");
    pH0->GetXaxis()->SetTitle("#font[12]{R'} [mg/cm^{2}]");
    pH0->GetYaxis()->SetTitle("#font[12]{N}");
    SetSize(pH0);
    pH0->GetXaxis()->SetRangeUser(3, 10);
    pH0->GetYaxis()->SetRangeUser(0, 13000);

    TCanvas *cR = new TCanvas("cR", "FF Range", 200, 10, 700, 500);
    gPad->SetTicks(1,1);
    pH0->Draw("hist");
    pH1->Draw("same hist"); pH2->Draw("same hist"); pH3->Draw("same hist");
    l->Draw("same");
    cR->Update();
}

void DrawGEFTree(string spectrum = "15MeV", Int_t isotope = 235)
{
    string xName = "TKEPost";
    string yName = "MassPost";
    char name[64] = "";
    sprintf(name, "~/Programme/GEF/results/GEF_U%i_nf_%s_1e6.root", isotope, spectrum.c_str());
    TFile *f = TFile::Open(name); if (!f) cout << "Could not open " << name << endl;
    TTree *t = (TTree*) f->Get("FissionFragmentTree"); if (!t) cout << "Could not get " << "FissionFragmentTree" << endl;
    Int_t N = t->GetEntries();
    Float_t x;
    Int_t y;
    t->SetBranchAddress(xName.c_str(), &x);
    t->SetBranchAddress(yName.c_str(), &y);
    for (Int_t event = 0; event < N; event++)
    {
        t->GetEvent(event);

    }

//    new TCanvas();
//    t->Draw("TKEPost:MassPost");
}

TGraphErrors* GetCarlson(TFile *fAna, string var_name = "eta", Int_t isotope = 235)
{
    char name[64] = "";
    sprintf(name, "Carlson/U%i/U%i_%s", isotope, isotope, var_name.c_str());
    TGraphErrors *g = (TGraphErrors*) fAna->Get(name); if (!g) cout << "Could not get " << name << endl;
    sprintf(name, "U%i_%s", isotope, var_name.c_str());
    g->SetName(name);
    return g;
}

void DrawA2(TFile *fAna)
{
    TGraphErrors *ge235a2 = GetCarlson(fAna, "a2", 235);
    TGraphErrors *ge238a2 = GetCarlson(fAna, "a2", 238);
    SetSize(ge235a2);
    ge235a2->SetLineColor(kBlue);
    ge238a2->SetLineColor(kRed);

    new TCanvas("cA2");
    gPad->SetTicks(1, 1);
    ge235a2->Draw("ap");
    ge238a2->Draw("same p");
    ge235a2->SetTitle("; #font[12]{E}_{n}; #font[12]{a}_{2}");
    ge235a2->GetXaxis()->SetNdivisions();
    ge235a2->GetYaxis()->SetNdivisions();
    TLegend *l = new TLegend(0.6, 0.6, 0.8, 0.8);
    l->AddEntry(ge235a2, "{}^{235}U(n,f)");
    l->AddEntry(ge238a2, "{}^{238}U(n,f)");
    l->Draw();
}

void DrawEta(TFile *fAna)
{
    TGraphErrors *ge235eta = GetCarlson(fAna, "eta", 235);
//    TGraphErrors *ge238eta = GetCarlson(fAna, "eta", 238);
    ge235eta->SetTitle("; #font[12]{E}_{n}; #eta");
    SetSize(ge235eta);

    new TCanvas("cEta");
    gPad->SetTicks(1, 1);
    ge235eta->Draw("AP");
    ge235eta->GetXaxis()->SetNdivisions();
    ge235eta->GetYaxis()->SetNdivisions();
    TLegend *l = new TLegend(0.2, 0.7, 0.5, 0.8);
    l->AddEntry(ge235eta, "#font[12]{n}ELBE, {}^{235}U(n,f)");
    l->Draw();
}

void DrawInefficiency(TFile *fAna, string spectrum = "nELBE", Int_t ch = 0)
{
    char name[64] = "";
    sprintf(name, "Carlson/U235/%s/I_U235_%s_%i", spectrum.c_str(), spectrum.c_str(), ch+1);
    TGraphErrors *ge235I = (TGraphErrors*) fAna->Get(name); if (!ge235I) cout << "Could not get " << name << endl;
    sprintf(name, "Carlson/U238/%s/I_U238_%s_%i", spectrum.c_str(), spectrum.c_str(), ch+1);
    TGraphErrors *ge238I = (TGraphErrors*) fAna->Get(name); if (!ge238I) cout << "Could not get " << name << endl;
    SetSize(ge235I);
    ge235I->SetLineColor(kBlue);
    ge238I->SetLineColor(kRed);

    new TCanvas("cI");
    gPad->SetTicks(1, 1);
    gPad->SetLogx(1);
    ge235I->Draw("ap");
    ge238I->Draw("same p");
    ge235I->SetTitle("; #font[12]{E}_{n} / MeV; #font[12]{I}");
    ge235I->GetXaxis()->SetRangeUser(0.1, 20);
    ge235I->GetYaxis()->SetRangeUser(0.03, 0.06);
    ge235I->GetXaxis()->SetNdivisions();
    ge235I->GetYaxis()->SetNdivisions();
    TLegend *l = new TLegend(0.4, 0.7, 0.6, 0.9);
    l->AddEntry(ge235I, "{}^{235}U(n,f)");
    l->AddEntry(ge238I, "{}^{238}U(n,f)");
    l->Draw();
}

void DrawQDC(TFile *fAna, string FC, Int_t ch, TFile *f, string Setup = "FG")
{
    char name[64] = "";

    sprintf(name, "Histograms/Raw/QDC/low/H1RawQDCl_%i", ch+1);
    TH1I *pH = (TH1I*)f->Get(name);
    pH->SetStats(0);
    SetSize(pH);
    pH->SetLineColor(kBlack);
    pH->SetTitle("; #font[12]{Q} / Channel; Counts");
    Double_t xP = pH->GetBinCenter(pH->GetMaximumBin()),
             x0 = QDCcut(ch, FC);
    Double_t x1 = 1520;
    pH->GetXaxis()->SetRangeUser(0, x1);
    pH->GetYaxis()->SetRangeUser(0.5, 2. * pH->GetMaximum());

    TH1I *pHff = (TH1I*)pH->Clone();
    pHff->GetXaxis()->SetRange(x0, x1);
    pHff->SetFillColorAlpha(kBlue, 0.5);

    sprintf(name, "Histograms/Analysis/FC/QDC/low/trig/H1AnaQDCl_trig_%i", ch+1);
    TH1I *pHt = (TH1I*)f->Get(name);
    pHt->SetLineColor(kBlue);

    sprintf(name, "%s/QDC/Fit/%s/pol4_%s_%i", FC.c_str(), Setup.c_str(), Setup.c_str(), ch+1);
    TF1 *pF = (TF1*)fAna->Get(name);
    pF->SetLineWidth(lw);
    pF->SetLineColor(kRed);

    sprintf(name, "QDCfit_%i", ch+1);
    TCanvas *pC = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);
    pH->Draw("hist");
    pHff->Draw("same hist");
    pHt->Draw("same hist");
    pF->Draw("same");
    pC->Update();

    Double_t y0 = pow(10, pC->GetUymin()),
             y1 = pow(10, pC->GetUymax());
    cout << x0 << " " << y0 << " " << x1 << " " << y1 << endl;
    TLine *lCut = new TLine(x0, y0, x0, y1);
    lCut->Draw("same");

    Double_t y = pF->GetMinimum();
    TGraph *g = new TGraph(4);
    g->SetPoint(0, xP, 0);
    g->SetPoint(1, xP, y);
    g->SetPoint(2, x0, y);
    g->SetPoint(3, x0, 0);
    g->SetLineColorAlpha(kRed, 0.5);
    g->SetFillColorAlpha(kRed, 0.5);
    g->Draw("same LF");

    pHt->GetXaxis()->SetRange(x0, x1);
    y0 = 2.0 * pHt->GetMaximum();
    pHt->GetXaxis()->SetRange();
    y1 = y0;

    TArrow *aFF = new TArrow(x0, y0, x1, y1, 0.02, "<|>");
    aFF->SetAngle(30);
    aFF->Draw();

    TLegend *lQDC = new TLegend(0.6, 0.7, 0.85, 0.8);
    lQDC->SetTextFont(132);
    sprintf(name, "%s Kanal %i", FC.c_str(), ch+1);
    lQDC->AddEntry(pHt, name);
//    lQDC->Draw("same");

    TLatex *t1 = new TLatex();
//    t1->SetTextColor();
    t1->DrawLatexNDC(0.17, 0.55, "#alpha's");
    t1->Draw();

    TLatex *t2 = new TLatex();
//    t2->SetTextColor();
    t2->DrawLatexNDC(0.5, 0.45, "fission fragments");
    t2->Draw();

    TLatex *t3 = new TLatex();
//    t3->SetTextColor();
    sprintf(name, "#varepsilon #approx %.0f%%", 100*pH->Integral(pH->FindBin(x0), 4096)/((x0-xP)*y + pH->Integral(pH->FindBin(x0), 4096)));
    t3->DrawLatexNDC(0.16, 0.18, name);
    t3->Draw();
}

void DrawQDCsum(TFile *fAna, string FC, Int_t ch)
{
//    char name[64] = "";
//    sprintf(name, "Histograms/Raw/QDC/low/H1RawQDCl_%i", ch+1);
//    TH1I *pH;
//    if (!strcpm(FC.c_str(), "PuFC"))
//    {
//        TFile* fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
//        pH = (TH1I*)fNIF->Get(name);
//        TFile* fSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
//        TH1I *pSB = (TH1I*)fSB->Get(name);
//        TFile* fSF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root");
//        TH1I *pSF = (TH1I*)fSF->Get(name);
//        pH->Add(pSB);
//        pH->Add(pSF);
//        fNIF->Close();
//        fSB->Close();
//        fSF->Close();
//    } else {
//        TFile* fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");
//        pH = (TH1I*)fNIF->Get(name);
//        TFile* fSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root");
//        TH1I *pSB = (TH1I*)fSB->Get(name);
//        pH->Add(pSB);
//        fNIF->Close();
//        fSB->Close();
//    }

//    sprintf(name, "%s/QDC/Fit/%s/pol4_%s_%i", FC.c_str(), Setup.c_str(), Setup.c_str(), ch+1);
//    TF1 *pF = (TF1*)fAna->Get(name);
//    pH->SetStats(0);
//    pH->GetXaxis()->SetTitleSize(0.05);
//    pH->GetYaxis()->SetTitleSize(0.05);
//    pH->SetLineColor(kBlue);
//    pH->Draw();
//    pF->SetLineWidth(lw);
//    pF->SetLineColor(kRed);
//    pF->Draw("same");
}

void DrawDtInt(TFile *fAna, string FC, Int_t ch, Int_t Left = 7, Int_t Right = 7, TGraphErrors *g = 0, Bool_t Save = 0)
{
    char name[64] = "";
    TFile *f;
    TH1I *pH;
    TF1 *fT;
    Double_t x = 0, y = 0;
    if (!strcmp(FC.c_str(), "PuFC"))
    {
        sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
        f = TFile::Open(name);
        sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
        pH = (TH1I*)f->Get(name);
        sprintf(name, "PuFC/ToF/Background/NIF/Total/NIF_fT_%i", ch+1);
        fT = (TF1*)fAna->Get(name);
        x = 0.5, y = 0.25;
    } else {
        sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");
        f = TFile::Open(name);
        sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
        pH = (TH1I*)f->Get(name);
        sprintf(name, "UFC/ToF/Background/UFC_NIF/Total/UFC_NIF_fT_%i", ch+1);
        fT = (TF1*)fAna->Get(name);
        x = 0.3, y = 0.3;
    }
    TH1I *pHfill = (TH1I*)pH->Clone();
    Double_t Bg = fT->GetParameter(0);
    Double_t DBg = fT->GetParError(0);
    Double_t chi2 = fT->GetChisquare();
    Int_t dof = fT->GetNDF();
    Double_t Peak = pH->Integral(Gate_1(ch, FC, Left), Gate_2(ch, FC, Right));
    Double_t DPeak = sqrt(Peak);
    Double_t C = Peak - (Left+Right+1) * Bg;
    Double_t DC = sqrt( pow(DPeak, 2) + pow((Left+Right+1) * DBg, 2) );

    if (g)
    {
        g->SetPoint(ch, ch+1, C);
        g->SetPointError(ch, 0, DC);
        cout << ch << " " << C << " +- " << DC << endl;
        return;
    }

    TGraph *gL = new TGraph(2);
    gL->SetPoint(0, pH->GetBinLowEdge(Gate_0(ch, FC)), Bg);
    gL->SetPoint(1, pH->GetBinLowEdge(Gate_a(ch, FC)), Bg);
    TGraph *gR = new TGraph(2);
    gR->SetPoint(0, pH->GetBinLowEdge(Gate_b(ch, FC)), Bg);
    gR->SetPoint(1, pH->GetBinLowEdge(Gate_3(ch, FC)), Bg);

    sprintf(name, "cDt_%s_%i", FC.c_str(), ch+1);
    TCanvas *cT = new TCanvas(name, name, 200, 10, 1400, 1000);
    gPad->SetTicks(1, 1);
    SetSize(pH);
    pH->GetXaxis()->SetNdivisions(810);
    pH->GetYaxis()->SetTitleOffset(1.1);
    gL->SetLineWidth(2);
    gR->SetLineWidth(2);
    gL->SetLineColor(kRed);
    gR->SetLineColor(kRed);
    pHfill->SetFillColorAlpha(kBlue, 0.5);
    pHfill->GetXaxis()->SetRange(Gate_1(ch, FC, Left), Gate_2(ch, FC, Right));
    pH->Draw("hist");
    pHfill->Draw("same hist");
    gL->Draw("same");
    gR->Draw("same");
    TLatex *t[5];
    t[0] = new TLatex();
    t[0]->SetTextColor(kRed);
    sprintf(name, "#font[12]{<B>} = %.2f #pm %.2f ns^{-1}", Bg, DBg);
    t[0]->DrawLatexNDC(x, y + 0.32, name);
    t[0]->Draw();
    t[1] = new TLatex();
    t[1]->SetTextColor(kBlue);
    sprintf(name, "#font[12]{C}_{P} = %.0f #pm %.0f", Peak, DPeak);
    t[1]->DrawLatexNDC(x, y + 0.16, name);
    t[1]->Draw();
    t[2] = new TLatex();
    t[2]->SetTextColor(kBlue);
    sprintf(name, "#font[12]{t}_{P} = %i ns", Left+Right+1);
    t[2]->DrawLatexNDC(x, y + 0.08, name);
    t[2]->Draw();
    t[3] = new TLatex();
    t[3]->SetTextColor(kBlack);
    sprintf(name, "#scale[1.25]{#font[12]{C}_{(n,f)} = %.f #pm %.f}", C, DC);
    t[3]->DrawLatexNDC(x, y - 0.01, name);
    t[3]->Draw();
    t[4] = new TLatex();
    t[4]->SetTextColor(kRed);
    sprintf(name, "#chi^{2} / #nu = %.2f, #nu = %i", chi2 / dof, dof);
    t[4]->DrawLatexNDC(x, y + 0.24, name);
    t[4]->Draw();
    sprintf(name, "/home/hoffma93/Pictures/TimeDiff/%s_Integral_%i.pdf", FC.c_str(), ch+1);
    if (Save)
        cT->SaveAs(name);
}

void DrawDtChange(TFile *fAna, string FC = "PuFC")
{
    char name[32] = "";
    TGraphErrors *ge40 = new TGraphErrors(8);
    ge40->SetNameTitle("ge40", "Flugzeitpeak ohne Untergrund; Deposit; #font[12]{C}");
    TGraphErrors *ge15 = new TGraphErrors(8);
    ge15->SetNameTitle("ge15", "Flugzeitpeak ohne Untergrund; Deposit; #font[12]{C}");
    Double_t x, y1, y0, fis0 = 0, fis1 = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        DrawDtInt(fAna, FC, i, 15, 15, ge40);
        DrawDtInt(fAna, FC, i, 7, 7, ge15);
        ge40->GetPoint(i, x, y0);
        fis0 += y0;
        ge15->GetPoint(i, x, y1);
        fis1 += y1;
        cout << i+1 << "   " << y1 / y0 << endl;
    }
    new TCanvas("cPG", "Peak-Inhalt");
    gPad->SetTicks(1, 1);
    SetSize(ge40);
    ge15->SetLineColor(kRed);
    ge40->Draw();
    ge15->Draw("same");
    TLegend *l = new TLegend(0.4, 0.65, 0.6, 0.85, "Gate-Breite");
    l->AddEntry(ge40, "31 ns");
    l->AddEntry(ge15, "15 ns");
    l->Draw();
    TLatex *t = new TLatex();
    sprintf(name, "Gate: %f", fis1 / fis0);
    t->DrawLatexNDC(0.4, 0.6, name);
}

void DrawDtGate(TFile *fAna, string FC, Int_t ch, TFile *f)
{
    char name[64] = "";
    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1I *pH0 = (TH1I*)f->Get(name);
    SetSize(pH0);
    pH0->SetTitle("");
    TH1I *pH1 = (TH1I*)pH0->Clone();
    pH1->GetXaxis()->SetRange(Gate_1(ch, FC), Gate_2(ch, FC));
    pH1->SetFillColorAlpha(kRed, 0.5);
    TH1I *pH2 = (TH1I*)pH0->Clone();
    pH2->GetXaxis()->SetRange(Gate_0(ch, FC), Gate_a(ch, FC)-1);
    pH2->SetFillColorAlpha(kGray, 0.5);
    TH1I *pH3 = (TH1I*)pH0->Clone();
    pH3->GetXaxis()->SetRange(Gate_b(ch, FC)+1, Gate_3(ch, FC)-1);
    pH3->SetFillColorAlpha(kGray, 0.5);

    sprintf(name, "cDt_%s_%i", FC.c_str(), ch+1);
    TCanvas *cT = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);

    pH0->Draw();
    pH1->Draw("same");
    pH2->Draw("same");
    pH3->Draw("same");
    cT->Update();

    TLegend *l;
    if (!strcmp(FC.c_str(), "PuFC")) l = new TLegend(0.5, 0.4, 0.8, 0.5);
    else l = new TLegend(0.2, 0.4, 0.5, 0.5);
    sprintf(name, "%s Kanal %i", FC.c_str(), ch+1);
    l->AddEntry(pH0, name);
    l->SetTextFont(132);
//    l->Draw();

    Int_t bin[] = {Gate_0(ch, FC), Gate_a(ch, FC), Gate_1(ch, FC), Gate_2(ch, FC)+1, Gate_b(ch, FC)+1, Gate_3(ch, FC)};
    Double_t y0 = cT->GetUymin(), y1 = cT->GetUymax();
    for (Int_t j = 0; j < 6; j++)
    {
        Double_t xval = pH0->GetBinLowEdge(bin[j]);
        TLine *line = new TLine(xval, y0, xval, y1);
        line->SetLineStyle(3);
        line->Draw("same");
    }
}

void DrawDtBackground(TFile *fAna, string Setup, Int_t ch, TFile *f)
{
    char name[64] = "";
    string FC = (Setup[0] == 'U') ? "UFC" : "PuFC";
    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1I *pH = (TH1I*)f->Get(name);
    SetSize(pH);
    sprintf(name, "%s/ToF/Background/%s/Left/%s_fL_%i", FC.c_str(), Setup.c_str(), Setup.c_str(), ch+1);
    TF1 *fL = (TF1*)fAna->Get(name);
    if (!fL) cout << "Could not get " << name << endl;
    sprintf(name, "%s/ToF/Background/%s/Right/%s_fR_%i", FC.c_str(), Setup.c_str(), Setup.c_str(), ch+1);
    TF1 *fR = (TF1*)fAna->Get(name);
    if (!fR) cout << "Could not get " << name << endl;

    Double_t x1, x2, y, Dy;
    fL->GetRange(x1, x2);
    x1 = -80;
    y = fL->GetParameter(0);
    Dy = fL->GetParError(0);
    sprintf(name, "%.1f#pm%.1f", y, Dy);
    TLatex *t1 = new TLatex(x1, y + pH->GetMaximum()/5, name);
    sprintf(name, "#chi^{2}/#font[12]{dof} = %.2f", fL->GetChisquare() / fL->GetNDF());
    TLatex *t2 = new TLatex(x1, y + pH->GetMaximum()/10, name);

    fR->GetRange(x1, x2);
    x1 = 80;
    y = fR->GetParameter(0);
    Dy = fR->GetParError(0);
    sprintf(name, "%.1f#pm%.1f", y, Dy);
    TLatex *t3 = new TLatex(x1, y + pH->GetMaximum()/5, name);
    sprintf(name, "#chi^{2}/#font[12]{dof} = %.2f", fR->GetChisquare() / fR->GetNDF());
    TLatex *t4 = new TLatex(x1, y + pH->GetMaximum()/10, name);

    fL->SetLineColor(kRed);
    t1->SetTextColor(kRed);
    t2->SetTextColor(kRed);
    fR->SetLineColor(kRed);
    t3->SetTextColor(kRed);
    t4->SetTextColor(kRed);

    TLegend *l = new TLegend(0.2, 0.4, 0.5, 0.5);
    sprintf(name, "%s Kanal %i", FC.c_str(), ch+1);
    l->AddEntry(pH, name);

    sprintf(name, "cDtBg_%s_%i", FC.c_str(), ch+1);
    TCanvas *cT = new TCanvas(name, name, 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pH->Draw("hist");
    fL->Draw("same");
    fR->Draw("same");
    l->Draw();
    t1->Draw();
    t2->Draw();
    t3->Draw();
    t4->Draw();
}

void DrawSponFisStability(TFile *fAna, string FC)
{
    char name[64] = "";
    TH1D *hC[8];
    TF1 *fC[8];
    Int_t color[] = {600, 632, 418, 867, 887, 820, 616, 807, 1, 921};
    Double_t t0 = 0, t1 = 0;
    Double_t y0 = 0, y1 = 0;
    if (FC == "UFC") {
        t0 = 1.40058E9; t1 = 1.40072E9;
    } else {
        t0 = 1.4011E9; t1 = 1.4017E9;
    }

    TLatex LaText, ChiText, NuText;
    LaText.SetTextSize(0.20);
    ChiText.SetTextSize(0.20); NuText.SetTextSize(0.20);
    ChiText.SetTextColor(kBlack); NuText.SetTextColor(kBlack);
    Double_t chi2 = 0, x, y;
    Int_t dof = 1;

    sprintf(name, "cST_%s", FC.c_str());
    TCanvas *cST = new TCanvas(name);
    cST->Divide(2,4);

    // Load and Draw 8 channels' constance
    for (Int_t i = 0; i < 8; i++)
    {
        // Switch pad
        cST->cd(i+1)->SetTicks(1, 1);

        // Histogram: Background vs Time
        sprintf(name, "%s/Stability/SF/%s_%i", FC.c_str(), FC.c_str(), i+1);
        hC[i] = (TH1D*) fAna->Get(name); if (!hC[i]) cout << "Could not get " << name << endl;
        hC[i]->SetStats(0);
        hC[i]->SetMarkerStyle(20);
        hC[i]->SetMarkerSize(0.5);
        hC[i]->SetMarkerColorAlpha(color[i], 0.5);
        hC[i]->SetLineColorAlpha(color[i], 0.5);
        hC[i]->GetXaxis()->SetTitle("");
        hC[i]->GetXaxis()->SetRangeUser(t0, t1);
        hC[i]->GetXaxis()->SetTimeDisplay(1);
        hC[i]->GetXaxis()->SetTimeFormat("%d.%m.");
        hC[i]->GetXaxis()->SetNdivisions(606);
        hC[i]->GetXaxis()->SetTitleSize(0.20);
        hC[i]->GetXaxis()->SetLabelSize(0.15);
        hC[i]->GetXaxis()->SetLabelOffset(0.03);
//        if (i % 2 == 0)
            hC[i]->GetYaxis()->SetTitle("#it{#dot{C}}_{BG} / s^{-1}");
//        else
//            hC[i]->GetYaxis()->SetTitle("");
        hC[i]->GetYaxis()->SetTitleSize(0.20);
        hC[i]->GetYaxis()->SetLabelSize(0.15);
        hC[i]->GetYaxis()->SetTitleOffset(0.3);

        // Constant fit: Background vs Time
        sprintf(name, "%s/Stability/SF/%s_%i_fit", FC.c_str(), FC.c_str(), i+1);
        fC[i] = (TF1*) fAna->Get(name); if (!fC[i]) cout << "Could not get " << name << endl;
        fC[i]->SetLineColor(color[i]);
        chi2 = fC[i]->GetChisquare();
        dof = fC[i]->GetNDF();
        y0 = fC[i]->GetParameter(0) - fC[i]->GetParError(0);
        y1 = fC[i]->GetParameter(0) + fC[i]->GetParError(0);
        hC[i]->GetYaxis()->SetRangeUser(y0, y1);

        hC[i]->Draw("P");
        fC[i]->Draw("same");

        //add text to plot
        sprintf(name, "%s Ch.%i", FC.c_str(), i+1);
        LaText.SetTextColor(color[i]);
        LaText.DrawLatexNDC(.7,.8, name);

        sprintf(name, "#chi^{2}/dof = %.2f", chi2 / dof);
        ChiText.DrawLatexNDC(.2,.235, name);

        sprintf(name, "dof = %i", dof);
        NuText.DrawLatexNDC(.75,.235, name);
    }

    cST->Update();

    // Drawing properties

//    SetSize(hC[0]);

        // print reduced chi^2
//        chi2 = fC[i]->GetChisquare();
//        dof = fC[i]->GetNDF();
//        cout << "Channel " << i+1 << ": chi^2 = " << chi2 << ", dof = " << dof << endl;
//        t[i] = new TLatex();
//        x = hC[i]->GetBinLowEdge(10);
//        y = fC[i]->GetParameter(0);
//        t[i]->SetTextColor(color[i]);
//        sprintf(name, "#chi^{2} / #nu = %.2f", chi2 / dof);
//        t[i]->DrawLatex(x, y - 0.5, name);
}

void IndFisSum(TFile *fAna, string FC)
{
    Int_t color = kBlack;
    char name[64] = "";

    /// sum up plots: Signal vs Run nr
    TGraphErrors *g[8];
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/Stability/%s_nf_%i", FC.c_str(), FC.c_str(), i+1);
        g[i] = (TGraphErrors*) fAna->Get(name); if (!g[i]) cout << "Could not get " << name << endl;
    }
    sprintf(name, "%s_nf_sum", FC.c_str());
    TGraphErrors *gCnf = (TGraphErrors*)g[0]->Clone(name);

    for (Int_t j = 0; j < 7; j++) // Iterate over graph points
    {
        Double_t x, y, ex, ey, sum = 0, sume2 = 0;
        for (Int_t i = 0; i < 8; i++)
        {
            g[i]->GetPoint(j, x, y);
            sum += y; // sum up y
//            sume2 += pow(g[i]->GetErrorY(j), 2); // sum up squared y error
            sume2 += g[i]->GetErrorY(j); // sum up y error
        }
        gCnf->SetPoint(j, x, sum);
//        gCnf->SetPointError(j, 0, sqrt(sume2));
        gCnf->SetPointError(j, 0, sume2);
    }

    /// Do constant fit
    sprintf(name, "%s_nf_FG_fit", FC.c_str());
    TF1 *fFG = new TF1(name, "pol0", 0, 5.5);
//    sprintf(name, "%s_nf_BG_fit", FC.c_str());
//    TF1 *fFG = new TF1(name, "pol0", 5.5, 8);
    gCnf->Fit(name, "LR0Q");

    /// Drawing properties
    sprintf(name, "%s all ch; Run number; #font[12]{C}_{(n,f)} / #font[12]{C}_{NM}", FC.c_str());
    gCnf->SetTitle(name);
    SetSize(gCnf);
    gCnf->SetLineColor(color);
    gCnf->SetMarkerColor(color);
    gCnf->SetMarkerStyle(20);
    gCnf->SetMarkerSize(0.5);
    fFG->SetLineColor(color);
//    fBG->SetLineColor(color);

    /// Draw
    TCanvas *cNFS = new TCanvas("cNFS");
    gPad->SetTicks(1, 1);
    gPad->SetTopMargin(0.06);
    gCnf->Draw("ap");
//    gCnf->GetYaxis()->SetRangeUser(0.0, 0.0001); //-0.00002, 0.00004);//
    gCnf->GetXaxis()->SetNdivisions(110);
    gCnf->GetXaxis()->SetLabelOffset(0.0);
    fFG->Draw("same");
    TLatex *tFG = new TLatex();
    tFG->SetNDC();
    sprintf(name, "#chi^{2} / #font[12]{dof} = %.2f", fFG->GetChisquare() / fFG->GetNDF());
    tFG->SetTextColor(color);
    tFG->SetTextSize(0.08);
    tFG->DrawLatex(0.3, 0.5, name);

//    fBG->Draw("same");
//    TLatex *tBG = new TLatex();
//    tBG->SetNDC();
//    sprintf(name, "#chi^{2} / #font[12]{dof} = %.2f", fBG->GetChisquare() / fBG->GetNDF());
//    tBG->SetTextColor(color[ch]);
//    tBG->SetTextSize(0.05);
//    tBG->DrawLatex(0.75, 0.6, name);

    TLine *line = new TLine(0.72, gPad->GetBottomMargin(), 0.72, 1 - gPad->GetTopMargin());//4.5, -0.00002, 4.5, 0.00004);//
    line->SetLineStyle(3);
    line->SetNDC();
    line->Draw();
}

void IndFisStability(TFile *fAna, string FC, Int_t ch)
{
    Int_t color[]   = {600,632,418,867,887,820,616,807,432,1};
    char name[64] = "";
    // Get plot: Signal vs Run nr
    sprintf(name, "%s/Stability/%s_nf_%i", FC.c_str(), FC.c_str(), ch+1);
    TGraphErrors *gCnf = (TGraphErrors*) fAna->Get(name); if (!gCnf) cout << "Could not get " << name << endl;
    sprintf(name, "%s/Stability/FG/%s_nf_%i_fit", FC.c_str(), FC.c_str(), ch+1);
    TF1 *fFG = (TF1*) fAna->Get(name); if (!fFG) cout << "Could not get " << name << endl;
    sprintf(name, "%s/Stability/BG/%s_nf_%i_fit", FC.c_str(), FC.c_str(), ch+1);
    TF1 *fBG = (TF1*) fAna->Get(name); if (!fBG) cout << "Could not get " << name << endl;

    // Drawing properties
    sprintf(name, "%s Ch. %i; Run number; #font[12]{C}_{(n,f)} / #font[12]{C}_{NM}", FC.c_str(), ch+1);
    gCnf->SetTitle(name);
    SetSize(gCnf);
    gCnf->SetLineColor(1);//Color(ch));
    gCnf->SetMarkerColor(1);//Color(ch));
    gCnf->SetMarkerStyle(20);
    gCnf->SetMarkerSize(1.5);
    fFG->SetLineColor(1);//Color(ch));
    fBG->SetLineColor(1);//Color(ch));

    // Draw
    gPad->SetTicks(1, 1);
    gPad->SetTopMargin(0.06);
    gCnf->Draw("ap");
    if (FC == "PuFC")
    {
        gCnf->GetYaxis()->SetRangeUser(-0.00002, 0.000045);
//        Double_t x, y;
//        for (uint i = 0; i < gCnf->GetN(); i++)
//        {
//            gCnf->GetPoint(i, x, y);
//            gCnf->SetPoint(i, x + 7, y);
//        }
    }
    else
        gCnf->GetYaxis()->SetRangeUser(0.0, 0.00012);
    gCnf->GetXaxis()->SetNdivisions(110);
    gCnf->GetXaxis()->SetLabelOffset(0.0);
    fFG->Draw("same");

    TLatex LaText;
    LaText.SetTextSize(0.08);
    sprintf(name, "%s Ch.%i", FC.c_str(), ch+1);
//    LaText.SetTextColor(color[ch]);
    LaText.DrawLatexNDC(.45,.85, name);

    TLatex *tFG = new TLatex();
    tFG->SetNDC();
    sprintf(name, "#chi^{2} / #font[12]{dof} = %.2f", fFG->GetChisquare() / fFG->GetNDF());
//    tFG->SetTextColor(color[ch]);
    tFG->SetTextSize(0.07);
    tFG->DrawLatex(0.3, 0.6, name);

    fBG->Draw("same");
//    TLatex *tBG = new TLatex();
//    tBG->SetNDC();
//    sprintf(name, "#chi^{2} / #font[12]{dof} = %.2f", fBG->GetChisquare() / fBG->GetNDF());
//    tBG->SetTextColor(color[ch]);
//    tBG->SetTextSize(0.05);
//    tBG->DrawLatex(0.75, 0.6, name);

    TLine *line = new TLine(0.72, gPad->GetBottomMargin(), 0.72, 1 - gPad->GetTopMargin());//0.605, gPad->GetBottomMargin(), 0.605, 1 - gPad->GetTopMargin());//
    line->SetLineStyle(3);
    line->SetNDC();
    line->Draw();
}

void DrawIndFisStability(TFile *fAna, string FC)
{
    char name[64] = "";
    sprintf(name, "cST_%s", FC.c_str());
//    TCanvas *cST = new TCanvas(name);
//    cST->Divide(2, 2);
    for (Int_t i = 0; i < 8; i++)
//    Int_t i = 7;
    {
//        cST->cd(i+1);
        sprintf(name, "cST_%s_%i", FC.c_str(), i+1);
        TCanvas *cST = new TCanvas(name);
        gPad->SetTopMargin(0.06);
        IndFisStability(fAna, FC, i);
        cST->Update();
    }
//    cST->Update();
}

TH1D* SignalStability(string Run, Int_t ch)
{
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *f = TFile::Open(name); if (!f) cout << "Could not open " << name << endl;
    sprintf(name, "Histograms/Raw/Scaler/Rates/H1RawRate_%i", ch+1);
    TH1D *h = (TH1D*) f->Get(name); if (!h) cout << "Could not get " << name << endl;
    sprintf(name, "%s_%i", Run.c_str(), ch+1);
    h->SetName(name);
    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
    sprintf(name, "%s Ch. %i; ; scaler rate / s^{-1}", FC.c_str(), ch+1);
    h->SetTitle(name);
    return h;
}

void DrawSignalStability()
{
    // Pu
    TH1D *hPu[8];
    TLegend *lPu = new TLegend(0.91, 0.1, 1.0, 0.9, "Ch");
    for (Int_t i = 0; i < 8; i++)
    {
        hPu[i] = SignalStability("NIF", i);
        TH1D *hPuSB = SignalStability("SB", i);
        hPu[i]->Add(hPuSB);
        TH1D *hPuSF = SignalStability("SF", i);
        hPu[i]->Add(hPuSF);
        hPu[i]->Rebin(4); // Minute binning
        hPu[i]->Scale(1.0 / 4);
        hPu[i]->SetLineColor(Color(i));
        lPu->AddEntry(hPu[i], to_string(i+1).c_str(), "PE");
    }
    new TCanvas("cSigPu");
    gPad->SetTicks(1, 1);
    hPu[0]->SetTitle("");
    hPu[0]->GetXaxis()->SetTimeDisplay(1);
    hPu[0]->GetXaxis()->SetTimeFormat("%d.%m. %H:%M");
    hPu[0]->GetXaxis()->SetNdivisions(605);
    hPu[0]->GetYaxis()->SetRangeUser(0, 100);
    SetSize(hPu[0]);
    hPu[0]->Draw("hist");
    for (Int_t i = 1; i < 8; i++)
        hPu[i]->Draw("same hist");
    lPu->Draw();
    // U
    TH1D *hU[8];
    TLegend *lU = new TLegend(0.91, 0.1, 1.0, 0.9, "Ch");
    for (Int_t i = 0; i < 8; i++)
    {
        hU[i] = SignalStability("UFC_NIF", i);
        TH1D *hUSB = SignalStability("UFC_SB", i);
        hU[i]->Add(hUSB);
        hU[i]->Rebin(4); // Minute binning
        hU[i]->Scale(1.0 / 4);
        hU[i]->SetLineColor(Color(i));
        lU->AddEntry(hU[i], to_string(i+1).c_str(), "PE");
    }
    new TCanvas("cSigU");
    gPad->SetTicks(1, 1);
    hU[0]->SetTitle("");
    hU[0]->GetXaxis()->SetTimeDisplay(1);
    hU[0]->GetXaxis()->SetTimeFormat("%d.%m. %H:%M");
    hU[0]->GetXaxis()->SetNdivisions(604);
    hU[0]->GetYaxis()->SetRangeUser(0, 240);
    SetSize(hU[0]);
    hU[0]->Draw("hist");
    for (Int_t i = 1; i < 8; i++)
        hU[i]->Draw("same hist");
    lU->Draw();
}
TH2I* QDCvsTime(string Run, Int_t ch)
{
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *f = TFile::Open(name); if (!f) cout << "Could not open " << name << endl;
    sprintf(name, "Histograms/Raw/QDC/low/QDCvsTime/H2RawQDCvsTime_%i", ch+1);
    TH2I *h = (TH2I*) f->Get(name); if (!h) cout << "Could not get " << name << endl;
    sprintf(name, "%s_QDCvsTime_%i", Run.c_str(), ch+1);
    h->SetName(name);
    return h;
}
void DrawQDCvsTime(Int_t ch)
{
    gROOT->SetStyle("SinglePadColZStyle");
    gROOT->ForceStyle(kTRUE);
    char name[128] = "";
    // Pu
    TH2I *hPu = QDCvsTime("NIF", ch);
    TH2I *hPuSB = QDCvsTime("SB", ch);
    hPu->Add(hPuSB);
    TH2I *hPuSF = QDCvsTime("SF", ch);
    hPu->Add(hPuSF);
    hPu->RebinY(4);
    SetSize(hPu);
    hPu->GetYaxis()->SetTimeDisplay(1);
    hPu->GetYaxis()->SetTimeFormat("%d.%m.");
    hPu->GetYaxis()->SetNdivisions(610);
    hPu->GetXaxis()->SetNdivisions(508);
    new TCanvas("cQDCtPu");
    gPad->SetRightMargin(0.15);
    gPad->SetTicks(1,1);
    gPad->SetLogz(1);
    hPu->Draw("colz");

    // U
    TH2I *hU = QDCvsTime("UFC_NIF", ch);
    TH2I *hUSB = QDCvsTime("UFC_SB", ch);
    hU->Add(hUSB);
    hU->RebinY(4);
//    hU->Scale(1.0 / 10 / 4);
    SetSize(hU);
    hU->GetYaxis()->SetTimeDisplay(1);
    hU->GetYaxis()->SetTimeFormat("%d.%m. %H:%M");
    hU->GetYaxis()->SetNdivisions(610);
    hU->GetXaxis()->SetNdivisions(508);
    new TCanvas("cQDCtU");
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.15);
    gPad->SetTicks(1,1);
    gPad->SetLogz(1);
    hU->Draw("colz");
}

void DrawTvsE(TFile *fAna, string Simulation = "Geant4", string FC = "PuFC", string key = "real", Int_t ch = 0, Bool_t Save = 0)
{
    gROOT->SetStyle("SinglePadColZStyle");
    gROOT->ForceStyle(kTRUE);
    char name[128] = "";
    sprintf(name, "Simulation/%s/%s_%s/ToFvsEkin/%s_ToFvsEkin_%s_Ch.%i", Simulation.c_str(), FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), ch+1);
    TH2D *h2TvsE = (TH2D*) fAna->Get(name); if (!h2TvsE) cout << "Could not get " << name << endl;

    h2TvsE->SetTitle("");
    h2TvsE->SetStats(0);
    SetSize(h2TvsE);
    h2TvsE->GetZaxis()->SetLabelSize(0.06);
    h2TvsE->GetXaxis()->SetRangeUser(25, 75);
    h2TvsE->GetZaxis()->SetRangeUser(1.E-12, strcmp(key.c_str(), "SB") ? 1.E-08 : 1.E-08);

    sprintf(name, "TvsE_%s_%s_%s_%i", Simulation.c_str(), FC.c_str(), key.c_str(), ch+1);
    TCanvas *c = new TCanvas(name, name, 1, 1, 1915, 1120);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetTicks(1,1);
    gPad->SetLogz(1);
    h2TvsE->Draw("colz");
    c->Update();
    sprintf(name, "~/Pictures/Simulation/ToFvsEkin/%s_%s_%s_%i.png", Simulation.c_str(), FC.c_str(), key.c_str(), ch+1);
    if (Save)
        c->SaveAs(name);
}

void DrawT(TFile *fAna, string Simulation = "Geant4", string FC = "PuFC", Int_t ch = 0, Bool_t Save = 0)
{
    char name[128] = "";
    sprintf(name, "Simulation/%s/%s_real/EffToF/%s_ProjT_real_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), ch+1);
    TH1D *h1real = (TH1D*) fAna->Get(name); if (!h1real) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/%s/%s_ideal/EffToF/%s_ProjT_ideal_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), ch+1);
    TH1D *h1ideal = (TH1D*) fAna->Get(name); if (!h1ideal) cout << "Could not get " << name << endl;

    // drawing properties
    h1real->GetYaxis()->SetRangeUser(1.E-12, 1.E-05);
    h1real->SetTitle("; #font[12]{t} / ns; #font[12]{N}_{(n,f)}/#font[12]{Y}");
    h1real->SetStats(0);
    SetSize(h1real);
    h1real->SetLineColor(4);
    h1ideal->SetLineColor(2);
    sprintf(name, "%s Ch.%i", FC.c_str(), ch+1);
    TLegend *l = new TLegend(0.6, 0.6, 0.85, 0.85, name);
    l->AddEntry(h1real, "Geometrie");
    l->AddEntry(h1ideal, "Vakuum");

    // Depict summation
    Double_t x, y; Int_t l1, l2;
    sprintf(name, "Simulation/Geant4/%s_Geant4_lim1", FC.c_str());
    TGraph *g1 = (TGraph*) fAna->Get(name); if (!g1) cout << "Could not get " << name << endl;
    g1->GetPoint(ch, x, y); l1 = y;
    sprintf(name, "Simulation/Geant4/%s_Geant4_lim2", FC.c_str());
    TGraph *g2 = (TGraph*) fAna->Get(name); if (!g2) cout << "Could not get " << name << endl;
    g2->GetPoint(ch, x, y); l2 = y;
    TH1D *h1realInt = (TH1D*) h1real->Clone();
    h1realInt->GetXaxis()->SetRange(l1, l2);
    h1realInt->SetFillColorAlpha(4, 0.2);
    TH1D *h1idealInt = (TH1D*) h1ideal->Clone();
    h1idealInt->GetXaxis()->SetRange(l1, l2);
    h1idealInt->SetFillColorAlpha(kRed - 9, 1);

    sprintf(name, "ToF_%s_%s_%i", Simulation.c_str(), FC.c_str(), ch+1);
    TCanvas *cT = new TCanvas(name, name, 1, 1, 1400, 1000);
    gPad->SetTopMargin(0.05);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    gPad->SetTicks(1,1);
    gPad->SetLogy(1);
    h1real->Draw("hist");
    h1idealInt->Draw("same hist");
    h1realInt->Draw("same hist");
    h1ideal->Draw("same hist");
    h1real->Draw("same hist");
    l->Draw();
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    sprintf(name, "%.7f", h1real->Integral(l1, l2));
    t1->SetTextColor(4);
    t1->DrawLatex(0.3, 0.7, name);
    TLatex *t2 = new TLatex();
    t2->SetNDC();
    sprintf(name, "%.7f", h1ideal->Integral(l1, l2));
    t2->SetTextColor(2);
    t2->DrawLatex(0.3, 0.6, name);
    sprintf(name, "~/Pictures/Simulation/FissionToF/ToF_%s_%s_%i.pdf", Simulation.c_str(), FC.c_str(), ch+1);
    if (Save)
        cT->SaveAs(name);
}

void DrawSimNotebook(TFile *fAna, Bool_t Expand = 0, Bool_t Save = 0)
{
    char name[64] = "";
    TH1F *h1Exp[3][8]; // {PuFC-FG, UFC-FG, UFC-SB} x ch
    TF1 *f1Sim[5][8]; // {G4-PuFC-FG, MCNP-PuFC-FG, G4-UFC-FG, MCNP-UFC-FG, G4-UFC-SB} x ch
    TCanvas *cFitToF[3][8]; // {PuFC-FG, UFC-FG, UFC-SB} x ch

    if (!Expand) {
        cFitToF[0][0] = new TCanvas("PuFC, FG");
        cFitToF[0][0]->Divide(4, 2);
        cFitToF[1][0] = new TCanvas("UFC, FG");
        cFitToF[1][0]->Divide(4, 2);
        cFitToF[2][0] = new TCanvas("UFC, SB");
        cFitToF[2][0]->Divide(4, 2);
    }
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "PuFC/ToF/Total/NIF/H1AnaHZDRDtG_%i", i+1);
        h1Exp[0][i] = (TH1F*) fAna->Get(name); if (!h1Exp[0][i]) cout << "Could not get " << name << endl;
        sprintf(name, "Dt_PuFC_FG_%i", i+1); h1Exp[0][i]->SetName(name);
        sprintf(name, "UFC/ToF/Total/UFC_NIF/H1AnaHZDRDtG_%i", i+1);
        h1Exp[1][i] = (TH1F*) fAna->Get(name); if (!h1Exp[1][i]) cout << "Could not get " << name << endl;
        sprintf(name, "Dt_UFC_FG_%i", i+1); h1Exp[1][i]->SetName(name);
        sprintf(name, "UFC/ToF/Total/UFC_SB/H1AnaHZDRDtG_%i", i+1);
        h1Exp[2][i] = (TH1F*) fAna->Get(name); if (!h1Exp[2][i]) cout << "Could not get " << name << endl;
        sprintf(name, "Dt_UFC_SB_%i", i+1); h1Exp[1][i]->SetName(name);

        sprintf(name, "Simulation/Geant4/PuFC_real/FitToF/fNIF_%i", i+1);
        f1Sim[0][i] = (TF1*) fAna->Get(name); if (!f1Sim[0][i]) cout << "Could not get " << name << endl;
        sprintf(name, "fDt_Geant4_PuFC_FG_%i", i+1);
        sprintf(name, "Simulation/MCNP/PuFC_real/FitToF/fNIF_%i", i+1);
        f1Sim[1][i] = (TF1*) fAna->Get(name); if (!f1Sim[1][i]) cout << "Could not get " << name << endl;
        sprintf(name, "fDt_MCNP_PuFC_FG_%i", i+1);
        sprintf(name, "Simulation/Geant4/UFC_real/FitToF/fUFC_NIF_%i", i+1);
        f1Sim[2][i] = (TF1*) fAna->Get(name); if (!f1Sim[2][i]) cout << "Could not get " << name << endl;
        sprintf(name, "fDt_Geant4_UFC_FG_%i", i+1);
        sprintf(name, "Simulation/MCNP/UFC_real/FitToF/fUFC_NIF_%i", i+1);
        f1Sim[3][i] = (TF1*) fAna->Get(name); if (!f1Sim[3][i]) cout << "Could not get " << name << endl;
        sprintf(name, "fDt_MCNP_UFC_FG_%i", i+1);
        sprintf(name, "Simulation/Geant4/UFC_SB/FitToF/fUFC_SB_%i", i+1);
        f1Sim[4][i] = (TF1*) fAna->Get(name); if (!f1Sim[4][i]) cout << "Could not get " << name << endl;
        sprintf(name, "fDt_Geant4_UFC_SB_%i", i+1);

        if (Expand) {
            sprintf(name, "PuFC, FG, %i", i+1);
            cFitToF[0][i] = new TCanvas(name);
        } else
            cFitToF[0][0]->cd(i+1);
        gPad->SetLeftMargin(0.125);
        gPad->SetTicks(1, 1);
        gPad->SetBottomMargin(0.12);
        h1Exp[0][i]->SetTitle("; #font[12]{t} / ns; Ereignisse / ns");
        h1Exp[0][i]->SetStats(0);
        SetSize(h1Exp[0][i]);
        h1Exp[0][i]->GetYaxis()->SetTitleOffset(1.1);
        Double_t xmin, xmax;
        f1Sim[0][i]->GetRange(xmin, xmax);
        h1Exp[0][i]->GetXaxis()->SetRangeUser(xmin, xmax);
        h1Exp[0][i]->Draw("p");
        f1Sim[1][i]->SetLineColor(800);
        f1Sim[1][i]->Draw("same");
        f1Sim[0][i]->SetLineColor(416);
        f1Sim[0][i]->Draw("same");
        sprintf(name, "PuFC Deposit %i", i+1);
        TLegend *l = new TLegend(0.5, 0.5, 0.8, 0.8, name);
        sprintf(name, "Geant4 #chi^{2}/#nu = %.2f", f1Sim[0][i]->GetChisquare() / f1Sim[0][i]->GetNDF());
        l->AddEntry(f1Sim[0][i], name);
        sprintf(name, "MCNP #chi^{2}/#nu = %.2f", f1Sim[1][i]->GetChisquare() / f1Sim[1][i]->GetNDF());
        l->AddEntry(f1Sim[1][i], name);
        l->Draw("same");
        cout << i+1 << " Geant4: " << f1Sim[0][i]->GetParameter(2) << "\\pm" << f1Sim[0][i]->GetParError(2) << " " << f1Sim[0][i]->GetNDF() << endl;
//        cout <<       "Pu  MCNP: " << f1Sim[1][i]->GetParameter(2) << "\\pm" << f1Sim[1][i]->GetParError(2) << " " << f1Sim[1][i]->GetNDF() << endl;

        if (Expand) {
            sprintf(name, "/home/hoffma93/Pictures/PeakForm/PuFC_FG_%i.pdf", i+1);
            if (Save)
                cFitToF[0][i]->SaveAs(name);
            sprintf(name, "UFC, FG, %i", i+1);
            cFitToF[1][i] = new TCanvas(name);
        } else
            cFitToF[1][0]->cd(i+1);
        gPad->SetLeftMargin(0.125);
        gPad->SetTicks(1, 1);
        gPad->SetBottomMargin(0.12);
        h1Exp[1][i]->SetTitle("; #font[12]{t} / ns; Ereignisse / ns");
        h1Exp[1][i]->SetTitle("");
        h1Exp[1][i]->SetStats(0);
        SetSize(h1Exp[1][i]);
        h1Exp[1][i]->GetYaxis()->SetTitleOffset(1.1);
        f1Sim[2][i]->GetRange(xmin, xmax);
        h1Exp[1][i]->GetXaxis()->SetRangeUser(xmin, xmax);
        h1Exp[1][i]->Draw("p");
        f1Sim[3][i]->SetLineColor(800);
        f1Sim[3][i]->Draw("Csame");
        f1Sim[2][i]->SetLineColor(416);
        f1Sim[2][i]->Draw("Csame");
        sprintf(name, "UFC Deposit %i", i+1);
        l = new TLegend(0.5, 0.5, 0.8, 0.8, name);
        sprintf(name, "Geant4 #chi^{2}/#nu = %.2f", f1Sim[2][i]->GetChisquare() / f1Sim[2][i]->GetNDF());
        l->AddEntry(f1Sim[2][i], name);
        sprintf(name, "MCNP #chi^{2}/#nu = %.2f", f1Sim[3][i]->GetChisquare() / f1Sim[3][i]->GetNDF());
        l->AddEntry(f1Sim[3][i], name);
        l->Draw("same");
//        cout << i+1 << " Geant4: " << f1Sim[2][i]->GetParameter(2) << "\\pm" << f1Sim[2][i]->GetParError(2) << " " << f1Sim[2][i]->GetNDF() << endl;
//        cout <<       "U   MCNP: " << f1Sim[3][i]->GetParameter(2) << "\\pm" << f1Sim[3][i]->GetParError(2) << " " << f1Sim[3][i]->GetNDF() << endl;

        if (Expand) {
            sprintf(name, "/home/hoffma93/Pictures/PeakForm/UFC_FG_%i.pdf", i+1);
            if (Save)
                cFitToF[1][i]->SaveAs(name);
            sprintf(name, "UFC, SB, %i", i+1);
            cFitToF[2][i] = new TCanvas(name);
        } else
            cFitToF[2][0]->cd(i+1);
        gPad->SetLeftMargin(0.125);
        gPad->SetTicks(1, 1);
        gPad->SetBottomMargin(0.12);
        h1Exp[2][i]->SetTitle("; #font[12]{t} / ns; Ereignisse / ns");
        h1Exp[2][i]->SetTitle("");
        h1Exp[2][i]->SetStats(0);
        SetSize(h1Exp[2][i]);
        h1Exp[2][i]->GetYaxis()->SetTitleOffset(1.1);
        f1Sim[4][i]->GetRange(xmin, xmax);
        h1Exp[2][i]->GetXaxis()->SetRangeUser(xmin, xmax);
        h1Exp[2][i]->Draw("p");
        f1Sim[4][i]->SetLineColor(416);
        f1Sim[4][i]->Draw("Csame");
        sprintf(name, "UFC Deposit %i", i+1);
        l = new TLegend(0.5, 0.5, 0.8, 0.7, name);
        sprintf(name, "Geant4 #chi^{2}/#nu = %.2f", f1Sim[4][i]->GetChisquare() / f1Sim[4][i]->GetNDF());
        l->AddEntry(f1Sim[4][i], name);
        l->Draw("same");
        if (Expand)
        {
            sprintf(name, "/home/hoffma93/Pictures/PeakForm/UFC_BG_%i.pdf", i+1);
            if (Save)
                cFitToF[0][i]->SaveAs(name);
        }
    }
}

void DrawG4vsMCNP(TFile *fAna, string FC, Int_t ch, Bool_t Save = 0)
{
    string key = "real";
    char name[128] = "";
    sprintf(name, "Simulation/Geant4/%s_%s/EffToF/%s_ProjT_%s_%i", FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), ch+1);
    TH1D *hT_G4 = (TH1D*) fAna->Get(name); if (!hT_G4) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/MCNP/%s_%s/EffToF/%s_ProjT_%s_%i", FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), ch+1);
    TH1D *hT_MCNP = (TH1D*) fAna->Get(name); if (!hT_MCNP) cout << "Could not get " << name << endl;
    hT_G4->SetTitle("; #font[12]{t} / ns;  #font[12]{N}_{(n,f)} / #font[12]{Y}");
    hT_G4->SetStats(0);
    SetSize(hT_G4);
    hT_G4->GetXaxis()->SetRangeUser(25, 50);
    hT_G4->GetYaxis()->SetRangeUser(0, 5E-8);
    hT_G4->SetLineColor(416);
    hT_MCNP->SetLineColor(800);
    sprintf(name, "%s Deposit %i", FC.c_str(), ch+1);
    TLegend *l = new TLegend(0.6, 0.6, 0.85, 0.85, name);
    l->AddEntry(hT_G4, "Geant4");
    l->AddEntry(hT_MCNP, "MCNP");
    sprintf(name, "cT_%s_%i", FC.c_str(), ch+1);
    TCanvas *cT = new TCanvas(name, name, 1, 1, 1400, 1000);
    gPad->SetTicks(1,1);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    hT_G4->Draw("hist");
    hT_MCNP->Draw("same hist");
    l->Draw();
    sprintf(name, "~/Pictures/Simulation/G4vsMCNP/T_%s_%s_%i.pdf", FC.c_str(), key.c_str(), ch+1);
    if (Save)
        cT->SaveAs(name);

    // E proj. 35 < t < 45
    Double_t tMin = 35, tMax = 45;
    sprintf(name, "Simulation/Geant4/%s_%s/ToFvsEkin/%s_ToFvsEkin_%s_Ch.%i", FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), ch+1);
    TH2D *hTE_G4 = (TH2D*) fAna->Get(name); if (!hT_G4) cout << "Could not get " << name << endl;
    sprintf(name, "G4_%s_%s_%i", FC.c_str(), key.c_str(), ch+1);
    TH1D *hE_G4 = (TH1D*) hTE_G4->ProjectionY(name, hTE_G4->GetXaxis()->FindBin(tMin), hTE_G4->GetXaxis()->FindBin(tMax)-1);
    sprintf(name, "Simulation/MCNP/%s_%s/ToFvsEkin/%s_ToFvsEkin_%s_Ch.%i", FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), ch+1);
    TH2D *hTE_MCNP = (TH2D*) fAna->Get(name); if (!hT_MCNP) cout << "Could not get " << name << endl;
    sprintf(name, "MCNP_%s_%s_%i", FC.c_str(), key.c_str(), ch+1);
    TH1D *hE_MCNP = (TH1D*) hTE_MCNP->ProjectionY(name, hTE_MCNP->GetXaxis()->FindBin(tMin), hTE_MCNP->GetXaxis()->FindBin(tMax)-1);
    cout << hE_G4->Integral() / hE_MCNP->Integral() << endl;
    hE_G4->SetTitle("; #font[12]{E} / MeV;  #font[12]{N}_{(n,f)} / #font[12]{Y}");
    hE_G4->SetStats(0);
    SetSize(hE_G4);
    hE_G4->SetLineColor(416);
    hE_MCNP->SetLineColor(800);
    sprintf(name, "cE_%s_%i", FC.c_str(), ch+1);
    TCanvas *cE = new TCanvas(name, name, 1, 1, 1400, 1000);
    gPad->SetTicks(1,1);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    hE_G4->Draw("hist");
    hE_MCNP->Draw("same hist");
    l->Draw();
    sprintf(name, "~/Pictures/Simulation/G4vsMCNP/E_%s_%s_%i.pdf", FC.c_str(), key.c_str(), ch+1);
    if (Save)
        cE->SaveAs(name);
}

TGraphErrors *GetResult(TFile *fAna, string FC, string Simulation = "Geant4")
{
    char name[64] = "";
    sprintf(name, "%s/Correction/%s_CS_%s_corrected", FC.c_str(), FC.c_str(), Simulation.c_str());
    TGraphErrors *g = (TGraphErrors*)fAna->Get(name); if (!g) cout << "Could not get " << name << endl;
    return g;
}

void DrawResult(TFile *fAna)
{
    TGraphErrors *gePu = GetResult(fAna, "PuFC", "Geant4");
    TF1 *fPu = new TF1("fPu", "pol0", 0.5, 8.5);
    fPu->SetRange(0.5, 8.5);
    gePu->Fit("fPu", "R0Q");
    cout << "Pu: " << fPu->GetParameter(0) << "+-" << fPu->GetParError(0) << endl;
    gePu->SetTitle("; Deposit; #sigma_{(n,f)} / b");
    gePu->SetMarkerStyle(20);
    gePu->SetMarkerColor(kBlue);
    gePu->SetLineColor(kBlue);
    fPu->SetLineColor(kBlue);

    TGraphErrors *geU = GetResult(fAna, "UFC", "Geant4");
    TF1 *fU = new TF1("fU", "pol0", 0.5, 8.5);
    fU->SetRange(0.5, 8.5);
    geU->Fit("fU", "R0Q");
    cout << "U: " << fU->GetParameter(0) << "+-" << fU->GetParError(0) << endl;
    geU->SetMarkerStyle(21);
    geU->SetMarkerColor(kRed);
    geU->SetLineColor(kRed);
    fU->SetLineColor(kRed);

    TLegend *l = new TLegend(0.7, 0.7, 0.85, 0.85);
    l->AddEntry(gePu, "{}^{242}Pu", "PE");
    l->AddEntry(geU, "{}^{235}U", "PE");
    BiasX(gePu, -0.05);
    BiasX(geU, +0.05);

    new TCanvas();
//    gPad->SetLeftMargin(0.11);
//    gPad->SetBottomMargin(0.11);
    gPad->SetTicks(1,1);
    gePu->Draw("AP");
    SetSize(gePu);
    fPu->Draw("same");
    geU->Draw("same P");
    fU->Draw("same");
    l->Draw();
}

void DrawCalN(TFile *fAna)
{
    TGraphErrors *geUo = (TGraphErrors*)fAna->Get("UFC/nAtoms/UFC_eN"); if (!geUo) cout << "Could not get " << "UFC/nAtoms/UFC_eN" << endl;
    TGraphErrors *geUn = (TGraphErrors*)fAna->Get("UFC/nAtoms/UFC_eN_TL"); if (!geUn) cout << "Could not get " << "UFC/nAtoms/UFC_eN_TL" << endl;
    TGraphErrors *gePTB = (TGraphErrors*)fAna->Get("UFC/nAtoms/UFC_eN_cal_Geant4"); if (!gePTB) cout << "Could not get " << "UFC/nAtoms/UFC_eN_cal_Geant4" << endl;
    geUo->SetTitle("UFC atom number calibration; Deposit; #varepsilon#font[12]{N}_{U} / 10^{19}");
    BiasX(geUo, -0.1, 1.E-19);
    BiasX(geUn, 0.0, 1.E-19);
    BiasX(gePTB, +0.1, 1.E-19);
    geUo->SetLineWidth(2);
    geUo->SetMarkerStyle(20);
    geUo->SetMarkerSize(2);
    geUo->SetMarkerColor(kBlack);
    geUo->SetLineColor(kBlack);
    geUn->SetLineWidth(2);
    geUn->SetMarkerStyle(22);
    geUn->SetMarkerSize(2);
    geUn->SetMarkerColor(kBlue);
    geUn->SetLineColor(kBlue);
    gePTB->SetLineWidth(2);
    gePTB->SetMarkerStyle(21);
    gePTB->SetMarkerSize(2);
    gePTB->SetMarkerColor(kRed);
    gePTB->SetLineColor(kRed);
    TLegend *lU = new TLegend(0.15, 0.65, 0.4, 0.85, "Calibration");
    lU->AddEntry(geUo, "UFC vs H19 @ #font[12]{n}ELBE, 2016 result", "pe");
    lU->AddEntry(geUn, "UFC vs H19 @ #font[12]{n}ELBE, with Track Length", "pe");
    lU->AddEntry(gePTB, "UFC @ PTB", "pe");

    new TCanvas("UFC_calN");
    gPad->SetBottomMargin(0.12);
    gPad->SetTicks(1,1);
    geUo->Draw("AP");
    SetSize(geUo);
    geUo->GetXaxis()->SetTitleOffset(0.8);
    geUo->GetXaxis()->SetNdivisions(110);
    geUo->GetYaxis()->SetRangeUser(3.5, 5);
    geUn->Draw("same P");
    gePTB->Draw("same P");
    lU->Draw();

    TGraphErrors *gePu_o = (TGraphErrors*)fAna->Get("PuFC/nAtoms/PuFC_eN"); if (!gePu_o) cout << "Could not get " << "PuFC/nAtoms/PuFC_eN" << endl;
    TGraphErrors *gePu_n = (TGraphErrors*)fAna->Get("PuFC/nAtoms/PuFC_eN_cal_Geant4"); if (!gePu_n) cout << "Could not get " << "PuFC/nAtoms/PuFC_eN_cal_Geant4" << endl;
    gePu_o->SetTitle("PuFC atom number calibration; Deposit; #varepsilon#font[12]{N}_{Pu} / 10^{19}");
    BiasX(gePu_o, -0.1, 1.E-19);
    BiasX(gePu_n, +0.1, 1.E-19);
    gePu_o->SetLineWidth(2);
    gePu_o->SetMarkerStyle(22);
    gePu_o->SetMarkerSize(2);
    gePu_o->SetMarkerColor(kBlue);
    gePu_o->SetLineColor(kBlue);
    gePu_n->SetLineWidth(2);
    gePu_n->SetMarkerStyle(21);
    gePu_n->SetMarkerSize(2);
    gePu_n->SetMarkerColor(kRed);
    gePu_n->SetLineColor(kRed);
    TLegend *lPu = new TLegend(0.15, 0.65, 0.4, 0.85, "Calibration");
    lPu->AddEntry(gePu_o, "#font[12]{T}_{1/2,SF}", "PE");
    lPu->AddEntry(gePu_n, "PTB, LC1", "PE");

    new TCanvas("PuFC_calN");
    gPad->SetBottomMargin(0.12);
    gPad->SetTicks(1,1);
    gePu_o->Draw("AP");
    SetSize(gePu_o);
    gePu_o->GetXaxis()->SetTitleOffset(0.8);
    gePu_o->GetYaxis()->SetRangeUser(0.75, 1.75);
    gePu_n->Draw("same P");
    lPu->Draw();

    Double_t x, y, yerr;
    for (Int_t i = 0; i < 8; i++)
    {
        gePu_n->GetPoint(i, x, y);
        yerr = gePu_n->GetErrorY(i);
        cout << x+1 << "   " << y << " +- " << yerr << endl;
    }
}

void DrawConstantBackground(TFile *fAna, string Run = "NIF", Int_t ch = 0)
{
    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
//    cout << FC << endl;
    char name[64] = "";
    // Dt-Bg
    TF1 *fL;
    TF1 *fR;
    Double_t lUg;
    Double_t DlUg;
    Double_t rUg;
    Double_t DrUg;
    Int_t nL, nR;
    Double_t Ug, DUg;
    Double_t lChi2dof;
    Double_t rChi2dof;
    TH1F *h1Dt;
    TLatex *tL1;
    TLatex *tL2;
    TLatex *tR1;
    TLatex *tR2;
    TCanvas *cLR;

    sprintf(name, "%s/ToF/Background/%s/Left/%s_fL_%i", FC.c_str(), Run.c_str(), Run.c_str(), ch+1);
    fL = (TF1*)fAna->Get(name); if (!fL) cout << "Could not get " << name << endl;
    sprintf(name, "%s/ToF/Background/%s/Right/%s_fR_%i", FC.c_str(), Run.c_str(), Run.c_str(), ch+1);
    fR = (TF1*)fAna->Get(name); if (!fR) cout << "Could not get " << name << endl;

    lUg = fL->GetParameter(0);
    DlUg = fL->GetParError(0);
    lChi2dof = fL->GetChisquare() / fL->GetNDF();
    nL = fL->GetNDF() + 1;
    rUg = fR->GetParameter(0);
    DrUg = fR->GetParError(0);
    rChi2dof = fR->GetChisquare() / fR->GetNDF();
    nR = fR->GetNDF() + 1;
    Ug = (nL * lUg + nR * rUg) / (nL + nR);
    DUg = sqrt(pow(nL * DlUg, 2) + pow(nR * DrUg, 2)) / (nL + nR);

    // Format: UFC 3 & 3.60(13) & 1.63 & 3.71(18) & 1.89 & 3.26(13) \\
//    sprintf(name, "HALLO!");
    sprintf(name, "\\num{%.2f(%.f)} & %.2f & \\num{%.2f(%.f)} & %.2f & \\num{%.2f(%.f)}", lUg, 100*DlUg, lChi2dof, rUg, 100*DrUg, rChi2dof, Ug, 100*DUg);
    cout << FC << " & " << ch+1 << " & " << name << " \\\\" << endl;
//    cout << name << endl;

    sprintf(name, "%s/ToF/Total/%s/H1AnaHZDRDtG_%i", FC.c_str(), Run.c_str(), ch+1);
    h1Dt = (TH1F*) fAna->Get(name); if (!h1Dt) cout << "Could not get " << name << endl;
    sprintf(name, "Dt_%s_%i", Run.c_str(), ch+1);
    h1Dt->SetName(name);
    h1Dt->SetTitle("; #Delta#font[12]{t} / ns; Counts");
    h1Dt->SetStats(0);
    sprintf(name, "cLR_%s_%i", Run.c_str(), ch+1);
    cLR = new TCanvas(name, "Konstanter Untergrund");
    gPad->SetTicks(1, 1);
    h1Dt->Draw();
    fL->SetLineColor(kRed);
    fL->Draw("same");
    fR->SetLineColor(kGreen);
    fR->Draw("same");

    tL1 = new TLatex();
    sprintf(name, "%.1f #pm %.1f", lUg, DlUg);
    tL1->DrawLatexNDC(0.2, 0.5, name);

    tL2 = new TLatex();
    sprintf(name, "#chi^{2}/#nu = %.2f", lChi2dof);
    tL2->DrawLatexNDC(0.2, 0.4, name);

    tR1 = new TLatex();
    sprintf(name, "%.1f #pm %.1f", rUg, DrUg);
    tR1->DrawLatexNDC(0.6, 0.5, name);

    tR2 = new TLatex();
    sprintf(name, "#chi^{2}/#nu = %.2f", rChi2dof);
    tR2->DrawLatexNDC(0.6, 0.4, name);
}

void DrawEvaluationPu()
{
    TGraphErrors *geB = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_B-VIII.dat", "%lg %lg %lg");
    TGraph *geC = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_CENDL.dat");
    TGraph *geF = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_JEFF.dat");
    TGraph *geJ = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_JENDL.dat");
    TGraph *geR = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_RUSFOND.dat");
    geB->SetLineColor(kBlack);
    geC->SetLineColor(kRed);
    geF->SetLineColor(kGreen);
    geJ->SetLineColor(kBlack);
    geR->SetLineColor(kBlue);
    geB->SetLineWidth(lw);
    geC->SetLineWidth(lw);
    geF->SetLineWidth(lw);
    geJ->SetLineWidth(lw);
    geR->SetLineWidth(lw);
    geB->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");
    TLegend *l = new TLegend(0.45, 0.2, 0.9, 0.5);
    l->SetTextFont(132);
//    l->AddEntry(geB, "ENDF/B-VIII.0");
    l->AddEntry(geC, "CENDL-3.1");
    l->AddEntry(geF, "JEFF-3.3");
    l->AddEntry(geB, "JENDL-4.0 / ENDF/B-VIII.0");
    l->AddEntry(geR, "RUSFOND-2010");
    TCanvas *c = new TCanvas();
    gPad->SetLogx(1);
    geB->SetFillColorAlpha(kGray, 0.5);
    geB->Draw("ALX");
    geB->Draw("3");
    geB->GetXaxis()->SetRangeUser(0.3, 20);
    SetSize(geB);
    geC->Draw("same");
    geF->Draw("same");
//    geJ->Draw("same");
    geR->Draw("same");
    l->Draw();
}

void DrawEvaluationU()
{
//    TGraph *geB = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_standard.dat");
//    TGraph *geB = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_B-VIII.dat");
    TGraph *geC = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_CENDL.dat");
    TGraph *geF = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_JEFF.dat");
    TGraph *geJ = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_JENDL.dat");
    TGraph *geR = new TGraph("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_BROND.dat");

//    geB->SetLineColor(kGray);
    geC->SetLineColor(kRed);
    geF->SetLineColor(kGreen);
    geJ->SetLineColor(kBlack);
    geR->SetLineColor(kBlue);
//    geB->SetLineWidth(lw);
    geC->SetLineWidth(lw);
    geF->SetLineWidth(lw);
    geJ->SetLineWidth(lw);
    geR->SetLineWidth(lw);
    geC->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");
    TLegend *l = new TLegend(0.45, 0.2, 0.9, 0.5);
    l->SetTextFont(132);
//    l->AddEntry(geB, "");
    l->AddEntry(geC, "CENDL-3.1");
    l->AddEntry(geF, "JEFF-3.3 / ENDF/B-VIII.0");
    l->AddEntry(geJ, "JENDL-4.0");
    l->AddEntry(geR, "BROND-2016");
    TCanvas *c = new TCanvas();
    gPad->SetLogx(1);
    geC->Draw();
    geC->GetXaxis()->SetRangeUser(0.3, 20);
    geC->GetYaxis()->SetRangeUser(0.0, 2.3);
    SetSize(geC);
//    geB->Draw("same");
    geF->Draw("same");
    geJ->Draw("same");
    geR->Draw("same");
    l->Draw();
}

void DrawExpData()
{
    TGraphErrors *ge[10];
    TLegend *l1 = new TLegend(0.45, 0.2, 0.9, 0.5);
    l1->SetTextFont(132);
    TLegend *l2 = new TLegend(0.45, 0.2, 0.9, 0.5);
    l2->SetTextFont(132);

    /// Data normed to U235

    ge[0] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Koegler.dat", "%lg %lg %*lg %lg %lg %*lg");
    ge[0]->SetLineColor(kBlue);
    ge[0]->SetMarkerColor(kBlue);
    ge[0]->SetMarkerStyle(20);
    l1->AddEntry(ge[0], "K#ddot{o}gler et al., 2016", "PE");

    ge[1] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Staples.dat", "%lg %lg %lg");
    ge[1]->SetLineColor(kMagenta+2);
    ge[1]->SetMarkerColor(kMagenta+2);
    ge[1]->SetMarkerStyle(24);
    l1->AddEntry(ge[1], "Staples et al., 1998", "PE");
    ge[1]->SetTitle("; #font[12]{E}_{n} / MeV; #sigma_{(n,f)} / #sigma_{Ref}");

    TCanvas *c1 = new TCanvas("c1");
    gPad->SetLogx(1);
    ge[1]->Draw("AP");
//    ge[1]->GetXaxis()->SetRangeUser(0.3, 20);
    SetSize(ge[1]);
    ge[0]->Draw("same P");
    l1->Draw();

    /// absolute data

    ge[5] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Koegler.dat", "%lg %*lg %lg %lg %*lg %lg");
    ge[5]->SetLineColor(kBlue);
    ge[5]->SetMarkerColor(kBlue);
    ge[5]->SetMarkerStyle(20);
    l2->AddEntry(ge[5], "K#ddot{o}gler et al., 2016", "PE");

    ge[6] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_S-C.dat", "%lg %lg %lg %lg");
    ge[6]->SetLineColor(kOrange);
    ge[6]->SetMarkerColor(kOrange);
    ge[6]->SetMarkerStyle(28);
    l2->AddEntry(ge[6], "Salvador-Casi#tilde{n}eira et al., 2016", "PE");

    ge[7] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Tovesson.dat", "%lg %lg %lg");
    ge[7]->SetLineColor(kRed);
    ge[7]->SetMarkerColor(kRed);
    ge[7]->SetMarkerStyle(25);
    l2->AddEntry(ge[7], "Tovesson et al., 2009", "PE");
    ge[7]->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");

    ge[8] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Weigmann.dat", "%lg %lg %lg");
    ge[8]->SetLineColor(kGreen);
    ge[8]->SetMarkerColor(kGreen);
    ge[8]->SetMarkerStyle(22);
    l2->AddEntry(ge[8], "Weigmann et al., 1984", "PE");

    TCanvas *c2 = new TCanvas("c2");
    gPad->SetLogx(1);
    ge[7]->Draw("AP");
//    ge[7]->GetXaxis()->SetRangeUser(0.3, 20);
    SetSize(ge[7]);
    ge[5]->Draw("same P");
    ge[6]->Draw("same P");
    ge[8]->Draw("same P");
    l2->Draw();
}

void DrawExpAbs()
{
    TGraphErrors *ge[10];
    TLegend *l = new TLegend(0.16, 0.6, 0.5, 0.94);
    l->SetTextFont(132);

    ge[0] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Tovesson.dat", "%lg %lg %lg");
    ge[0]->SetLineColor(kOrange);
    ge[0]->SetMarkerColor(kOrange);
    ge[0]->SetMarkerStyle(24);
    ge[0]->SetMarkerSize(2);
    l->AddEntry(ge[0], "Tovesson et al., 2009", "PE");
    ge[0]->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");

    ge[1] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Manabe.dat", "%lg %lg %*lg %lg %lg");
    ge[1]->SetLineColor(kBlue);
    ge[1]->SetMarkerColor(kBlue);
    ge[1]->SetMarkerStyle(27);
    ge[1]->SetMarkerSize(3);
    l->AddEntry(ge[1], "Manabe et al., 1988", "PE");
    ge[1]->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");

    ge[2] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Maedows.dat", "%lg %lg %lg %lg");
    ge[2]->SetLineColor(kBlack);
    ge[2]->SetMarkerColor(kBlack);
    ge[2]->SetMarkerStyle(22);
    ge[2]->SetMarkerSize(3);
    l->AddEntry(ge[2], "Maedows et al., 1988", "PE");

    ge[3] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Gul.dat", "%lg %lg %lg");
    ge[3]->SetLineColor(8);
    ge[3]->SetMarkerColor(8);
    ge[3]->SetMarkerStyle(23);
    ge[3]->SetMarkerSize(3);
    l->AddEntry(ge[3], "Gul et al., 1986", "PE");

    ge[4] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_Alkhazov.dat", "%lg %lg %lg");
    ge[4]->SetLineColor(6);
    ge[4]->SetMarkerColor(6);
    ge[4]->SetMarkerStyle(29);
    ge[4]->SetMarkerSize(3);
    l->AddEntry(ge[4], "Alkhazov et al., 1983", "PE");

    ge[5] = new TGraphErrors(1);
    ge[5]->SetPoint(0, 14.97, 2.21);
    ge[5]->SetPointError(0, 0.3, 0.08);
    ge[5]->SetLineColor(kGreen);
    ge[5]->SetMarkerColor(kGreen);
    ge[5]->SetMarkerStyle(20);
    ge[5]->SetMarkerSize(3);
    l->AddEntry(ge[5], "Diese Arbeit", "PE");

    ge[6] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/Pu242_nf_B-VIII.dat", "%lg %lg %lg");
    ge[6]->SetLineColor(kGray);
//    ge[6]->SetMarkerColor(6);
//    ge[6]->SetMarkerStyle(29);
//    ge[6]->SetMarkerSize(3);
    l->AddEntry(ge[6], "ENDF/B-VIII.0");

    TCanvas *c2 = new TCanvas("c2");
//    gPad->SetLogx(1);
    ge[0]->Draw("AP");
    ge[0]->GetXaxis()->SetRangeUser(13.5, 15.5);
    ge[0]->GetYaxis()->SetRangeUser(1.9, 2.3);
//    ge[0]->GetYaxis()->SetNdivisions(206);
    SetSize(ge[0]);
    ge[1]->Draw("same P");
    ge[2]->Draw("same P");
    ge[3]->Draw("same P");
    ge[4]->Draw("same P");
    ge[5]->Draw("same P");
    ge[6]->SetFillColorAlpha(kGray, 0.3);
    ge[6]->Draw("same LX");
    ge[6]->Draw("same 3");
    l->Draw();
}

void DrawExpU()
{
    TGraphErrors *ge[10];
    TLegend *l = new TLegend(0.6, 0.2, 0.94, 0.55);
    l->SetTextFont(132);

    ge[0] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_Carlson.dat", "%lg %lg %lg");
    ge[0]->SetLineColor(kRed);
    ge[0]->SetMarkerColor(kRed);
    ge[0]->SetMarkerStyle(28);
    ge[0]->SetMarkerSize(2);
    l->AddEntry(ge[0], "Carlson et al., 1991", "PE");
    ge[0]->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");

    ge[1] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_Jingwen88.dat", "%lg %lg %lg");
    ge[1]->SetLineColor(4);
    ge[1]->SetMarkerColor(4);
    ge[1]->SetMarkerStyle(22);
    ge[1]->SetMarkerSize(3);
    l->AddEntry(ge[1], "Jingwen et al., 1988", "PE");
    ge[1]->SetTitle("; #font[12]{E}_{n} / MeV; #sigma / b");

    ge[2] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_Jingwen83.dat", "%lg %lg %lg");
    ge[2]->SetLineColor(9);
    ge[2]->SetMarkerColor(9);
    ge[2]->SetMarkerStyle(23);
    ge[2]->SetMarkerSize(3);
    l->AddEntry(ge[2], "Jingwen et al., 1983", "PE");

    ge[3] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_Mahdavi.dat", "%lg %lg %lg %lg");
    ge[3]->SetLineColor(kOrange);
    ge[3]->SetMarkerColor(kOrange);
    ge[3]->SetMarkerStyle(33);
    ge[3]->SetMarkerSize(3);
    l->AddEntry(ge[3], "Mahdavi, 1982", "PE");

    ge[4] = new TGraphErrors(1);
    ge[4]->SetPoint(0, 14.99, 1.8);
    ge[4]->SetPointError(0, 0.26, 0.04);
    ge[4]->SetLineColor(8);
    ge[4]->SetMarkerColor(8);
    ge[4]->SetMarkerStyle(24);
    ge[4]->SetMarkerSize(3);
    l->AddEntry(ge[4], "#font[12]{n}ELBE", "PE");

    ge[5] = new TGraphErrors(1);
    ge[5]->SetPoint(0, 14.97, 2.20);
    ge[5]->SetPointError(0, 0.3, 0.017);
    ge[5]->SetLineColor(3);
    ge[5]->SetMarkerColor(3);
    ge[5]->SetMarkerStyle(20);
    ge[5]->SetMarkerSize(3);
    l->AddEntry(ge[5], "PTB", "PE");

    ge[6] = new TGraphErrors("/gpfs/home/hoffma93/Programme/ROOT/Data/U235_nf_B-VIII.dat", "%lg %lg %lg");
    ge[6]->SetLineColor(kBlack);
//    ge[6]->SetMarkerColor(6);
//    ge[6]->SetMarkerStyle(29);
//    ge[6]->SetMarkerSize(3);
    l->AddEntry(ge[6], "ENDF/B-VIII.0");

    TCanvas *c2 = new TCanvas("c2");
//    gPad->SetLogx(1);
    ge[0]->Draw("AP");
    ge[0]->GetXaxis()->SetRangeUser(13.0, 17.0);
    ge[0]->GetYaxis()->SetRangeUser(1.6, 2.4);
    ge[0]->GetXaxis()->SetNdivisions(205);
    SetSize(ge[0]);
    ge[1]->Draw("same P");
    ge[2]->Draw("same P");
    ge[3]->Draw("same P");
    ge[4]->Draw("same P");
    ge[5]->Draw("same P");
//    ge[6]->SetFillColorAlpha(kGray, 0.3);
    ge[6]->Draw("same LX");
//    ge[6]->Draw("same 3");
    l->Draw();
}

int DrawPics()
{
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle(kTRUE);
    gStyle->SetLegendFont(132);
//    gStyle->SetCanvasPreferGL(kTRUE);
//    gStyle->SetPalette(kRainBow);
//    gStyle->SetPalette(kSunset);
//    TColor::InvertPalette();

    TFile* fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");
    TFile* fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    TFile* fUNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");

//    DrawSigma(242);
//    DrawCrossSectionRuns(fAna);
//    DrawUscattering(fAna);
//    DrawPuCorrection(fAna);
//    DrawMonitorRate();
//    DrawPeakWidth(fAna, "UFC");
//    DrawPeakWidth(fAna, "UFC");
//    DrawSimPeak(fAna, "UFC_NIF", "real", 1, "Geant4");
//        DrawSimPeak(fAna, "NIF", "real", i, "Geant4");
//    DrawMonitorRatio(fAna, 1);
//    DrawMonitorRatio(fAna, "UFC");
//    DrawPuSponFis();
//    DrawTimeDiff(fNIF, 0);
//    DrawExfor();
//    DrawSourceE();
//    DrawEnELBE();
//    DrawTargetE(fAna, "", 100000000);
//    DrawTargetToF(fAna, "");
//    DrawTargetAng(fAna);
//    DrawTargetAir(fAna);
//    DrawGEF(235, "nELBE");
//    DrawGEF("Mass");
//    DrawGEF("TKE", "#font[12]{TKE} / MeV");
//    DrawFFRange(fAna);
//    DrawGEFTree("15MeV", 235);
//    DrawA2(fAna);
//    DrawEta(fAna);
//    DrawInefficiency(fAna);
//    DrawQDC(fAna, "UFC", 0, fUNIF);
//        DrawDtInt(fAna, "UFC", 0, 15, 40, 0, 0);
//        DrawDtInt(fAna, "PuFC", 0, 15, 35, 0, 0);
//    DrawDtChange(fAna, "PuFC");
//    for (Int_t i = 0; i < 8; i++)
//        DrawDtGate(fAna, "UFC", 7, fUNIF);
//    DrawDtGate(fAna, "PuFC", 7, fNIF);
//    DrawDtBackground(fAna, "NIF", 0, fNIF);
//        DrawIndFisStability(fAna, "UFC", i);
//    DrawIndFisStability(fAna, "PuFC");
//    DrawIndFisStability(fAna, "UFC");
//    IndFisSum(fAna, "UFC");
//    DrawSponFisStability(fAna, "UFC");
//    DrawSponFisStability(fAna, "PuFC");
//    DrawSignalStability();
//    DrawQDCvsTime(0);
//    DrawTvsE(fAna);
//        DrawTvsE(fAna, "Geant4", "PuFC", "real", i, 0);
//        DrawTvsE(fAna, "Geant4", "UFC", "real", i, 0);
//        DrawTvsE(fAna, "MCNP", "PuFC", "real", i, 0);
//        DrawTvsE(fAna, "MCNP", "UFC", "real", i, 0);
//        DrawTvsE(fAna, "Geant4", "UFC", "SB", i, 1);
//        DrawT(fAna, "Geant4", "UFC", 0, 1);
//    DrawSimNotebook(fAna, 1, 1);
//    DrawG4vsMCNP(fAna, "PuFC", 0, 0);
//    DrawResult(fAna);
    DrawCalN(fAna);
//    DrawConstantBackground(fAna, "UFC_NIF", 0);
//    DrawConstantBackground(fAna, "UFC_SB", 0);
//    DrawConstantBackground(fAna, "NIF", i);
//    DrawEvaluationPu();
//    DrawEvaluationU();
//    DrawExpData();
//    DrawExpAbs();
//    DrawExpU();

    return 1;
}
#endif
