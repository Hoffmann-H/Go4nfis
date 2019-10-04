#include "TMathBase.h"
#include "FC.C"
#include "Runs.C"
#include "MCNPtoROOT.C"
#define lw 2
//#define a 0.25 // alpha
using namespace std;

void DrawSigma()
{
    TGraphErrors *pSigma = new TGraphErrors("/home/hoffma93/Programme/ROOT/Data/Pu242.dat", "%lg %lg %lg");

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

void BiasX(TGraphErrors* pG, Double_t dx)
{
    Double_t x = 0, y = 0;
    for (int i = 0; i < 8; i++)
    {
        pG->GetPoint(i, x, y);
        pG->SetPoint(i, x + dx, y);
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
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/CrossSection/%s_CS_raw_%i", FC.c_str(), FC.c_str(), i+1);
        geCS[i] = (TGraphErrors*)f->Get(name);
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
    TGraphErrors *geSim = (TGraphErrors*)f->Get("Simulation/Geant4/UFC_Open/Correction/S");
    TGraphErrors *geSimSB = (TGraphErrors*)f->Get("UFC/Correction/SimS");
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

void DrawMonitorRate()
{
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "READ");
    TGraphErrors *gU = new TGraphErrors("~/Programme/Go4nfis/FC-Analysis/results/NeutronField_UFC.txt", "%lg %*s %*lg %*lg %lg %lg");
    TGraphErrors *gPu = new TGraphErrors("~/Programme/Go4nfis/FC-Analysis/results/NeutronField_PuFC.txt", "%lg %*s %*lg %*lg %lg %lg");
    Double_t x, y, yerr, tStart, tLength;
    for (Int_t j = 0; j < 7; j++)
    {
        gU->GetPoint(j, x, y);
        yerr = gU->GetErrorY(j);
        string Run = GetRunName("UFC", j);
        tStart = GetStartTime(Run);
        tLength = GetRealTime(Run);
        gU->SetPoint(j, tStart + 0.5*tLength, y);
        gU->SetPointError(j, 0.5*tLength, yerr);

        gPu->GetPoint(j, x, y);
        yerr = gPu->GetErrorY(j);
        Run = GetRunName("PuFC", j);
        tStart = GetStartTime(Run);
        tLength = GetRealTime(Run);
        gPu->SetPoint(j, tStart + 0.5*tLength, y);
        gPu->SetPointError(j, 0.5*tLength, yerr);
    }
    gU->SetLineWidth(2);
    gPu->SetLineWidth(2);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gU);
    mg->Add(gPu);
    mg->GetXaxis()->SetTimeDisplay(1);
    mg->GetXaxis()->SetTimeFormat("%d.%m. %H:%M %F1970-01-01");
    mg->GetXaxis()->SetTitle("#font[12]{t}");
    mg->GetYaxis()->SetTitle("Monitor-Rate [1/s]");
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetTitleSize(0.07);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleSize(0.07);
    mg->GetYaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetRangeUser(0, 1800);
    TLegend *l = new TLegend(0.4, 0.2, 0.6, 0.3);
//    l->AddEntry(gU, "KW 21");
    l->AddEntry(gPu, "KW22");
    TCanvas *c1 = new TCanvas("name", "title", 1, 1, 1200, 500);
    gPad->SetTicks(1, 1);
    mg->Draw("AP");
    l->Draw();
    c1->Update();
}

void DrawPeakWidth(TFile *f, string FC = "PuFC")
{
    char name[64] = "";
    sprintf(name, "%s/ToF/Gate/Peak/gPeak_av", FC.c_str());
    TGraphErrors *geAv = (TGraphErrors*)f->Get(name);
    if (!geAv)
        cout << "Could not open " << name << endl;
    TGraphErrors *gePeak[8];
    TMultiGraph *mg = new TMultiGraph();
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.5, 0.15, 0.85, 0.4, name);
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/ToF/Gate/Peak/gPeak_%i", FC.c_str(), i+1);
        gePeak[i] = (TGraphErrors*)f->Get(name);
        BiasX(gePeak[i], -0.7 + 0.2*i);
        gePeak[i]->SetLineColor(i+2);
        gePeak[i]->SetMarkerColor(i+2);
        gePeak[i]->SetMarkerStyle(20);
        gePeak[i]->SetMarkerSize(2);
        mg->Add(gePeak[i]);
        sprintf(name, "%i", i+1);
        l->AddEntry(gePeak[i], name);
    }
    geAv->SetLineColor(1);
    geAv->SetMarkerColor(1);
    geAv->SetMarkerStyle(20);
    geAv->SetMarkerSize(2);
    mg->Add(geAv);
    l->AddEntry(geAv, "Mittel");
    mg->GetXaxis()->SetRangeUser(0, 120);
    mg->SetTitle("Integrationsgrenzen");
    mg->GetXaxis()->SetTitle("Gate [ns]");
    mg->GetYaxis()->SetTitle("Peak");
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

void DrawSimPeak(TFile *f, string FC, string key, Int_t ch, string SimulationPath)
{
    char name[64] = "";
    if (!strcmp(FC.c_str(), "PuFC")) // if PuFC
        sprintf(name, "PuFC/ToF/Signal/NIF/H1AnaHZDRDtG_%i", ch+1);
    else
        sprintf(name, "UFC/ToF/Signal/UFC_NIF/H1AnaHZDRDtG_%i", ch+1);
    TH1F *pH1Exp = (TH1F*)f->Get(name);
    if (!pH1Exp)
        cout << "Could not open " << name << endl;
    sprintf(name, "%s/SimToF/%s%s_FoldT_%s%i", SimulationPath.c_str(), PathTag(key).c_str(), FC.c_str(), NameTag(key).c_str(), ch+1);
    TH1D *pH1Sim = (TH1D*)f->Get(name);
    if (!pH1Sim)
        cout << "Could not open " << name << endl;
    Double_t max_pos = PeakCenter(ch, FC);
    Double_t Gate[] = {max_pos - 15.0, max_pos + 15.0};
    Double_t IntExp = pH1Exp->Integral(pH1Exp->GetXaxis()->FindBin(Gate[0]), pH1Exp->GetXaxis()->FindBin(Gate[1]));
    Double_t IntSim = pH1Sim->Integral(pH1Sim->GetXaxis()->FindBin(Gate[0]), pH1Sim->GetXaxis()->FindBin(Gate[1]));
    Double_t IntSimTotal = pH1Sim->Integral();
    cout << Gate[0] << " " << Gate[1] << " " << IntExp << " " << IntSim << " " << IntSimTotal << endl;
    sprintf(name, "%s Deposit %i", FC.c_str(), ch+1);
    TLegend *l = new TLegend(0.5, 0.15, 0.85, 0.4, name);
    pH1Exp->SetLineWidth(2);
    pH1Exp->SetStats(0);
    pH1Exp->GetXaxis()->SetLabelSize(0.06);
    pH1Exp->GetXaxis()->SetTitleSize(0.08);
    pH1Exp->GetXaxis()->SetTitleOffset(0.8);
    pH1Exp->GetYaxis()->SetLabelSize(0.06);
    pH1Exp->GetYaxis()->SetTitleSize(0.08);
    pH1Exp->GetYaxis()->SetTitleOffset(0.8);
    l->AddEntry(pH1Exp, "Messung");

    Double_t Scale = IntExp / IntSim * pH1Exp->GetXaxis()->GetBinWidth(1) / pH1Sim->GetXaxis()->GetBinWidth(1);
    pH1Sim->Scale(Scale);
    pH1Sim->SetLineColor(kRed);
    pH1Sim->SetLineWidth(2);
    l->AddEntry(pH1Sim, "Simulation");

//    Double_t min = pH1Exp->GetMinimum();
    Double_t max = pH1Exp->GetMaximum();
    TLine *ll = new TLine(Gate[0], 0, Gate[0], max);
    TLine *lr = new TLine(Gate[1], 0, Gate[1], max);
    Double_t Efficiency = IntSim / IntSimTotal;
    sprintf(name, "%f im Zeitfenster", Efficiency);
    TLatex *t = new TLatex(max_pos, 0.5*max, name);
    t->SetTextColor(kRed);

    TCanvas *c = new TCanvas("name", "title", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pH1Exp->Draw("hist");
    pH1Sim->Draw("same hist");
    l->Draw();
    ll->Draw();
    lr->Draw();
    t->Draw();
    c->Update();
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
    TGraphErrors *geNIF = (TGraphErrors*)fNIF->Get("Analysis/results/effSponFis");
    TGraphErrors *geSB = (TGraphErrors*)fSB->Get("Analysis/results/effSponFis");
    TGraphErrors *geSF = new TGraphErrors(8);
    TGraphErrors *geToni = new TGraphErrors(8);

    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
    Double_t DrSF_Toni[] = {0.0013, 0.0012, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0013};
    TH1D *pHt = (TH1D*)fSF->Get("Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t t_live = pHt->Integral();
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I *pHDt = (TH1I*)fSF->Get(name);
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

void DrawTargetE(TFile *f, string Subfolder, Int_t Population)
{
    char name[64] = "";
    sprintf(name, "Simulation/Target/%sDirect", Subfolder.c_str());
    TH1D *pDir = (TH1D*)f->Get(name);
    if (!pDir)
        cout << "Could not open " << name << endl;
    sprintf(name, "Simulation/Target/%sScattered", Subfolder.c_str());
    TH1D *pSc = (TH1D*)f->Get(name);
    sprintf(name, "Simulation/Target/%sTotal", Subfolder.c_str());
    TH1D *pTot = (TH1D*)f->Get(name);
    pDir->SetLineColor(kRed);
    pSc->SetLineColor(kBlue);
    pTot->SetLineColor(0);
    sprintf(name, "%i simulierte Ionen", Population);
    TLegend *l = new TLegend(0.4, 0.2, 0.6, 0.3, name);
    sprintf(name, "%.2f%% direkte Neutronen", 100 * pDir->Integral());
    l->AddEntry(pDir, name);
    sprintf(name, "%.2f%% indirekte Neutronen", 100 * pSc->Integral());
    l->AddEntry(pSc, name);
    pTot->GetXaxis()->SetTitleSize(0.08);
    pTot->GetXaxis()->SetLabelSize(0.06);
    pTot->GetXaxis()->SetTitleOffset(0.8);
    pTot->GetYaxis()->SetTitleSize(0.08);
    pTot->GetYaxis()->SetLabelSize(0.06);
    pTot->GetYaxis()->SetTitleOffset(0.7);
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
    sprintf(name, "Simulation/Target/%sToF_Target", Subfolder.c_str());
    TH1D *pTarget = (TH1D*)f->Get(name);
    if (!pTarget)
        cout << "Could not open " << name << endl;
    sprintf(name, "Simulation/Target/%sToF_E", Subfolder.c_str());
    TH1D *pE = (TH1D*)f->Get(name);

    pE->SetLineColor(kBlue);
    pE->SetStats(0);
    pTarget->SetLineColor(kRed);
    pE->GetXaxis()->SetLabelSize(0.06);
    pE->GetXaxis()->SetTitleSize(0.07);
    pE->GetXaxis()->SetTitleOffset(0.8);
    pE->GetYaxis()->SetLabelSize(0.06);
    pE->GetYaxis()->SetTitleSize(0.07);
    pE->GetYaxis()->SetTitleOffset(0.7);

    TLegend *l = new TLegend(0.6, 0.7, 0.85, 0.85);
    l->AddEntry(pTarget, "TARGET-Ausgabe");
    l->AddEntry(pE, "ToF(E)");
    TCanvas *cT = new TCanvas("cT", "TARGET simulated ToF", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPad->SetLogy(1);
    pE->Draw("hist");
    pTarget->Draw("same hist");
    l->Draw();
    cT->Update();
}

int DrawPics()
{
//    gStyle->SetCanvasPreferGL();

    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");
//    TFile* fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");


//    DrawSigma();
//    DrawCrossSectionRuns(fAna);
//    DrawUscattering(fAna);
//    DrawPuCorrection(fAna);
//    DrawMonitorRate();
//    DrawPeakWidth(fAna, "UFC");
//    DrawSimPeak(fAna, "UFC", "2", 0, "Simulation/Geant4/UFC_Open");
//    DrawSimPeak(fAna, "PuFC", "2", 0, "Simulation/MCNP");
//    DrawMonitorRatio(fAna, 1);
//    DrawMonitorRatio(fAna, "UFC");
//    DrawPuSponFis();
//    DrawTimeDiff(fNIF, 0);
//    DrawExfor();
//    DrawSourceE();
//    DrawTargetE(fAna, "", 100000000);
    DrawTargetToF(fAna, "");

    return 1;
}
