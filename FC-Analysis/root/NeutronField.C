#ifndef NEUTRON_FIELD_H
#define NEUTRON_FIELD_H
#include "SaveToFile.C"
#include "Runs.C"
#include "FC.C"
#include "/home/hoffma93/StyleSheets/StyleSheet.C"
#include "TGraphErrors.h"
#include <fstream>

#define neutron_field_data "/home/hoffma93/Programme/Go4nfis/FC-Analysis/data/runs.txt"

void Freifeld(Bool_t Draw = 1)
{
    Int_t nFreifeld = 11;
    Double_t Q[] = {9.502182, 5.283354, 1.756110, 1.747356, 46.947156, 57.481554, 53.542836, 96.288870, 6.377172, 50.632404, 8.310576};
    Double_t PLC[] = {1205191, 670297, 222819, 221708, 5837042, 7104705, 6559310, 11616244, 767556, 6071205, 1008025};
    Double_t He3[] = {1661687, 924619, 307917, 306382, 8109248, 9959170, 9238458, 16413936, 1089187, 8592403, 1427724};
    Double_t GM[] = {57278, 31877, 10750, 10697, 283758, 349482, 324393, 575540, 37931, 301346, 49700};
    Double_t NM[] = {5671053, 3154900, 1051251, 1046011, 27492079, 33478370, 30916955, 54792679, 3623069, 28646614, 4757385};
    Double_t sQ = 0, sHe3 = 0, sGM = 0, sNM = 0;
    for (Int_t j = 0; j < nFreifeld; j++)
    {
        sQ += Q[j] / PLC[j];
        sHe3 += He3[j] / PLC[j];
        sGM += GM[j] / PLC[j];
        sNM += NM[j] / PLC[j];
    }
    TGraphErrors* gQ = new TGraphErrors(nFreifeld);
    TGraphErrors* gHe3 = new TGraphErrors(nFreifeld);
    TGraphErrors* gGM = new TGraphErrors(nFreifeld);
    TGraphErrors* gNM = new TGraphErrors(nFreifeld);
    for (Int_t j = 0; j < nFreifeld; j++)
    {
        gQ->SetPoint(j, j + 1, Q[j] / PLC[j] * nFreifeld / sQ);
        gHe3->SetPoint(j, j + 1, He3[j] / PLC[j] * nFreifeld / sHe3);
        gHe3->SetPointError(j, 0, sqrt(He3[j]) / PLC[j] * nFreifeld / sHe3);
        gGM->SetPoint(j, j + 1, GM[j] / PLC[j] * nFreifeld / sGM);
        gGM->SetPointError(j, 0, sqrt(GM[j]) / PLC[j] * nFreifeld / sGM);
        gNM->SetPoint(j, j + 1, NM[j] / PLC[j] * nFreifeld / sNM);
        gNM->SetPointError(j, 0, sqrt(NM[j]) / PLC[j] * nFreifeld / sNM);
        cout << j << " " << Q[j] / PLC[j] * nFreifeld / sQ << " " << He3[j] / PLC[j] * nFreifeld / sHe3 << endl;
    }
    if (Draw)
    {
        LoadStyles();
        gROOT->SetStyle("SinglePadStyle");
        gROOT->ForceStyle(kTRUE);

        gQ->SetMarkerColor(kRed);
        gQ->SetLineColor(kRed);
        gQ->SetLineWidth(2);
        gQ->SetMarkerSize(2);
        gQ->SetMarkerStyle(21);
        gHe3->SetMarkerColor(kGreen);
        gHe3->SetLineColor(kGreen);
        gHe3->SetLineWidth(2);
        gHe3->SetMarkerSize(2);
        gHe3->SetMarkerStyle(21);
        gNM->SetMarkerColor(kBlue);
        gNM->SetLineColor(kBlue);
        gNM->SetLineWidth(2);
        gNM->SetMarkerSize(2);
        gNM->SetMarkerStyle(21);
        TLegend *l = new TLegend(0.65, 0.2, 0.9, 0.45);
        l->SetTextFont(132);
        l->AddEntry(gQ, "Q / PLC", "PE");
        l->AddEntry(gHe3, "{}^{3}He / PLC", "PE");
        l->AddEntry(gNM, "NM / PLC", "PE");
        TCanvas *c = new TCanvas();
        c->SetGrid();
        gPad->SetTicks(1,1);
        gQ->Draw("AP");
        gQ->SetTitle("; Freifeld-Messung; Z#ddot{a}hlraten-Verh#ddot{a}ltnis");
        gQ->GetXaxis()->SetRangeUser(0.5, 11.5);
        gQ->GetXaxis()->SetNdivisions(112);
        gQ->GetYaxis()->SetNdivisions(207);
        gHe3->Draw("sameP");
        gNM->Draw("sameP");
        l->Draw();
        c->Update();
    }
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Save(fAna, "Freifeld", gQ, "Q");
    Save(fAna, "Freifeld", gHe3, "He3");
    Save(fAna, "Freifeld", gGM, "GM");
    Save(fAna, "Freifeld", gNM, "NM");
    fAna->Save();
    fAna->Close();
}

Bool_t GetRun(string RunName, Double_t &monitor, Double_t &delta_rel, Double_t &t_real)
{
    // Open txt tabular, break if not successful
    std::ifstream input(neutron_field_data);
    if (!input.is_open())
    {
        cout << "Could not open " << neutron_field_data << endl;
        return 0;
    }
    string Run = "", lastRun;
    Double_t Monitor, DeltaRel, tReal;
    do {
        lastRun = Run; // keep last run name
        // Read line
        input >> Run >> Monitor >> DeltaRel >> tReal;
        // Break if last line was reached
        if (!strcmp(Run.c_str(), lastRun.c_str()))
            return 0;
        // repeat until names match
    } while (strcmp(Run.c_str(), RunName.c_str()));
    // Save numbers to reference
    monitor = Monitor;
    delta_rel = DeltaRel;
    t_real = tReal;
    return 1;
}

string NeutronFieldRun(string Run)
{ // Find 1 run's neutron flux. Return output string
    TFile* fAna = TFile::Open(results_file, "UPDATE");
    Double_t Yield = 2.217E4;
    Double_t DYield = 240;// 0.003 * Yield;
    char name[64] = "";

    // choose fission chamber
    string FC = (Run[0] == 'U') ? "UFC" : "PuFC";
    // Get monitor data
    Double_t Monitor, DeltaRel, tReal;
    if (!GetRun(Run, Monitor, DeltaRel, tReal))
        cout << "Could not find run " << Run << endl;

    // Create graphs
    TGraphErrors *geFluence = new TGraphErrors(8);
    TGraphErrors *geFlux = new TGraphErrors(8);
    TGraphErrors *geDensity = new TGraphErrors(8);

    // Initialize returned string
    std::stringstream line;
    line << Run << " " << Monitor << " " << DeltaRel * Monitor << " " << Monitor / tReal << " " << DeltaRel * Monitor / tReal;

    // loop over deposits
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t w = SolidAngle(i, FC);
        Double_t area = DepositArea(i, FC);
        Double_t nFluence = Yield * Monitor * w;
        Double_t DnFluence = sqrt( /*pow(DYield * Monitor * w, 2) +*/
                                   pow(Yield * DeltaRel * Monitor * w, 2) +
                                   pow(Yield * Monitor * 0, 2) );
        Double_t nFlux = nFluence / area;
        Double_t DnFlux = sqrt( pow(DnFluence / area, 2) +
                                pow(0 * nFluence / (area*area), 2) );
        Double_t nDensity = nFlux / tReal;
        Double_t DnDensity = DnFlux / tReal;

//        cout << i+1 << " " << Monitor << " " << Yield << " " << w << " " << nFluence << " " << area << " " << nFlux << " " << tReal << " " << nDensity << endl;

        geFluence->SetPoint(i, i+1, nFluence);
        geFluence->SetPointError(i, 0, DnFluence);
        geFlux->SetPoint(i, i+1, nFlux);
        geFlux->SetPointError(i, 0, DnFlux);
        geDensity->SetPoint(i, i+1, nDensity);
        geDensity->SetPointError(i, 0, DnDensity);

        line << " " << nDensity << " " << DnDensity;
    }
    // Save
    Save(fAna, FC+"/NeutronField/"+Run, geFluence, Run+"_Fluence");
    Save(fAna, FC+"/NeutronField/"+Run, geFlux, Run+"_Flux");
    Save(fAna, FC+"/NeutronField/"+Run, geDensity, Run+"_Density");
    fAna->Save();
    fAna->Close();
    return line.str();
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

void DoNeutronField(string FC)
{
    string FileName = "NeutronField_"+FC+".txt";
    ofstream output("../results/"+FileName);
    output << "# nr run monitor unc monitor[1/s] unc[1/s] {flux[1/(s*mm^2)] unc[1/(s*mm^2)]}" << endl;
    Int_t nRuns = GetnRuns(FC);
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(FC, j);
        // Analyze monitor data
        output << j+1 << " " << NeutronFieldRun(Run) << endl;
    }
    output.close();

    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gMon = new TGraphErrors(("../results/"+FileName).c_str(), "%lg %*s %*lg %*lg %lg %lg");
    if (!gMon) cout << "Could not create TGraphErrors " << FileName << endl;
    gMon->SetName((FC+"_MonitorRate").c_str());
    Double_t x, y, yerr, tStart, tLength;
    for (Int_t j = 0; j < 7; j++)
    {
        gMon->GetPoint(j, x, y);
        yerr = gMon->GetErrorY(j);
        string Run = GetRunName(FC, j);
        tStart = GetStartTime(Run);
        tLength = GetRealTime(Run);
        gMon->SetPoint(j, tStart + 0.5*tLength, y);
        gMon->SetPointError(j, 0.5*tLength, yerr);
    }
    SaveToFile(fAna, FC+"/NeutronField", gMon);
    fAna->Save();
    fAna->Close();
    new TCanvas();
    gMon->Draw();
}

void DoNeutronField()
{
    string FileName = "NeutronField.txt";
    ofstream output("../results/"+FileName);
    output << "# nr run monitor unc monitor[1/s] unc[1/s] {flux[1/(s*mm^2)] unc[1/(s*mm^2)]}" << endl;
    Int_t nRuns = GetnRuns();
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(j);
        // Analyze monitor data
        output << j+1 << " " << NeutronFieldRun(Run) << endl;
    }
    output.close();
}

void NeutronField()
{
    Freifeld();

    DoNeutronField("UFC");
//    DoNeutronField("UFC_FG");
//    DoNeutronField("UFC_BG");
    NeutronFieldRun("UFC_NIF");
    NeutronFieldRun("UFC_SB");

    DoNeutronField("PuFC");
//    DoNeutronField("PuFC_FG");
//    DoNeutronField("PuFC_BG");
    NeutronFieldRun("NIF");
    NeutronFieldRun("SB");
}

#endif
