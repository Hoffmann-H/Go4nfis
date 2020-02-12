#ifndef SHADOW_CONE_H
#define SHADOW_CONE_H
#include "FC.C"
#include "SaveToFile.C"
#include "Target.C"
#include "DrawPics.C"
#include "TLegend.h"

TGraphErrors* ExpShadowCone(TFile *fAna, string RunFG = "NIF", string RunBG = "SB")
{
    string FC = (RunFG[0] == 'U') ? "UFC" : "PuFC";
    char name[64] = "";
    sprintf(name, "%s/ToF/Signal/%s/InducedFission", FC.c_str(), RunFG.c_str());
    TGraphErrors *gFisFG = (TGraphErrors*)fAna->Get(name); if (!gFisFG) cout << "Could not get " << name << endl;
    sprintf(name, "%s/ToF/Signal/%s/InducedFission", FC.c_str(), RunBG.c_str());
    TGraphErrors *gFisBG = (TGraphErrors*)fAna->Get(name); if (!gFisBG) cout << "Could not get " << name << endl;
    sprintf(name, "%s/NeutronField/%s/%s_Fluence", FC.c_str(), RunFG.c_str(), RunFG.c_str());
    TGraphErrors *gFluFG = (TGraphErrors*)fAna->Get(name); if (!gFluFG) cout << "Could not get " << name << endl;
    sprintf(name, "%s/NeutronField/%s/%s_Fluence", FC.c_str(), RunBG.c_str(), RunBG.c_str());
    TGraphErrors *gFluBG = (TGraphErrors*)fAna->Get(name); if (!gFluBG) cout << "Could not get " << name << endl;

    TGraphErrors *gS = new TGraphErrors(8);
    sprintf(name, "%s_SB_Exp", FC.c_str());
    gS->SetName(name);
    gS->SetTitle("Shadow Cone scattered portion, experimental; Deposit; (n,f) #font[12]{C}_{sc} / #font[12]{C}_{tot}");

    Double_t x, fisFG, fisBG, fluFG, fluBG, dfisFG, dfisBG, dfluFG, dfluBG, S, dS;
    cout << "ch\t fisFG\t fisBG\t fluFG\t fluBG" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        gFisFG->GetPoint(i, x, fisFG);
        dfisFG = gFisFG->GetErrorY(i) / fisFG;
        gFisBG->GetPoint(i, x, fisBG);
        dfisBG = gFisBG->GetErrorY(i) / fisBG;
        gFluFG->GetPoint(i, x, fluFG);
        dfluFG = gFluFG->GetErrorY(i) / fluFG;
        gFluBG->GetPoint(i, x, fluBG);
        dfluBG = gFluBG->GetErrorY(i) / fluBG;
        S = fisBG / fisFG * fluFG / fluBG;
        dS = sqrt(pow(dfisFG, 2) + pow(dfisBG, 2) + pow(dfluFG, 2) + pow(dfluBG, 2));
        cout << i+1 << "  \t" << fisFG << "  \t" << fisBG << "  \t" << fluFG << "  \t" << fluBG << "  \t" << S << endl;
        gS->SetPoint(i, i+1, S);
        gS->SetPointError(i, 0, dS * S);
    }
    return gS;
}

TGraphErrors* ExpShadowConeHardCode()
{
    char name[64] = "";
    TGraphErrors *gS = new TGraphErrors(8);
    gS->SetName("UFC_SB_Exp");
    gS->SetTitle("Shadow Cone scattered portion, experimental; Deposit; (n,f)_{sc} / (n,f)_{tot}");
    Double_t S[] = {0.0974296259, 0.1051633466, 0.1202697731, 0.1647621571, 0.142541863, 0.1067283051, 0.0952544876, 0.1264492991};
    Double_t DS[] = {0.0109502724, 0.0131367478, 0.0138942925, 0.014743786, 0.0135914887, 0.0125611275, 0.0128016384, 0.0135880239};
    for (Int_t i = 0; i < 8; i++)
    {
        gS->SetPoint(i, i+1, S[i]);
        gS->SetPointError(i, 0, DS[i]);
    }
    return gS;
}

TGraphErrors* SimShadowCone(TFile *fAna, string Simulation = "Geant4", string FC = "UFC")
{
    char name[128] = "";
    sprintf(name, "%s/Correction/%s_%s_real_C", FC.c_str(), Simulation.c_str(), FC.c_str());
    TGraphErrors *gFG = (TGraphErrors*)fAna->Get(name); if (!gFG) cout << "Could not get " << name << endl;
    sprintf(name, "%s/Correction/%s_%s_SB_C", FC.c_str(), Simulation.c_str(), FC.c_str());
    TGraphErrors *gBG = (TGraphErrors*)fAna->Get(name); if (!gBG) cout << "Could not get " << name << endl;

    TGraphErrors *gS = new TGraphErrors(8);
    sprintf(name, "%s_SB_Sim", FC.c_str());
    gS->SetName(name);
    sprintf(name, "Shadow Cone scattered portion, %s; Deposit; (n,f)_{sc} / (n,f)_{tot}", Simulation.c_str());
    gS->SetTitle(name);

    Double_t x, cFG, cBG, dcFG, dcBG, S, dS;
    for (Int_t i = 0; i < 8; i++)
    {
        gFG->GetPoint(i, x, cFG);
        dcFG = gFG->GetErrorY(i) / cFG;
        gBG->GetPoint(i, x, cBG);
        dcBG = gBG->GetErrorY(i) / cBG;
        S = cFG / cBG;
        dS = sqrt(pow(dcFG, 2) + pow(dcBG, 2));
        cout << i+1 << "  \t" << cFG << "  \t" << cBG << "  \t" << S << endl;
        gS->SetPoint(i, i+1, S);
        gS->SetPointError(i, 0, dS * S);
    }
    return gS;
}

TGraphErrors* SimShadowCone2(TFile *fAna, string Simulation = "Geant4", string FC = "UFC")
{
    char name[128] = "";
    TGraphErrors *gS = new TGraphErrors(8);
    sprintf(name, "%s_SB_Sim2", FC.c_str());
    gS->SetName(name);
    sprintf(name, "Shadow Cone scattered portion, %s; Deposit; (n,f)_{sc} / (n,f)_{tot}", Simulation.c_str());
    gS->SetTitle(name);

    Double_t nFG, nBG, fisFG, dfisFG, fisBG, dfisBG, S, dS;
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "Simulation/%s/%s_real/ToFvsEkin/%s_ToFvsEkin_real_Ch.%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH2D *hFG = (TH2D*)fAna->Get(name); if (!hFG) cout << "Could not get " << name << endl;
        sprintf(name, "Simulation/%s/%s_SB/ToFvsEkin/%s_ToFvsEkin_SB_Ch.%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH2D *hBG = (TH2D*)fAna->Get(name); if (!hBG) cout << "Could not get " << name << endl;
        nFG = hFG->GetEntries();
        fisFG = hFG->Integral();
        dfisFG = 1.0 / sqrt(nFG);
        nBG = hBG->GetEntries();
        fisBG = hBG->Integral();
        dfisBG = 1.0 / sqrt(nBG);
        S = fisBG / fisFG;
        dS = sqrt(pow(dfisFG, 2) + pow(dfisBG, 2));
        cout << i+1 << "  \t" << fisFG << "  \t" << fisBG << "  \t" << S << endl;
        gS->SetPoint(i, i+1, S);
        gS->SetPointError(i, 0, dS * S);
    }
    return gS;
}

TGraphErrors* SimScatteredPart(TFile *fAna, string Simulation = "Geant4", string FC = "UFC")
{
    char name[128] = "";
    // Get correction factor
    sprintf(name, "%s/Correction/%s_%s_real_C", FC.c_str(), Simulation.c_str(), FC.c_str());
    TGraphErrors *gCf = (TGraphErrors*)fAna->Get(name); if (!gCf) cout << "Could not get " << name << endl;

    // Create scattered portion graph
    TGraphErrors *gSc = new TGraphErrors(8);
    sprintf(name, "%s_SB_SimSc", FC.c_str());
    gSc->SetName(name);
    sprintf(name, "%s Simulated Scattered Portion; Deposit; (n,f)_{sc} / (n,f)_{tot}", FC.c_str());
    gSc->SetTitle(name);

    Double_t x, cf, Dcf;
    // Loop over Deposits
    for (Int_t i = 0; i < 8; i++)
    {
        // Read correction factor
        gCf->GetPoint(i, x, cf);
        Dcf = gCf->GetErrorY(i);

        // Write scattered portion to graph
        gSc->SetPoint(i, i+1, 1.0 - cf);
        gSc->SetPointError(i, 0, Dcf);
    }
    return gSc;
}

TGraphErrors* AnalyticTransmission(string FC = "PuFC", string file = "TransMis.dat")
{
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Experiment/Calculation/%s", file.c_str());
    string format = strcmp(FC.c_str(), "PuFC") ? "%lg %lg %lg" : "%lg %*lg %*lg %lg %lg";
    TGraphErrors *ge = new TGraphErrors(name, format.c_str());
    return ge;
}

TGraphErrors* SimTransmission(TFile* fAna, string FC = "PuFC", string Simulation = "Geant4")
{ // Only works for Geant4
    char name[128] = "";
    Double_t fTarget = TargetFactor(fAna);
    TGraphErrors *ge = new TGraphErrors(8);
    ge->SetNameTitle((FC+"_SimT").c_str(), "Geant4 Simulated Transmission; Deposit; T");

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "Simulation/%s/%s_real/ToFvsEkin/%s_ToFvsEkin_real_Ch.%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH2D *h2Tot = (TH2D*)fAna->Get(name); if (!h2Tot) cout << "Could not get " << name << endl;
        sprintf(name, "Simulation/%s/%s_real/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH2D *h2Sc = (TH2D*)fAna->Get(name); if (!h2Sc) cout << "Could not get " << name << endl;
        sprintf(name, "Simulation/%s/%s_ideal/ToFvsEkin/%s_ToFvsEkin_ideal_Ch.%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH2D *h2Ideal = (TH2D*)fAna->Get(name); if (!h2Ideal) cout << "Could not get " << name << endl;
        Double_t T1 = (h2Tot->GetEntries() - h2Sc->GetEntries()) / h2Ideal->GetEntries() / fTarget;
        Double_t DT1 = T1 * sqrt(1.0 / (h2Tot->GetEntries() - h2Sc->GetEntries()) + 1.0 / h2Ideal->GetEntries());
        Int_t b0 = h2Ideal->GetYaxis()->FindBin(14.0);
        Int_t b1 = h2Ideal->GetYaxis()->FindBin(16.0);
        Double_t T2 = (h2Tot->Integral(0, -1, b0, b1) - h2Sc->Integral(0, -1, b0, b1)) / h2Ideal->Integral(0, -1, b0, b1) / fTarget;
        Double_t DT2 = T2 * sqrt(1.0 / (h2Tot->Integral(0, -1, b0, b1) - h2Sc->Integral(0, -1, b0, b1)) * (h2Tot->Integral() - h2Sc->Integral()) / (h2Tot->GetEntries() - h2Sc->GetEntries()) + 1.0 / h2Ideal->Integral(0, -1, b0, b1) * h2Ideal->Integral() / h2Ideal->GetEntries());
//        cout << i+1 << "\t " << T1 << " +/- " << DT1 << endl;

        ge->SetPoint(i, i+1, T2);
        ge->SetPointError(i, 0, DT2);
    }
    return ge;
}

//void BiasX(TGraphErrors *pG, Double_t dx)
//{
//    Double_t x, y;
//    for (Int_t i = 0; i < pG->GetN(); i++)
//    {
//        pG->GetPoint(i, x, y);
//        pG->SetPoint(i, x + dx, y);
//    }
//}

void Transmission(string file = "TransMis.dat")
{
    TFile* fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");

    // UFC
    TGraphErrors *geAnaU = AnalyticTransmission("UFC", file);
    TGraphErrors *geSimU = SimTransmission(fAna, "UFC");
    // PuFC
    TGraphErrors *geAnaPu = AnalyticTransmission("PuFC", file);
    TGraphErrors *geSimPu = SimTransmission(fAna, "PuFC");
    geAnaU->SetTitle("; Deposit; #font[12]{T}");
//    SetSize(geAnaU);
    geAnaU->GetXaxis()->SetNdivisions(109);
    geAnaU->SetLineColor(kRed);
//    geAnaU->SetLineStyle(0);
    geAnaU->SetLineWidth(2);
    geAnaU->SetMarkerSize(2);
    geAnaU->SetMarkerStyle(20);
    geAnaU->SetMarkerColor(kRed);
    geSimU->SetLineColor(kRed);
    geSimU->SetLineWidth(2);
    geSimU->SetMarkerSize(2);
    geSimU->SetMarkerStyle(21);
    geSimU->SetMarkerColor(kRed);
    geAnaPu->SetLineColor(kBlue);
    geAnaPu->SetLineWidth(2);
    geAnaPu->SetMarkerSize(2);
    geAnaPu->SetMarkerStyle(20);
    geAnaPu->SetMarkerColor(kBlue);
    geSimPu->SetLineWidth(2);
    geSimPu->SetLineColor(kBlue);
    geSimPu->SetMarkerSize(2);
    geSimPu->SetMarkerStyle(21);
    geSimPu->SetMarkerColor(kBlue);
    BiasX(geAnaU, -0.15);
    BiasX(geSimU, -0.05);
    BiasX(geAnaPu, +0.05);
    BiasX(geSimPu, +0.15);
    TLegend *l = new TLegend(0.6, 0.2, 0.85, 0.4);
    l->AddEntry(geAnaU, "UFC, Analytisch, #Delta_{sys}");
    l->AddEntry(geSimU, "UFC, Geant4, #Delta_{stat}");
    l->AddEntry(geAnaPu, "PuFC, Analytisch, #Delta_{sys}");
    l->AddEntry(geSimPu, "PuFC, Geant4, #Delta_{stat}");

    new TCanvas();
    gPad->SetTicks(1, 1);
    geAnaU->Draw("AP");
    geSimU->Draw("sameP");
    geAnaPu->Draw("sameP");
    geSimPu->Draw("sameP");
    l->Draw();
}

void ShadowCone()
{
    LoadStyles();
    gROOT->SetStyle("SinglePadStyle");
    gROOT->ForceStyle(kTRUE);
    TFile* fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
//    TGraphErrors *gExp = ExpShadowCone(fAna, "UFC_NIF", "UFC_BG_MS20_5");
//    TGraphErrors *gExp = ExpShadowCone(fAna, "UFC_NIF", "UFC_BG_MS21_4");
    TGraphErrors *gExp = ExpShadowCone(fAna, "UFC_NIF", "UFC_SB");
//    TGraphErrors *gExp = ExpShadowConeHardCode();
    TGraphErrors *gSim = SimShadowCone(fAna, "Geant4", "UFC");
//    TGraphErrors *gSim2 = SimShadowCone2(fAna, "Geant4", "UFC");
    TGraphErrors *gSimSc = SimScatteredPart(fAna, "Geant4", "UFC");

    BiasX(gExp, -0.05);
    BiasX(gSim, +0.05);
    SetSize(gExp);
    gExp->SetLineWidth(2);
    gExp->SetMarkerSize(2);
    gExp->SetMarkerStyle(22);
    gSim->SetLineWidth(2);
    gSim->SetMarkerSize(2);
    gSim->SetMarkerStyle(22);
    gSim->SetMarkerColor(kRed);
    gSim->SetLineColor(kRed);
    gSimSc->SetLineWidth(2);
    gSimSc->SetMarkerSize(2);
    gSimSc->SetMarkerStyle(20);
    gSimSc->SetMarkerColor(kBlue);
    gSimSc->SetLineColor(kBlue);
    TLegend *l = new TLegend(0.3, 0.2, 0.8, 0.4);
    l->AddEntry(gExp, "Messung: SB / offen", "PE");
    l->AddEntry(gSim, "Simulation: SB / offen", "P");
    l->AddEntry(gSimSc, "Simulation: offen / Vakuum", "P");
    l->SetTextFont(132);
    new TCanvas();
    gExp->Draw("AP");
    gExp->GetYaxis()->SetRangeUser(0, 0.14);
    gExp->GetYaxis()->SetTitle("#font[12]{C}_{sc} / #font[12]{C}_{tot}");
    gExp->GetYaxis()->SetNdivisions(505);
    gExp->GetXaxis()->SetNdivisions(110);
    gSim->Draw("same P");
    gSimSc->Draw("same P");
    l->Draw();

//    Transmission("TransMis.dat");
}
#endif
