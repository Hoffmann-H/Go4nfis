#ifndef SHADOW_CONE_H
#define SHADOW_CONE_H
#include "FC.C"
#include "SaveToFile.C"
#include "AnaSim.C"
//#include "DrawPics.C"

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
    gS->SetTitle("Shadow Cone scattered portion, experimental; Deposit; (n,f)_{sc} / (n,f)_{tot}");

    Double_t x, fisFG, fisBG, fluFG, fluBG, dfisFG, dfisBG, dfluFG, dfluBG, S, dS;
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
        Int_t b0 = h2Ideal->GetYaxis()->FindBin(14.7);
        Int_t b1 = h2Ideal->GetYaxis()->FindBin(15.1);
        Double_t T2 = (h2Tot->Integral(0, -1, b0, b1) - h2Sc->Integral(0, -1, b0, b1)) / h2Ideal->Integral(0, -1, b0, b1) / fTarget;
        Double_t DT2 = T2 * sqrt(1.0 / (h2Tot->GetEntries() - h2Sc->GetEntries()) + 1.0 / h2Ideal->GetEntries());
//        cout << i+1 << "\t " << T1 << " +/- " << DT1 << endl;

        ge->SetPoint(i, i+1, T1);
        ge->SetPointError(i, 0, DT1);
    }
    return ge;
}

void BiasX(TGraphErrors *pG, Double_t dx)
{
    Double_t x, y;
    for (Int_t i = 0; i < pG->GetN(); i++)
    {
        pG->GetPoint(i, x, y);
        pG->SetPoint(i, x + dx, y);
    }
}

void Transmission()
{
    TFile* fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");

    // UFC
    TGraphErrors *geAnaU = AnalyticTransmission("UFC", "TransMis.dat");
    TGraphErrors *geSimU = SimTransmission(fAna, "UFC");
    // PuFC
    TGraphErrors *geAnaPu = AnalyticTransmission("PuFC", "TransMis.dat");
    TGraphErrors *geSimPu = SimTransmission(fAna, "PuFC");

//    SetSize(geAnaU);
    geAnaU->GetXaxis()->SetNdivisions(109);
    geAnaU->SetLineColor(kRed);
//    geAnaU->SetLineStyle(0);
    geAnaU->SetMarkerStyle(20);
    geAnaU->SetMarkerColor(kRed);
    geSimU->SetLineColor(kRed);
    geSimU->SetMarkerStyle(21);
    geSimU->SetMarkerColor(kRed);
    geAnaPu->SetLineColor(kBlue);
    geAnaPu->SetMarkerStyle(20);
    geAnaPu->SetMarkerColor(kBlue);
    geSimPu->SetLineColor(kBlue);
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
    TFile* fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
//    TGraphErrors *gExp = ExpShadowCone(fAna, "UFC_NIF", "UFC_SB");
//    TGraphErrors *gSim = SimShadowCone(fAna, "Geant4", "UFC");
//    TGraphErrors *gSim2 = SimShadowCone2(fAna, "Geant4", "UFC");

//    gSim2->SetLineColor(kRed);
//    TLegend *l = new TLegend(0.2, 0.4, 0.2, 0.4);
//    l->AddEntry(gExp, "exp.");
//    l->AddEntry(gSim2, "sim.");
//    new TCanvas();
//    gExp->Draw();
////    gSim->Draw("same");
//    gSim2->Draw("same");
//    l->Draw();

    Transmission();
}
#endif
