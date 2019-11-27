#ifndef VGLSIM_H
#define VGLSIM_H
#include "TMath.h"
#include "TLatex.h"
#include "DrawPics.C"
#include "FC.C"

Double_t Hits(TFile *f, string FC, Int_t ch)
{
    char name[128] = "";
    sprintf(name, "Simulation/Geant4/%s_real/ToFvsEkin/%s_ToFvsEkin_real_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2D *pH2 = (TH2D*)f->Get(name); if (!pH2) cout << "Could not get " << name << endl;
    return pH2->GetEntries();
}

Double_t ScatteredHits(TFile *f, string FC, Int_t ch)
{
    char name[128] = "";
    sprintf(name, "Simulation/Geant4/%s_ideal/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2D *pH2 = (TH2D*)f->Get(name); if (!pH2) cout << "Could not get " << name << endl;
    return pH2->GetEntries();
}

Double_t DirectHits(TFile *f, string FC, Int_t ch)
{
    return Hits(f, FC, ch) - ScatteredHits(f, FC, ch);
}

Double_t F()
{ // Deposit area
    return TMath::Pi() * 3.7 * 3.7; // in cm^2
}

TGraphErrors *GetSigma(Bool_t G4 = 0)
{
    if (!G4)
    {
        TGraphErrors *gSigma = new TGraphErrors("/home/hoffma93/Programme/ROOT/Data/Pu242.dat", "%lg %lg %lg");
        if (!gSigma) cout << "Could not open " << "/home/hoffma93/Programme/ROOT/Data/Pu242.dat" << endl;
        gSigma->SetNameTitle("gSigma", "ENDF WQ; E/MeV; n");
        return gSigma;
    } else {
        TGraphErrors *gSigma = new TGraphErrors("/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/Pu_ENDF-VIII.dat", "%lg %lg");
        if (!gSigma) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/Pu_ENDF-VIII.dat" << endl;
        gSigma->SetNameTitle("gSigma", "G4 WQ; E/MeV; n");
        return gSigma;
    }
}

Double_t AvSigma(TGraphErrors *gSigma, TH1F *hE)
{
    Int_t nBins = hE->GetNbinsX();
    Double_t sum_val = 0;
    for (Int_t j = 1; j < nBins+1; j++)
    {
        Double_t SigmaBinCenter = gSigma->Eval(hE->GetXaxis()->GetBinCenter(j));
        sum_val += SigmaBinCenter * hE->GetBinContent(j);
    }
    return sum_val / hE->Integral();
}

Double_t oldCorrection(string FC, Int_t ch)
{ // hits, FG, trackID, solid angle, 15MeV
    char name[128] = "";
    Double_t DepositArea = F();
    // AnaSim.C
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");
    if (!fAna) cout << "/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root" << endl;
    sprintf(name, "Simulation/Target/Direct"); TH1D *hTargetDir = (TH1D*)fAna->Get(name); if (!hTargetDir) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/Target/Total"); TH1D *hTargetTot = (TH1D*)fAna->Get(name); if (!hTargetTot) cout << "Could not get " << name << endl;
    Double_t fTarget = hTargetDir->Integral() / hTargetTot->Integral();
    cout << "FG/(FG+BG) = " << fTarget << endl;
    // 1/cos(theta), sigma(E), real, FG
    TFile *f1 = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_1E8.root");
    if (!f1) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_1E8.root" << endl;
    TH1F *hT1 = (TH1F*)f1->Get("Source/Source_Theta"); if (!hT1) cout << "Coudl not get " << "Source/Source_Theta" << endl;
    Double_t N1 = hT1->Integral();
    // 1/cos(theta, sigma(E), ideal, FG+BG
    TFile *f2 = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_ideal_5E7.root");
    if (!f2) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_ideal_5E7.root" << endl;
    TH1F *hT2 = (TH1F*)f2->Get("Source/Source_Theta"); if (!hT2) cout << "Coudl not get " << "Source/Source_Theta" << endl;
    Double_t N2 = hT2->Integral();
    TH1F *hEkin = (TH1F*)f2->Get("Source/nEnergy/Source_Ekin"); if (!hEkin) cout << "Could not get " << "Source/nEnergy/Source_Ekin" << endl;
    Double_t SigmaSource = AvSigma(GetSigma(1), hEkin);
    // FG+BG-Spektrum
    TFile *f3 = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_1E7.root");
    if (!f3) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_1E7.root" << endl;
    TH1F *hT3 = (TH1F*)f3->Get("Source/Source_Theta"); if (!hT3) cout << "Coudl not get " << "Source/Source_Theta" << endl;
    Double_t N3 = hT3->Integral();
    // Hit, sigma(bc), real, FG
    TFile *f = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_Hit_5E7.root");
    if (!f) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_Hit_5E7.root" << endl;
    TH1F *hT = (TH1F*)f->Get("Source/Source_Theta"); if (!hT) cout << "Coudl not get " << "Source/Source_Theta" << endl;
    Double_t N = hT->Integral();

    TGraphErrors *gSigma = GetSigma();
    Double_t Sigma15MeV = gSigma->Eval(15.0);
    TH1F *hE[2][8];
    TGraph *gF = new TGraph(8);
    gF->SetNameTitle("gF", "Korrekturfaktor; Deposit; Faktor");
    TGraph *gT = new TGraph(8);
    gT->SetNameTitle("gT", "Korrekturfaktor; Deposit; T");
    TGraph *gS = new TGraph(8);
    gS->SetNameTitle("gS", "Korrekturfaktor; Deposit; S");
    TGraph *gF2 = new TGraph(8);
    gF2->SetNameTitle("gF2", "Korrekturfaktor; Deposit; Faktor");
    TGraph *gSigmaDir = new TGraph(8);
    gSigmaDir->SetNameTitle("gSigmaDir", "Mittlerer WQ direkter Neutronen; Deposit; #sigma / b");
    sprintf(name, "%s/Correction/%s_Geant4_C", FC.c_str(), FC.c_str());
    TGraph *gTL = (TGraph*)fAna->Get(name); if (!gTL) cout << "Could not get " << name << endl;
    TGraph *g3 = new TGraph(8);
    g3->SetNameTitle("gFis_TL_E", "Track length, E; Deposit; #font[12]{p}(n,f)");
    TGraph *gSource = new TGraph(8);
    gSource->SetNameTitle("gSource", "Neutronen; Deposit; #font[12]{N}_{n}");
    TGraph *gFis = new TGraph(8);
    gFis->SetNameTitle("gFis", "Neutronen; Deposit; #font[12]{N}_{n}");
    Double_t n0 = 0, n1 = 0;

    Double_t fis0 = 0, fis1 = 0, sig0 = 0;
    cout << "ch \tF       \tdFG      \tdHit      \tdVoid       \tAnaSim2 / F \tAnaSim2   \tSigma      \tnIdeal / nQu2" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get hit neutron spectra
        sprintf(name, "%s/nEnergy/%s_Ekin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH1F *hTot = (TH1F*)f->Get(name); if (!hTot) cout << "Could not get " << name << endl;
        sprintf(name, "%s/nEnergy/Scattered/%s_Ekin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH1F *hSc = (TH1F*)f->Get(name); if (!hSc) cout << "Could not get " << name << endl;
        TH1F *hDir = (TH1F*)hTot->Clone();
        hDir->Add(hSc, -1);

        // Scale with sigma(E)
        Int_t nBinsE = hTot->GetNbinsX();
        sprintf(name, "%s_ProjT_Tot_%i", FC.c_str(), i+1);
        hE[0][i] = new TH1F(name, "Spektrum spaltender Neutronen; #font[12]{t} / ns; #font[12]{N}_{eff}",
                            nBinsE, hTot->GetXaxis()->GetBinLowEdge(1), hTot->GetXaxis()->GetBinLowEdge(nBinsE + 1));
        sprintf(name, "%s_ProjT_Dir_%i", FC.c_str(), i+1);
        hE[1][i] = new TH1F(name, "Spektrum spaltender Neutronen; #font[12]{t} / ns; #font[12]{N}_{eff}",
                            nBinsE, hTot->GetXaxis()->GetBinLowEdge(1), hTot->GetXaxis()->GetBinLowEdge(nBinsE + 1));
        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
        TH1F *hQu = (TH1F*)f->Get(name); if (!hQu) cout << "Could not get " << name << endl;
        for (Int_t j = 1; j < nBinsE+1; j++)
        {
            Double_t SigmaBinCenter = gSigma->Eval(hTot->GetXaxis()->GetBinCenter(j));
            Double_t pTot = hTot->GetBinContent(j);
            Double_t pSc = hSc->GetBinContent(j);
            Double_t pDir = hDir->GetBinContent(j);
//            cout << hProjTot->Integral() << " " << hProjSc->Integral() << "\t";
            hE[0][i]->SetBinContent(j, pTot * SigmaBinCenter / Sigma15MeV);
            hE[1][i]->SetBinContent(j, pDir * SigmaBinCenter / Sigma15MeV);
        }
        // correction factors
        Double_t cTot = hE[0][i]->Integral();
        Double_t cDir = hE[1][i]->Integral();
        Double_t nDir = hTot->Integral() - hSc->Integral();
        Double_t nQu = hQu->Integral();
        Double_t F = cDir / cTot * nQu / nDir;
        Double_t T = nDir / nQu;
        Double_t S = cDir / cTot;
        Double_t SigmaDir = AvSigma(gSigma, hDir);
        Double_t F2 = SigmaDir / cTot * nQu / Sigma15MeV;
//        cout << i+1 << "\t" << F << "\t" << F2 << "\t" << T << "\t" << S << endl;

        // fill graphs
        gF->SetPoint(i, i+1, F);
        gT->SetPoint(i, i+1, T);
        gS->SetPoint(i, i+1, S);
        gSigmaDir->SetPoint(i, i+1, SigmaDir);
        gF2->SetPoint(i, i+1, F2);

        // sum up
        fis0 += F;
        Double_t x, AnaSim;
        gTL->GetPoint(i, x, AnaSim);
        fis1 += AnaSim;
        sig0 += SigmaDir;

        // recalculate new correction factor
        sprintf(name, "%s/nEnergy/%s_Ekin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH1F *hEreal = (TH1F*)f1->Get(name);
        TH1F *hEideal = (TH1F*)f2->Get(name); if (!hEreal || !hEideal) cout << "Could not get " << name << endl;
        Double_t pReal = hEreal->Integral() / N1;
        Double_t pIdeal = hEideal->Integral() / N2;
        Double_t nIdeal = hEideal->GetEntries();
        Double_t AnaSim2 = pIdeal / pReal * fTarget;
//        cout << i+1 << "\t" << AnaSim << "\t" << AnaSim2 << "\t" << AnaSim / AnaSim2 << endl;

        // TARGET FG / BG
        sprintf(name, "%s/nEnergy/%s_Ekin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH1F *hE_FG = (TH1F*)f3->Get(name); if (!hE_FG) cout << "Could not get " << name << endl;
        Double_t pFG = hE_FG->Integral() / N3;
        Double_t dFG = pFG / pReal * fTarget;

        // Hit / cos(theta)
        Double_t pHit = cTot / N * Sigma15MeV / DepositArea;
        Double_t dHit = pHit / pFG;
//        cout << i+1 << "\t" << dFG << "\t" << dHit << "\t" << AnaSim2 / F << "\t" << endl;

        // Void sim. / dir. n // Void / source angular
///        d = AnaSim2 / F = pIdeal / pReal * fTarget / cDir * cTot / nQu * nDir = dFG * dHit * dVoid
///          = pFG / pReal * fTarget * pHit / pFG * dVoid
///        dVoid = pIdeal / pReal * fTarget / cDir * cTot / nQu * nDir / (pFG / pReal * fTarget * pHit / pFG)
///              = pIdeal / cDir * cTot / nQu * nDir / pHit
        Double_t dVoid = pIdeal / pHit * (cTot / nQu) / (cDir / nDir);
        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
        TH1F *hAng = (TH1F*)f2->Get(name); if (!hAng) cout << "Could not get " << name << endl;
        Double_t nQu2 = hAng->Integral();
        gSource->SetPoint(i, i+1, nQu2);
        gFis->SetPoint(i, i+1, nIdeal);
        n1 += nIdeal;
        n0 += nQu2;

        cout << i+1 << " \t" << F << " \t" << dFG << " \t" << dHit << " \t" << dVoid << " \t" << AnaSim2 / F << " \t" << AnaSim2 << " \t" << SigmaSource / SigmaDir << " \t" << nIdeal / nQu2 << endl;

    }
    // Draw one channel's neutron spectra
    sprintf(name, "cDir_%i", ch+1);
    new TCanvas(name, "Effektives n-Spektrum");
    gPad->SetTicks(1, 1);
    SetSize(hE[0][ch]);
    hE[1][ch]->SetLineColor(kGreen);
    sprintf(name, "%s Kanal %i", FC.c_str(), ch+1);
    TLegend *l = new TLegend(0.15, 0.7, 0.4, 0.85, name);
    l->AddEntry(hE[0][ch], "Total");
    l->AddEntry(hE[1][ch], "Direkt");
    hE[0][ch]->Draw("hist");
    hE[1][ch]->Draw("same hist");
    l->Draw();
    TLatex *t[6];
    t[0] = new TLatex();
    sprintf(name, "#Sigma = %.0f", hE[0][ch]->Integral());
    t[0]->DrawLatexNDC(0.4, 0.5, name);
    t[1] = new TLatex();
    sprintf(name, "#Sigma = %.0f", hE[1][ch]->Integral());
    t[1]->SetTextColor(kGreen);
    t[1]->DrawLatexNDC(0.4, 0.4, name);

    // Draw correction factors
    new TCanvas("cAll", "Correction correction");
    gPad->SetTicks(1, 1);
    SetSize(gF);
    gF->GetYaxis()->SetRangeUser(0, 1);
    gF->SetLineColor(kBlack);
    gF2->SetLineColor(kBlue);
    gTL->SetLineColor(kRed);
    l = new TLegend(0.15, 0.15, 0.7, 0.3);
    l->AddEntry(gF, "Hit, (n,f) / n");
    l->AddEntry(gF2, "Hit, <#sigma>");
    l->AddEntry(gTL, "TL, vac, FG+BG, #sigma(E)");
    gF->Draw();
    gF2->Draw("same");
    gTL->Draw("same");
    l->Draw();
    t[2] = new TLatex();
    sprintf(name, "#delta = %.3f", fis1 / fis0);
    t[2]->DrawLatexNDC(0.3, 0.5, name);

    // Draw cross ssection
    new TCanvas("cSig", "Cross section");
    gPad->SetTicks(1, 1);
    SetSize(gSigmaDir);
    TGraph *gSigmaSource = new TGraph(2);
    gSigmaSource->SetPoint(0, 1, SigmaSource);
    gSigmaSource->SetPoint(1, 8, SigmaSource);
    gSigmaSource->SetLineColor(kRed);
    gSigmaDir->Draw();
    gSigmaSource->Draw("same");
    l = new TLegend(0.15, 0.15, 0.7, 0.3, "Spektrum");
    l->AddEntry(gSigmaDir, "Direkte n");
    sprintf(name, "%.3f b", SigmaSource);
    l->AddEntry(gSigmaSource, name);
    l->Draw();
    t[3] = new TLatex();
    sprintf(name, "Normierung #sigma: %.5f", 8 * SigmaSource / sig0);
    t[3]->DrawLatexNDC(0.3, 0.5, name);

    // Draw void simulation influence
    new TCanvas("cAng", "Winkelverteilung");
    gPad->SetTicks(1, 1);
    SetSize(gFis);
    gFis->SetLineColor(kRed);
    l = new TLegend(0.15, 0.7, 0.4, 0.85);
    l->AddEntry(gFis, "Void-Simulation");
    l->AddEntry(gSource, "Emissionswinkel");
    gFis->Draw();
    gSource->Draw("same");
    l->Draw();
    t[4] = new TLatex();
    sprintf(name, "Raumwinkel: %f", n1 / n0);
    t[4]->DrawLatexNDC(0.5, 0.4, name);
    return 0;
}

void VglSigmaPu()
{
    TGraphErrors *gOld = GetSigma();
    TGraphErrors *gNew = new TGraphErrors("/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/Pu_ENDF-VIII.dat", "%lg %lg %lg");
    if (!gNew) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/Pu_ENDF-VIII.dat" << endl;
    new TCanvas();
    gPad->SetTicks(1, 1);
    gOld->SetTitle("ENDF WQ; #font[12]{E} / MeV; #sigma / b");
    SetSize(gOld);
    gOld->GetXaxis()->SetNdivisions();
    gNew->SetLineColor(kRed);
    TLegend *l = new TLegend(0.15, 0.7, 0.4, 0.85);
    l->AddEntry(gOld, "ENDF");
    l->AddEntry(gNew, "Geant4");
    gOld->Draw();
    gNew->Draw("same");
    l->Draw();
}

void VglSigmaU()
{
    TGraphErrors *gOld = new TGraphErrors("/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/U_ENDF-VIII.dat", "%lg %lg");
    if (!gOld) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/U_ENDF-VIII.dat" << endl;
    gOld->SetName("gOld");
    TGraphErrors *gNew = new TGraphErrors("/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/U_G4NeutronHP.dat", "%lg %lg");
    if (!gNew) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/G4PuFCvsH19/Data/U_G4NeutronHP.dat" << endl;
    gNew->SetName("gNew");
    new TCanvas();
    gPad->SetTicks(1, 1);
    gOld->SetTitle("ENDF WQ; #font[12]{E} / MeV; #sigma / b");
    SetSize(gOld);
    gOld->GetXaxis()->SetNdivisions();
    gOld->GetYaxis()->SetRangeUser(0, 2.5);
    gNew->SetLineColor(kRed);
    TLegend *l = new TLegend(0.15, 0.7, 0.4, 0.85);
    l->AddEntry(gOld, "ENDF");
    l->AddEntry(gNew, "G4NeutronHP");
    gOld->Draw();
    gNew->Draw("same");
    l->Draw();
}

TH1F* FissionFold(TFile *fHits, TGraphErrors *gSigma, string method, Int_t i = 0)
{
    char name[128] = "";
    // Get simulation population
    TH1F *hTheta = (TH1F*)fHits->Get("Source/Source_Theta");
    if (!hTheta) cout << "Could not get " << "Source/Source_Theta" << endl;
    Double_t nSim = hTheta->Integral();
    // Get neutron energy histogram
    sprintf(name, "PuFC/nEnergy/PuFC_Ekin_Ch.%i", i+1);
    TH1F *hHits = (TH1F*)fHits->Get(name);
    if (!hHits) cout << "Could not get " << name << endl;
    // ready folded: return
    if (!gSigma)
    {
        hHits->Scale(1.0 / nSim);
        return hHits;
    }

    Int_t N = hHits->GetNbinsX();
    sprintf(name, "fis%s_%i", method.c_str(), i+1);
    TH1F *hFis = new TH1F(name, "(n,f)-Wahrscheinlichkeit vs Ekin; #font[12]{E} / MeV; d#font[12]{p}(n,f)", N, 0, hHits->GetBinLowEdge(N+1));
    for (Int_t j = 1; j < N+1; j++)
    {
        Double_t hits = hHits->GetBinContent(j);
        Double_t E = hHits->GetBinCenter(j);
        Double_t CS = gSigma->Eval(E);
        hFis->SetBinContent(j, hits / F() * CS / nSim);
    }
    return hFis;
}

void VglTrackLength(Int_t ch = 0)
{
    char name[64] = "";
    TGraphErrors *gSigma = GetSigma();

    TFile *f1 = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_Hit_5E7.root");
    if (!f1) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_Hit_5E7.root" << endl;
    TGraph *g1 = new TGraph(8);
    g1->SetNameTitle("gFis_Hit_BC", "Hits, bin center; Deposit; #font[12]{p}(n,f)");

    TFile *f2 = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_TL_2E7.root");
    if (!f2) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_TL_2E7.root" << endl;
    TGraph *g2 = new TGraph(8);
    g2->SetNameTitle("gFis_TL_BC", "Track length, bin center; Deposit; #font[12]{p}(n,f)");

    TFile *f3 = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_1E7.root");
    if (!f3) cout << "Could not open " << "/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_1E7.root" << endl;
    TGraph *g3 = new TGraph(8);
    g3->SetNameTitle("gFis_TL_E", "Track length, E; Deposit; #font[12]{p}(n,f)");

    TH1F *h[3][8];
    Double_t fis0 = 0, fis1 = 0, fis2 = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        h[0][i] = FissionFold(f1, gSigma, "Hit_BC", i);
        g1->SetPoint(i, i+1, h[0][i]->Integral());
        fis0 += h[0][i]->Integral();
        h[1][i] = FissionFold(f2, gSigma, "TL_BC", i);
        g2->SetPoint(i, i+1, h[1][i]->Integral());
        fis1 += h[1][i]->Integral();
        h[2][i] = FissionFold(f3, 0, "TL_E", i);
        g3->SetPoint(i, i+1, h[2][i]->Integral());
        fis2 += h[2][i]->Integral();
        cout << fis0 << " \t" << fis2 << " \t" << h[0][i]->Integral() / h[2][i]->Integral() << endl;
    }

    TLatex *t[5];
    new TCanvas("cTL", "Track length");
    gPad->SetTicks(1, 1);
    SetSize(g1);
    g1->GetYaxis()->SetRangeUser(0, 0.005);
    g2->SetLineColorAlpha(kGreen, 0.5);
    g3->SetLineColorAlpha(kRed, 0.5);
    TLegend *l = new TLegend(0.15, 0.2, 0.6, 0.4);
    l->AddEntry(g1, "Hits, sigma(bin center)");
    l->AddEntry(g2, "Track length, sigma(bin center)");
    l->AddEntry(g3, "Track length, sigma(E)");
    g1->Draw();
    g2->Draw("same");
    g3->Draw("same");
    l->Draw();
    t[3] = new TLatex();
    sprintf(name, "1/cos(#theta): %f", fis0 / fis1);
    t[3]->DrawLatex(4, 0.003, name);
    t[4] = new TLatex();
    sprintf(name, "#sigma(E): %f", fis1 / fis2);
    t[4]->DrawLatex(4, 0.0025, name);

    // Example hostograms channel ch
    sprintf(name, "cTL_%i", ch+1);
    new TCanvas(name, "Track length");
    gPad->SetTicks(1, 1);
    gPad->SetLogy();
    SetSize(h[2][ch]);
    h[2][ch]->GetYaxis()->SetRangeUser(1E-9, 2E-5);
    h[0][ch]->SetLineColorAlpha(kBlack, 0.8);
    h[1][ch]->SetLineColorAlpha(kGreen, 0.8);
    h[2][ch]->SetLineColorAlpha(kRed, 0.8);
    h[2][ch]->Draw("hist");
    h[1][ch]->Draw("same");
    h[0][ch]->Draw("same hist");
    sprintf(name, "Kanal %i", ch+1);
    l = new TLegend(0.15, 0.65, 0.65, 0.85, name);
    l->AddEntry(g1, "Hits, sigma(bin center)");
    l->AddEntry(g2, "Track length, sigma(bin center)");
    l->AddEntry(g3, "Track length, sigma(E)");
    l->Draw();
    t[0] = new TLatex();
    sprintf(name, "#Sigma = %f", h[0][ch]->Integral());
    t[0]->DrawLatexNDC(0.75, 0.8, name);
    t[1] = new TLatex();
    t[1]->SetTextColor(kGreen);
    sprintf(name, "#Sigma = %f", h[1][ch]->Integral());
    t[1]->DrawLatexNDC(0.75, 0.7, name);
    t[2] = new TLatex();
    t[2]->SetTextColor(kRed);
    sprintf(name, "#Sigma = %f", h[2][ch]->Integral());
    t[2]->DrawLatexNDC(0.75, 0.6, name);
}

void VglSolidAngle(TFile *f, string FC = "PuFC")
{ // Simulate ideal
    char name[128] = "";
    // source spectrum
    TH1F *hE = (TH1F*)f->Get("Source/nEnergy/Source_Ekin"); if (!hE) cout << "Could not get " << "Source/nEnergy/Source_Ekin" << endl;
    Double_t SigmaSource = AvSigma(GetSigma(1), hE);
    Double_t w = SigmaSource / F();
    cout << SigmaSource << "\t" << w << endl;
    // create graphs
    TGraph *gSource = new TGraph(8);
    gSource->SetNameTitle("gSource", "Neutronen; Deposit; #font[12]{N}_{n}");
    TGraph *gFis = new TGraph(8);
    gFis->SetNameTitle("gFis", "Neutronen; Deposit; #font[12]{N}_{n}");
    Double_t n0 = 0, n1 = 0;

    for (Int_t i = 0; i < 8; i++)
    {
        // source angular
        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
        TH1F *hSource = (TH1F*)f->Get(name); if (!hSource) cout << "Could not get " << name << endl;
        // ideal fission
        sprintf(name, "%s/nEnergy/%s_Ekin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH1F *hFis = (TH1F*)f->Get(name); if (!hFis) cout << "Could not get " << name << endl;
        Double_t nQu = hSource->Integral();
        Double_t nVac = hFis->Integral() / w;
        cout << i+1 << "\t" << nQu << "\t" << nVac << "\t" << hFis->Integral() / hFis->GetEntries() << endl;
        gSource->SetPoint(i, i+1, nQu);
        gFis->SetPoint(i, i+1, nVac);
        n1 += nVac;
        n0 += nQu;
    }
    new TCanvas("cVac", "Void-Simulation");
    gPad->SetTicks(1, 1);
    SetSize(gFis);
//    gFis->GetYaxis()->SetRangeUser(0, 1);
    gFis->SetLineColor(kRed);
    TLegend *l = new TLegend(0.15, 0.7, 0.4, 0.85);
    l->AddEntry(gFis, "Void-Simulation");
    l->AddEntry(gSource, "Emissionswinkel");
    gFis->Draw();
    gSource->Draw("same");
    l->Draw();
    TLatex *t = new TLatex();
    sprintf(name, "Void Simulation: %f", n1 / n0);
    t->DrawLatexNDC(0.5, 0.4, name);
}

void VglSourceSpectrum(TFile *fFG, TFile *fTot, Int_t ch = 0)
{
    char name[64] = "";
    // Simulation population
    Double_t ScPart = 0.016; // TARGET scattered neutrons portion
    TH1F *hE1 = (TH1F*)fFG->Get("Source/nEnergy/Source_Ekin");
    if (!hE1) cout << "Could not get " << "Source/nEnergy/Source_Ekin" << endl;
    Double_t nSim1 = hE1->Integral();
    hE1->Scale(1.0 / nSim1);
    TH1F *hE2 = (TH1F*)fTot->Get("Source/nEnergy/Source_Ekin");
    if (!hE2) cout << "Could not get " << "Source/nEnergy/Source_Ekin" << endl;
    Double_t nSim2 = hE2->Integral();
    hE2->Scale(1.0 / (1.0 - ScPart) / nSim2);

    // Histograms and graphs
    TH1F *h[2][8];
    TGraph *g1 = new TGraph(8);
    g1->SetNameTitle("gFG", "(n,f)-Wahrscheinlichkeit; Deposit; d#font[12]{p}(n,f)");
    TGraph *g2 = new TGraph(8);
    g2->SetNameTitle("gTot", "(n,f)-Wahrscheinlichkeit; Deposit; d#font[12]{p}(n,f)");

    Double_t fis0 = 0, fis1 = 0;
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "PuFC/nEnergy/PuFC_Ekin_Ch.%i", i+1);
        h[0][i] = (TH1F*)fFG->Get(name);
        h[1][i] = (TH1F*)fTot->Get(name);
        if (!h[0][i] || !h[1][i]) cout << "Could not get " << name << endl;
        h[0][i]->Scale(1.0 / nSim1);
        h[1][i]->Scale(1.0 / (1.0 - ScPart) / nSim2);
        fis0 += h[0][i]->Integral();
        fis1 += h[1][i]->Integral();
        g1->SetPoint(i, i+1, h[0][i]->Integral());
        g2->SetPoint(i, i+1, h[1][i]->Integral());
    }
    new TCanvas("cS", "Source spectrum 1");
    gPad->SetTicks(1, 1);
    g1->GetYaxis()->SetRangeUser(0, 0.005);
    SetSize(g1);
    g2->SetLineColor(kRed);
    TLegend *l = new TLegend(0.15, 0.7, 0.4, 0.85, "Quellspektrum");
    l->AddEntry(g1, "FG");
    l->AddEntry(g2, "FG+BG");
    g1->Draw();
    g2->Draw("same");
    l->Draw();
    TLatex *t[3];
    t[0] = new TLatex();
    sprintf(name, "Target-Streuung: %f", fis0 / fis1);
    t[0]->DrawLatexNDC(0.5, 0.5, name);

    sprintf(name, "cS_%i", ch+1);
    new TCanvas(name, "fission neutron spectrum");
    gPad->SetTicks(1, 1);
    sprintf(name, "Spaltung Neutronenspektrum Deposit %i", ch+1);
    h[0][ch]->SetTitle(name);
    h[0][ch]->Rebin(10);
    h[1][ch]->Rebin(10);
    SetSize(h[0][ch]);
    h[1][ch]->SetLineColor(kRed);
    h[0][ch]->Draw("hist");
    h[1][ch]->Draw("same hist");
    l->Draw();
    t[1] = new TLatex();
    sprintf(name, "#Sigma = %f", h[0][ch]->Integral());
    t[1]->DrawLatexNDC(0.5, 0.6, name);
    t[2] = new TLatex();
    t[2]->SetTextColor(kRed);
    sprintf(name, "#Sigma = %f", h[1][ch]->Integral());
    t[2]->DrawLatexNDC(0.5, 0.5, name);

    new TCanvas("cE", "Source spectrum 3");
    gPad->SetTicks(1, 1);
    SetSize(hE1);
    hE2->SetLineColor(kRed);
//    h[0][ch]->Rebin(10);
//    h[1][ch]->Rebin(10);
    hE1->Draw("hist");
    hE2->Draw("same hist");
    l->Draw();
}

void DrawAngular()
{
    Double_t maxTheta = 180.0 / TMath::Pi() * atan(12.7 / 150.0);
    TFile *f = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_1E8.root"); if (!f) cout << "Could not open " << "PuFC_real_1E8.root" << endl;
    TH1F *hAng = (TH1F*)f->Get("Source/Source_Theta"); if (!hAng) cout << "Could not get " << "Source/Source_Theta" << endl;
    hAng->SetNameTitle("hAng", "Winkelverteilung Quelle; #Theta / deg; #font[12]{N}");
    Double_t N = hAng->Integral();
    TH1F *hTheo = (TH1F*)hAng->Clone();
    hTheo->SetNameTitle("hTheo", "Winkelverteilung Quelle; #Theta / deg; #font[12]{N}");
    for (Int_t j = 1; j < hTheo->GetNbinsX()+1; j++)
    {
        Double_t low = hTheo->GetBinLowEdge(j);
        Double_t up = hTheo->GetBinLowEdge(j) + hTheo->GetBinWidth(j);
        if (low > maxTheta)
            up = low;
        if (low < maxTheta && up > maxTheta)
            up = maxTheta;
        Double_t c = N * (pow(up, 2) - pow(low, 2)) / pow(maxTheta, 2);
        hTheo->SetBinContent(j, c);
    }

    new TCanvas("cTheta", "Source Angular Distribution");
    gPad->SetTicks(1, 1);
    SetSize(hAng);
    hTheo->SetLineColor(kRed);
    hAng->Draw("hist");
    hTheo->Draw("same hist");
    TLegend *l = new TLegend(0.15, 0.7, 0.4, 0.85);
    l->AddEntry(hAng, "Geant4");
    l->AddEntry(hTheo, "Erwartung");
    l->Draw();
}

void G4vsMCNP(string Setup = "real", Int_t ch = 0)
{
    char name[64] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root");
    sprintf(name, "Simulation/Geant4/PuFC_%s/ToFvsEkin/PuFC_ToFvsEkin_%s_Ch.%i", Setup.c_str(), Setup.c_str(), ch+1);
    TH2D *hG4 = (TH2D*)fAna->Get(name); if (!hG4) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/MCNP/PuFC_%s/ToFvsEkin/PuFC_ToFvsEkin_%s_Ch.%i", Setup.c_str(), Setup.c_str(), ch+1);
    TH2D *hMCNP = (TH2D*)fAna->Get(name); if (!hMCNP) cout << "Could not get " << name << endl;
    cout << "G4 " << hG4->Integral() << "  MCNP " << hMCNP->Integral() << endl;

    TH1D *hTG4 = (TH1D*)hG4->ProjectionX("tG4");
    TH1D *hTMCNP = (TH1D*)hMCNP->ProjectionX("tMCNP");

    sprintf(name, "cT_%s_%i", Setup.c_str(), ch+1);
    new TCanvas(name, "Time projection");
    gPad->SetLogy();
    hTG4->SetStats(0);
    SetSize(hTG4);
    hTMCNP->SetLineColor(kGreen);
    TLegend *l = new TLegend(0.4, 0.7, 0.6, 0.85);
    l->AddEntry(hTG4, "Geant4");
    l->AddEntry(hTMCNP, "MCNP");
    hTG4->Draw("hist");
    hTMCNP->Draw("same hist");
    l->Draw();

    TH1D *hEG4 = (TH1D*)hG4->ProjectionY("EG4");
    TH1D *hEMCNP = (TH1D*)hMCNP->ProjectionY("EMCNP");

    sprintf(name, "cE_%s_%i", Setup.c_str(), ch+1);
    new TCanvas(name, "Energy projection");
    gPad->SetLogy();
    hEG4->SetStats(0);
    SetSize(hEG4);
    hEMCNP->SetLineColor(kGreen);
    hEG4->Draw("hist");
    hEMCNP->Draw("same hist");
    l->Draw();
}

void VglSim()
{
//    oldCorrection("PuFC", 7);
//    VglSigmaU();
//    VglTrackLength();
//    TFile *fVac = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_ideal_5E7.root"); if (!fVac) cout << "Could not open " << "fVac" << endl;
//    VglSolidAngle(fVac);
//    TFile *fFG = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_FG_1E7.root"); if (!fFG) cout << "Could not open " << "fVac" << endl;
//    TFile *fTot = TFile::Open("/home/hoffma93/Programme/Geant4-Work/results/PuFC_real_1E8.root"); if (!fTot) cout << "Could not open " << "fVac" << endl;
//    VglSourceSpectrum(fFG, fTot);
//    DrawAngular();
    G4vsMCNP("real");
}

#endif
