#include "SaveToFile.C"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "FC.C"

void GetLimits(Int_t i, string FC, Int_t left, Int_t right, Int_t* l0, Int_t* l1, Int_t* l2, Int_t* l3)
{
    *l0 = 42;
    *l3 = 402;
    if (!strcmp(FC.c_str(), "PuFC")) // if PuFC
    {
        Int_t m[] = {140, 128, 126, 129, 129, 128, 127, 127};
        *l1 = m[i] - left;
        *l2 = m[i] + right;
    } else {
        Int_t m[] = {273, 269, 260, 265, 265, 264, 263, 263};
        *l1 = m[i] - left;
        *l2 = m[i] + right;
    }
}

void GetPeak(Int_t i, TH1I* pH, string FC, Int_t left, Int_t right, Int_t tail, Double_t* nif, Double_t* Dnif)
{
    Int_t l0, la, l1, l2, lb, l3;
//    GetLimits(i, FC, left, right, &l0, &l1, &l2, &l3);
    l0 = Gate_0(i, FC);
    la = Gate_a(i, FC, left);
    l1 = Gate_1(i, FC, left);
    l2 = Gate_2(i, FC, right);
    lb = Gate_b(i, FC, tail);
    l3 = Gate_3(i, FC);
    Double_t PeakInt = pH->Integral(l1, l2);
    Double_t UgInt = pH->Integral(l0, la - 1) + pH->Integral(lb + 1, l3);
    *nif = PeakInt - (l2-l1+1.0) / (l3-lb+la-l0) * UgInt;
    *Dnif = sqrt( PeakInt + pow((l2-l1+1) / (l3-lb+la-l0), 2) * UgInt );
}

void PeakWidth(string file_name, string subfolder, string FC, Int_t lStart, Int_t lStop, Int_t rStart, Int_t rStop, Int_t n)
{
    char name[64] = "";
    Double_t avNIF[n];
    Double_t D2avNIF[n];

    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s", file_name.c_str());
    TFile *f = TFile::Open(name, "READ");
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *ge[8];
    TGraphErrors *geAv = new TGraphErrors(n);
    TMultiGraph *mg = new TMultiGraph("mgPeak", "ToF Peak; Gate[ns]; Peak");
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.9, 0.5, 1.0, 1.0, name);

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(name);
        ge[i] = new TGraphErrors(n);
        for (Int_t j = 0; j < n; j++)
        {
            if (i == 0)
            {
                avNIF[j] = 0;
                D2avNIF[j] = 0;
            }
            Int_t r = rStart + j * (rStop - rStart) / (n - 1.0);
            Int_t l = lStart + j * (lStop - lStart) / (n - 1.0);
            Double_t NIF, DNIF;
            GetPeak(i, pH, FC, l, r, rStop, &NIF, &DNIF);
            ge[i]->SetPoint(j, r, NIF);
            ge[i]->SetPointError(j, 0, DNIF);
//            if (i != 4 && i != 7 && i != 3) {
                avNIF[j] += NIF / 8.0;
                D2avNIF[j] += pow(DNIF / 8.0, 2);
//            }
            if (i == 7)
            { // After 8th deposit, fill average
                geAv->SetPoint(j, r, avNIF[j]);
                geAv->SetPointError(j, 0, sqrt(D2avNIF[j]));
            }
        }
        Save(fAna, FC+"/ToF/Gate/Peak", ge[i], "gPeak_"+to_string(i+1));
        ge[i]->SetLineColorAlpha(i+2, 0.5);
        ge[i]->SetLineWidth(1);
        ge[i]->SetMarkerColorAlpha(i+2, 0.5);
        ge[i]->SetMarkerStyle(20);
        ge[i]->SetMarkerSize(2);
        sprintf(name, "%i", i+1);
        l->AddEntry(ge[i], name, "ep");
        mg->Add(ge[i]);
    }
    Save(fAna, FC+"/ToF/Gate/Peak/"+subfolder, geAv, "gPeak_av");
    geAv->SetLineColorAlpha(1, 0.5);
    geAv->SetLineWidth(2);
    geAv->SetMarkerColorAlpha(1, 0.5);
    geAv->SetMarkerStyle(20);
    geAv->SetMarkerSize(2);
    l->AddEntry(geAv, "av", "ep");
    mg->Add(geAv);
//    new TCanvas();
//    mg->Draw("AP");
//    l->Draw();
    fAna->Save();
    fAna->Close();
    f->Close();
}

void SimWidth(string Simulation, string FC, Int_t left, Int_t rStart, Int_t rStop, Int_t rStep)
{
    char name[64] = "";
    Int_t n = (rStop - rStart) / rStep;
    Double_t avNIF[n];
    Double_t D2avNIF[n];

//    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s", file_name.c_str());
//    TFile *f = TFile::Open(name, "READ");
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *ge[8];
    TGraphErrors *geAv = new TGraphErrors(n);
    TMultiGraph *mg = new TMultiGraph("mgSim", "ToF Peak; Gate[ns]; Peak");
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.9, 0.5, 1.0, 1.0, name);

    for (Int_t i = 0; i < 8; i++)
    {
//        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
//        TH1I* pExp = (TH1I*)f->Get(name);
        sprintf(name, "Simulation/%s/%s_real/FitToF/%s_FitT_real_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH1D *pH = (TH1D*)fAna->Get(name); if (!pH) cout << "Could not get " << name << endl;
        ge[i] = new TGraphErrors(n);
        for (Int_t j = 0; j < n; j++)
        {
            if (i == 0)
            {
                avNIF[j] = 0;
                D2avNIF[j] = 0;
            }
            Int_t r = rStart + j * rStep;
            Double_t NIF, DNIF = 0;
            Double_t bg = pH->GetBinContent(pH->GetNbinsX());
            Double_t xc = pH->GetBinCenter(pH->GetMaximumBin());
            Double_t x0 = xc - left; // - width;
            Double_t x1 = xc + r; // + width;
            Int_t bin0 = pH->FindBin(x0);
            Int_t bin1 = pH->FindBin(x1);
            Double_t w0 = (pH->GetBinLowEdge(bin0 + 1) - x0) / pH->GetBinWidth(bin0); // Constant bin interpolation
            Double_t w1 = (x1 - pH->GetBinLowEdge(bin1)) / pH->GetBinWidth(bin1);
            Double_t GatedIntegral = w0 * pH->GetBinContent(bin0) + pH->Integral(bin0 + 1, bin1 - 1) + w1 * pH->GetBinContent(bin1) - (w0 + bin1 - bin0 + 1 + w1) * bg;
            NIF = GatedIntegral * pH->GetBinWidth(1);

            ge[i]->SetPoint(j, left + r, NIF);
            ge[i]->SetPointError(j, 0, DNIF);
            avNIF[j] += NIF / 8.0;
            D2avNIF[j] += pow(DNIF / 8.0, 2);
            if (i == 7)
            { // After 8th deposit, fill average
                geAv->SetPoint(j, left + r, avNIF[j]);
                geAv->SetPointError(j, 0, sqrt(D2avNIF[j]));
            }
        }
        Save(fAna, FC+"/ToF/Gate/Peak", ge[i], "gSim_"+to_string(i+1));
        ge[i]->SetLineColorAlpha(i+2, 0.5);
        ge[i]->SetLineWidth(1);
        ge[i]->SetMarkerColorAlpha(i+2, 0.5);
        ge[i]->SetMarkerStyle(20);
        ge[i]->SetMarkerSize(2);
        sprintf(name, "%i", i+1);
        l->AddEntry(ge[i], name, "ep");
        mg->Add(ge[i]);
    }
    Save(fAna, FC+"/ToF/Gate/Peak", geAv, "gSim_av");
    geAv->SetLineColorAlpha(1, 0.5);
    geAv->SetLineWidth(2);
    geAv->SetMarkerColorAlpha(1, 0.5);
    geAv->SetMarkerStyle(20);
    geAv->SetMarkerSize(2);
    l->AddEntry(geAv, "av", "ep");
    mg->Add(geAv);
//    new TCanvas();
//    mg->Draw("AP");
//    l->Draw();
    fAna->Save();
    fAna->Close();
//    f->Close();
}

void GetBackground(Int_t i, TH1I* pH, string FC, Int_t left, Int_t right, Double_t* L, Double_t* DL, Double_t* R, Double_t* DR)
{
    Int_t l0, l1, l2, l3;
    GetLimits(i, FC, left, right, &l0, &l1, &l2, &l3);
    char name[64] = "";
    sprintf(name, "fL_%s_%i", FC.c_str(), i+1);
    TF1 *fL = new TF1(name, "[0]", pH->GetBinLowEdge(l0), pH->GetBinLowEdge(l1));
    pH->Fit(fL, "LR0Q");
    *L = fL->GetParameter(0);
    *DL = fL->GetParError(0);
    sprintf(name, "fR_%s_%i", FC.c_str(), i+1);
    TF1 *fR = new TF1(name, "[0]", pH->GetBinLowEdge(l2+1), pH->GetBinLowEdge(l3+1));
    pH->Fit(fR, "LR0Q");
    *R = fR->GetParameter(0);
    *DR = fR->GetParError(0);
    cout << TMath::Prob(fR->GetChisquare(), fR->GetNDF()) << "\t";
//    *L = pH->Integral(l0, l1 - 1) / (l1 - l0);
//    *DL = sqrt(pH->Integral(l0, l1 - 1)) / (l1 - l0);
//    *R = pH->Integral(l2 + 1, l3) / (l3 - l2);
//    *DR = sqrt(pH->Integral(l2 + 1, l3)) / (l3 - l2);
}

void Background(string file_name, string FC, Int_t left, Int_t rStart, Int_t rStop, Int_t rStep)
{
    char name[64] = "";
    Int_t n = (rStop - rStart) / rStep;
    Double_t avNIF[n];
    Double_t D2avNIF[n];

    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s", file_name.c_str());
    TFile *f = TFile::Open(name, "READ");
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gL[8];
    TGraphErrors *gR[8];
    TMultiGraph *mg = new TMultiGraph("mgBg", "ToF Background; Gate[ns]; Background/ns[1/s]");
    sprintf(name, "%s Deposit", FC.c_str());
    TLegend *l = new TLegend(0.9, 0.5, 1.0, 1.0, name);

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*)f->Get(name);
        gL[i] = new TGraphErrors(n);
        gR[i] = new TGraphErrors(n);
        cout << i+1 << "\t";
        for (Int_t j = 0; j < n; j++)
        {
            Int_t r = rStart + j * rStep;
            Double_t L, DL, R, DR;
            GetBackground(i, pH, FC, left, r, &L, &DL, &R, &DR);
            gL[i]->SetPoint(j, left + r, L);
            gL[i]->SetPointError(j, 0, DL);
            gR[i]->SetPoint(j, left + r, R);
            gR[i]->SetPointError(j, 0, DR);
//            cout << r << "  " << L << "+-" << DL << " " << R << "+-" << DR << endl;
        }
        cout << endl;
        Save(fAna, FC+"/ToF/Gate/Background/Left", gL[i], "gBgL_"+to_string(i+1));
        Save(fAna, FC+"/ToF/Gate/Background/Right", gR[i], "gBgR_"+to_string(i+1));
        gR[i]->SetLineColorAlpha(i+2, 0.5);
        gR[i]->SetLineWidth(1);
        gR[i]->SetMarkerColorAlpha(i+2, 0.5);
        gR[i]->SetMarkerStyle(20);
        gR[i]->SetMarkerSize(2);
        sprintf(name, "%i", i+1);
        l->AddEntry(gR[i], name, "ep");
        mg->Add(gR[i]);
    }
    new TCanvas();
    mg->Draw("AP");
    l->Draw();
    fAna->Save();
    fAna->Close();
    f->Close();
}

void PeakWidth()
{
//    PeakWidth("NIF.root", "Left", "PuFC", 0, 45, 35, 35, 16);
    PeakWidth("NIF.root", "Right", "PuFC", 15, 15, 5, 100, 20);
//    SimWidth("Geant4", "PuFC", 15, 5, 200, 5);
//    Background("NIF.root", "PuFC", 15, 0, 50, 2);
//    PeakWidth("UFC_NIF.root", "Left", "UFC", 0, 30, 40, 40, 16);
    PeakWidth("UFC_NIF.root", "Right", "UFC", 15, 15, 5, 100, 20);
//    SimWidth("Geant4", "UFC", 15, 5, 100, 5);
//    Background("UFC_SB.root", "UFC", 15, 0, 100, 5);
//    Background("NIF.root", "PuFC", 40, 0, 30, 1);
//    PeakWidth("UFC_NIF.root", "UFC", 40, 0, 30, 1);
}
