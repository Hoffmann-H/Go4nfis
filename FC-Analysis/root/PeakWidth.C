#include "SaveToFile.C"
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

void GetPeak(Int_t i, TH1I* pH, string FC, Int_t left, Int_t right, Double_t* nif, Double_t* Dnif)
{
    Int_t l0, l1, l2, l3;
    GetLimits(i, FC, left, right, &l0, &l1, &l2, &l3);
    Double_t PeakInt = pH->Integral(l1, l2);
    Double_t UgInt = pH->Integral(l0, l3) - PeakInt;
    *nif = PeakInt - (l2-l1+1.0) / (l3-l2+l1-l0) * UgInt;
    *Dnif = sqrt( PeakInt + pow((l2-l1+1) / (l3-l2+l1-l0), 2) * UgInt );
}

void PeakWidth(string file_name, string FC, Int_t left, Int_t rStart, Int_t rStop, Int_t rStep)
{
    char name[64] = "";
    Int_t n = (rStop - rStart) / rStep;
    Double_t avNIF[n];
    Double_t D2avNIF[n];

    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s", file_name.c_str());
    TFile *f = TFile::Open(name, "UPDATE");
    TDirectory *pDir = Prepare(f, "Analysis/ToF/Gate");
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
            Int_t r = rStart + j * rStep;
            Double_t NIF, DNIF;
            GetPeak(i, pH, FC, left, r, &NIF, &DNIF);
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
        sprintf(name, "gPeak_%i", i+1);
        Save(pDir, ge[i], name);
        ge[i]->SetLineColorAlpha(i+2, 0.5);
        ge[i]->SetLineWidth(1);
        ge[i]->SetMarkerColorAlpha(i+2, 0.5);
        ge[i]->SetMarkerStyle(20);
        ge[i]->SetMarkerSize(2);
        sprintf(name, "%i", i+1);
        l->AddEntry(ge[i], name, "ep");
        mg->Add(ge[i]);
    }
    pDir = Prepare(f, "Analysis/ToF/Gate");
    Save(pDir, geAv, "gPeak_av");
    geAv->SetLineColorAlpha(1, 0.5);
    geAv->SetLineWidth(2);
    geAv->SetMarkerColorAlpha(1, 0.5);
    geAv->SetMarkerStyle(20);
    geAv->SetMarkerSize(2);
    l->AddEntry(geAv, "av", "ep");
    mg->Add(geAv);
    new TCanvas();
    mg->Draw("AP");
    l->Draw();
    f->Save();
    f->Close();
}

void GetBackground(Int_t i, TH1I* pH, string FC, Int_t left, Int_t right, Double_t* L, Double_t* DL, Double_t* R, Double_t* DR)
{
    Int_t l0, l1, l2, l3;
    GetLimits(i, FC, left, right, &l0, &l1, &l2, &l3);
    *L = pH->Integral(l0, l1 - 1) / (l1 - l0);
    *DL = sqrt(pH->Integral(l0, l1 - 1)) / (l1 - l0);
    *R = pH->Integral(l2 + 1, l3) / (l3 - l2);
    *DR = sqrt(pH->Integral(l2 + 1, l3)) / (l3 - l2);
}

void Background(string file_name, string FC, Int_t left, Int_t rStart, Int_t rStop, Int_t rStep)
{
    char name[64] = "";
    Int_t n = (rStop - rStart) / rStep;
    Double_t avNIF[n];
    Double_t D2avNIF[n];

    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s", file_name.c_str());
    TFile *f = TFile::Open(name, "UPDATE");
    TDirectory *pDir = Prepare(f, "Analysis/ToF/BG");
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
        sprintf(name, "gBgL_%i", i+1);
        Save(pDir, gL[i], name);
        sprintf(name, "gBgR_%i", i+1);
        Save(pDir, gR[i], name);
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
    f->Save();
    f->Close();
}

void PeakWidth()
{
    PeakWidth("NIF.root", "PuFC", 15, 0, 100, 5);
    Background("NIF.root", "PuFC", 15, 0, 50, 2);
    PeakWidth("UFC_NIF.root", "UFC", 15, 0, 100, 5);
    Background("UFC_NIF.root", "UFC", 15, 0, 100, 5);
}
