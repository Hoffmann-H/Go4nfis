#ifndef STABILITY_H
#define STABILITY_H
#include "FC.C"
#include "Runs.C"
#include "NeutronField.C"

TH1D* RunBackground(string FC_Setup, Int_t nRun, Int_t ch, Int_t r = 1)
{
    string FC = FC_Setup[0] == 'U' ? "UFC" : "PuFC";
    char name[64] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", GetRunName(FC_Setup, nRun).c_str());
    TFile *f = TFile::Open(name); if (!f) cout << "Could not open " << name << endl;

    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", ch+1);
    TH2I *h2 = (TH2I*)f->Get(name); if (!h2) cout << "Could not get " << name << endl;
    Int_t l[] = {Gate_0(ch, FC), Gate_a(ch, FC), Gate_b(ch, FC), Gate_3(ch, FC)};
    TH1D *h1T = (TH1D*) h2->ProjectionY("_pTotal");
    TH1D *h1C = (TH1D*) h2->ProjectionY("_pConst", l[0], l[3]);
    TH1D *h1P = (TH1D*) h2->ProjectionY("_pPeak", l[1], l[2]);
    h1T->Add(h1C, -1);
    h1T->Sumw2();
    // h1T: Flanken außerhalb l0, l3
    h1C->Add(h1P, -1);
    h1C->Sumw2();
    // h1C: Untergrund-Bereich
    h1T->Add(h1C, (l[3] - l[0] + 1.0) / (l[1] - l[0] + l[3] - l[2]));
    h1T->Rebin(r);

    // reject incomplete time bins
    TH1D *h1treal = (TH1D*) f->Get("Histograms/Raw/Scaler/Rates/H1RawRate_48");
    if (!h1treal) cout << "Could not get " << "Real time" << endl;
    Double_t fRes = h1T->GetBinWidth(1) / h1treal->GetBinWidth(1); // time resolution factor
    for (Int_t bin = 1; bin < h1T->GetNbinsX() + 1; bin++)
    {
        Int_t tBin = h1treal->FindBin(h1T->GetBinLowEdge(bin));
        if (h1treal->Integral(tBin, tBin + fRes - 1) < h1T->GetBinWidth(1) /*&& h1T->GetBinContent(bin) > 0*/)
        {
//            cout << "reject " << bin << ", " << tBin << ": " << h1treal->Integral(tBin, tBin + fRes - 1) << " < " << h1treal->GetBinWidth(1) * fRes << endl;
            h1T->SetBinContent(bin, 0);
            h1T->SetBinError(bin, 0);
        }
    }
    return h1T;
}

TH1D* ConstantBackground(string FC_Setup, Int_t ch, Int_t r = 1)
{
    char name[128] = "";
    // Get first Run's yield
    TH1D *h = RunBackground(FC_Setup, 0, ch, r);
    sprintf(name, "%s_%i", FC_Setup.c_str(), ch+1);
    h->SetNameTitle(name, "Background rate; #font[12]{t}; #font[12]{#dot{C}}_{(SF)} / s^{-1}");

    for (Int_t j = 1; j < GetnRuns(FC_Setup); j++)
    {
        // Get Run's yield
        TH1D *h1R = RunBackground(FC_Setup, j, ch, r);
        h->Add(h1R);
    }
    if (FC_Setup[0] != 'U')
    { // PuFC: Add Spontaneous fission run
        // Open file, get 2D histogram
        sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/SF.root");
        TFile *f = TFile::Open(name); if (!f) cout << "Could not open " << name << endl;
        sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", ch+1);
        TH2I *h2SF = (TH2I*) f->Get(name); if (!h2SF) cout << "Could not get " << name << endl;

        // time projection
        TH1D *h1T = /*(TH1D*) */h2SF->ProjectionY();
        h1T->Sumw2();
        h1T->Rebin(r);

        // reject incomplete time bins
        TH1D *h1treal = (TH1D*) f->Get("Histograms/Raw/Scaler/Rates/H1RawRate_48");
        if (!h1treal) cout << "Could not get " << "Real time" << endl;
        Double_t fRes = h1T->GetBinWidth(1) / h1treal->GetBinWidth(1); // time resolution factor
        for (Int_t bin = 1; bin < h1T->GetNbinsX() + 1; bin++)
        {
            Int_t tBin = h1treal->FindBin(h1T->GetBinLowEdge(bin));
            if (h1treal->Integral(tBin, tBin + fRes - 1) < h1T->GetBinWidth(1) /*&& h1T->GetBinContent(bin) > 0*/)
            {
                h1T->SetBinContent(bin, 0);
                h1T->SetBinError(bin, 0);
            }
        }
        h->Add(h1T);
    }
    h->Scale(1.0 / 60.0 / r); // r minutes -> 1 second
    return h;
}

TH1D* RunSignal(string FC_Setup, Int_t nRun, Int_t ch, Int_t r = 1)
{
    string FC = FC_Setup[0] == 'U' ? "UFC" : "PuFC";
    char name[64] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", GetRunName(FC_Setup, nRun).c_str());
    TFile *f = TFile::Open(name); if (!f) cout << "Could not open " << name << endl;

    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", ch+1);
    TH2I *h2 = (TH2I*)f->Get(name); if (!h2) cout << "Could not get " << name << endl;
    Int_t l[] = {Gate_0(ch, FC), Gate_a(ch, FC), Gate_b(ch, FC), Gate_3(ch, FC)};
    TH1D *h1C = (TH1D*) h2->ProjectionY("_pConst", l[0], l[3]);
    TH1D *h1P = (TH1D*) h2->ProjectionY("_pPeak", l[1], l[2]);
    h1P->Sumw2();
    // h1P: Peak-Bereich
    h1C->Add(h1P, -1);
    h1C->Sumw2();
    // h1C: Untergrund-Bereich
    h1P->Add(h1C, - (l[2] - l[1] + 1.0) / (l[1] - l[0] + l[3] - l[2]));
    h1P->Rebin(r);

    // reject incomplete time bins
    TH1D *h1treal = (TH1D*) f->Get("Histograms/Raw/Scaler/Rates/H1RawRate_48");
    if (!h1treal) cout << "Could not get " << "Real time" << endl;
    Double_t fRes = h1P->GetBinWidth(1) / h1treal->GetBinWidth(1); // time resolution factor
    for (Int_t bin = 1; bin < h1P->GetNbinsX() + 1; bin++)
    {
        Int_t tBin = h1treal->FindBin(h1P->GetBinLowEdge(bin));
        if (h1treal->Integral(tBin, tBin + fRes - 1) < h1P->GetBinWidth(1) /*&& h1P->GetBinContent(bin) > 0*/)
        {
            h1P->SetBinContent(bin, 0);
            h1P->SetBinError(bin, 0);
        }
    }
    // divide by Monitor counts to get Yield
    Double_t Monitor, DeltaRel, tReal;
    GetRun(GetRunName(FC_Setup, nRun), Monitor, DeltaRel, tReal);
//    cout << "Run " << GetRunName(FC_Setup, nRun) << " " << h1P->GetMaximum() << " " << Monitor / tReal * h1P->GetBinWidth(1) << endl;
    h1P->Scale(1.0 / Monitor * tReal / h1P->GetBinWidth(1));
    return h1P;
}

TH1D* ConstantSignal(string FC_Setup, Int_t ch, Int_t r = 1)
{ // r: Rebin, in min
    char name[128] = "";
    // Get first Run's yield
    TH1D *h = RunSignal(FC_Setup, 0, ch, r);
    sprintf(name, "%s_%i", FC_Setup.c_str(), ch+1);
    h->SetNameTitle(name, "; #font[12]{t}; #font[12]{C}_{(n,f)} / #font[12]{NM}");

    for (Int_t j = 1; j < GetnRuns(FC_Setup); j++)
    {
        // Get Run's yield
        TH1D *h1R = RunSignal(FC_Setup, j, ch, r);
        h->Add(h1R);
    }
    return h;
}

TF1* Statistics(TGraphErrors *g, Double_t start = 0, Double_t stop = 0)
{
    char name[64] = "";
    sprintf(name, "%s_fit", g->GetName());
    TF1 *fFit = new TF1(name, "pol0", start, stop ? stop : g->GetN() + 0.5);
    g->Fit(name, "R0Q");
    cout << fFit->GetParameter(0) << " +- " << fFit->GetParError(0) << endl << fFit->GetChisquare() << " / " << fFit->GetNDF() << " = " << fFit->GetChisquare() / fFit->GetNDF() << endl;
    return fFit;
}

void IndFisYield(TFile *fAna, string RunName, Int_t ch, Double_t &Y, Double_t &DY)
{
    string FC = RunName[0] == 'U' ? "UFC" : "PuFC";
    char name[128] = "";
    Double_t monitor, delta_rel, t_real;
    Double_t x, IndFis, DIndFis;

    // Get Run's monitor counts
    GetRun(RunName, monitor, delta_rel, t_real);
    // Get Run's (n,f)-Signal
    sprintf(name, "%s/ToF/Signal/%s/InducedFission", FC.c_str(), RunName.c_str());
    TGraphErrors *geIndFis = (TGraphErrors*) fAna->Get(name); if (!geIndFis) cout << "Could not get " << name << endl;
    geIndFis->GetPoint(ch, x, IndFis);
    DIndFis = geIndFis->GetErrorY(ch);
    // Divide (n,f) / monitor
    Y = IndFis / monitor;
    DY = TMath::Abs(Y) * sqrt( pow(DIndFis / IndFis, 2) + pow(delta_rel, 2) );
    cout << Y << " +- " << DY << endl;
    return;
}

void IndFisStabilityPu(TFile *fAna)
{
    char name[64] = "";
//    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    for (Int_t i = 0; i < 8; i++)
    {
        TGraphErrors *geCnf = new TGraphErrors(7);
        sprintf(name, "PuFC_nf_%i", i+1);
        geCnf->SetName(name);
        sprintf(name, "PuFC Ch. %i; Run nr.; #font[12]{C}_{(n,f)} / #font[12]{NM}", i+1);
        geCnf->SetTitle(name);
        Double_t Y, DY;
        IndFisYield(fAna, "PuFC_FG_MS4", i, Y, DY);
        geCnf->SetPoint(0, 1, Y);
        geCnf->SetPointError(0, 0, DY);
        IndFisYield(fAna, "PuFC_FG_MS5", i, Y, DY);
        geCnf->SetPoint(1, 2, Y);
        geCnf->SetPointError(1, 0, DY);
        IndFisYield(fAna, "PuFC_FG_MS6", i, Y, DY);
        geCnf->SetPoint(2, 3, Y);
        geCnf->SetPointError(2, 0, DY);
        IndFisYield(fAna, "PuFC_FG_MS7", i, Y, DY);
        geCnf->SetPoint(3, 4, Y);
        geCnf->SetPointError(3, 0, DY);
        IndFisYield(fAna, "PuFC_BG_MS9", i, Y, DY);
        geCnf->SetPoint(4, 5, Y);
        geCnf->SetPointError(4, 0, DY);
        IndFisYield(fAna, "PuFC_BG_MS10", i, Y, DY);
        geCnf->SetPoint(5, 6, Y);
        geCnf->SetPointError(5, 0, DY);
        IndFisYield(fAna, "PuFC_BG_MS11", i, Y, DY);
        geCnf->SetPoint(6, 7, Y);
        geCnf->SetPointError(6, 0, DY);
        Save(fAna, "PuFC/Stability", geCnf);

        TF1 *fFG = Statistics(geCnf, 0.0, 4.5);
        Save(fAna, "PuFC/Stability/FG", fFG);
        TF1 *fBG = Statistics(geCnf, 4.5);
        Save(fAna, "PuFC/Stability/BG", fBG);
    }
}

void IndFisStabilityU(TFile *fAna)
{
    char name[64] = "";
//    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    for (Int_t i = 0; i < 8; i++)
    {
        TGraphErrors *geCnf = new TGraphErrors(7);
        sprintf(name, "UFC_nf_%i", i+1);
        geCnf->SetName(name);
        sprintf(name, "UFC Ch. %i; Run nr.; #font[12]{C}_{(n,f)} / #font[12]{NM}", i+1);
        geCnf->SetTitle(name);
        Double_t Y, DY;
        IndFisYield(fAna, "UFC_FG_MS20_2", i, Y, DY);
        geCnf->SetPoint(0, 1, Y);
        geCnf->SetPointError(0, 0, DY);
        IndFisYield(fAna, "UFC_FG_MS20_3", i, Y, DY);
        geCnf->SetPoint(1, 2, Y);
        geCnf->SetPointError(1, 0, DY);
        IndFisYield(fAna, "UFC_FG_MS20_4", i, Y, DY);
        geCnf->SetPoint(2, 3, Y);
        geCnf->SetPointError(2, 0, DY);
        IndFisYield(fAna, "UFC_FG_MS21_2", i, Y, DY);
        geCnf->SetPoint(3, 4, Y);
        geCnf->SetPointError(3, 0, DY);
        IndFisYield(fAna, "UFC_FG_MS21_3", i, Y, DY);
        geCnf->SetPoint(4, 5, Y);
        geCnf->SetPointError(4, 0, DY);
        IndFisYield(fAna, "UFC_BG_MS20_5", i, Y, DY);
        geCnf->SetPoint(5, 6, Y);
        geCnf->SetPointError(5, 0, DY);
        IndFisYield(fAna, "UFC_BG_MS21_4", i, Y, DY);
        geCnf->SetPoint(6, 7, Y);
        geCnf->SetPointError(6, 0, DY);
        Save(fAna, "UFC/Stability", geCnf);

        TF1 *fFG = Statistics(geCnf, 0.0, 5.5);
        Save(fAna, "UFC/Stability/FG", fFG);
        TF1 *fBG = Statistics(geCnf, 5.5);
        Save(fAna, "UFC/Stability/BG", fBG);
    }
}

TH1D *hFit;
Bool_t reject;
Double_t fconst(Double_t *x, Double_t *par)
{
    if (reject)
    {
        if (hFit->GetBinContent(hFit->FindBin(x[0])) != 0)
            return par[0];
        else {
            TF1::RejectPoint();
            return 0;
        }
    } else
        return par[0];
}

TF1* Statistics(TH1D* h)
{ // Test hypothesis: Constant mean
    Double_t N = 0, Mean = 0;
    for (Int_t j = 1; j < h->GetNbinsX(); j++)
    {
        if (h->GetBinContent(j) == 0)
            continue;
        N++; Mean += h->GetBinContent(j);
    }
    Mean /= N;
    Double_t Variance = 0, Chi2 = 0;
    for (Int_t j = 1; j < h->GetNbinsX(); j++)
    {
        if (h->GetBinContent(j) == 0)
            continue;
        Variance += pow(h->GetBinContent(j) - Mean, 2) / (N - 1.0);
        Chi2 += pow(h->GetBinContent(j) - Mean, 2) / h->GetBinError(j);
    }
    cout << "N " << N << endl << "Mean " << Mean << endl << "Std.dev " << sqrt(Variance) << endl << "red.chi^2 " << Chi2 / (N - 1.0) << endl;

    reject = kTRUE;
    hFit = (TH1D*)h->Clone();
    char name[64] = "";
    sprintf(name, "%s_fit", h->GetName());
    TF1 *fFit = new TF1(name, fconst, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 1);
    h->Fit(name, "LR0Q");
    fFit->SetNpx(1000);
    reject = kFALSE;

//    h->Draw();
//    fFit->Draw("same");
    cout << fFit->GetParameter(0) << " +- " << fFit->GetParError(0) << endl << fFit->GetChisquare() << " / " << fFit->GetNDF() << " = " << fFit->GetChisquare() / fFit->GetNDF() << endl;
    return fFit;
}

void Stability()
{
//    TFile *fNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
//    if (!fNIF) cout << "Could not open " << "/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root" << endl;
//    TFile *fUNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");
//    if (!fUNIF) cout << "Could not open " << "/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root" << endl;
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    for (Int_t i = 0; i < 8; i++)
    {
/*        TGraphErrors* gUFC_FG = SignalStability(fAna, "UFC_FG", i);
        TF1 *fUFC_FG = Statistics(gUFC_FG);
        Save(fAna, "UFC/Stability/FG", gUFC_FG);
        Save(fAna, "UFC/Stability/FG", fUFC_FG);
        TGraphErrors* gUFC_BG = SignalStability(fAna, "UFC_BG", i);
        TF1 *fUFC_BG = Statistics(gUFC_BG);
        Save(fAna, "UFC/Stability/BG", gUFC_BG);
        Save(fAna, "UFC/Stability/BG", fUFC_BG);
        TGraphErrors* gPuFC_FG = SignalStability(fAna, "PuFC_FG", i);
        TF1 *fPuFC_FG = Statistics(gPuFC_FG);
        Save(fAna, "PuFC/Stability/FG", gPuFC_FG);
        Save(fAna, "PuFC/Stability/FG", fPuFC_FG);
        TGraphErrors* gPuFC_BG = SignalStability(fAna, "PuFC_BG", i);
        TF1 *fPuFC_BG = Statistics(gPuFC_BG);
        Save(fAna, "PuFC/Stability/BG", gPuFC_BG);
        Save(fAna, "PuFC/Stability/BG", fPuFC_BG);
        TGraphErrors* gPuFC = SignalStability(fAna, "PuFC", i);
        Save(fAna, "PuFC/Stability/FG+BG", gPuFC);*/

        TH1D *hUFC_C = ConstantBackground("UFC", i, 60);
        TF1 *fUFC_C = Statistics(hUFC_C);
        Save(fAna, "UFC/Stability/SF", hUFC_C);
        Save(fAna, "UFC/Stability/SF", fUFC_C);
        TH1D *hPuFC_C = ConstantBackground("PuFC", i, 60);
        TF1 *fPuFC_C = Statistics(hPuFC_C);
        Save(fAna, "PuFC/Stability/SF", hPuFC_C);
        Save(fAna, "PuFC/Stability/SF", fPuFC_C);
    }
//    fAna->Save();
//    fAna->Close();
    IndFisStabilityPu(fAna);
    IndFisStabilityU(fAna);
    fAna->Save();
    fAna->Close();
}

#endif
