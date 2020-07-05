//// See also: ~/FC-Analysis/2014.05/src/SponFis.cpp
//// Use results in: TnfisFCAnalysis::MakeConditions(),
////                 PuFC::HardCodedThresholds(), UFC::HardCodedThresholds()

#include "SaveToFile.C"
#include "TDirectory.h"
#include "TFile.h"

void FitLimits(Int_t i, Bool_t PuFC = 1)
{
    Int_t xmin = 800, xmax = 1300;
    if (PuFC)
    {
        xmin = 750, xmax = 1300;
    }
    else
    {
        switch (i) {
        case 7:
            xmin = 50, xmax = 250;
            break;
        case 3:
        case 2:
        case 0:
            xmin = 150, xmax = 350;
            break;
        case 1:
            xmin = 350, xmax = 800;
            break;
        default:
            xmin = 100, xmax = 300;
            break;
        }
    }
    cout << "Limits " << xmin << ", " << xmax << endl;
}

TF1* FindMin(TH1I *pH, Int_t xmin, Int_t xmax)
{
    TF1 *pPoly = new TF1("pol4", "pol4", 0, 4096);
    pPoly->SetRange(xmin, xmax);
//    TCanvas *c1 = new TCanvas();
//    c1->cd(1); c1->cd(1)->SetLogy();
//    pH->Draw();
    pH->Fit("pol4", "RFQ0");
//    Double_t xlow = pPoly->GetMinimumX(xmin, xmax);
//    cout << " The local minimum is at: " << xlow << endl;
    return pPoly;
}

void FindPuMinima()
{
    char name[128] = "";
    sprintf(name, "%sNIF.root", hist_data_path);
    TFile *fPuNIF = TFile::Open(name, "READ");

    sprintf(name, "%sSB.root", hist_data_path);
    TFile *fPuSB = TFile::Open(name, "READ");

    sprintf(name, "%sSF.root", hist_data_path);
    TFile *fPuSF = TFile::Open(name, "READ");

    TFile *fAna = TFile::Open(results_file, "UPDATE");

    Int_t x[] = {950, 900, 950, 900, 1050, 950, 950, 900}; // rough minimum positions
    Int_t xmin, xmax;
    Double_t mFG, mBG, mSF, mAv, mSum;
    TGraph *gFG = new TGraph(8);
    TGraph *gBG = new TGraph(8);
    TGraph *gSF = new TGraph(8);
    TGraph *gAv = new TGraph(8);
    TGraph *gSum = new TGraph(8);
    gFG->SetName("gFG");
    gBG->SetName("gBG");
    gSF->SetName("gSF");
    gAv->SetName("gAv");
    gSum->SetName("gSum");

    cout << endl << "PuFC QDC min" << endl
         << "Ch   FG fit   BG fit   SF fit   average   sum fit" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Open histograms and fit minima
        sprintf(name, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i+1);
        xmin = x[i] - 200;
        xmax = x[i] + 300;
        TH1I *hFG = (TH1I*) fPuNIF->Get(name);
        TF1 *fFG = FindMin(hFG, xmin, xmax);
        mFG = fFG->GetMinimumX(xmin, xmax);
        TH1I *hBG = (TH1I*) fPuSB->Get(name);
        TF1 *fBG = FindMin(hBG, xmin, xmax);
        mBG = fBG->GetMinimumX(xmin, xmax);
        TH1I *hSF = (TH1I*) fPuSF->Get(name);
        TF1 *fSF = FindMin(hSF, xmin, xmax);
        mSF = fSF->GetMinimumX(xmin, xmax);

        // Save fit functions
        Save(fAna, "PuFC/QDC/Fit/FG", fFG, "pol4_FG_"+to_string(i+1));
        Save(fAna, "PuFC/QDC/Fit/BG", fBG, "pol4_BG_"+to_string(i+1));
        Save(fAna, "PuFC/QDC/Fit/SF", fSF, "pol4_SF_"+to_string(i+1));

        // Calculate average
        Int_t intFG = hFG->Integral(hFG->FindBin(mFG), 4096);
        Int_t intBG = hBG->Integral(hBG->FindBin(mBG), 4096);
        Int_t intSF = hSF->Integral(hSF->FindBin(mSF), 4096);
        mAv = (intFG * mFG + intBG * mBG + intSF * mSF) / (intFG + intBG + intSF);

        // Fit sum of all runs
        hSF->Add(hFG, 1);
        hSF->Add(hBG, 1);
        TF1 *fSum = FindMin(hSF, xmin, xmax);
        mSum = fSum->GetMinimumX(xmin, xmax);

        // Save as graphs
        gFG->SetPoint(i, i+1, mFG);
        gBG->SetPoint(i, i+1, mBG);
        gSF->SetPoint(i, i+1, mSF);
        gAv->SetPoint(i, i+1, mAv);
        gSum->SetPoint(i, i+1, mSum);

        cout << " " << i+1 << "   " << mFG << "   " << mBG << "   " << mSF << "   "
             << mAv << "   " << mSum << endl;
    }
//    TMultiGraph *mg = new TMultiGraph("mg", "QDC pol4 fit minima; Deposit; QDC channel");
//    gFG->SetTitle("QDC pol4 fit minima");
//    gBG->SetTitle("QDC pol4 fit minima");
//    gSF->SetTitle("QDC pol4 fit minima");
//    gAv->SetTitle("QDC average fit minima");
//    gSum->SetTitle("QDC sum fit minima");
    Save(fAna, "PuFC/QDC", gFG, "Minima_FG");
    Save(fAna, "PuFC/QDC", gBG, "Minima_BG");
    Save(fAna, "PuFC/QDC", gSF, "Minima_SF");
    Save(fAna, "PuFC/QDC", gAv, "Average");
    Save(fAna, "PuFC/QDC", gSum, "Sum");

    fAna->Save();
    fAna->Close();
    fPuNIF->Close();
    fPuSB->Close();
    fPuSF->Close();
}

void FindUMinima()
{
    char name[128] = "";
    sprintf(name, "%sUFC_NIF.root", hist_data_path);
    TFile *fUNIF = TFile::Open(name, "READ");

    sprintf(name, "%sUFC_SB.root", hist_data_path);
    TFile *fUSB = TFile::Open(name, "READ");

    TFile *fAna = TFile::Open(results_file, "UPDATE");

    Int_t xmin[] = {150, 350, 150, 150, 100, 100, 100, 50}; // rough minimum ranges
    Int_t xmax[] = {400, 900, 400, 400, 400, 350, 350, 300};
    Double_t mFG, mBG, mAv, mSum;
    TGraph *gFG = new TGraph(8);
    TGraph *gBG = new TGraph(8);
    TGraph *gAv = new TGraph(8);
    TGraph *gSum = new TGraph(8);
    gFG->SetName("gFG");
    gBG->SetName("gBG");
    gAv->SetName("gAv");
    gSum->SetName("gSum");

    cout << endl << "UFC QDC min" << endl
         << "Ch   FG fit   BG fit   average   sum fit" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Open histograms and fit minima
        sprintf(name, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i+1);
        TH1I *hFG = (TH1I*) fUNIF->Get(name);
        TF1 *fFG = FindMin(hFG, xmin[i], xmax[i]);
        mFG = fFG->GetMinimumX(xmin[i], xmax[i]);
        TH1I *hBG = (TH1I*) fUSB->Get(name);
        TF1 *fBG = FindMin(hBG, xmin[i], xmax[i]);
        mBG = fBG->GetMinimumX(xmin[i], xmax[i]);

        // Save fit functions
        Save(fAna, "UFC/QDC/Fit/FG", fFG, "pol4_FG_"+to_string(i+1));
        Save(fAna, "UFC/QDC/Fit/BG", fBG, "pol4_BG_"+to_string(i+1));

        // Calculate average
        Int_t intFG = hFG->Integral(hFG->FindBin(mFG), 4096);
        Int_t intBG = hBG->Integral(hBG->FindBin(mBG), 4096);
        mAv = (intFG * mFG + intBG * mBG) / (intFG + intBG);

        // Fit sum of all runs
        hFG->Add(hBG, 1);
        TF1 *fSum = FindMin(hFG, xmin[i], xmax[i]);
        mSum = fSum->GetMinimumX(xmin[i], xmax[i]);

        // Save as graphs
        gFG->SetPoint(i, i+1, mFG);
        gBG->SetPoint(i, i+1, mBG);
        gAv->SetPoint(i, i+1, mAv);
        gSum->SetPoint(i, i+1, mSum);

        cout << " " << i+1 << "   " << mFG << "   " << mBG << "   "
             << mAv << "   " << mSum << endl;
    }
    Save(fAna, "UFC/QDC", gFG, "Minima_FG");
    Save(fAna, "UFC/QDC", gBG, "Minima_BG");
    Save(fAna, "UFC/QDC", gAv, "Average");
    Save(fAna, "UFC/QDC", gSum, "Sum");

    fAna->Save();
    fAna->Close();
    fUNIF->Close();
    fUSB->Close();
}

void CheckIntegral(TFile *f, string FC, Int_t ch)
{
    char name[64] = "";
//    sprintf(name, "Histograms/Raw/QDC/low/H1RawQDCl_%i", ch+1);
    sprintf(name, "Histograms/Analysis/FC/QDC/low/trig/H1AnaQDCl_trig_%i", ch+1);
    TH1I *hQDC = (TH1I*)f->Get(name); if (!hQDC) cerr << "Could not get " << name << endl;
    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1I *hTDC = (TH1I*)f->Get(name); if (!hTDC) cerr << "Could not get " << name << endl;
    sprintf(name, "Histograms/Analysis/FC/QDC/low/H1AnaQDCl_%i", ch+1);
    TH1I *hQDCg = (TH1I*)f->Get(name); if (!hQDCg) cerr << "Could not get " << name << endl;

    Double_t qdc_min_Pu[] = {899.24,  853.668, 895.652, 849.393, 1046.41, 891.396, 906.123, 837.486}; // PuFC
    Double_t qdc_min_U[]  = {220.396, 500.964, 219.016, 214.101, 199.19,  171.929, 169.046, 121.56}; // UFC
    Double_t qdc_min = 0;
    if (FC[0] == 'U')
        qdc_min = qdc_min_U[ch];
    else
        qdc_min = qdc_min_Pu[ch];
    Int_t bin_min = hQDC->GetXaxis()->FindBin(qdc_min);

    cout << f->GetName() << ", " << FC << ", " << ch+1 << endl;
    cout << "\tQDC\t" << hQDC->Integral(bin_min+1, hQDC->GetNbinsX()) << endl;
    cout << "\tQDCg\t" << hQDCg->Integral(bin_min+1, hQDCg->GetNbinsX()) << endl;
    cout << "\tTDC\t" << hTDC->Integral(0, hTDC->GetNbinsX()+1) << endl;

//    sprintf(name, "Diff_%s_%i", FC.c_str(), ch+1);
//    TH1I *hDiff = (TH1I*) hQDC->Clone(name);
//    hDiff->Add(hQDCg, -1.0);
    cout << "\tDiff\t" << hQDC->GetBinContent(hQDC->GetNbinsX()+1) << " - " << hQDCg->GetBinContent(hQDCg->GetNbinsX()+1) << endl;
//    new TCanvas();
//    hDiff->Draw();
//    hQDC->Draw();
//    hQDCg->Draw("same");
}

void CheckIntegrals(string file)
{
    TFile *f = TFile::Open(file.c_str());
    for (Int_t i = 0; i < 8; i++)
        CheckIntegral(f, "PuFC", i);
}

void QDCmin()
{
    FindPuMinima();
    FindUMinima();

//    CheckIntegrals("~/Programme/Go4nfis/offline/results/UFC_NIF.root");
}
