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
    TFile *fPuNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root", "UPDATE");
    TFile *fPuSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root", "UPDATE");
    TFile *fPuSF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root", "UPDATE");
    TDirectory *pDir;
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
    char name[64] = "";
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
        pDir = Prepare(fPuNIF, "Analysis/QDC/Fit");
        Save(pDir, fFG, "pol4_FG_" + to_string(i+1));
        pDir = Prepare(fPuSB, "Analysis/QDC/Fit");
        Save(pDir, fBG, "pol4_BG_" + to_string(i+1));
        pDir = Prepare(fPuSF, "Analysis/QDC/Fit");
        Save(pDir, fSF, "pol4_SF_" + to_string(i+1));

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
    pDir = Prepare(fPuNIF, "Analysis/QDC");
    Save(pDir, gFG, "Minima_FG");
    Save(pDir, gAv, "Average");
    Save(pDir, gSum, "Sum");
    pDir = Prepare(fPuSB, "Analysis/QDC");
    Save(pDir, gBG, "Minima_BG");
//    Save(pDir, gAv, "Average");
//    Save(pDir, gSum, "Sum");
    pDir = Prepare(fPuSF, "Analysis/QDC");
    Save(pDir, gSF, "Minima_SF");
//    Save(pDir, gAv, "Average");
//    Save(pDir, gSum, "Sum");

    fPuNIF->Save();
    fPuSB->Save();
    fPuSF->Save();
    fPuNIF->Close();
    fPuSB->Close();
    fPuSF->Close();
}

void FindUMinima()
{
    TFile *fUNIF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root", "UPDATE");
    TFile *fUSB = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root", "UPDATE");
    TDirectory *pDir;
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
    char name[64] = "";
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
        pDir = Prepare(fUNIF, "Analysis/QDC/Fit");
        Save(pDir, fFG, "pol4_FG_" + to_string(i+1));
        pDir = Prepare(fUSB, "Analysis/QDC/Fit");
        Save(pDir, fBG, "pol4_BG_" + to_string(i+1));

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
    pDir = Prepare(fUNIF, "Analysis/QDC");
    Save(pDir, gFG, "Minima_FG");
    Save(pDir, gAv, "Average");
    Save(pDir, gSum, "Sum");
    pDir = Prepare(fUSB, "Analysis/QDC");
    Save(pDir, gBG, "Minima_BG");
//    Save(pDir, gAv, "Average");
//    Save(pDir, gSum, "Sum");

    fUNIF->Save();
    fUSB->Save();
    fUNIF->Close();
    fUSB->Close();
}

void QDCmin()
{
//    FindPuMinima();
    FindUMinima();
}
