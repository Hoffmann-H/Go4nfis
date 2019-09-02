#include "PuFC.h"
using namespace std;

PuFC::PuFC(Bool_t two_runs, Bool_t draw)
{
    Name = "PuFC";
    cout << endl << "Creating fission chamber " << Name << endl;
    CommentFlag = kTRUE;
    InitVar(draw);
    InitPuVar();

    if (two_runs)
    {
        ///////////////////////////////////////////////////////////////////////////////////////
        Run("Open", "NIF", 27492078.81+33478369.95+30916955.27+54792678.94, 0.0014, 16794.37+20286.24+18736.55+34490.15);
        Run("SB", "SB", 3623068.621401+28646613.7419+4757385.125897, 0.0015, 2287.18+18331.64+3098.89);
        ///////////////////////////////////////////////////////////////////////////////////////
    } else {
        ///////////////////////////////////////////////////////////////////////////////////////
        Run("Open", "PuFC_FG_MS4", 27492078.81, 0.0014, 16794.37);
        Run("Open", "PuFC_FG_MS5", 33478369.95, 0.0014, 20286.24);
        Run("Open", "PuFC_FG_MS6", 30916955.27, 0.0014, 18736.55);
        Run("Open", "PuFC_FG_MS7", 54792678.94, 0.0014, 34490.15);
        Run("SB", "PuFC_BG_MS9", 3623068.621401, 0.0015, 2287.18);
        Run("SB", "PuFC_BG_MS10", 28646613.7419, 0.0014, 18331.64);
        Run("SB", "PuFC_BG_MS11", 4757385.125897, 0.0014, 3098.89);
        ///////////////////////////////////////////////////////////////////////////////////////
    }
    cout << endl << "Created " << Name << endl;
}


//void PuFC::plt()
//{// So funktioniert die Übergabe an Plot!
//    Double_t x1[NumCh];
//    Double_t y1[NumCh];
//    Double_t x2[NumCh];
//    Double_t y2[NumCh];
//    for (int i = 0; i < NumCh; i++)
//    {
//        x1[i] = 4*i;
//        y1[i] = i;
//        x2[i] = 4-i;
//        y2[i] = 0.5;
//    }
//    plot->ExpT(x1, y1, x2, y2);
//}


PuFC::~PuFC()
{

}


void PuFC::InitPuVar()
{
    cout << endl << "Initializing Pu variables.." << endl;
    SetLimits();
    PuSFT2 = 6.76E10 * 365.24*24*60*60; // Pu-242 spontaneaus fission half-life period in s
    DPuSFT2 = 7E8 * 365.24*24*60*60;
//    MonitorFG = 146680082.972163; // MS#4 - MS#7
//    DMonitorFG = 0.0014;
//    MonitorBG = 3623069+28646614+4757385;
//    DMonitorBG = 0.0015;

    /// Distances..
    /// Doppelspaltkammer_Sobiella.pdf
    ///
    /// Hinterstes Deposit richtig bemasst (bis Mitte Si):
    /// Quelle - Deposit =
    ///   (1500+-1) mm Quelle - Cover
    /// + 2 mm Cover
    /// + 207.400 mm front Flansch - center of last Si wafer
    /// - 0.2 mm Front surface of last Si wafer
    /// -> Last Deposit (nr 1, index 0): 209.2 mm after cover
    ///
    /// Falsche Abstaende von 10.4 mm. Richtige:
    /// 100 um Tantal
    /// 10 mm spacer
    /// 400 um Si
    /// 10 mm spacer
    /// -> Distance btw Deposits: 20.5 mm

    cout << endl << "Distances... " << endl << " Deposit   S-D   solid angle" << endl;
    Double_t L = 1500, DL = 1;
    for (int i = 0; i < NumHist; i++)
    {// Distances source-deposit
        sd[i] = L + 209.2 - 20.5 * i;
        Dsd[i] = DL;
        cout << " " << i+1 << "   " << sd[i] << " mm   " << Area / pow(sd[i], 2) << endl;
    }
    cout << "Done: Pu variables" << endl;
}


void PuFC::HardCodedThresholds()
{
    Double_t avCut[] = {907, 853, 896, 849, 979, 892, 904, 837};

    cout << endl << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " <<
                       avCut[3] << " " << avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
}


//void PuFC::AnalyzeQDC()
//{// Calculate the average QDC minimum position
//    if (CommentFlag)
//        cout << endl << "Calculate QDC thesholds" << endl;
//    Double_t avCut[NumHist];
//    Double_t DavCut[NumHist];

//    for (Int_t j = 0; j < FgRuns; j++)
//        pHBG[j]->DoAnalyzeQDC();
//    for (Int_t j = 0; j < BgRuns; j++)
//        pHFG[j]->DoAnalyzeQDC();
//    pHSF->DoAnalyzeQDC();

//    for (int i = 0; i < NumHist; i++)
//    {
//        Double_t w[FgRuns + BgRuns + 1]; // weights
//        Double_t wSum = 0;
//        avCut[i] = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//        {
//            w[j] = pHFG[j]->GetNevents(i);
//            avCut[i] += w[j] * pHFG[j]->CutQDC[i];
//            wSum += w[j];
//        }
//        for (Int_t j = 0; j < BgRuns; j++)
//        {
//            w[FgRuns + j] = pHBG[j]->GetNevents(i);
//            avCut[i] += w[FgRuns + j] * pHBG[j]->CutQDC[i];
//            wSum += w[FgRuns + j];
//        }
//        w[FgRuns + BgRuns] = pHSF->GetNevents(i);
//        avCut[i] += w[FgRuns + BgRuns] * pHSF->CutQDC[i];
//        wSum += w[FgRuns + BgRuns];
//        avCut[i] /= wSum;
//        DavCut[i] = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//            DavCut[i] += w[j] * pow(pHFG[j]->CutQDC[i] - avCut[i], 2);
//        for (Int_t j = 0; j < BgRuns; j++)
//            DavCut[i] += w[FgRuns + j] * pow(pHBG[j]->CutQDC[i] - avCut[i], 2);
//        DavCut[i] += w[FgRuns + BgRuns] * pow(pHSF->CutQDC[i] - avCut[i], 2);
//        DavCut[i] = sqrt(DavCut[i] / wSum);
//        if (CommentFlag)
//            cout << " ch " << i+1 << ": " << avCut[i] << "+-" << DavCut[i] << endl;
//    }
//    cout << "QDC thresholds: " << avCut[0] << " " << avCut[1] << " " << avCut[2] << " " << avCut[3] << " " <<
//                                  avCut[4] << " " << avCut[5] << " " << avCut[6] << " " << avCut[7] << endl;
//    cout << "Done: QDC thresholds" << endl;
//    DoneQDC = kTRUE;
//}


//void PuFC::AnalyzeDtBG()
//{   // Calculate average fission background per livetime and Dt bin
//    if (!DoneLimits)
//        GetLimits();
//    cout << endl << "Analyzing ToF background..." << endl;
//    if (CommentFlag)
//        cout << " Ch   FG   BG   SF   avBg" << endl;
//    for (int i = 0; i < NumHist; i++)
//    {
//        // Integrate background
//        Double_t counts[4];
//        counts[0] = 0; counts[1] = 0; counts[2] = 0;
//        Double_t tBin[4];
//        tBin[0] = 0; tBin[1] = 0; tBin[2] = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//        {
//            counts[0] += pHFG[j]->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHFG[j]->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
//            tBin[0] += pHFG[j]->t_live * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]);
//        }
//        for (Int_t j = 0; j < BgRuns; j++)
//        {
//            counts[1] += pHBG[j]->pHDtG[i]->Integral(lim[0][i], lim[3][i]) - pHBG[j]->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
//            tBin[1] += pHBG[j]->t_live * (lim[1][i] - lim[0][i] + lim[3][i] - lim[2][i]);
//        }
//        counts[2] += pHSF->pHDtG[i]->Integral(lim[0][i], lim[3][i]);
//        tBin[2] += pHSF->t_live * (lim[3][i] - lim[0][i] + 1); // Counted bins * live time
//        counts[3] = counts[0] + counts[1] + counts[2];
//        tBin[3] = tBin[0] + tBin[1] + tBin[2];
//        avBg[i] = counts[3] / tBin[3];
//        DavBg[i] = sqrt(counts[3]) / tBin[3];
//        if (CommentFlag)
//            cout << " " << i+1 << "  " << counts[0]/tBin[0] << "+-" << sqrt(counts[0])/tBin[0] << "  " << counts[1]/tBin[1] << "+-" << sqrt(counts[1])/tBin[1] << "  " << counts[2]/tBin[2] << "  " << avBg[i] << "+-" << DavBg[i] << endl;
//    }
//    DoneDtBG = kTRUE;
//    cout << "Done: ToF background" << endl;
//}


void PuFC::GetNatoms()
{
//    if (!DoneDt)
//        AnalyzeDt();
    cout << endl << "Calculating effective number of Pu atoms..." << endl;
    cout << " T2SF(Pu-242) = " << PuSFT2 << "+-" << DPuSFT2 << endl;
    char name[64] = "";
    TFile *f;

    // Intergation limits in bin numbers (1ns binning)
    Int_t left = 10; Int_t right = 50;
    Int_t l0 = 42;
    Int_t l1[NumCh];
    Int_t m[] = {140, 128, 126, 129, 129, 128, 127, 127};
    Int_t l2[NumCh];
    Int_t l3 = 402;

    Double_t nSF[NumCh];
    Double_t D2nSF[NumCh];
    Double_t t_live = 0;

    cout << " Run   ch   t_live   constUg   factor   sf" << endl;
    f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root");
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    TH1I* pT = (TH1I*) f->Get(name);
    t_live += pT->Integral();
    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*) f->Get(name);
        Double_t sf = pH->Integral();
        Double_t D2sf = pH->Integral();
        cout << " SF   " << i+1 << "   " << pT->Integral() << "   -   -   " << sf << " +- " << sqrt(D2sf) << endl;
        nSF[i] = sf;
        D2nSF[i] = D2sf;
    }
    f->Close();

    f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pT = (TH1I*) f->Get(name);
    t_live += pT->Integral();
    for (Int_t i = 0; i < NumCh; i++)
    {
        l1[i] = m[i] - left;
        l2[i] = m[i] + right;

        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*) f->Get(name);
        Double_t constUg = (Double_t)pH->Integral(l0, l1[i]-1) + (Double_t)pH->Integral(l2[i], l3-1);
        Double_t factor = (Double_t)(l3 - l0) / (Double_t)(l1[i] - l0 + l3 - l2[i]);
        Double_t sf = pH->Integral()
                    - pH->Integral(l0, l3-1)
                    + factor * constUg;
        Double_t D2sf = pH->Integral()
                      - pH->Integral(l0, l3-1)
                      + factor*factor * constUg;
        cout << " FG   " << i+1 << "   " << pT->Integral() << "   " << constUg << "   " << factor << "   " << sf << " +- " << sqrt(D2sf) << endl;
        nSF[i] += sf;
        D2nSF[i] += D2sf;
    }
    f->Close();
    f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pT = (TH1I*) f->Get(name);
    t_live += pT->Integral();
    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I* pH = (TH1I*) f->Get(name);
        Double_t constUg = (Double_t)pH->Integral(l0, l1[i]-1) + (Double_t)pH->Integral(l2[i], l3-1);
        Double_t factor = (Double_t)(l3 - l0) / (Double_t)(l1[i] - l0 + l3 - l2[i]);
        Double_t sf = pH->Integral()
                    - pH->Integral(l0, l3-1)
                    + factor * constUg;
        Double_t D2sf = pH->Integral()
                      - pH->Integral(l0, l3-1)
                      + factor*factor * constUg;
        cout << " SB   " << i+1 << "   " << pT->Integral() << "   " << constUg << "   " << factor << "   " << sf << " +- " << sqrt(D2sf) << endl;
        nSF[i] += sf;
        D2nSF[i] += D2sf;
    }
    f->Close();
    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
    Double_t DrSF_Toni[] = {0.0013, 0.0012, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0013};
    cout << " ch   sf-rate   vgl.Toni   nAtoms" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
//        nAtoms[i] = rSF_Toni[i] * PuSFT2 / log(2.0);
//        DnAtoms[i] = DrSF_Toni[i] * PuSFT2 / log(2.0);
        nAtoms[i] = nSF[i] / t_live * PuSFT2 / log(2.0);
        DnAtoms[i] = sqrt(D2nSF[i]) / t_live * PuSFT2 / log(2.0);
        cout << " " << i+1 << "   " << nSF[i] / t_live << " +- " << sqrt(D2nSF[i]) / t_live << "   " << rSF_Toni[i] << " +- " << DrSF_Toni[i] << "   " << nAtoms[i] / 0.986 << " -+ " << DnAtoms[i] / 0.986 << endl;
    }

    cout << "Done: effective number of Pu atoms" << endl;
    DoneNatoms = kTRUE;
}


void PuFC::DrawStability()
{
    cout << endl << "Drawing PuFC sf stability..." << endl;

    char name[64] = "";
    char title[128] = "";
    TFile *f;

    Int_t tlive_width;
    Int_t tlive_nbins;
    Int_t stab_width = 60 * 30;
    Int_t stab_nbins;

    TH1F* pH1SF[NumCh];

//    cout << " " << endl;
    f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root");
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    TH1D* pT = (TH1D*) f->Get(name);
    Double_t tmin = pT->GetXaxis()->GetXmin();
    Double_t tmax = pT->GetXaxis()->GetXmax();
    tlive_width = pT->GetBinWidth(1);
    tlive_nbins = pT->GetNbinsX();
    stab_nbins = tlive_nbins * tlive_width / stab_width;
    cout << tlive_width << " " << tlive_nbins << " " << stab_width << " " << stab_nbins << endl;

    pHtLive = new TH1D("PuFC_tLive", "PuFC Live time; #font[12]{t}; t_live", tlive_nbins, tmin, tmax);
    pHtLive->Add(pT, 1);

    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "Stab_%i", i+1);
        sprintf(title, "Effektive Spontanspaltrate Deposit %i; #font[12]{t}; Eff. rate [1/s]", i+1);
        pH1SF[i] = new TH1F(name, title, stab_nbins, tmin, tmax);
//        pH1SF[i]->Fill(tmin + 1);
//        pH1SF[i]->Fill(tmax - 1);
    }

    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", i+1);
        TH2I* pH = (TH2I*) f->Get(name);
        Int_t h2y_nbins = pH->GetNbinsY();
        for (Int_t sbin = 1; sbin < stab_nbins+1; sbin++)
        {
            Int_t tbin0 = (sbin - 1) * tlive_nbins / stab_nbins + 1;
            Int_t tbin1 = sbin * tlive_nbins / stab_nbins;
            Double_t t_live = pT->Integral(tbin0, tbin1);
            if (t_live > 0.5 * stab_width)
            {
                tbin0 = (sbin - 1) * h2y_nbins / stab_nbins + 1;
                tbin1 = sbin * h2y_nbins / stab_nbins;
                Double_t sf = pH->Integral(1, 441, tbin0, tbin1);
                pH1SF[i]->SetBinContent(sbin, sf / t_live);
                pH1SF[i]->SetBinError(sbin, sqrt(sf) / t_live);
//                cout << " " << sbin << "   " << tbin0 << "-" << tbin1 << "   " << sf << endl;
            }
        }
    }

    f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pT = (TH1D*) f->Get(name);
    pHtLive->Add(pT, 1);

    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", i+1);
        TH2I* pH = (TH2I*) f->Get(name);
        Int_t h2y_nbins = pH->GetNbinsY();
        for (Int_t sbin = 1; sbin < stab_nbins+1; sbin++)
        {
            Int_t tbin0 = (sbin - 1) * tlive_nbins / stab_nbins + 1;
            Int_t tbin1 = sbin * tlive_nbins / stab_nbins;
            Double_t t_live = pT->Integral(tbin0, tbin1);
            if (t_live > 0.5 * stab_width)
            {
                tbin0 = (sbin - 1) * h2y_nbins / stab_nbins + 1;
                tbin1 = sbin * h2y_nbins / stab_nbins;
                Double_t constUg = (Double_t)pH->Integral(l0, l1[i]-1, tbin0, tbin1) + (Double_t)pH->Integral(l2[i], l3-1, tbin0, tbin1);
                Double_t factor = (Double_t)(l3 - l0) / (Double_t)(l1[i] - l0 + l3 - l2[i]);
                Double_t sf = pH->Integral(0, l0-1, tbin0, tbin1)
                        + pH->Integral(l3, 441, tbin0, tbin1)
                        + factor * constUg;
                Double_t D2sf = pH->Integral(0, l0-1, tbin0, tbin1)
                        + pH->Integral(l3, 441, tbin0, tbin1)
                        + factor*factor * constUg;
                pH1SF[i]->SetBinContent(sbin, sf / t_live);
                pH1SF[i]->SetBinError(sbin, sqrt(D2sf) / t_live);
                cout << " " << sbin << "   " << tbin0 << "-" << tbin1 << "   " << sf << endl;
            }
        }
    }

    f = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
    sprintf(name, "/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pT = (TH1D*) f->Get(name);
    pHtLive->Add(pT, 1);

    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H2DtGvsTime_%i", i+1);
        TH2I* pH = (TH2I*) f->Get(name);
        Int_t h2y_nbins = pH->GetNbinsY();
        for (Int_t sbin = 1; sbin < stab_nbins+1; sbin++)
        {
            Int_t tbin0 = (sbin - 1) * tlive_nbins / stab_nbins + 1;
            Int_t tbin1 = sbin * tlive_nbins / stab_nbins;
            Double_t t_live = pT->Integral(tbin0, tbin1);
            if (t_live > 0.5 * stab_width)
            {
                tbin0 = (sbin - 1) * h2y_nbins / stab_nbins + 1;
                tbin1 = sbin * h2y_nbins / stab_nbins;
                Double_t constUg = (Double_t)pH->Integral(l0, l1[i]-1, tbin0, tbin1) + (Double_t)pH->Integral(l2[i], l3-1, tbin0, tbin1);
                Double_t factor = (Double_t)(l3 - l0) / (Double_t)(l1[i] - l0 + l3 - l2[i]);
                Double_t sf = pH->Integral(0, l0-1, tbin0, tbin1)
                        + pH->Integral(l3, 441, tbin0, tbin1)
                        + factor * constUg;
                Double_t D2sf = pH->Integral(0, l0-1, tbin0, tbin1)
                        + pH->Integral(l3, 441, tbin0, tbin1)
                        + factor*factor * constUg;
                pH1SF[i]->SetBinContent(sbin, sf / t_live);
                pH1SF[i]->SetBinError(sbin, sqrt(D2sf) / t_live);
                cout << " " << sbin << "   " << tbin0 << "-" << tbin1 << "   " << sf << endl;
            }
        }
    }

    TCanvas* c1 = new TCanvas("cStab", "Stability", 200, 10, 700, 500);
    pH1SF[0]->Draw();
    c1->Modified();
    c1->Update();
}

//void PuFC::GetNatoms()
//{
////    if (!DoneDt)
////        AnalyzeDt();
//    cout << endl << "Calculating effective number of Pu atoms..." << endl;
//    cout << " T2SF(Pu-242) = " << PuSFT2 << "+-" << DPuSFT2 << endl;
//    Double_t eMinimum = 0.986;
//    Double_t DeMinimum = 0.010;
//    Double_t cAtoms[NumCh];
//    Double_t DcAtoms[NumCh];
//    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
//    cout << "Ch   SF-rate   eff.Atoms   Atoms" << endl;
//    for (int i = 0; i < NumCh; i++)
//    {
//        Double_t nSF = 0;
//        Double_t t_live = 0;
//        for (Int_t j = 0; j < FgRuns + BgRuns; j++)
//        {
//            nSF += pH[j]->nSF[i];
//            t_live += pH[j]->t_live;
//        }
//        nSF += pHSF->nSF[i];
//        t_live += pHSF->t_live;
//        nAtoms[i] = rSF_Toni[i] * PuSFT2 / log(2.0);
////        nAtoms[i] = nSF / t_live * PuSFT2 / log(2.0);
//        DnAtoms[i] = sqrt(nSF) / t_live * PuSFT2 / log(2.0);
//        cAtoms[i] = nAtoms[i] / eMinimum;
//        DcAtoms[i] = cAtoms[i] * sqrt( pow(DnAtoms[i] / nAtoms[i], 2) + pow(DeMinimum / eMinimum, 2) );
//        cout << " " << i+1 << "   " << nSF / t_live << "+-" << sqrt(nSF) / t_live << "   " << nAtoms[i] << "+-" << DnAtoms[i] << "   " << cAtoms[i] << "+-" << DcAtoms[i] << endl;
//    }
//    for (Int_t k = 0; k < FgRuns + BgRuns; k++)
//        pR[k]->SetNatoms(nAtoms, DnAtoms);
//    cout << "Done: effective number of Pu atoms" << endl;
//    DoneNatoms = kTRUE;
//}


void PuFC::IsoVec()
{
    cout << endl << "PuFC::IsoVec() defining isotope vector correction factor..." << endl;
    // PuFC: Do nothing.
    fIsoVec = 1;
    DfIsoVec = 0;
    sIsoVec = 0;
    DsIsoVec = 0;
    cout << "Done: IsoVec" << endl;
    DoneIso = kTRUE;
}


//void PuFC::GetExpT()
//{
//    cout << endl << "PuFC::GetExpT() experimental transmission" << endl;
//    for (int i = 0; i < NumCh; i++)
//    {
//        Double_t nSF = 0;
//        for (Int_t j = 0; j < FgRuns; j++)
//            nSF += pHFG[j]->nSF[i];
//        uT[i] = nFG[i] / nSF; // ExpT := N(n,f) / N(sf)
//        DuT[i] = sqrt( pow(DnFG[i] / nFG[i], 2) +
//                   1.0 / nSF ) * uT[i];
//    }
//}
