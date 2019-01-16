// Class for loading and accessing histograms from a root file. 
// This is meant to be created once for each exp. setup: NIF, SB. 
// Open histograms: Hist *pHist = new Hist("path/to/file.root")

#include "Hist.h"

#define NumHist 8

using namespace std;

Hist::Hist(string file_path)   // Constructor
{
    CommentFlag = kFALSE;

    //define surrounding peak limits
    ToF_low = 66500;
    ToF_up  = 68500;
    //define time window for underground integration
    Dt_min = 63300;
    Dt_max = 78500;

    FilePath = file_path.c_str();
    file = TFile::Open(FilePath);//, "UPDATE");
    if (file == 0) //check if root file exists
    {   cerr << "File not found!" << endl;
        return;
    }
    pHtLive = GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pHtReal = GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_48");
    char hist_name[64] = "";
    for(int i_ch = 0; i_ch < 8; i_ch++) {
        sprintf(hist_name, "/Histograms/Analysis/PuFC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i_ch + 1);
        pHDtG[i_ch] = GetTH1I(hist_name);
        sprintf(hist_name, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i_ch + 1);
        pHRawQDCl[i_ch] = GetTH1I(hist_name);
        sprintf(hist_name, "/Histograms/Analysis/PuFC/QDC/low/trig/H1AnaQDCl_trig_%i", i_ch + 1);
        pHAnaQDCl[i_ch] = GetTH1I(hist_name);
    }

    if(CommentFlag)
        cout << "Created instance Hist using " << file_path << endl;

//    DoAnalyzeDt(kTRUE);
//    file->Close();
}

Hist::~Hist()    // Destructor
{
    if(CommentFlag)
        cout << "Destroy" << endl;
}

TH1I* Hist::GetTH1I(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TH1I* histo = (TH1I*) file->Get(hname);
    if(CommentFlag)
        cout << "Opening histogram " << hname << endl;
    if (histo==0)
    {   cerr << "Histogram " << hname << " not found!" << endl;
        return 0;
    }
    else
        return histo;
}

TH1D* Hist::GetTH1D(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TH1D* histo = (TH1D*) file->Get(hname);
    if(CommentFlag)
        cout << "Opening histogram " << hname << endl;
    if (histo==0)
    {   cerr << "Histogram " << hname << " not found!" << endl;
        return 0;
    }
    else
        return histo;
}

void Hist::SaveToFile(string fname, string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    TFile *f = new TFile(fname.c_str(), "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
        //check if folder "Analysis" already exists, otherwise create it
        if (f->Get(path.c_str())!=0)
            EvalDir = (TDirectory*) f->Get(path.c_str());
        else
        {   f->mkdir(path.c_str(), "Folder containing offline Analysis objects");
            EvalDir = f->GetDirectory(path.c_str());
        }
        EvalDir->cd();
        if (EvalDir->Get(GraphName.c_str())!=0)
            EvalDir->Delete((GraphName+";*").c_str());
        pGraph->Clone()->Write();
    f->Save(); f->Close();
}

void Hist::SaveTo(string path, string name)
{   // Prepare saving a root object.
    // Check if folder exists. If not, create.
    TDirectory *SaveDir;
    if (file->Get(path.c_str())!=0)
        SaveDir = (TDirectory*) file->Get(path.c_str());
    else
    {
        if(CommentFlag)
            cout << "Create dir " << path << endl;
        SaveDir = file->mkdir(path.c_str(), "subdir");
    }
    string dir_name = "/"+path;
    SaveDir->cd(dir_name.c_str());
    if(CommentFlag)
        SaveDir->pwd();
    // Delete old data
    if (SaveDir->Get(name.c_str())!=0) {
        string del_name = name+";*";
        SaveDir->Delete(del_name.c_str());
        if (CommentFlag)
            cout << "Delete old " << name << endl;
    }
    else if (CommentFlag)
        cout << "Object " << name << " not existing yet." << endl;
    return;
    // Write object now.
}

void Hist::SaveTH1I(TH1I* pHtoSave, string hpath, string hname)
{
    // prepare
    SaveTo(hpath, hname);
    // save
    pHtoSave->Write(hname.c_str());
    file->Save();
    return;
}

void Hist::SaveTF1(TF1* pFtoSave, string fpath, string fname)
{
    SaveTo(fpath, fname);
    pFtoSave->Write(fname.c_str());
    file->Save();
    return;
}

void Hist::SavePad(TPad* pPtoSave, string ppath, string pname)
{
    SaveTo(ppath, pname);
    pPtoSave->Write(pname.c_str());
    file->Save();
    return;
}

void Hist::SaveCanvas(TCanvas* pCtoSave, string cpath, string cname)
{
    SaveTo(cpath, cname);
    pCtoSave->Write(cname.c_str());
    file->Save();
    return;
}

void Hist::DoAnalyzeDt(Bool_t peak)
{   // Analyze TimeDiff spectra for all channels: Determine underground and, if peak is set to kTRUE, peak content. Set peak = kFALSE for measurements without beam.
    char hname[64] = "";
    for (int i_ch = 0; i_ch < NumHist; i_ch++)
    {
        sprintf(hname, "H1AnaDtFit_%i", i_ch + 1);
        TPad* pC = new TPad();//(hname, "Pulse-height gated Time difference fit", 600, 400);
        if(peak)
            AnalyzeDtPeak(pHDtG[i_ch], (TCanvas*)pC, &nNIF_raw[i_ch], &DnNIF_raw[i_ch], &nSF_raw[i_ch], &DnSF_raw[i_ch]);
        else
            AnalyzeDtUnderground(pHDtG[i_ch], (TCanvas*)pC, &nSF_raw[i_ch], &DnSF_raw[i_ch]);
//        SaveTH1I((TH1I*)pHDtG[i_ch], "Analysis/TimeDiff", hname);
//        SavePad(pC, "Analysis/TimeDiff", hname);
    }

    DeadTimeCorrection(peak);

//    if (CommentFlag)
        cout << "Times: " << t_real << ", " << t_live << endl;
    for (int i = 0; i < NumHist; i++)
    {
        SFRate[i] = nSF[i] / t_real;
        DSFRate[i] = DnSF[i] / t_real;
    }

    TCanvas* pCSF = new TCanvas("SFRate", "Title", 600, 400);
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    TGraphErrors* g1 = new TGraphErrors(NumHist, x, SFRate, xerr, DSFRate);
    g1->SetMarkerStyle(20);
    g1->SetTitle("SF detection rate vs channel");
    g1->Draw("ap");
//    SaveCanvas(pCSF, "Analysis/TimeDiff", "SFRate");

    if(peak)
    {   // If measured with beam, draw NIF rates
        for (int i = 0; i < NumHist; i++)
        {
            NIFRate[i] = nNIF[i] / t_real;
            DNIFRate[i] = DnNIF[i] / t_real;
        }
        TCanvas* pCNIF = new TCanvas("NIFRate", "Neutron-induced fission rate without in-scattering correction", 600, 400);
        TGraphErrors* g2 = new TGraphErrors(NumHist, x, NIFRate, xerr, DNIFRate);
        g2->Draw("ap");
//        SaveCanvas(pCNIF, "Analysis/TimeDiff", "NIFRate");
    }
    return;
}

void Hist::AnalyzeDtPeak(TH1I *pH, TCanvas* pCresult, Double_t *pNIF, Double_t *pDNIF, Double_t *pSF, Double_t *pDSF)
{
    Double_t ChPerBin = pH->GetBinWidth(0);
    Double_t BinOffset = pH->GetBinCenter(0);
    // Calculate integration limit bins
    Double_t lim_0 = (Dt_min - BinOffset) / ChPerBin;
    Double_t lim_1 = (ToF_low - BinOffset) / ChPerBin;
    Double_t lim_2 = (ToF_up - BinOffset) / ChPerBin;
    Double_t lim_3 = (Dt_max - BinOffset) / ChPerBin;
    if(CommentFlag)
        cout << "Analyzing peak. Ch per bin: " << ChPerBin << ", Offset: " << BinOffset << endl <<
                "Integration limits: Bin " << lim_0 << " " << lim_1 << " " << lim_2 << " " << lim_3 << endl;
    // Integrations
    Int_t PeakCount = pH->Integral(lim_1, lim_2); // Sum within peak region
    Double_t UgCount = pH->Integral(lim_0, lim_3) - PeakCount; // Sum within underground region; follows Poisson-stat
    Double_t DUgCount = sqrt(UgCount); // Underground deviation
    Double_t UgPerBin = 1.0 * UgCount / (lim_1 - lim_0 + lim_3 - lim_2);
    Double_t DUgPerBin = 1.0 * DUgCount / (lim_1 - lim_0 + lim_3 - lim_2);
    *pNIF = PeakCount - (lim_2 - lim_1) * UgPerBin;
    *pDNIF = (lim_2 - lim_1) * DUgPerBin;
    *pSF = pH->Integral() - *pNIF; // Calculate total underground
    *pDSF = sqrt((*pDNIF)*(*pDNIF) + *pSF); // Error according to Gauss
    if(CommentFlag)
        cout << "NIF: " << *pNIF << " +- " << *pDNIF << endl <<
                "SF:  " << *pSF  << " +- " << *pDSF  << endl;
    // Draw
    pH->GetXaxis()->SetRangeUser(62000, 80000);
    pH->Draw();
    TLine *lH = new TLine(Dt_min, UgPerBin, Dt_max, UgPerBin);
    lH->SetLineColor(kRed);
    lH -> Draw("same");
    Double_t max = pH->GetMaximum();
    TLine *lV0 = new TLine(ToF_low, 0, ToF_low, max);
    lV0->Draw();
    TLine *lV1 = new TLine(ToF_up , 0, ToF_up , max);
    lV1->Draw();
    char message[32] = "";
    sprintf(message, "NIF counts: %i", (int)*pNIF);
    TText *t = new TText();
    t->SetNDC();
    t->DrawText(0.4, 0.8, message);
//    pH->GetListOfFunctions()->Add(lH);
    pH->GetListOfFunctions()->Add(lV0);
//    pH->GetListOfFunctions()->Add(lV1);
}

void Hist::AnalyzeDtUnderground(TH1I *pH, TCanvas* pCresult, Double_t *pSF, Double_t *pDSF)
{
    Double_t ChPerBin = pH->GetBinWidth(0);
    Double_t BinOffset = pH->GetBinCenter(0);
    // Calculate integration limit bins
    Int_t lim_0 = (Dt_min - BinOffset) / ChPerBin;
    Int_t lim_3 = (Dt_max - BinOffset) / ChPerBin;
    // Analyze underground analogue to AnalyzePeak(), but without peak area exclusion.
    Double_t SF = pH->Integral();
    Double_t DSF = sqrt(SF);
    *pSF = SF;
    *pDSF = DSF;
    return;
}

void Hist::DoAnalyzeQDC()
{
    for (int i_ch = 0; i_ch < NumHist; i_ch++)
    {
        AnalyzeQDC(pHRawQDCl[i_ch], i_ch, 0);
    }
    TCanvas *pC = new TCanvas("Name", "Title", 600, 400);
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
//    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    TGraph* g1 = new TGraph(NumHist, x, CutQDC);
    g1->SetMarkerStyle(20);
    g1->SetTitle("Relative minimum position vs channel");
    g1->Draw("ap");
//    SaveCanvas(pC, "Analysis/QDC", "relCut");
    return;
}

void Hist::AnalyzeQDC(TH1I *pH, Int_t channel, Double_t pedestal)
{   // Analyze single QDC spectrum. Use hard-coded fit positions acc. to channel. For self-triggered QDC spectra, use pedestal from before
    // rough peak positions   channel 0,    1,    2,    3,    4,    5,    6,    7
    Int_t RoughPeakPosition[3][8] = {{125,  75,   155,  90,   80,   80,   100,  60},
                                     {1000, 900,  950,  900,  1000, 950,  950,  900},
                                     {1600, 1500, 1550, 1500, 1600, 1550, 1550, 1600} };
    // 1. Fit pedestal
//    TF1 *Fit = new TF1("fit", func2, 0.0, 1.0, 3);
    Double_t Pedestal;
    Double_t DPedestal;
    char fname[64] = "";
//    char hname[64] = "";
    Int_t low, up;
    if (pedestal == 0)
    {
        Int_t PBin = pH->GetMaximumBin(); // bin number of maximum
        Int_t Max = pH->GetBinContent(PBin); // maximum value
        sprintf(fname, "fPed_%i", channel + 1);
        // Bounds: FWHM
        low = pH->GetBinCenter(PBin)-1;
        up = pH->GetBinCenter(PBin)+1;
        while(2 * pH->GetBinContent(low) > Max) {
            low--;
        }
        while(2 * pH->GetBinContent(up) > Max) {
            up++;
        }
        TF1* fP = new TF1(fname, func2, low, up+1, 3);//(TF1*)Fit->Clone();
        fP->SetRange(low, up+1);
        fP->SetParameters(-100000.0, (Double_t)pH->GetBinCenter(PBin), (Double_t)Max);
        pH->Fit(fname, "RB");
        Pedestal = (Double_t)fP->GetParameter(1);
        DPedestal = (Double_t)fP->GetParError(1);
//        SaveTF1(fP, "Histograms/Raw/QDC/low/fit", fname);
    }
    else
        Pedestal = pedestal;

    // 2. Find (FF,alpha)-cut and underground level
    low = RoughPeakPosition[1][channel] - 200;
    up = RoughPeakPosition[1][channel] + 300;
    sprintf(fname, "fCut_%i", channel + 1);
    TF1* fC = new TF1(fname, "pol4", low, up);
    fC->SetRange(low, up);
    fC->SetParameters(1.E3, -10.0, 0.01, -1.E-5, 1.E-8);
    pH->Fit(fname, "RB");
    Double_t Cut = fC->GetMinimumX();
//    SaveTF1(fC, "Histograms/Raw/QDC/low/fit", fname);
    cout << "**** QDC cut for channel " << channel+1 << ": " << Cut << endl;

/// Find underground level

    // 3. Find FF-maximum
    Int_t xMax = 0, yMax = 0;
    for(Int_t bin = 1300; bin < 2000; bin++)
    {   // search for bin with maximum content in range
        if(pH->GetBinContent(bin) > yMax)
            xMax = bin;
        yMax = pH->GetBinContent(xMax);
    }
    low = xMax - 200;
    up = xMax + 300;
    sprintf(fname, "fMax_%i", channel + 1);
    TF1* fM = new TF1(fname, "pol4", low, up);
    fM->SetRange(low, up);
    fM->SetParameters(-1.E4, 40.0, -0.01, 1.E-7, -1.E-10);
    pH->Fit(fname, "RB");
    Double_t Max = fM->GetMaximumX();
//    SaveTF1(fM, "Histograms/Raw/QDC/low/fit", fname);
    if(CommentFlag)
        cout << "Pedestal " << Pedestal << ", Cut " << Cut << ", Max " << Max << endl;

    // 4. Calculate relative cut position and efficiency
    Double_t relCutPos = (Cut - Pedestal) / (Max - Pedestal);
    PedQDC[channel] = Pedestal;
//    UgQDC[channel] = Ug;
//    DUgQDC[channel] = DUg;
    CutQDC[channel] = relCutPos;
    Double_t Integral = pH->Integral(pH->GetBin(Cut), pH->GetSize() + 1);
    Double_t DIntegral = sqrt(Integral);
//    efficiency[channel] = Integral / ( (pH->GetBin(Cut) - pH->GetBin(Pedestal)) * Ug + Integral );
    //*/
    return;
}

Double_t Hist::func2(Double_t *x, Double_t *par)
{   // polynom (deg 2) for fitting minima and maxima
    return par[0] * pow(x[0]-par[1], 2) + par[2];
}

void Hist::GetTimes()
{   // Sum up all entries in real-time and live-time histograms.
    // Include Overflow bin for the case of wrong binning.
    t_live = pHtLive->Integral(0, pHtLive->GetNbinsX()+1);
    t_real = pHtReal->Integral(0, pHtReal->GetNbinsX()+1);
}

void Hist::DeadTimeCorrection(Bool_t peak)
{
    GetTimes();
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {
        if(peak)
        {
            nNIF[i_ch] = nNIF_raw[i_ch] * t_real / t_live;
            DnNIF[i_ch] = DnNIF_raw[i_ch] * t_real / t_live;
        }
        nSF[i_ch] = nSF_raw[i_ch] * t_real / t_live;
        DnSF[i_ch] = DnSF_raw[i_ch] * t_real / t_live;
    }
}
/*
Double_t Hist::Fit2(TH1I *pH, Double_t xmin, Double_t xmax)
{
    TF1 *fit2 = new TF1("fit2", Hist::func2, xmin, xmax, 3);
    fit2->SetParNames("spread", "x", "y");
    fit2->SetParameter(0, -100000.0);
    fit2->SetParameter(1, (xmin+xmax)/2);
    fit2->SetParLimits(1, xmin, xmax);
    fit2->SetParameter(2, pH->GetBinContent(pH->GetBin((xmin+xmax)/2)));
    pH->Fit(fit2,"RB");
    return fit2->GetParameter("x");
}


Int_t Hist::FindMax(TH1I *pH, Int_t rough_pos, Int_t range)
{   // Find a maximum in (rough_pos-range, rough_pos+range).
    Int_t low = rough_pos-range;
    Int_t up = rough_pos+range;
    Double_t Max_fit = Fit2(pH, low, up);
    if(Max_fit < low || Max_fit > up)
    {
        cout << "Error: Maximum " << Max_fit << " not in (" << low << ";" << ")" << endl;
        return 0;
    }
    if(CommentFlag)
        cout << "Found maximum in (" << low << "," << up << ") at " << Max_fit << endl;
    return Max_fit;
}*/
/*
TH1F* Hist::GetTH1F(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile *f = TFile::Open(FilePath);
    if (f==0) //check if root file exists
    {   cerr << "File not found!" << endl;
        return 0;
    }
    else
    {   TH1F* histo = (TH1F*) f->Get(hname)->Clone();
        if (histo==0)
        {   cerr << "Histogram not found!" << endl;
            return 0;
        }
        else
            return histo;
    }
    f->Close();
}
*/
