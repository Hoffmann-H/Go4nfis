// Class for loading and accessing histograms from a root file. 
// This is meant to be created once for each exp. setup: NIF, SB, SF.
// Open histograms: Hist *pHist = new Hist("path/to/file.root")

#include "Hist.h"

#define NumHist 8

using namespace std;

Hist::Hist(string file_path, string setup, string FC)   // Constructor
{
    CommentFlag = kFALSE;
    Setup = setup;
    UFC = strcmp(FC.c_str(), "PuFC");

    if(CommentFlag)
        cout << endl << "Creating instance Hist  " << (UFC?"UFC":"PuFC") << "  " << Setup << "  " << file_path << endl;
    //define surrounding peak limits
    Dt_min = 63600;
    Dt_max = 78500;
    //define time window for underground integration
    if(!UFC)
    { // PuFC
        ToF_low = 66500;
        ToF_up  = 68500;
    } else { // UFC
        ToF_low = 72400;
        ToF_up = 73600;
    }

    FilePath = file_path;

    pHtLive = GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_47");
    pHtReal = GetTH1D("/Histograms/Raw/Scaler/Rates/H1RawRate_48");
    char hist_name[64] = "";
    for(int i_ch = 0; i_ch < 8; i_ch++) {
        sprintf(hist_name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i_ch + 1);
        pHDtG[i_ch] = GetTH1I(hist_name);
        sprintf(hist_name, "/Histograms/Raw/QDC/low/H1RawQDCl_%i", i_ch + 1);
        pHRawQDCl[i_ch] = GetTH1I(hist_name);
        sprintf(hist_name, "/Histograms/Analysis/FC/QDC/low/trig/H1AnaQDCl_trig_%i", i_ch + 1);
        pHAnaQDCl[i_ch] = GetTH1I(hist_name);
        sprintf(hist_name, "/Conditions/Analysis/QDC/QDCl_Thr_Ch%i", i_ch + 1);
        TGo4WinCond* con = GetWinCond(hist_name);
        CutUsed[i_ch] = con->GetXLow();
    }
    GetTimes();

}

Hist::~Hist()    // Destructor
{
    if(CommentFlag)
        cout << "Destroy" << endl;
}

TH1I* Hist::GetTH1I(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile* file = TFile::Open(FilePath.c_str());
    TH1I* histo = (TH1I*) file->Get(hname)->Clone();
    if(CommentFlag)
        cout << "Opening histogram " << hname << endl;
    if (histo==0)
    {   cout << "Histogram " << hname << " not found!" << endl;
        return 0;
    }
    else
        return histo;
    file->Close();
}

TH1D* Hist::GetTH1D(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile* file = TFile::Open(FilePath.c_str());
    TH1D* histo = (TH1D*) file->Get(hname)->Clone();
    if(CommentFlag)
        cout << "Opening histogram " << hname << endl;
    if (histo==0)
    {   cerr << "Histogram " << hname << " not found!" << endl;
        return 0;
    }
    else
        return histo;
    file->Close();
}

TGo4WinCond* Hist::GetWinCond(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile* file = TFile::Open(FilePath.c_str());
    TGo4WinCond* obj = (TGo4WinCond*) file->Get(hname)->Clone();
    if(CommentFlag)
        cout << "Opening Condition " << hname << endl;
    if (obj==0)
    {   cout << "Condition " << hname << " not found!" << endl;
        return 0;
    }
    else
        return obj;
    file->Close();
}

void Hist::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    TFile* file = TFile::Open(FilePath.c_str(), "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder "Analysis" already exists, otherwise create it
    if (file->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) file->Get(path.c_str());
    else
    {   file->mkdir(path.c_str(), "Folder containing offline Analysis objects");
        EvalDir = file->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    file->Save(); file->Close();
}


void Hist::SetNeutronField(Double_t yield, Double_t Dyield, Double_t monitor, Double_t Dmonitor, Double_t l, Double_t Dl)
{
    Yield = yield;
    DYield = Dyield;
    Monitor = monitor;
    DMonitor = Dmonitor;
    L = l;
    DL = Dl;
    Double_t r;
    if(CommentFlag)
        cout << "Neutron flux density" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        if (UFC)
            r = L + 45 + 70 - 10 * i;
        else
            r = L + 45 + 35 - 5 * i;
        NeutronFlux[i] = Yield * Monitor / (r*r * t_real);
        DNeutronFlux[i] = sqrt( pow(DYield * Monitor / (r*r * t_real), 2) +
                                pow(Yield * DMonitor / (r*r * t_real), 2) +
                                pow(Yield * Monitor * 2*DL / (r*r*r * t_real), 2) );
        if(CommentFlag)
            cout << "ch " << i+1 << ": " << NeutronFlux[i] << " +- " << DNeutronFlux[i] << endl;
    }
}


Double_t Hist::GetPeakLow(Int_t i_ch)
{
    char name[64] = "";
    sprintf(name, "/Conditions/Analysis/ToF/ToF_Thr_Ch%i", i_ch + 1);
    TGo4WinCond* con = GetWinCond(name);
    return con->GetXLow();
}


Double_t Hist::GetPeakUp(Int_t i_ch)
{
    char name[64] = "";
    sprintf(name, "/Conditions/Analysis/ToF/ToF_Thr_Ch%i", i_ch + 1);
    TGo4WinCond* con = GetWinCond(name);
    return con->GetXUp();
}

/*
void Hist::DoAnalyzeDt()
{   // Analyze TimeDiff spectra for all channels: Determine underground and, if the beam was on, peak content. Setup == "SF" for measurements without beam.
    Bool_t peak;
    if (strcmp(Setup.c_str(), "SF") == 0) // decide if to analyze a Dt peak
        peak = kFALSE;
    else
        peak = kTRUE;
    string FC;
    if(UFC)
        FC = "UFC";
    else
        FC = "PuFC";
    char hname[64] = "";
    if(CommentFlag)
        cout << endl;
    for (int i_ch = 0; i_ch < NumHist; i_ch++)
    { // Analyze one channel's TimeDiff spectrum
        if(CommentFlag)
            cout << "Channel " << i_ch+1 << endl;
        Double_t niveau;
        if(peak)
            niveau = AnalyzeDtPeak(i_ch, pHDtG[i_ch]);
        else
            niveau = AnalyzeDtUnderground(i_ch, pHDtG[i_ch]);
        Double_t x[] = {Dt_min, Dt_max};
        Double_t y[] = {niveau, niveau};
        TGraph* g0 = new TGraph(2, x, y);
        sprintf(hname, "f%s%sUg_%i", FC.c_str(), Setup.c_str(), i_ch+1);
        g0->SetNameTitle(hname, "Pulse-heigt gated TimeDiff Underground Fit");
        SaveToFile("Analysis/TimeDiff/Ug", g0);
        x[0] = ToF_low;
        x[1] = ToF_up;
        TGraph* g1 = new TGraph(2, x, y);
        sprintf(hname, "f%s%sUgPeak_%i", FC.c_str(), Setup.c_str(), i_ch+1);
        g1->SetNameTitle(hname, "TimeDiff Peak area Underground");
        SaveToFile("Analysis/TimeDiff/Ug", g1);
    }

    DeadTimeCorrection(peak);

    if (CommentFlag)
        cout << "Times: " << t_real << ", " << t_live << endl
             << "ch   NIF   SF   SF rate" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        SFRate[i] = nSF[i] / t_real;
        DSFRate[i] = DnSF[i] / t_real;
        if(CommentFlag)
            cout << " " << i+1 << "  " << nNIF[i] << "+-" << DnNIF[i] << "  " << nSF[i] << "+-" << DnSF[i] << "  " << SFRate[i] << "+-" << DSFRate[i] << endl;
    }

    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    TGraphErrors* g0 = new TGraphErrors(NumHist, x, SFRate, xerr, DSFRate);
    g0->SetNameTitle("SF_Rate", "SF detection rate");
    SaveToFile("Analysis/TimeDiff", g0);

    TGraphErrors* g5 = new TGraphErrors(NumHist, x, nSF, xerr, DnSF);
    g5->SetNameTitle("SF", "SF count");
    SaveToFile("Analysis/TimeDiff", g5);

    if(peak)
    {   // If measured with beam, draw NIF rates, NIF/SF rates ratio and NIF/Flux ratio
        Double_t Ratio_SF[NumHist];
        Double_t DRatio_SF[NumHist];
        Double_t Ratio_Flux[NumHist];
        Double_t DRatio_Flux[NumHist];
        for (int i = 0; i < NumHist; i++)
        {
            NIFRate[i] = nNIF[i] / t_real;
            DNIFRate[i] = DnNIF[i] / t_real;
            Ratio_SF[i] = nNIF[i] / nSF[i];
            DRatio_SF[i] = sqrt( pow(DnNIF[i] / nSF[i], 2) +
                              pow(nNIF[i] * DnSF[i] / pow(nSF[i], 2), 2) );
            Ratio_Flux[i] = NIFRate[i] / NeutronFlux[i];
            DRatio_Flux[i] = sqrt( pow(DNIFRate[i] / NeutronFlux[i], 2) +
                                   pow(NIFRate[i] * DNeutronFlux[i] / pow(NeutronFlux[i], 2), 2) );
        }
        TGraphErrors* g1 = new TGraphErrors(NumHist, x, nNIF, xerr, DnNIF);
        g1->SetNameTitle("NIF", "NIF count");
        SaveToFile("Analysis/TimeDiff", g1);

        TGraphErrors* g2 = new TGraphErrors(NumHist, x, NIFRate, xerr, DNIFRate);
        g2->SetNameTitle("NIF_Rate", "NIF detection rate");
        SaveToFile("Analysis/TimeDiff", g2);

        TGraphErrors* g3 = new TGraphErrors(NumHist, x, Ratio_SF, xerr, DRatio_SF);
        g3->SetNameTitle("NIF_SF_Ratio", "NIF to SF detection ratio");
        SaveToFile("Analysis/TimeDiff", g3);

        TGraphErrors* g4 = new TGraphErrors(NumHist, x, Ratio_Flux, xerr, DRatio_Flux);
        g4->SetNameTitle("NIF_Flux_Ratio", "NIF detection to neutron flux ratio");
        SaveToFile("Analysis/TimeDiff", g4);
    }
    return;
}

Double_t Hist::AnalyzeDtPeak(Int_t i_ch, TH1I *pH)
{
//    TH1I* pH = (TH1I*)pHDtG[i_ch]->Clone();
    Double_t ChPerBin = pH->GetBinWidth(0);
    Double_t BinOffset = pH->GetBinCenter(0);
    // Calculate integration limit bins
    Double_t lim_0 = (Dt_min - BinOffset) / ChPerBin;
    Double_t lim_1 = (ToF_low - BinOffset) / ChPerBin;
    Double_t lim_2 = (ToF_up - BinOffset) / ChPerBin;
    Double_t lim_3 = (Dt_max - BinOffset) / ChPerBin;
    if(CommentFlag)
        cout << "Analyzing Dt peak. Ch per bin: " << ChPerBin << ", Offset: " << BinOffset << endl <<
                " Integration limits: Bin nr. " << lim_0 << " " << lim_1 << " " << lim_2 << " " << lim_3 << endl;
    // Integrations
    Int_t PeakCount = pH->Integral(lim_1, lim_2); // Sum within peak region
    Double_t UgCount = pH->Integral(lim_0, lim_3) - PeakCount; // Sum within underground region; follows Poisson-stat
    Double_t DUgCount = sqrt(UgCount); // Underground deviation
    Double_t UgPerBin = 1.0 * UgCount / (lim_1 - lim_0 + lim_3 - lim_2);
    Double_t DUgPerBin = 1.0 * DUgCount / (lim_1 - lim_0 + lim_3 - lim_2);
    nNIF_raw[i_ch] = PeakCount - (lim_2 - lim_1 + 1) * UgPerBin;
    DnNIF_raw[i_ch] = (lim_2 - lim_1 + 1) * DUgPerBin;
    nSF_raw[i_ch] = pH->Integral() - nNIF_raw[i_ch]; // Calculate total underground
    DnSF_raw[i_ch] = sqrt(pow(DnNIF_raw[i_ch], 2) + pH->Integral()); // Error according to Gauss
    if(CommentFlag)
        cout << "NIF: " << nNIF_raw[i_ch] << " +- " << DnNIF_raw[i_ch] << endl <<
                "SF:  " << nSF_raw[i_ch]  << " +- " << DnSF_raw[i_ch]  << endl;
    return UgPerBin;
}

Double_t Hist::AnalyzeDtUnderground(Int_t i_ch, TH1I *pH)
{
//    TH1I* pH = (TH1I*)pHDtG[i_ch];
    Double_t ChPerBin = pH->GetBinWidth(0);
    Double_t BinOffset = pH->GetBinCenter(0);
    // Calculate integration limit bins
    Int_t lim_0 = (Dt_min - BinOffset) / ChPerBin;
    Int_t lim_3 = (Dt_max - BinOffset) / ChPerBin;
//    cout << lim_0 << " " << lim_3 << " " << ChPerBin << " " << pH->GetBinContent(lim_0) << " " << pH->GetBinContent(lim_0+1) << endl;
    // Analyze underground analogue to AnalyzePeak(), but without peak area exclusion.
    Double_t SF = pH->Integral();
    Double_t DSF = sqrt(SF);
    nSF_raw[i_ch] = SF;
    DnSF_raw[i_ch] = DSF;
    if(CommentFlag)
        cout << "SF raw " << nSF_raw[i_ch] << " +- " << DnSF_raw[i_ch] << endl;
    return pH->Integral(lim_0, lim_3) / (lim_3 - lim_0);
}
*/
void Hist::DoAnalyzeQDC()
{
    Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8};
    Double_t xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
    Double_t relCutPos[NumHist];

    ///Analyze raw QDC spectrum
    for (int i = 0; i < NumHist; i++)
    {
        AnalyzeQDC(pHRawQDCl[i], i);
        relCutPos[i] = (CutQDC[i] - PedQDC[i]) / (MaxQDC[i] - PedQDC[i]);
    }
    // print minima
    cout << (UFC?"UFC":"PuFC") << "  " << Setup << "  " << "Minima: ";
    for (int i = 0; i < 8; i++)
        cout << CutQDC[i] << "  ";
    cout << endl;

    TGraph* g4 = new TGraph(NumHist, x, CutQDC);
    g4->SetMarkerStyle(20);
    g4->SetNameTitle("gCut", "Absolute minimum position");
    SaveToFile("Analysis/QDC", g4);

    TGraph* g5 = new TGraph(NumHist, x, MaxQDC);
    g5->SetMarkerStyle(20);
    g5->SetNameTitle("gMax", "Absolute maximum position");
    SaveToFile("Analysis/QDC", g5);

    TGraph* g6 = new TGraph(NumHist, x, relCutPos);
    g6->SetMarkerStyle(20);
    g6->SetNameTitle("gRelCut", "Relative minimum position");
    SaveToFile("Analysis/QDC", g6);

    TGraph* g7 = new TGraph(NumHist, x, PedQDC);
    g6->SetMarkerStyle(20);
    g7->SetNameTitle("gPed", "Absolute pedestal position");
    SaveToFile("Analysis/QDC", g7);

    TGraphErrors* g8 = new TGraphErrors(NumHist, x, eInt, xerr, DeInt);
    g8->SetNameTitle("eInt", "Intrinsic efficiency");
    SaveToFile("Analysis/QDC", g8);

    return;
}

void Hist::AnalyzeQDC(TH1I *pH, Int_t channel)
{   // Analyze single QDC spectrum. Use hard-coded fit positions acc. to channel. For self-triggered QDC spectra, use pedestal from before
    // fit areas channel        0,    1,    2,    3,    4,    5,    6,    7
    Int_t FitRange[4][8] = {{ 800,  700,  750,  700,  800,  750,  750,  700}, // PuFC cut min
                            {1300, 1200, 1250, 1200, 1300, 1250, 1250, 1200}, // PuFC cut max
                            { 100,  350,  100,  100,  100,   50,  100,   50}, // UFC cut min
                            { 400,  900,  400,  400,  400,  350,  400,  350} }; // UFC cut max

    /// 1. Find pedestal
    Double_t Pedestal;
    Double_t DPedestal;
    char fname[64] = "";
    char hname[64] = "";
    sprintf(hname, "Analysis/QDC/Fit");
    string FC;
    if(UFC)
        FC = "UFC";
    else
        FC = "PuFC";
    Int_t low, up;
    if(UFC) { // if analyzing UFC, don't fit but find pedestal.
        Int_t PBin = pH->GetMaximumBin();
        Pedestal = pH->GetBinCenter(PBin);
        DPedestal = pH->GetBinWidth(PBin);
    } else { // If analyzing PuFC, fit pedestal position
        Int_t PBin = pH->GetMaximumBin(); // bin number of maximum
        Int_t Max = pH->GetBinContent(PBin); // maximum value
        sprintf(fname, "fPuFC%sPed_%i", Setup.c_str(), channel + 1);
        // Bounds (low, up): Find limits of FWHM
        low = PBin-1;
        up = PBin+1;
        while(2 * pH->GetBinContent(low-1) > Max) {
            low--;
        }
        while(2 * pH->GetBinContent(up+1) > Max) {
            up++;
        }
        TF1* fP = new TF1(fname, func2, pH->GetBinCenter(low), pH->GetBinCenter(up), 3);
        fP->SetRange(pH->GetBinCenter(low), pH->GetBinCenter(up));
        fP->SetParameters(-100000.0, (Double_t)pH->GetBinCenter(PBin), (Double_t)Max);
        pH->Fit(fname, "RBQ");
        Pedestal = (Double_t)fP->GetParameter(1);
        DPedestal = (Double_t)fP->GetParError(1);
        SaveToFile(hname, fP);
    }
    PedQDC[channel] = Pedestal;

    /// 2. Find (FF,alpha)-cut and underground level
    low = FitRange[UFC?2:0][channel];
    up = FitRange[UFC?3:1][channel];
    if(CommentFlag)
        cout << "low: " << low << ", up: " << up << endl;
    sprintf(fname, "f%s%sCut_%i", FC.c_str(), Setup.c_str(), channel + 1);
    TF1* fC = new TF1(fname, "pol4", low, up);
    fC->SetRange(low, up);
    fC->SetParameters(1.E3, -10.0, 0.01, -1.E-5, 1.E-8);
    pH->Fit(fname, "RBQ"); //"RBQ"
    Double_t Cut = fC->GetMinimumX();
    SaveToFile(hname, fC);

    // Find underground level
    low = pH->GetBin(Cut);
    while (fC->Eval(pH->GetBinCenter(low), 0, 0, 0) < 1.1 * fC->Eval(Cut))
        low--;
    up = pH->GetBin(Cut);
    while (fC->Eval(pH->GetBinCenter(up), 0, 0, 0) < 1.1 * fC->Eval(Cut))
        up++;
    Double_t UgLevel = pH->Integral(low, up) / (up - low + 1.0);
    Double_t DUgLevel = sqrt(pH->Integral(low, up)) / (up - low + 1);
    if(CommentFlag)
        cout << "low: " << low << ", up: " << up << ", level: " << UgLevel << " +- " << DUgLevel << endl;
    if (pH->GetBinCenter(low) > CutUsed[channel] || pH->GetBinCenter(up) < CutUsed[channel])
        cout << "*** " << FC << "," << Setup << ",ch" << channel+1 <<
                " Warning: Individual minimum region " << pH->GetBinCenter(low) << "," << pH->GetBinCenter(up) <<
                " does not fit into used cut " << CutUsed[channel] << endl;

    // Find internal efficiency
    Double_t Integral = pH->Integral(CutUsed[channel], 4096);
    eInt[channel] = Integral / (Integral + UgLevel * (CutUsed[channel] - Pedestal));
    DeInt[channel] = sqrt( Integral * pow(UgLevel * (CutUsed[channel] - Pedestal), 2) +
                           pow(DUgLevel * Integral * (CutUsed[channel] - Pedestal), 2)
                         ) / pow(Integral + UgLevel * (CutUsed[channel] - Pedestal), 2);

    Double_t x[] = {Pedestal, Pedestal, pH->GetBinCenter(up), pH->GetBinCenter(up)};
    Double_t y[] = {0, UgLevel, UgLevel, 0};
//    cout << "Check" << endl;
    TGraph* g1 = new TGraph(4, x, y);
    sprintf(fname, "f%s%sUg_%i", FC.c_str(), Setup.c_str(), channel+1);
    g1->SetNameTitle(fname, "QDC (low gain) constant level fit");
    SaveToFile("Analysis/QDC/Ug", g1);

    /// 3. Fit FF-maximum
    // find smeared FF-maximum
    Int_t xMax = 0, yMax = 0;
    Int_t smooth = 5;
    for(Int_t bin = FitRange[UFC?3:1][channel]; bin < 2000; bin++)
    {   // search for bin with maximum content in range
        Double_t y = pH->Integral(bin-smooth, bin+smooth) / (2*smooth+1);
        if(y > yMax)
        {
            xMax = bin;
            yMax = (Int_t)y;
        }
    }
    // Fit fission fragment maximum
    if(UFC)
    {
        low = xMax * 0.8;
        up = xMax * 1.4;
    } else {
        low = xMax * 0.9;
        up = xMax * 1.2;
    }
    sprintf(fname, "f%s%sMax_%i", FC.c_str(), Setup.c_str(), channel + 1);
    TF1* fM = new TF1(fname, "pol4", low, up);
    fM->SetRange(low, up);
    fM->SetParameters(-1.E4, 40.0, -0.01, 1.E-7, -1.E-10);
    pH->Fit(fname, "RBQ");
    Double_t Max = fM->GetMaximumX();
    SaveToFile(hname, fM);
    fC->Delete();
    if(CommentFlag)
        cout << "Cut " << Cut << ", Max " << Max << endl;

    // 4. Write values for graphs. Calculate efficiency
//    UgQDC[channel] = Ug;
//    DUgQDC[channel] = DUg;
    CutQDC[channel] = Cut;
    MaxQDC[channel] = Max;
//    Double_t Integral = pH->Integral(pH->GetBin(Cut), pH->GetSize() + 1);
//    Double_t DIntegral = sqrt(Integral);
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
//    GetTimes();
    for(int i_ch = 0; i_ch < NumHist; i_ch++)
    {
        if(peak) {
            nNIF[i_ch] = nNIF_raw[i_ch] * t_real / t_live;
            DnNIF[i_ch] = DnNIF_raw[i_ch] * t_real / t_live;
        } else {
            nNIF[i_ch] = 0;
            DnNIF[i_ch] = 0;
        }
        nSF[i_ch] = nSF_raw[i_ch] * t_real / t_live;
        DnSF[i_ch] = DnSF_raw[i_ch] * t_real / t_live;
    }
}

Double_t Hist::GetNumberEvents(int i)
{ // return total number of self-triggered events for uncertainty weighting purposes
    Double_t number = (Double_t)pHAnaQDCl[i]->Integral();
    return number;
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
    file = TFile::Open(FilePath, "UPDATE");
    SaveTo(hpath, hname);
    // save
    pHtoSave->Write(hname.c_str());
    file->Save();
    file->Close();
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
*/
