#ifndef HIST_H
#define HIST_H

#include <string>
#include <iostream>
#include <sstream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGo4WinCond.h"
#include "Plot.h"

#define NumHist 8

using namespace std;

class Hist
{
public:
    Hist(std::string file_path, string setup/* = "NIF"*/, string FC/* = "PuFC"*/, string name, Int_t run_nr);  // Constructor
//     Hist(); // Standard constructor, needed for object i/o
    ~Hist(); // Destructor
    void SetDraw(Plot *p);
    Bool_t IsForeground();

    void DoAnalyzeDt();
    void DoAnalyzeQDC();
    void SetNeutronField(Double_t monitorCounts, Double_t monitorRelUnc, Double_t monitorTreal, Double_t l = 1500, Double_t Dl = 1);
    Double_t GetPeakLow(Int_t i_ch);
    Double_t GetPeakUp(Int_t i_ch);
    void GetLimits(Double_t n = 3);
//    void SetAvBg(Double_t *avBg, Double_t *DavBg);
    void AnalyzeDtBg(Int_t i);
    void AnalyzeDtPeak(Int_t i);
    void Stability(Int_t i, TGraphErrors *ge1, TGraphErrors *ge2, Int_t tbins = 1890);
    Double_t SFvsTime(Int_t i, Int_t start, Int_t stop);
    Double_t NIFvsTime(Int_t i, Int_t start, Int_t stop);
    Double_t DNIFvsTime(Int_t i, Int_t start, Int_t stop);
    Double_t GetNevents(Int_t i);

//    Double_t SFRate[NumHist]; // SF rates assuming efficiency==100%
//    Double_t DSFRate[NumHist];
//    Double_t NIFRate[NumHist]; // NIF rates assuming efficiency==100%
//    Double_t DNIFRate[NumHist];
    string Name;
    Double_t t_live; // times
    Double_t t_real;
    Int_t t_start;
    Int_t t_stop;
    TH1I *pHDtG[NumHist];
    TH2I *pHDtGvsT[NumHist];
    Double_t t_mon;
    Double_t NeutronFlux[NumHist];
    Double_t DNeutronFlux[NumHist];
    Double_t DstatNeutronFlux[NumHist];
    Int_t RunNr;
    Double_t eInt[NumHist]; // intrinsic detection efficiency
    Double_t DeInt[NumHist];
    Double_t PedQDC[NumHist]; // extrema positions
    Double_t CutQDC[NumHist];
    Double_t MaxQDC[NumHist];
    Double_t Monitor, DMonitor;
    Double_t nNIF[NumHist];
    Double_t DnNIF[NumHist];
    Double_t nSF[NumHist];
    Double_t DnSF[NumHist];
private:
    // run parameters
    Bool_t CommentFlag;
    string Setup; // measurement setup. Possible values are "NIF", "SB", "SF".
    Bool_t UFC; // UFC==kFALSE: Analyze PuFC. UFC==kTRUE: Analyze UFC.
    string FilePath;
    TFile *file;
    Bool_t Draw;
    Plot *plot;
    // Neutron field
    Double_t Yield,   DYield;
    Double_t L,       DL;
    // analysis parameters
    Double_t DtPeakLow[NumCh];
    Double_t DtPeakUp[NumCh];
    Double_t DtBgLow;
    Double_t DtBgUp;
    Int_t lim[4][NumCh]; // integration limits' bin numbers
    // Variables
    Double_t avBg[NumHist]; // average constant background
    Double_t DavBg[NumHist];
    Double_t CutUsed[NumHist];
//    Double_t efficiency[NumHist]; // efficiency
    // Histograms
    TH1D *pHtLive;
    TH1D *pHtReal;
    TH1I *pHRawQDCl[NumHist];
    TH1I *pHAnaQDCl[NumHist];
//    TCanvas *pCQDClFit;
//    TCanvas *pCQDClCut;
    // Methods
    TH1I* GetTH1I(const char *hname);
    TH1D* GetTH1D(const char *hname);
    TH2I* GetTH2I(const char *hname);
    TGo4WinCond* GetWinCond(const char* hname);
    void SaveToFile(string path, TObject *pObj);
//    void SaveTo(string path, string name);
//    void SaveTH1I(TH1I* pHtoSave, string hpath, string hname);
//    void SaveTF1(TF1* pFtoSave, string fpath, string fname);
//    void SaveTN(TNamed* pNtoSave, string fpath, string fname);
//    void SavePad(TPad* pPtoSave, string ppath, string pname);
//    void SaveCanvas(TCanvas* pCtoSave, string hpath, string hname);
    void AnalyzeQDC(TH1I *pH, Int_t channel);
//    Double_t Fit2(TH1I *pH, Double_t xmin, Double_t xmax);
    static Double_t func2(Double_t *x, Double_t *p);
    void GetTimes();
//    void DeadTimeCorrection(Bool_t peak);
//    Int_t FindMax(TH1I *pH, Int_t rough_pos, Int_t range);
//    TH1F* GetTH1F(const char *hname);

//    void Save(TH1I* pHtoSave, string hpath, string hname);

//    void AnalyzePeaks();
//    void AnalyzeUndergrounds();
//    void AnalyzeUnderground(TH1I *pH, Double_t *pSF, Double_t *pDSF);
//    void GetTimes();
//    void DeadTimeCorrection();
//    Int_t GetMinBin(TH1I *pH, Int_t low, Int_t up);




//    ClassDef(Hist, 1);
};
#endif
