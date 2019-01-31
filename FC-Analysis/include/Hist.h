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

#define NumHist 8

using namespace std;

class Hist
{
public:
    Hist(std::string file_path, string setup = "NIF", string FC = "PuFC");  // Constructor
//     Hist(); // Standard constructor, needed for object i/o
    ~Hist(); // Destructor

    void DoAnalyzeDt();
    void DoAnalyzeQDC();
    void SetNeutronField(Double_t yield, Double_t Dyield, Double_t monitor, Double_t Dmonitor, Double_t l, Double_t Dl);
    Double_t GetNumberEvents(int i);

    Double_t SFRate[NumHist]; // SF rates assuming efficiency==100%
    Double_t DSFRate[NumHist];
    Double_t NIFRate[NumHist]; // NIF rates assuming efficiency==100%
    Double_t DNIFRate[NumHist];
    Double_t t_live; // times
    Double_t t_real;
    Double_t NeutronFlux[NumHist];
    Double_t DNeutronFlux[NumHist];
    Double_t eInt[NumHist]; // intrinsic detection efficiency
    Double_t DeInt[NumHist];
    Double_t PedQDC[NumHist];
    Double_t CutQDC[NumHist];
    Double_t MaxQDC[NumHist];
private:
    // run parameters
    Bool_t CommentFlag;
    string Setup; // measurement setup. Possible values are "NIF", "SB", "SF".
    Bool_t UFC; // UFC==kFALSE: Analyze PuFC. UFC==kTRUE: Analyze UFC.
//    string FC;
    string FilePath;
    TFile *file;
    // Neutron field
    Double_t Yield,   DYield;
    Double_t Monitor, DMonitor;
    Double_t L,       DL;
    // analysis parameters
    Double_t ToF_low;
    Double_t ToF_up;
    Double_t Dt_min;
    Double_t Dt_max;
    // Variables
    Double_t nNIF_raw[NumHist]; // not dead-time corrected
    Double_t DnNIF_raw[NumHist];
    Double_t nSF_raw[NumHist];
    Double_t DnSF_raw[NumHist];
    Double_t nNIF[NumHist]; // NIF count after dead-time correction
    Double_t DnNIF[NumHist];
    Double_t nSF[NumHist];
    Double_t DnSF[NumHist];
    Double_t UgQDC[NumHist];
    Double_t DUgQDC[NumHist];
    Double_t CutUsed[NumHist];
//    Double_t efficiency[NumHist]; // efficiency
    // Histograms
    TH1D *pHtLive;
    TH1D *pHtReal;
    TH1I *pHDtG[NumHist];
    TH1I *pHRawQDCl[NumHist];
    TH1I *pHAnaQDCl[NumHist];
    TCanvas *pCQDClFit;
    TCanvas *pCQDClCut;
    // Methods
    TH1I* GetTH1I(const char *hname);
    TH1D* GetTH1D(const char *hname);
    TGo4WinCond* GetWinCond(const char* hname);
    void SaveToFile(string path, TObject *pObj);
//    void SaveTo(string path, string name);
//    void SaveTH1I(TH1I* pHtoSave, string hpath, string hname);
//    void SaveTF1(TF1* pFtoSave, string fpath, string fname);
//    void SaveTN(TNamed* pNtoSave, string fpath, string fname);
    void SavePad(TPad* pPtoSave, string ppath, string pname);
    void SaveCanvas(TCanvas* pCtoSave, string hpath, string hname);
    Double_t AnalyzeDtPeak(Int_t i_ch, TH1I *pH);//, Double_t *pNIF, Double_t *pDNIF, Double_t *pSF, Double_t *pDSF);
    Double_t AnalyzeDtUnderground(Int_t i_ch, TH1I *pH);//, Double_t *pSF, Double_t *pDSF);
    void AnalyzeQDC(TH1I *pH, Int_t channel);
    Double_t Fit2(TH1I *pH, Double_t xmin, Double_t xmax);
    static Double_t func2(Double_t *x, Double_t *p);
    Int_t FindMax(TH1I *pH, Int_t rough_pos, Int_t range);
    void GetTimes();
    void DeadTimeCorrection(Bool_t peak);
/*    TH1F* GetTH1F(const char *hname);

    void Save(TH1I* pHtoSave, string hpath, string hname);

    void AnalyzePeaks();
    void AnalyzeUndergrounds();
    void AnalyzeUnderground(TH1I *pH, Double_t *pSF, Double_t *pDSF);
    void GetTimes();
    void DeadTimeCorrection();
    Int_t GetMinBin(TH1I *pH, Int_t low, Int_t up);



    */
//    ClassDef(Hist, 1);
};


#endif
