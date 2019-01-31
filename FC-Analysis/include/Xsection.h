#ifndef XSECTION_H
#define XSECTION_H

#include <string>
#include <iostream>
#include <sstream>

#include "Hist.h"
#include "TF1.h"
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
#include "TLegend.h"

class Xsection
{
public:
    Xsection();
    ~Xsection();
    
    void Evaluation1();
    void Evaluation2();
    void Evaluation3();
    void Evaluation4();

    Hist *pHNIF; // PuFC
    Hist *pHSB;
    Hist *pHSF;
    Hist *pHUNIF; // UFC
    Hist *pHUSB;
private:
    void CalculateThresholds();
    void PrintInScat();
    void ScatCorrNIF();
    void CalculateNPu();
    void CalculateCrossSection();
    void SaveToFile(string path, TObject *pObj);
    void CompareFiles(string path, Int_t start, Int_t stop);
    Bool_t CommentFlag;

    // physics parameters
    Double_t u;
    Double_t PuLit; // given Pu-242 (n,fis) cross section
    Double_t DPuLit;
    Double_t ULit;
    Double_t DULit;
    Double_t PuSFT2; // Spontaneaus fission halflife
    Double_t DPuSFT2;
    Double_t Yield;
    Double_t DYield;
    Double_t AreaPuFC;
    Double_t DAreaPuFC;
    Double_t AreaUFC;
    Double_t DAreaUFC;
    Double_t mPu; // given masses
    Double_t DmPu;
    Double_t mU;
    Double_t DmU;

    // physics variables
    Double_t MonitorNIF;
    Double_t DMonitorNIF;
    Double_t MonitorSB;
    Double_t DMonitorSB;
    Double_t MonitorUNIF;
    Double_t DMonitorUNIF;
    Double_t MonitorUSB;
    Double_t DMonitorUSB;
    Double_t SFRate[NumHist]; // SF rate for PuFC channels
    Double_t DSFRate[NumHist];
    Double_t NIFRate[NumHist]; // NIF rate for PuFC channels, NIF measurement, In-scattering corrected
    Double_t DNIFRate[NumHist];
    Double_t UNIFRate[NumHist]; // In-scattering corrected NIF rate for UFC channels
    Double_t DUNIFRate[NumHist];
    Double_t nPu[NumHist];
    Double_t DnPu[NumHist];
    Double_t XSec[NumHist];
    Double_t DXSec[NumHist];
    Double_t CrossSection;
    Double_t DCrossSection;
    Double_t eNIF[NumHist]; // efficiency
    Double_t DeNIF[NumHist];
    Double_t eSF;
    Double_t DeSF;
    Double_t eRel[NumHist]; // relative efficiency NIF:SF
    Double_t DeRel[NumHist];
};

#endif // XSECTION_H
