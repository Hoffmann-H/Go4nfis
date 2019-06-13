#ifndef XSECTION_H
#define XSECTION_H

#include <string>
#include <iostream>
#include <sstream>

#include "PuFC.h"
#include "UFC.h"
#include "FC.h"
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
    Xsection(Bool_t draw = kFALSE);
    ~Xsection();
    void RelativeCS();
    PuFC *Pu;
    UFC *U;
    AnaSim *simPu;
    AnaSim *simU;
private:
    Bool_t CommentFlag;
    Bool_t Draw;
};

#endif // XSECTION_H
/*
void Evaluation1(); // PuFC
void Evaluation2(); // PuFC additional
void Evaluation3(); // PuFC additional
void Evaluation4(); // UFC

Hist *pHNIF; // PuFC
Hist *pHSB;
Hist *pHUG;
Hist *pHUNIF; // UFC
Hist *pHUSB;
private:
void SetParam();
void CalculateThresholds();
//    void UgPuFC();
//    void UgUFC();
void DoAnalyzeDt(string UFC, Bool_t method = 0);
void CalculateNPu();
//    void CalculateNU();
//    void CalculateEfficiency(); // no success
Double_t GetCorrectionFactor(string FC, Int_t i);
void CalculateCrossSection(); // PuFC
void SaveToFile(string path, TObject *pObj);
static Double_t func_peak(Double_t *x, Double_t *p);
//    void CompareFiles(string path, Int_t start, Int_t stop);
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
//    Double_t eSimGayther; // simulated efficiencies
//    Double_t DeSimGayther;
Double_t eSimMinimum;
Double_t DeSimMinimum;

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
Double_t USFRate[NumHist];
Double_t DUSFRate[NumHist];
Double_t nU235eff[NumHist]; // number of U atoms
Double_t DnU235eff[NumHist];
Double_t nU238eff[NumHist];
Double_t DnU238eff[NumHist];
Double_t nU[NumHist];
Double_t DnU[NumHist];
Double_t UXSec[NumHist];
Double_t DUXSec[NumHist];
Double_t nPuSF[NumHist]; // number of Pu atoms
Double_t DnPuSF[NumHist];
Double_t nPuNIF[NumHist];
Double_t DnPuNIF[NumHist];
Double_t nPu[NumHist];
Double_t DnPu[NumHist];
Double_t XSec[NumHist];
Double_t DXSec[NumHist];
//    Double_t CrossSection;
//    Double_t DCrossSection;
Double_t eNIF[NumHist]; // efficiency
Double_t DeNIF[NumHist];
Double_t eSF[NumHist];
Double_t DeSF[NumHist];
Double_t eRel[NumHist]; // relative efficiency NIF:SF
Double_t DeRel[NumHist];
Double_t eU[NumHist]; // simplified UFC efficiency
Double_t DeU[NumHist];
*/
