#ifndef XSECTION_H
#define XSECTION_H

#include "Hist.h"
#include "TH1.h"
#include "TF1.h"

class Xsection
{
public:
    Xsection();
    ~Xsection();

//    void DoAnalyze();
    void CalculateNPu();
    void CalculateCrossSection();

    Hist *pHNIF;
    Hist *pHSB;
    Hist *pHSF;
    static Double_t xmin;
    static Double_t xmax;
private:
    Bool_t CommentFlag;

    // physics parameters
    Double_t PuSFT2; // Spontaneaus fission halflife
    Double_t DPuSFT2;
    Double_t Yield;
    Double_t DYield;
    Double_t Area;
    Double_t DArea;
//    Double_t NeutronFlux;

    // physics variables
    Int_t MonitorNIF;
    Double_t DMonitorNIF;
    Int_t MonitorSB;
    Double_t DMonitorSB;
    Double_t SFRate[NumHist];
    Double_t DSFRate[NumHist];
    Double_t NIFRate[NumHist];
    Double_t DNIFRate[NumHist];
    Double_t nPu[NumHist];
    Double_t DnPu[NumHist];
    Double_t XSec[NumHist];
    Double_t DXSec[NumHist];
    Double_t CrossSection;
    Double_t DCrossSection;
    void PrintInScat();
};

#endif // XSECTION_H
