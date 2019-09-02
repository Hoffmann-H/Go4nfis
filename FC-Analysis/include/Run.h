#ifndef RUN_H
#define RUN_H

#include <string>
#include "ToF.h"
#include "Plot.h"
#include "TTimeStamp.h"
#define NumCh 8
#define MaxHists 10

using namespace std;

class Run
{
public:
    Run(string FC, string setup, Int_t nr = 0, Bool_t draw = kFALSE);
    ~Run();
    Bool_t IsForeground();
    Double_t SetNeutronField(Double_t monitorCounts, Double_t monitorRelUnc, Double_t monitorTreal, Double_t d_start = 0, Double_t t_start = 0/*, Double_t l = 1500, Double_t Dl = 1*/);
    void SetToF(string name);
    void SetLimits(Int_t *pl0, Int_t *pl1, Int_t *pl2, Int_t *pl3);
//    void GetHist(Hist *hist);
//    void SetNatoms(Double_t *nAt, Double_t *DnAt);
//    void AnalyzeDt(Int_t i, Int_t lim0, Int_t lim1, Int_t lim2, Int_t lim3);
    Double_t GetMonitor() {return Monitor;}
//    void CrossSection(Int_t i);
    Double_t GetnfoverPhi(Int_t i);
    Double_t GetDnfoverPhi(Int_t i);
    string Name;
    string Setup;
//    Int_t nHist;
//    Hist *pH[MaxHists];
    ToF *pToF;
    Double_t Monitor, DMonitor;
    Double_t t_mon;
    Int_t tStart, tStop;
    Double_t NeutronFlux[NumCh];
    Double_t DNeutronFlux[NumCh];
    Double_t DstatNeutronFlux[NumCh];
    Double_t nNIF[NumHist];
    Double_t DnNIF[NumHist];
    Double_t nSF[NumHist];
    Double_t DnSF[NumHist];
    Double_t t_live;
//    Double_t uncCS[NumHist];
//    Double_t DuncCS[NumHist];

private:
    Bool_t CommentFlag;
    Bool_t DrawSingle;
    Bool_t DrawMulti;
    string FC;
    Bool_t UFC;
    Double_t t_real;
//    Double_t nAtoms[NumHist];
//    Double_t DnAtoms[NumHist];
};

#endif
