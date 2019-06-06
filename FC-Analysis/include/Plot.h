#ifndef PLOT_H
#define PLOT_H

#include <string>
#include <iostream>
#include <stdlib.h>
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
#define NumCh 8
using namespace std;

class Plot
{
public:
    Plot();
    Plot(string FC, string setup, TFile *file = 0);
    ~Plot();
    string FC, Setup;
    TFile *f;
    void simply(TH1F *pH, string Name, string xTitle, string yTitle);
    void Sigma(TGraph* gPu242, TGraph* gU235, TGraph* gU238);
    void Source_E(TGraph *g, Double_t E, Double_t DE, Int_t Variante);
    void Source_E_v0(TGraph *g, Double_t E, Double_t DE);
    void Source_E_v1(TGraph *g, Double_t E, Double_t DE);
    void QDCfit(Int_t ch, TH1I *pH, Double_t cut, TF1 *fCut);
    void QDCeff(Int_t ch, TH1I *pH, Double_t pedestal, Double_t cut, Double_t level, Double_t eInt, Double_t DeInt);
    void Dt(Int_t ch, TH1I *pH, Double_t nf, Double_t Dnf, Int_t l0, Int_t l1, Int_t l2, Int_t l3, Double_t level, string Setup);
    void ExpT(Double_t* expT, Double_t* DexpT, Double_t* simT, Double_t* DsimT);
    void SimF(Double_t* T, Double_t* DT, Double_t* S, Double_t* DS, Double_t* F, Double_t* DF);
    void SimF(Double_t* T, Double_t* DT, Double_t* S, Double_t* DS, Double_t* F, Double_t* DF, Double_t* T2, Double_t* DT2, Double_t* S2, Double_t* DS2, Double_t* F2, Double_t* DF2);
    void CompT(Double_t* pExpT, Double_t* DpExpT, Double_t* pSimT, Double_t* DpSimT);
    void CompSc(Double_t* pExp, Double_t* DpExp, Double_t* pSim1, Double_t* DpSim1, Double_t* pSim2, Double_t* DpSim2);
    void TvsE(Int_t ch, TH2F *pH2TvsE, Double_t binToFmin, Double_t binToFmax, Int_t directN);
    void TvsE(Int_t ch, TH2F *pH2TvsE, Double_t binToFmin, Double_t binToFmax, Double_t Emin, Double_t Emax);
    void DtPeakForm(Int_t ch, TH1F *pSc, TH1F *pUnsc, Double_t Emin, Int_t Projectiles);
    void Eproj(Int_t ch, TH1F *pSc, TH1F *pUnsc, Int_t Projectiles);
    void Eproj(Int_t ch, TH1F *pSc, TH1F *pUnsc, Int_t Projectiles, TGraph *gSigma);
    void Result(Double_t* uCS, Double_t* DuCS, Double_t* CS, Double_t* DCS);
//    void T();
private:
    TH1F* CopyBins(TH1I *pH, char *name, Int_t b0, Int_t b1, Double_t yoffset);
    TH1F* CopyRange(TH1I* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset);
    void SaveToFile(string RootPath, TObject *pObj);
    Int_t round_digits(Double_t x, Int_t digits);
    Int_t x[NumCh];
    Int_t xerr[NumCh];
};
#endif // PLOT_H
