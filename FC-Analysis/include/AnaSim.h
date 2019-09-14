#ifndef ANASIM_H
#define ANASIM_H

#include <string>
#include "Hist.h"
#include "Sim.h"
#include "TMath.h"
#define NumCh 8
#define m0 939.565 // MeV
#define c 299792458 // m/s

using namespace std;

class AnaSim
{
public:
    AnaSim(string fc, Bool_t use_track_id, Plot *p = 0);
    ~AnaSim();
    void Corrections();
    void ShadowCone();
    string FC, Setup;
    Bool_t PuFC;
    string FgPath, BgPath;
    Bool_t tID;
    Bool_t DrawSingle;
    Bool_t DrawMulti;
    Plot *plot;
    Bool_t CommentFlag;
    Sim *Fg;
    Sim *Bg;
    Double_t T[NumCh];
    Double_t DT[NumCh];
    Double_t S[NumCh];
    Double_t DS[NumCh];
    Double_t F[NumCh];
    Double_t DF[NumCh];
    Double_t SC[NumCh]; // Shadow Cone simulation
    Double_t DSC[NumCh];

private:
    void Uncertainties();
    Double_t dF_dnE(Double_t w, Int_t ch);
    Double_t dF_dwE(Double_t n, Int_t ch);
    Double_t dS_dnE(Double_t w, Int_t ch);
    Double_t dS_dwE(Double_t n, Int_t ch);
    void SaveToFile(string RootPath, TObject *pObj);
    TH1D* CopyRange(TH1D* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset);

    Bool_t DoneFg,
           DoneBg,
           DoneCorrections;

};

#endif
