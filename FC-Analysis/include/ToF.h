#ifndef TOF_H
#define TOF_H

#include <string>
#include "Hist.h"
#define NumCh 8

using namespace std;

class ToF
{
public:
    ToF(string file_name, string fc, string setup);
    ~ToF();
    string GetFCname();
    Bool_t IsPuFC();
    void SetSave(Bool_t save) {Save = save;}
    void SetLimits(Int_t *l0, Int_t *l1, Int_t *l2, Int_t *l3);
    Double_t GetnfRate(Int_t i) {return nf[i] / t_live;}
    Double_t GetDnfRate(Int_t i) {return Dnf[i] / t_live;}
    string GetName() {return Name;}
    void DrawDtCh(Int_t i);
    void DrawDtCh();
    void DrawDt(Int_t i);
    void DrawDt();
    void DrawPeakLim(Int_t l, Int_t r_start, Int_t r_stop, Int_t r_step);
    Double_t nf[NumCh]; // neutron-induced fissions
    Double_t Dnf[NumCh];
    Double_t t_live;

private:
    TH1F* TH1ItoTH1F(TH1I* pH);
    void SpontaneousFission();
    void OpenToF(string file_name);
    void MakeLimits(Double_t left, Double_t right);
    void FitCommonBackground();
    void FitBackground();
    void FitBackground(Int_t  i);
    void SaveToFile(string path, TObject *pObj);
    void SubtractBackground();
    void SubtractBackground(Int_t i);
    void InducedFission();
    void InducedFission(Int_t i);
    void Print();

    string FilePath;
    string FC;
    string Setup;
    string Name;
    Bool_t PuFC;
    Bool_t Save;
    Int_t l0, l3;
    Double_t m[NumCh]; // Peak center
    Int_t l1[NumCh], l2[NumCh];
    Double_t ns_bin;
    TFile* f;
    TH1F* pH1Dt[NumCh];
    TH1F* pH1DtSub[NumCh];
    TF1* fTotal[NumCh];
    TF1* fLeft[NumCh];
    TF1* fRight[NumCh];
    TF1* fZero[NumCh];
    Double_t ug[NumCh];
    Double_t Dug[NumCh];
    Double_t unf[NumCh]; // manually implemented neutron-induced fissions
    Double_t Dunf[NumCh];
    Double_t sf[NumCh];
    Double_t Dsf[NumCh];
};

#endif
