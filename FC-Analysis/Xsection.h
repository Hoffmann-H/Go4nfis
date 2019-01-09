#ifndef XSECTION_H
#define XSECTION_H
#include "Hist.h"
#include "TH1I.h"

class Xsection
{
public:
    Xsection();
    void DoAnalyze();
    void AnalyzePeak(TH1I *pH, Double_t *pNIF, Double_t *pDNIF, Double_t *pSF, Double_t *pDSF);
    void AnalyzeUnderground(TH1I *pH, Double_t *pSF, Double_t *pDSF);
    void GetTimes();
    void DeadTimeCorrection();
    void CalculateNPu();
    void CalculateCrossSection();
    void Fit2(TH1I *pH, Double_t xmin, Double_t xmax);
    static Double_t f2(Double_t *x, Double_t *p);
    Int_t GetMinBin(TH1I *pH, Int_t low, Int_t up);

    Hist *pHNIF;
    Hist *pHSB;
    TH1I *pNIF_DtG[NumHist];
    TH1D *pNIF_t_live;
    TH1D *pNIF_t_real;
    TH1I *pSB_DtG[NumHist];
    TH1D *pSB_t_live;
    TH1D *pSB_t_real;
    TH1D *pSF_t_live;
    TH1D *pSF_t_real;
    TH1I *pHRawQDCl[NumHist];
    TH1I *pHAnaQDCl[NumHist];
    TCanvas *pNIFNPu;
    TCanvas *pSBNPu;
    static Double_t xmin;
    static Double_t xmax;
private:
    Bool_t CommentFlag;
    // analysis parameters
    Double_t ToF_low;
    Double_t ToF_up;
    Double_t Dt_min;
    Double_t Dt_max;
    // physics parameters
    Double_t PuSFT2;
    Double_t DPuSFT2;
    Double_t Neutron_flux;
    // physics variables
    // '1': NIF measurement, '2': Shadow bar, '3': Spontaneaus fission.
    // 'SF': Spontaneaus fission, 'NIF': neutron-induced fission
    // 'raw': without, 'dtc': with Death-Time Correction
    Double_t N1SF_raw[NumHist]; // not dead-time corrected
    Double_t DN1SF_raw[NumHist];
    Double_t N1NIF_raw[NumHist];
    Double_t DN1NIF_raw[NumHist];
    Double_t N2SF_raw[NumHist];
    Double_t DN2SF_raw[NumHist];
    Double_t N2NIF_raw[NumHist];
    Double_t DN2NIF_raw[NumHist];
    Double_t N3SF_raw[NumHist];
    Double_t DN3SF_raw[NumHist];
    Double_t N1SF_dtc[NumHist]; // after dead-time correction
    Double_t DN1SF_dtc[NumHist];
    Double_t N1NIF_dtc[NumHist];
    Double_t DN1NIF_dtc[NumHist];
    Double_t N2SF_dtc[NumHist];
    Double_t DN2SF_dtc[NumHist];
    Double_t N2NIF_dtc[NumHist];
    Double_t DN2NIF_dtc[NumHist];
    Double_t N3SF_dtc[NumHist];
    Double_t DN3SF_dtc[NumHist];
    Double_t t_live_SB; // times
    Double_t t_real_SB;
    Double_t t_live_NIF;
    Double_t t_real_NIF;
    Double_t t_live_SF;
    Double_t t_real_SF;
    Double_t NPu[NumHist];
    Double_t DNPu[NumHist];
    Double_t NPuNIF[NumHist];
    Double_t DNPuNIF[NumHist];
    Double_t NPuSB[NumHist];
    Double_t DNPuSB[NumHist];
    Double_t CrossSection;
    Double_t DCrossSection;
};

#endif // XSECTION_H
