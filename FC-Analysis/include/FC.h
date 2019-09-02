#ifndef FC_H
#define FC_H

#include <string>
#include "Run.h"
#include "Hist.h"
#include "AnaSim.h"
#include "Plot.h"

#define NumCh 8
#define MaxRuns 5
#define MaxHist 10

using namespace std;

class FC
{
public:
    string Name;
//    string SimPath;
    Bool_t CommentFlag;
    Bool_t DrawSingle, DrawMulti;
    Plot *plot;
//    Int_t ScatterMethod;
    Int_t FgRuns;
    Int_t BgRuns;
    Run* pR[2 * MaxRuns];
    Int_t nHist;
    Hist *pH[2*MaxHist];
    TH1D *pHtLive;
    Double_t FgMon;
    Double_t BgMon;
    AnaSim *sim;
    FC(); //Double_t nFieldFG, Double_t DnFieldFG, Double_t nFieldBG, Double_t DnFieldBG);
    ~FC();

    // virtual methods: Differences between FC's
    virtual void HardCodedThresholds() = 0;//{} // hard-code parameters
//    virtual void AnalyzeQDC() = 0;//{}
//    virtual void GetDistance(Int_t i) = 0;//{} // hard-code parameters
//    virtual void GetSimTransmission(Int_t i) = 0;//{} // hard-code parameters
//    virtual void GetSimScattering(Int_t i) = 0;//{} // hard-code parameters
//    virtual void AnalyzeDtBG() = 0;//{cout << "virtual FC::AnalyzeDtBG. Don't call me" << endl;}
//    virtual void GetExpT() = 0;
    virtual void GetNatoms() = 0;
    virtual void DrawStability() = 0;
    virtual void IsoVec() = 0;
    void InitVar(Bool_t draw);
    void SetLimits(Int_t left = 12, Int_t right = 37);
//    void UseHists(Int_t start, Int_t stop, string setup, Int_t run);
//    void UseHist(string file_name, string setup, Int_t run);
//    void AnalyzeDt();
    void CrossSection();
//    void RegisterHists();
//    void Stability();
//    void ScatCorrDiff();
//    void ScatCorrFit();
    void ScatCorrSim();
//    void ExpTrans();
    void GetSimFg();
//    void GetSimBg();
    void Corrections();
//    void CompareShadowCone();
//    void CompareTransmission();

    // Workflow variables
    Bool_t DoneQDC,
           DoneThresholds,
           DoneLimits,
           DoneNatoms,
           DoneDtBG,
           DoneDt,
           DoneScatCorr,
           DoneSimFg,
           DoneSimBg,
           DoneIso,
           DoneRawCS,
           DoneCorrections,
           DoneTransmission;

    // Physics
    Double_t tFG, tBG;
    Double_t u;
    Double_t Area;
    Double_t DArea;
    Double_t nAtoms[NumCh];
    Double_t DnAtoms[NumCh];
    Double_t sd[NumCh]; // Distance source-detector
    Double_t Dsd[NumCh];

    // TimeDiff integration parameters
    Double_t DtPeakLow[NumCh];
    Double_t DtPeakUp[NumCh];
    Double_t DtBgLow;
    Double_t DtBgUp;

    // Analysis
    Int_t l0;
    Int_t l1[NumCh];
    Int_t l2[NumCh];
    Int_t l3;
    Double_t avBg[NumCh]; // average fission background per livetime and Dt bin
    Double_t DavBg[NumCh];
    Double_t nFG[NumCh]; // fission count correlated with Neutron pulse
    Double_t DnFG[NumCh];
    Double_t nBG[NumCh];
    Double_t DnBG[NumCh];
//    Double_t nFlux[NumCh];
//    Double_t DnFlux[NumCh];
//    Double_t cFG[NumCh]; // fission count Constant in time
//    Double_t DcFG[NumCh];
//    Double_t cBG[NumCh];
//    Double_t DcBG[NumCh];
//    Double_t nfDirect[NumCh]; // number of fissions induced by direct neutrons
//    Double_t DnfDirect[NumCh];
    Double_t sIsoVec; // correction subtrahend for isotope vector
    Double_t DsIsoVec;
    Double_t fIsoVec; // correction factor for isotope vector
    Double_t DfIsoVec;
    Double_t fTS[NumCh]; // Transmission and scattering correction factor
    Double_t DfTS[NumCh];
    Double_t uT[NumCh]; // NIF over SF rate
    Double_t DuT[NumCh];
    Double_t ExpT[NumCh]; // Experimental transmission factor
    Double_t DExpT[NumCh];
    Double_t SimT[NumCh]; // Simulated transmission factor
    Double_t DSimT[NumCh];
    Double_t pDirect[NumCh][3]; // Direct portion: ratio of direct to all fissions
    Double_t DpDirect[NumCh][3]; // 2nd index: 0 Experimental, 1 one sim, 2 shadow cone sim
    Double_t uCS[NumCh]; // un-corrected, raw cross section
    Double_t D2uCS[NumCh];
    Double_t CS[NumCh]; // corrected cross section
    Double_t DCS[NumCh];

private:
    void SaveToFile(string path, TObject *pObj);
    static Double_t func_peak(Double_t *x, Double_t *p);
    void ScatCorr(Int_t i);
};
#endif
