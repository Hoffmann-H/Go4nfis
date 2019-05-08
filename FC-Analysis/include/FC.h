#ifndef FC_H
#define FC_H

#include <string>
#include "Hist.h"
#include "AnaSim.h"
#include "Plot.h"

#define NumHist 8

using namespace std;

class FC
{
public:
    string Name;
    string SimPath;
    Bool_t CommentFlag;
    Bool_t Draw;
    Plot *plot;
    Hist *pHFG;
    Hist *pHBG;
    AnaSim *sim;
    FC(); //Double_t nFieldFG, Double_t DnFieldFG, Double_t nFieldBG, Double_t DnFieldBG);
    ~FC();

    // virtual methods: Differences between FC's
    virtual void AnalyzeQDC() = 0;//{}
    virtual void HardCodedThresholds() = 0;//{} // hard-code parameters
//    virtual void GetDistance(Int_t i) = 0;//{} // hard-code parameters
//    virtual void GetSimTransmission(Int_t i) = 0;//{} // hard-code parameters
//    virtual void GetSimScattering(Int_t i) = 0;//{} // hard-code parameters
    virtual void AnalyzeDtBG() = 0;//{cout << "virtual FC::AnalyzeDtBG. Don't call me" << endl;}
    virtual void GetNatoms() = 0;
    virtual void IsoVec() = 0;
    virtual void ExpTrans() = 0;
    virtual void SetDraw(Plot *p) = 0;
    void InitVar();
    void GetLimits(Double_t n = 3);
    void AnalyzeDt();
    void ScatCorrDiff();
    void ScatCorrFit();
    void ScatCorrSim();
    void CrossSection();
    void GetSimRes();
    void Corrections();
    void CompareShadowCone(string BgPath);

    // Workflow variables
    Bool_t DoneQDC,
           DoneThresholds,
           DoneLimits,
           DoneNatoms,
           DoneDtBG,
           DoneDt,
           DoneScatCorr,
           DoneSim,
           DoneIso,
           DoneRawCS,
           DoneCorrections;

    // Physics
    Double_t u;
    Double_t Area;
    Double_t DArea;
    Double_t nFluenceFG;
    Double_t DnFluenceFG;
    Double_t nFluenceBG;
    Double_t DnFluenceBG;
    Double_t nAtoms[NumHist];
    Double_t DnAtoms[NumHist];
    Double_t sd[NumHist]; // Distance source-detector
    Double_t Dsd[NumHist];
    Double_t Yield, DYield;
    Double_t MonitorFG, DMonitorFG, MonitorBG, DMonitorBG;

    // TimeDiff integration parameters
    Double_t DtPeakLow[NumHist];
    Double_t DtPeakUp[NumHist];
    Double_t DtBgLow;
    Double_t DtBgUp;

    // Analysis
    Int_t lim[4][NumHist]; // integration limits' bin numbers
    Double_t avBg[NumHist]; // average fission background per livetime and Dt bin
    Double_t DavBg[NumHist];
    Double_t nFG[NumHist]; // fission count correlated with Neutron pulse
    Double_t DnFG[NumHist];
    Double_t nBG[NumHist];
    Double_t DnBG[NumHist];
    Double_t cFG[NumHist]; // fission count Constant in time
    Double_t DcFG[NumHist];
    Double_t cBG[NumHist];
    Double_t DcBG[NumHist];
//    Double_t nfDirect[NumHist]; // number of fissions induced by direct neutrons
//    Double_t DnfDirect[NumHist];
    Double_t sIsoVec; // correction subtrahend for isotope vector
    Double_t DsIsoVec;
    Double_t fIsoVec; // correction factor for isotope vector
    Double_t DfIsoVec;
    Double_t fTS[NumCh];
    Double_t DfTS[NumCh];
//    Double_t Transm[NumHist]; // transmission factor
//    Double_t DTransm[NumHist];
    Double_t pDirect[NumHist][3]; // Direct portion: ratio of direct to all fissions
    Double_t DpDirect[NumHist][3]; // 2nd index: 0 Experimental, 1 one sim, 2 shadow cone sim
//    Double_t nInc[NumHist]; // Incident neutron fluence = #neutrons / area
//    Double_t DnInc[NumHist];
    Double_t uCS[NumHist]; // un-corrected, raw cross section
    Double_t DuCS[NumHist];
    Double_t CS[NumHist]; // corrected cross section
    Double_t DCS[NumHist];

private:
    void SaveToFile(string path, TObject *pObj);
    static Double_t func_peak(Double_t *x, Double_t *p);
    void ScatCorr(Int_t i);
};
#endif
