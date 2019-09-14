#ifndef SIM_H
#define SIM_H

#include <string>
#include "Plot.h"
#include "TGraph.h"
#include "TMath.h"
#define NumCh 8
#define m0 939.565 // MeV
#define c 299792458 // m/s

using namespace std;

class Sim
{
public:
    Sim(string file_name, string fc, string setup, Bool_t use_track_id, Bool_t draw_flag);
    ~Sim();
    void Calculate();
    void PrintEffN();
    Bool_t CommentFlag;
    string FileName;
    TFile *f;
    string FC, Setup;
    Bool_t PuFC;
    Bool_t tID;
    Bool_t DrawSingle;
    Bool_t DrawMulti;
    Plot *plot;

    // physics
    Double_t Distance[NumCh];
    TGraph *gPu242, *DgPu242, *gU235, *DgU235, *gU238, *DgU238;
    Double_t UisoVec[2], DUisoVec[2];
    Double_t En, DEn; // neutron source energy
    TH2F *pH2TvsE[NumCh][2];
    TH1F *pH1Eproj[NumCh][2]; // time-gated projection on energy axis
    TH1F *pH1Eeff[NumCh][2];
    TH1F *pH1Tproj[NumCh][2];
    TH1F *pH1ToF[NumCh][2];
    Double_t Emin, Emax;
    Int_t binEmin, binEmax;
    Double_t ToFmin[NumCh];
    Double_t ToFmax[NumCh];
    Int_t binToFmin[NumCh];
    Int_t binToFmax[NumCh];
    Int_t nProj[NumCh];

    Double_t nDirect[NumCh];
    Double_t DnDirect[NumCh];
    Double_t effDirect[NumCh];
    Double_t DeffDirect[NumCh];
    Double_t nScat[NumCh];
    Double_t DnScat[NumCh];
    Double_t effScat[NumCh];
    Double_t DeffScat[NumCh];
    Double_t CsDir[NumCh];
    Double_t DCsDir[NumCh];
    Double_t CsSc[NumCh];
    Double_t DCsSc[NumCh];

    Int_t nN(Int_t ch, Int_t binElow, Int_t binEup, Int_t Sc = 0);
    Double_t effN(Int_t ch, Int_t binElow, Int_t binEup, Int_t Sc = 0);
    Double_t DeffN(Int_t ch, Int_t binElow, Int_t binEup, Int_t Sc = 0);

    void relSigma(Double_t E, Double_t &w, Double_t &Dw);
    void SimToF();
private:
    void OpenHists();
    void SetDistances();
    void GetSigma(string path);
    void nProjectiles();
    void GetEwidth();
    void GetTwidth();
    void Projections();
    void DirectN();
    void DirectEff();
    void ScatN();
    void ScatEff();
    void SigmaEff();
    void SaveToFile(string path, TObject *pObj);

    Double_t ToF(Double_t Ekin, Double_t FlightPath);

    Bool_t  DoneHist,
            DoneDistance,
            DoneSigma,
            DoneProjectiles,
            DoneEwidth,
            DoneTwidth,
            DoneProjections,
            DoneDirect,
            DoneCalc;
};
#endif
