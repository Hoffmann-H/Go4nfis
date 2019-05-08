#include <string>
#include "Hist.h"
#include "TMath.h"
#define NumCh 8
#define m0 939.565 // MeV
#define c 299792458 // m/s

using namespace std;

class AnaSim
{
public:
    AnaSim(string file_name, string fc, string setup, Plot *p = 0);
    ~AnaSim();
    void OpenHists();
    void Projections();
    void DrawSigma();
//    void PeakForm(Bool_t scatter, Double_t *nFis, Double_t *DnFis);
    void nFissions();
    void nProjectiles();
    void Corrections();
    Bool_t CommentFlag;
    TFile *f;
    string FileName;
    string FC, Setup;
    Bool_t PuFC;
    Bool_t DrawSingle;
    Bool_t DrawMulti;
    Plot *plot;
    TH2F *pH2TvsE[NumCh];
    TH1F *pH1Eproj[NumCh]; // time-gated projection on energy axis
    Double_t Distance[NumCh];
    Double_t En, DEn; // neutron source eenrgy
    Int_t binEmin, binEmax;
    Double_t Emin, Emax;
    Double_t ToFmin[NumCh];
    Double_t ToFmax[NumCh];
    Int_t binToFmin[NumCh];
    Int_t binToFmax[NumCh];
    Double_t nProj[NumCh]; // number of projectiles
//    Double_t DnProj[NumCh];
    Double_t nDirect[NumCh]; // unscattered n
    Double_t DnDirect[NumCh];
    Double_t nFullE[NumCh]; // incident n in energy range
    Double_t DnFullE[NumCh];
    Double_t effFullE[NumCh];
    Double_t DeffFullE[NumCh];
    Double_t nFis[NumCh]; // effective n for fission
    Double_t DnFis[NumCh];
    Double_t T[NumCh];
    Double_t DT[NumCh];
    Double_t S[NumCh];
    Double_t DS[NumCh];
    Double_t F[NumCh];
    Double_t DF[NumCh];
//    Double_t S[NumCh];
//    Double_t DS[NumCh];

private:
    void SaveToFile(string RootPath, TObject *pObj);
    TH1D* CopyRange(TH1D* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset);
    Double_t Et(Double_t ToF, Double_t FlightPath);
    Double_t ToF(Double_t Ekin, Double_t FlightPath);
    void SetDistances();
    void SetQDCwidths();
    void SetToFwidths();
    void OpenSigma(string path);
    void relSigma(Double_t E, Double_t &w, Double_t &Dw);
//    Int_t nSource();
    void DirectN();
    void PeakForm(Int_t ch);
    void GetEwidth(TGraph *gE);

    TGraph *gPu242, *DgPu242, *gU235, *DgU235, *gU238, *DgU238;
    Double_t UisoVec[2], DUisoVec[2];

    Bool_t DoneHist,
           DoneProjection,
           DoneDistance,
           DoneToFwidth,
           DoneSigma,
           DoneProjectiles,
           DoneDirect,
           DoneFis,
           DoneCorrections;

};
