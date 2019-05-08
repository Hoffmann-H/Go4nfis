#ifndef PUFC_H
#define PUFC_H

#include <string>
#include "FC.h"

using namespace std;

class PuFC : public FC
{
public:
    Hist *pHSF;
    PuFC(Plot *p = 0);
    ~PuFC();
    void AnalyzeQDC() override;
    void HardCodedThresholds() override;
    void AnalyzeDtBG() override;
    void GetNatoms() override;
    void IsoVec() override;
    void ExpTrans() override;
    void SetDraw(Plot *p);
    void plt();

private:
    void InitPuVar();
    Double_t PuSFT2, DPuSFT2;
    Double_t cSF[NumHist];
    Double_t DcSF[NumHist];

};

#endif // PUFC_H
