#ifndef UFC_H
#define UFC_H

#include <string>
#include "FC.h"
using namespace std;

class UFC : public FC
{
public:
    UFC(Plot *p = 0);
    ~UFC();
    void AnalyzeQDC() override;
    void HardCodedThresholds() override;
    void AnalyzeDtBG() override;
    void GetNatoms() override;
    void IsoVec() override;
    void ExpTrans() override;
    void SetDraw(Plot *p) override;
//    void CalcCS(Int_t i);

private:
    void InitUVar();
    Double_t frac238, Dfrac238;
    Double_t frac235, Dfrac235;
    Double_t sigma238, Dsigma238;
};

#endif
