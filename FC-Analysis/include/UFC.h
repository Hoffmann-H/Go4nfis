#ifndef UFC_H
#define UFC_H

#include <string>
#include "FC.h"
using namespace std;

class UFC : public FC
{
public:
    UFC(Bool_t draw);
    ~UFC();
//    void AnalyzeQDC() override;
    void HardCodedThresholds() override;
//    void AnalyzeDtBG() override;
    void GetNatoms() override;
    void DrawStability() override;
    void IsoVec() override;
//    void GetExpT() override;
//    void Calibrate();
//    void CalcCS(Int_t i);


//private:
    void InitUVar();
    Double_t M;
    Double_t DM;
    Double_t frac238, Dfrac238;
    Double_t frac235, Dfrac235;
    Double_t n238[NumCh], Dn238[NumCh];
    Double_t n235[NumCh], Dn235[NumCh];
    Double_t emA[NumCh];
    Double_t DemA[NumCh];
    Double_t sigma238, Dsigma238;
};

#endif
