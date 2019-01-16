#ifndef NTOTHISTOGRAMS_H
#define NTOTHISTOGRAMS_H

#include "nfisGlobals.h"
#include "TnfisParam.h"

#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <TTimeStamp.h>

#include "TGo4EventProcessor.h"

using namespace std;

class nfisHistograms : public TGo4EventProcessor {
public:
    nfisHistograms();
    virtual ~nfisHistograms();

private:

    vector<string> *ScalerChannelNames();

    //parameter
    TnfisParamGlobal  *pParGlobal;

    //raw step
    void DefineRawScaler(const char *step_name, const char *dir_name,
                        unsigned int NumHist = NumScaler);
    void DefineRawQDC(const char *step_name, const char *dir_name,
                        unsigned int NumHist = NumQDC);
    void DefineRawTDC(const char *step_name, const char *dir_name,
                        unsigned int NumHist = NumTDC);

    //analysis step
    void DefineAnaCommon(const char *step_name, const char *dir_name);
    void DefineAnaTrigAcc(const char *step_name, const char *dir_name);
    void DefineAnaHZDR(const char *step_name, const char *dir_name,
                        unsigned int NumHist = NumHZDRFC);
    void DefineCalHZDR(const char *step_name, const char *dir_name,
                        unsigned int NumHist = NumHZDRFC);

    //calibration step

    //physics step

public:

    void DefineHistograms(const char *step_name, const char *det = "");

    //declaration of used histograms

    //raw step
    TH1I    *pH1RawVeto[NumVeto], *pH1RawVetoTot,
            *pH1RawTRGTime;
    TH1D    *pH1RawRate[NumScaler];
    TH1I    *pH1RawTDC[NumTDC];
    TH1I    *pH1RawQDCh[NumQDC], *pH1RawQDCl[NumQDC];

    TH2I    *pH2RawVeto, *pH2RawTDC,
    *pH2RawQDCl, *pH2RawQDCh, *pH2QDCTime[NumHZDRFC];
    TH2D    *pH2RawRate;

    //analysis step
    TH1I    *pH1AnaHitHZDR[NumHZDRFC],                              //Hit
            *pH1AnaHitAcc, *pH1AnaHitTrig;
    TH1I    *pH1AnaDtHZDR[NumHZDRFC],                               //TimeDiff
            *pH1AnaDtHZDR_g[NumHZDRFC];                             //PH-gated

    TH1I    *pH1AnaQDCl[NumQDC], *pH1AnaQDCl_SF[NumQDC], *pH1AnaQDCl_NIF[NumQDC], *pH1AnaQDCl_trig[NumQDC];

    TH2I    *pH2AnaHit, *pH2AnaDt, *pH2AnaDt_g;

    /*/calibration step
    TH1I    *pH1CalDtHZDR[NumHZDRFC],  *pH1CalDtHZDRSum;
    TH1I    *pH1CalDtHZDR_g[NumHZDRFC], *pH1CalDtHZDRSum_g;
    TH1I    *pH1CalQDCh[NumQDC], *pH1CalQDCl[NumQDC];
    TH1I    *pH1CalQDClSum, *pH1CalQDChSum;

    TH1I    *pH1CalEkin[NumHZDRFC], *pH1CalEkinSum;

    TH2I    *pH2CalDt, *pH2CalDt_g;
    TH2I    *pH2CalQDCl, *pH2CalQDCh;

    TH2I    *pH2CalToFvsQDC;
//*/
    //physics step

    // ClassDef makes the class able to return its own class type.
    ClassDef(nfisHistograms,1)
};

#endif // NTOTHISTOGRAMS_H
