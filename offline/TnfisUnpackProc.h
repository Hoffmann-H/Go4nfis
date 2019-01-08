#ifndef TUNPACKPROCESSOR_H
#define TUNPACKPROCESSOR_H

#include "TGo4EventProcessor.h"
#include "TGo4MbsEvent.h"
#include "TGo4Log.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TTimeStamp.h>
#include <string>

#include "nfisHistograms.h"
#include "TnfisUnpackEvent.h"
#include "TnfisParam.h"
#include "TGo4Picture.h"
#include "TnfisAnalysis.h"
#include "TLine.h"
#include "TText.h"

#include "TRandom2.h"

class TnfisUnpackEvent;
class nfisHistograms;
class TnfisTrigAccClk;


class TnfisMakeUnp : public TGo4EventProcessor {
   public:
      TnfisMakeUnp() ;
      TnfisMakeUnp(const char* name);
      virtual ~TnfisMakeUnp() ;

      // event processing function, default name
      Bool_t BuildEvent(TGo4EventElement * target);
      Bool_t FillHistograms(TnfisUnpackEvent *event);

      nfisHistograms    *GetnfisHistograms()    {return pHist;}
      TGo4MbsEvent      *GetRaw()               {return pRaw;}
      string            ExtractRun (const string& str);

 private:
      nfisHistograms    *pHist;
      TGo4MbsEvent      *pRaw;

      //parameter
      TnfisParamGlobal  *pParGlobal;

      //pictures
      TGo4Picture       *pPicQDCvsTime[NumHZDRFC], *pPicQDCvsTimeAll,
                        *pPicQDCl, *pPicQDCh, *pPicRate;

      Bool_t            CommentFlag[10];

      Int_t             PrintScalerFlag;
      Int_t             channel;      // current channel
      Int_t             channelok;    // used for identify known TDC channels while readout;
      ULong_t           dat[7];       // data read from one data word
      Long_t            l_scaler, GlobalScaler, count;
      Long_t            triggertimediff;
      Long_t            lasttriggertime;      // value of last trigger time tag
      Long_t            lastAccTrigger;

      char              date[MAXSTRINGLENGTH];

      Double_t          LNE_time;

      TTimeStamp        *t0;

      vector<string>    *pChannelNames;

      string            LastFile, CurrentRun;

      Double_t          t_live, t_real, t_dead;

      Double_t          randtriang();

      TnfisTrigAccClk   *pTrig, *pAcc, *pClock, *pVeto, *pNotVeto;
      TnfisVetoLength   *pVetoLength;
      TnfisTimer        *pTimer;
      TnfisEvtCounter   *pEvtCounter;
      TnfisAbsorber     *pAbs;
      TnfisUnpFC        *pHZDRFC;
      TnfisUnpPreAmp    *pHZDRPreAmp[NumHZDRFC];

      TTimeStamp        *AbsTime;

      void              AnalyzeMeasurementTimes(nELBE_data *u_data);
      void              AnalyzeTimeFlag(nELBE_data *u_data);
      void              AnalyzeScaler(nELBE_data *u_data, Int_t **pdata,
                                   Int_t *lwords,
                                   vector<string> *ChannelNames, TnfisTimer *Timer);
      void              AnalyzeVetoScaler(nELBE_data *u_data, Int_t **pdata,
                                   Int_t *lwords,
                                   TnfisUnpackEvent *pTarget,
                                   TnfisEvtCounter *EvtCounter);
      void              AnalyzeTriggerTimeTag(nELBE_data *u_data);
      void              AnalyzeAbsorberFlag(nELBE_data *u_data);

      void              AnalyzeADC(nELBE_data *u_data, TnfisUnpackEvent *pTarget,
                                   TnfisEvtCounter *EvtCounter);
      void              AnalyzeADCHeader(nELBE_data *u_data,
                                         TnfisEvtCounter *EvtCounter);
      void              AnalyzeADCDataWord(nELBE_data *u_data,
                                           TnfisUnpackEvent *pTarget,
                                           TnfisEvtCounter *EvtCounter);


      void              AnalyzeQDC(nELBE_data *u_data, TnfisUnpackEvent *pTarget,
                                   TnfisEvtCounter *EvtCounter);
      void              AnalyzeQDCDataWord(nELBE_data *u_data,
                                           TnfisUnpackEvent *pTarget,
                                           TnfisEvtCounter *EvtCounter);
      void              AnalyzeQDCRange(ULong_t *data, TnfisUnpackEvent *pTarget,
                                        TnfisEvtCounter *EvtCounter);
      void              AnalyzeQDCHeader(nELBE_data *u_data,
                                         TnfisEvtCounter *EvtCounter);

      void              AnalyzeTDC(nELBE_data *u_data, TnfisUnpackEvent *pTarget,
                                   TnfisEvtCounter *EvtCounter);
      void              AnalyzeTDCGlobalTrailer(TnfisEvtCounter *EvtCounter);
      void              AnalyzeTDCTrigAccClkVeto(ULong_t *data,
                                                 TnfisEvtCounter *EvtCounter,
                                                 TnfisTrigAccClk *pTrigAccClk);
      void              AnalyzeAbsoluteTime(nELBE_data *u_data,
                                            TnfisEvtCounter *EvtCounter);
      void              InitEventStructure(TnfisUnpackEvent *pTarget);
      void              ResetEventStructure(Option_t *t="");

      void              DrawLine(TGo4Picture *pPic, TH2I *pHist, Double_t time,
                                 const char *kind, UInt_t value=1);

      void              MakePictures();

      vector<string>    *ScalerChannelNames();


   ClassDef(TnfisMakeUnp,1)
};




#endif //TUNPACKPROCESSOR_H


//----------------------------END OF GO4 SOURCE FILE ---------------------
