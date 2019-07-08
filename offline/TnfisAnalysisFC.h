#ifndef TNFISANALYSISFC_H
#define TNFISANALYSISFC_H

#include "TGo4EventProcessor.h"
#include "TGo4WinCond.h"

#include "nfisHistograms.h"
#include "TnfisParam.h"
#include "TnfisUnpackProc.h"
#include "TTimeStamp.h"

#include <iostream>


class nfisHistograms;
class TnfisAnaFCEvent;
class TnfisAnaFC;
class TnfisAnaPreAmp;


using namespace std;

class TnfisFCAnalysis : public TGo4EventProcessor {
   public:
      TnfisFCAnalysis() ;
      TnfisFCAnalysis(const char* name);
      virtual ~TnfisFCAnalysis() ;

      // event processing function, default name
      Bool_t BuildEvent(TGo4EventElement * target);
      Bool_t FillHistograms();

 private:

      Bool_t CommentFlag;

      void InitEventStructure(TnfisAnaFCEvent *pTarget);
      void ResetEventStructure(Option_t *t);

      //histograms
      nfisHistograms    *pHistFC;

      //conditions
      TGo4WinCond       *pConQDC[NumHZDRFC];
      TGo4WinCond       *pConToF[NumHZDRFC];
      TGo4WinCond       *pConDt[NumHZDRFC];

      //pictures
      TGo4Picture       *pPicToF, *pPicToF_g;

      TnfisAnaFC        *pHZDRFCOut;
      TnfisAnaPreAmp    *pHZDRPreAmpOut[NumHZDRFC];

      void              MakeConditions();
      void              MakePictures();


   ClassDef(TnfisFCAnalysis,1)
};


////////////////////////////////////////////////////////////////////////////////
///                         Fission Chamber Event                             //
////////////////////////////////////////////////////////////////////////////////

#include "TGo4EventElement.h"
#include "TnfisUnpackEvent.h"

class TnfisAnaFCEvent: public TGo4CompositeEvent{
    public:
    TnfisAnaFCEvent();
    TnfisAnaFCEvent(TnfisAnaFCEvent*);
    TnfisAnaFCEvent(const char* name);
    virtual ~TnfisAnaFCEvent();

    void                SetNumTrigEvents(UInt_t value)  {NumTrigEvents = value;}
    UInt_t              GetNumTrigEvents()              {return NumTrigEvents;}

    virtual void Clear(Option_t *t="");

    protected:
    Long_t              NumTrigEvents;

   ClassDef(TnfisAnaFCEvent,1)
};

//______________________________________________________________________________

// Fission Chamber
class TnfisAnaFC: public TGo4CompositeEvent {
    public:
    TnfisAnaFC();
    TnfisAnaFC(const char* name, Short_t id);
    virtual ~TnfisAnaFC();

    virtual void Clear(Option_t *t="");

    TnfisAnaPreAmp      *GetPreAmp(UInt_t index)    {return PreAmp[index];}


    protected:
    TnfisAnaPreAmp      *PreAmp[NumHZDRFC];          //!

    ClassDef(TnfisAnaFC,1)
};

//______________________________________________________________________________

// pre-amplifier
class TnfisAnaPreAmp: public TGo4EventElement {
    public:
    TnfisAnaPreAmp();
    TnfisAnaPreAmp(const char* name, Short_t id);
    virtual ~TnfisAnaPreAmp();

    virtual void Clear(Option_t *t="");

    void SetQDCl(Long_t value)          {pVecQDCl->push_back(value);}
    void SetQDCh(Long_t value)          {pVecQDCh->push_back(value);}
    void SetTime(Long_t value)          {pVecTime->push_back(value);}

    vector<Long_t> *GetQDCl()           {return pVecQDCl;}
    vector<Long_t> *GetQDCh()           {return pVecQDCh;}
    vector<Long_t> *GetTime()           {return pVecTime;}

    Long_t GetQDCl(UInt_t index)        {return pVecQDCl->at(index);}
    Long_t GetQDCh(UInt_t index)        {return pVecQDCh->at(index);}
    Long_t GetTime(UInt_t index)        {return pVecTime->at(index);}

    protected:
    vector<Long_t> *pVecQDCl;
    vector<Long_t> *pVecQDCh;
    vector<Long_t> *pVecTime;


    ClassDef(TnfisAnaPreAmp,1)
};

//----------------------------END OF GO4 SOURCE FILE ---------------------


#endif // TNFISANALYSISFC_H
//*/
