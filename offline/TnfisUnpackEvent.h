#ifndef TEVENTUNP_H
#define TEVENTUNP_H

#include "nfisGlobals.h"

#include "TGo4EventElement.h"
#include "TGo4CompositeEvent.h"
#include "TTimeStamp.h"

class TnfisUnpackEvent;
class TnfisAbsorber;
class TnfisTrigAccClk;
class TnfisVetoLength;
class TnfisPMT;
class TnfisEvtCounter;
class TnfisCounter;
class TnfisTimer;
class TnfisHZDRFC;
class TnfisUnpPreAmp;

using namespace std;

class TnfisUnpackEvent : public TGo4CompositeEvent {
   public:
      TnfisUnpackEvent();
      TnfisUnpackEvent(TnfisUnpackEvent*);
      TnfisUnpackEvent(const char* name);
      virtual ~TnfisUnpackEvent();

      virtual void Clear(Option_t *t="");

   ClassDef(TnfisUnpackEvent,1)
};


//______________________________________________________________________________

// HZDR Fission Chamber
class TnfisUnpFC: public TGo4CompositeEvent {
    public:
    TnfisUnpFC();
    TnfisUnpFC(const char* name, Short_t id);
    virtual ~TnfisUnpFC();

    virtual void Clear(Option_t *t="");

    TnfisUnpPreAmp    *GetPreAmp(UInt_t index)     {return PreAmp[index];}

    protected:
    TnfisUnpPreAmp    *PreAmp[NumHZDRFC];          //!

    ClassDef(TnfisUnpFC,1)
};



//______________________________________________________________________________

// pre-amplifier
class TnfisUnpPreAmp: public TGo4EventElement {
    public:
    TnfisUnpPreAmp();
    TnfisUnpPreAmp(const char* name, Short_t id);
    virtual ~TnfisUnpPreAmp();

    virtual void Clear(Option_t *t="");

    void SetHit(UInt_t index1, UInt_t index2,
                Long_t value)                       {Hit[index1][index2]= value;}
    void SetHitCounter(UInt_t index, Long_t value)  {HitCounter[index]  = value;}
    void SetFirstHit(UInt_t index, Long_t value)    {FirstHit[index]    = value;}
    void SetTime(UInt_t index, Long_t value)        {Time[index]        = value;}
    void SetQDCl(UInt_t index, Long_t value)        {QDCl[index]        = value;}
    void SetQDCh(UInt_t index, Long_t value)        {QDCh[index]        = value;}

    Long_t GetHit(UInt_t index1, UInt_t index2)     {return Hit[index1][index2];}
    Long_t GetHitCounter(UInt_t index)              {return HitCounter[index];}
    Long_t GetFirstHit(UInt_t index)                {return FirstHit[index];}
    Long_t GetTime(UInt_t index)                    {return Time[index];}
    Long_t GetQDCl(UInt_t index)                    {return QDCl[index];}
    Long_t GetQDCh(UInt_t index)                    {return QDCh[index];}


    void IncreaseHitCounter(UInt_t index)           {HitCounter[index]++;}
    void IncreaseHit(UInt_t index1, UInt_t index2)  {Hit[index1][index2]++;}


    protected:
    Long_t HitCounter[MAXEVENTS];
    Long_t Hit[MAXEVENTS][MAXHITS];
    Long_t FirstHit[MAXEVENTS];
    Long_t Time[MAXEVENTS];
    Long_t QDCl[MAXEVENTS];
    Long_t QDCh[MAXEVENTS];

    ClassDef(TnfisUnpPreAmp,1)
};

//______________________________________________________________________________

//absorber
class TnfisAbsorber: public TGo4EventElement{
    public:
    TnfisAbsorber();
    TnfisAbsorber(const char* name, Short_t id);
    TnfisAbsorber(TnfisAbsorber *copy);
    virtual ~TnfisAbsorber();

    virtual void Clear(Option_t *t="");

    void    SetAbsorber(UInt_t value)   {Absorber = value;}
    UInt_t  GetAbsorber()               {return Absorber;}

    private:
    UInt_t  Absorber;                   //Absorber which was currently in beam

    ClassDef(TnfisAbsorber,1)

};

//______________________________________________________________________________

//Veto length
class TnfisVetoLength: public TGo4EventElement{
    public:
    TnfisVetoLength();
    TnfisVetoLength(const char* name, Short_t id);
    TnfisVetoLength(TnfisVetoLength *copy);
    virtual ~TnfisVetoLength();

    virtual void Clear(Option_t *t="");

    void        SetVetoLength(UInt_t index,Double_t value)
                                            {VetoLength[index] = value;}
    Double_t    GetVetoLength(UInt_t index) {return VetoLength[index];}

    private:
    Double_t          VetoLength[MAXEVENTS];

    ClassDef(TnfisVetoLength,1)

};

//______________________________________________________________________________

// trigger, acc and veto
class TnfisTrigAccClk: public TGo4EventElement{
    public:
    TnfisTrigAccClk();
    TnfisTrigAccClk(const char* name, Short_t id);
    virtual ~TnfisTrigAccClk();

    virtual void Clear(Option_t *t="");

    void IncreaseHit(UInt_t index1, UInt_t index2)  {Hit[index1][index2]++;}
    void IncreaseHitCounter(UInt_t index)           {HitCounter[index]++;}

    void SetHit(UInt_t index1, UInt_t index2,
                Long_t value)                       {Hit[index1][index2]= value;}
    void SetHitCounter(UInt_t index, Long_t value)  {HitCounter[index]  = value;}
    void SetFirstHit(UInt_t index, Long_t value)    {FirstHit[index]    = value;}


    Long_t GetHit(UInt_t index1, UInt_t index2)     {return Hit[index1][index2];}
    Long_t GetHitCounter(UInt_t index)              {return HitCounter[index];}
    Long_t GetFirstHit(UInt_t index)                {return FirstHit[index];}

    protected:

    Long_t Hit[MAXEVENTS][MAXHITS];
    Long_t HitCounter[MAXEVENTS];
    Long_t FirstHit[MAXEVENTS];


    ClassDef(TnfisTrigAccClk,1)
};




// event counter
class TnfisEvtCounter: public TGo4EventElement
//class TnfisEvtCounter
{
    public:
    TnfisEvtCounter();
    TnfisEvtCounter(const char* name, Short_t id);
    virtual ~TnfisEvtCounter();

    virtual void Clear(Option_t *t="");

    void IncreaseTDC()                      {TDC++;}
    void IncreaseQDC()                      {QDC++;}
    void IncreaseVeto()                     {Veto++;}
    void IncreaseTotal()                    {Total++;}

    void SetTDC(UInt_t value)               {TDC    = value;}
    void SetQDC(UInt_t value)               {QDC    = value;}
    void SetVeto(UInt_t value)              {Veto   = value;}
    void SetTotal(UInt_t value)             {Total  = value;}
    void SetAbsTime(TTimeStamp value)       {AbsTime= value;}

    UInt_t GetTDC()                         {return TDC;}
    UInt_t GetQDC()                         {return QDC;}
    UInt_t GetVeto()                        {return Veto;}
    UInt_t GetTotal()                       {return Total;}

    TTimeStamp GetAbsTime()                 {return AbsTime;}


    protected:
    UInt_t TDC;
    UInt_t QDC;
    UInt_t Veto;
    UInt_t Total;
    TTimeStamp AbsTime;

    ClassDef(TnfisEvtCounter,1)
};



//TotalCounter
class TnfisCounter: public TGo4EventElement{
    public:
    TnfisCounter();
    TnfisCounter(const char* name, Short_t id);
    virtual ~TnfisCounter();

    virtual void Clear(Option_t *t="");

    void IncreaseTDC()                      {TDC++;}
    void IncreaseQDC()                      {QDC++;}
    void IncreaseVeto()                     {Veto++;}
    void IncreaseTotal()                    {Total++;}

    void SetTDC(UInt_t value)               {TDC    = value;}
    void SetQDC(UInt_t value)               {QDC    = value;}
    void SetVeto(UInt_t value)              {Veto   = value;}
    void SetTotal(UInt_t value)             {Total  = value;}

    UInt_t GetTDC()                         {return TDC;}
    UInt_t GetQDC()                         {return QDC;}
    UInt_t GetVeto()                        {return Veto;}
    UInt_t GetTotal()                       {return Total;}

    protected:
    UInt_t TDC;
    UInt_t QDC;
    UInt_t Veto;
    UInt_t Total;


    ClassDef(TnfisCounter,1)

};

//TotalTimer
class TnfisTimer: public TGo4EventElement{
    public:
    TnfisTimer();
    TnfisTimer(const char* name, Short_t id);
    virtual ~TnfisTimer();

    virtual void Clear(Option_t *t="");

    void IncreaseDeadTime(Double_t value)   {t_dead+= value;}
    void IncreaseLiveTime(Double_t value)   {t_live+= value;}
    void IncreaseRealTime(Double_t value)   {t_real+= value;}

    void SetDeadTime(Double_t value)        {t_dead = value;}
    void SetLiveTime(Double_t value)        {t_live = value;}
    void SetRealTime(Double_t value)        {t_real = value;}

    Double_t GetDeadTime()                  {return t_dead;}
    Double_t GetLiveTime()                  {return t_live;}
    Double_t GetRealTime()                  {return t_real;}

    protected:
    Double_t t_dead, t_live, t_real;

    ClassDef(TnfisTimer,1)

};



typedef union {
  unsigned long value;

    //----------- common ------------------
    struct {
        unsigned long data :27;
        unsigned long geo  : 5;
    } common;

    //----------- Time structure -----------
    struct {
        unsigned long data :26;
        unsigned long flag : 1;
        unsigned long geo  : 5;
    } time;

    struct {
        unsigned long data :27;
        unsigned long geo  : 5;
    } time_flag;

    //----------- Scaler--------------------
    struct {
      unsigned long words :27;
      unsigned long geo   : 5;
    } scaler;

    //----------- Veto scaler---------------
    struct {
      unsigned long words :27;
      unsigned long geo   : 5;
    } veto;

    //----------- ADC ----------------------
    struct {
        unsigned long      : 15;
        unsigned long flag :  1;
        unsigned long      : 11;
        unsigned long geo  :  5;
    } adc;
    struct {
        unsigned long cnt  : 12;
        unsigned long mult :  3;
        unsigned long flag :  1;
        unsigned long      : 11;
        unsigned long geo  :  5;
    } adc_header;
    struct {
        unsigned long data : 12;
        unsigned long ch   :  3;
        unsigned long flag :  1;
        unsigned long      : 11;
        unsigned long geo  :  5;
    } adc_data;

  //----------- QDC structure -----------
    struct {
      unsigned long      :24;
      unsigned long flag : 3;
      unsigned long geo  : 5;
    } qdc;
    struct {
      unsigned long      : 8;
      unsigned long cnt  : 6;
      unsigned long      : 2;
      unsigned long crt  : 8;
      unsigned long flag : 3;
      unsigned long geo  : 5;
    } qdc_header;
    struct {
      unsigned long data :12;
      unsigned long ov   : 1;
      unsigned long un   : 1;
      unsigned long      : 3;
      unsigned long r    : 1;
      unsigned long ch   : 3;
      unsigned long      : 3;
      unsigned long flag : 3;
      unsigned long geo  : 5;
    } qdc_data;
    struct {
      unsigned long cnt  :24;
      unsigned long flag : 3;
      unsigned long geo  : 5;
    } qdc_trailer;

    //----------- TDC structure -----------
    struct {
      unsigned long      :26;
      unsigned long flag : 1;
      unsigned long geo  : 5;
    } tdc;
    struct {
      unsigned long data :21;
      unsigned long ch   : 5;
      unsigned long flag : 1;
      unsigned long geo  : 5;
    } tdc_data;
    struct {
      unsigned long cnt  :16;
      unsigned long      : 3;
      unsigned long stat : 3;
      unsigned long      : 4;
      unsigned long flag : 1;
      unsigned long geo  : 5;
    } tdc_trailer;
    struct {
      unsigned long data :27;
      unsigned long geo  : 5;
    } tdc_trg_time;

    //----------- OPC ----------------------
    struct {
      unsigned long words :27;
      unsigned long geo   : 5;
    } opc;

    //----------- Absorber Flag ------------
    struct {
      unsigned long abs :  3;
      unsigned long     : 24;
      unsigned long geo :  5;
    } absorber;

    //----------- Test data ----------------
    struct {
      unsigned long data :27;
      unsigned long geo  : 5;
    } test;

    //----------- Absolute Time ------------
    struct {
      unsigned long data :  16;
      unsigned long      :  9;
      unsigned long id   :  2;
      unsigned long geo  :  5;
    } abstime;

  } nELBE_data;


// event data
typedef struct {
  // trigger, acc and veto
  struct {
    long hit_cnt;
    long hit[MAXHITS];
    long firsthit;
  } trig,acc;

  struct {
    long hit_cnt;
    long hit[MAXHITS];
    long length;
  } veto, nveto;

  // HZDR fission chamber
  struct {
    long hit_cnt;
    long hit[MAXHITS];
    long firsthit;
    long time;
    long qdcl;
    long qdch;
  } fchzdr[NumHZDRFC];

  struct {
    struct {
      long hit_cnt;
      long hit[MAXHITS];
      long firsthit;
      long time;
    } pmt[2];
    long diff;
    long sum;
  } pl;

  double vetolength;

} event_struct;

// multplicity
typedef struct {
  struct {
    int cnt;
    int ch[NumHZDRFC];
  } fchzdr;
  struct {
    int cnt;
    int ch;
  } pl;
} multiplicity_struct;

// eventcounter
typedef struct {
  long tdc;
  long qdc;
  long adc;
  long veto;
  long total;
} eventcounter_struct;


#endif //TEVENTUNP_H




//----------------------------END OF GO4 SOURCE FILE ---------------------
