#include "TnfisUnpackEvent.h"
#include "TGo4EventElement.h"
#include "TGo4CompositeEvent.h"
#include "Riostream.h"

////////////////////////////////////////////////////////////////////////////////

TnfisAbsorber::TnfisAbsorber():TGo4EventElement()
{ //void constructor to enable libary in TTree
}

TnfisAbsorber::TnfisAbsorber(const char* name, Short_t id) :
    TGo4EventElement(name,name,id)
{
  cout << "**** TnfisAbsorber: Create instance " << name << endl;
  //init array with 0
  Clear();
}

//______________________________________________________________________________
TnfisAbsorber::~TnfisAbsorber()
{
    cout << "**** Tnfis:" << this->GetName()
         << " Delete TnfisAbsorber instance " << endl;
}

//______________________________________________________________________________
void  TnfisAbsorber::Clear(Option_t *t)
{
    // all members should be cleared.
    Absorber = 0;

}


////////////////////////////////////////////////////////////////////////////////

TnfisVetoLength::TnfisVetoLength():TGo4EventElement()
{ //void constructor to enable libary in TTree
}

TnfisVetoLength::TnfisVetoLength(const char* name, Short_t id) :
    TGo4EventElement(name,name,id)
{
  cout << "**** TnfisVetoLength: Create instance " << name << endl;
  //init array with 0
  Clear();
}

//______________________________________________________________________________
TnfisVetoLength::~TnfisVetoLength()
{
    cout << "**** TnfisVetoLength:" << this->GetName()
         << " Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisVetoLength::Clear(Option_t *t)
{
  // all members should be cleared.
    memset(VetoLength,0, sizeof(VetoLength));

}


////////////////////////////////////////////////////////////////////////////////

TnfisTrigAccClk::TnfisTrigAccClk():TGo4EventElement()
{ //void constructor to enable libary in TTree
}

TnfisTrigAccClk::TnfisTrigAccClk(const char* name, Short_t id) :
    TGo4EventElement(name,name,id)
{
  cout << "**** TnfisTrigAccVeto: Create instance " << name << endl;
  //init array with 0
  Clear();
}

//______________________________________________________________________________
TnfisTrigAccClk::~TnfisTrigAccClk()
{
    cout << "**** TnfisTrigAccVeto:" << this->GetName()
         << " Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisTrigAccClk::Clear(Option_t *t)
{
  // all members should be cleared.
    if (strcmp(t, "init")==0)
    {   memset(Hit,0, sizeof(Hit));
        memset(FirstHit,0, sizeof(FirstHit));
    }

    memset(HitCounter,0, sizeof(HitCounter));

}

////////////////////////////////////////////////////////////////////////////////

TnfisUnpFC::TnfisUnpFC():TGo4CompositeEvent()
{ //void constructor to enable libary in TTree
}

TnfisUnpFC::TnfisUnpFC(const char* name, Short_t id) :
    TGo4CompositeEvent(Form("%s_Common", name),name,id)
{
  cout << "**** TnfisUnpFC: Create instance " << name << endl;

  TString chname;

  //loop over all PreAmps in FC
  for(Int_t i = 0; i < NumHZDRFC; i++)
  {   chname.Form("%s_PreAmpCh%i",name, i);
      PreAmp[i] = new TnfisUnpPreAmp(chname.Data(), i);
      addEventElement(PreAmp[i]);
  }
}
//______________________________________________________________________________
TnfisUnpFC::~TnfisUnpFC()
{
  cout << "**** TnfisUnpFC: Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisUnpFC::Clear(Option_t *t)
{
  // all members should be cleared.

}

////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TnfisUnpPreAmp::TnfisUnpPreAmp() : TGo4EventElement()
{ //void constructor to enable libary in TTree

}

//______________________________________________________________________________
TnfisUnpPreAmp::TnfisUnpPreAmp(const char* name, Short_t id) :
    TGo4EventElement(name, name, id)
{
    cout << "\t**** TnfisHZDRPreAmp: Create instance " << name << endl;
    //init arrays with 0
    Clear();
}

//______________________________________________________________________________
TnfisUnpPreAmp::~TnfisUnpPreAmp()
{
    cout << "\t**** TnfisHZDRPreAmp:" << this->GetName()
         << " Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisUnpPreAmp::Clear(Option_t *t)
{
  // all members should be cleared.
    if (strcmp(t, "init")==0)
    {   memset(QDCh,0, sizeof(QDCh));
        memset(QDCl,0, sizeof(QDCl));
        memset(Hit,0, sizeof(Hit));
        memset(FirstHit,0, sizeof(FirstHit));
        memset(Time,0, sizeof(Time));
    }

    memset(HitCounter,0, sizeof(HitCounter));

}


////////////////////////////////////////////////////////////////////////////////

TnfisEvtCounter::TnfisEvtCounter():TGo4EventElement()
{ //void constructor to enable libary in TTree
    Clear();
}

//______________________________________________________________________________
TnfisEvtCounter::TnfisEvtCounter(const char* name, Short_t id) :
    TGo4EventElement(name, name, id)
{
    cout << "**** TnfisEvtCounter: Create instance " << name << endl;
    //init global paramters with 0
    Clear();
}
//______________________________________________________________________________
TnfisEvtCounter::~TnfisEvtCounter()
{
    cout << "**** TnfisEvtCounter:" << this->GetName()
         << " Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisEvtCounter::Clear(Option_t *t)
{
  // all members should be cleared.
  TDC = Veto = Total = 0;
  QDC = -1;

}

////////////////////////////////////////////////////////////////////////////////

TnfisCounter::TnfisCounter():TGo4EventElement()
{ //void constructor to enable libary in TTree
}

TnfisCounter::TnfisCounter(const char* name, Short_t id) :
    TGo4EventElement(name,name,id)
{
  cout << "**** TnfisCounter: Create instance " << name << endl;
  //init array with 0
  Clear();
}

//______________________________________________________________________________
TnfisCounter::~TnfisCounter()
{
    cout << "**** TnfisCounter:" << this->GetName()
         << " Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisCounter::Clear(Option_t *t)
{
    // all members should be cleared.
    TDC = Veto = Total = 0;
    QDC = -1;
}

////////////////////////////////////////////////////////////////////////////////

TnfisTimer::TnfisTimer():TGo4EventElement()
{ //void constructor to enable libary in TTree
}

TnfisTimer::TnfisTimer(const char* name, Short_t id) :
    TGo4EventElement(name,name,id)
{
  cout << "**** TnfisTimer: Create instance " << name << endl;
  //init array with 0
  Clear();
}

//______________________________________________________________________________
TnfisTimer::~TnfisTimer()
{
    cout << "**** TnfisTimer:" << this->GetName()
         << " Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisTimer::Clear(Option_t *t)
{
    // all members should be cleared.
    t_dead = t_live = t_real = 0.0;

}


////////////////////////////////////////////////////////////////////////////////

TnfisUnpackEvent::TnfisUnpackEvent():TGo4CompositeEvent()
{ //void constructor to enable libary in TTree
}

TnfisUnpackEvent::TnfisUnpackEvent(const char* name) : TGo4CompositeEvent()
{
  cout << "**** TnfisUnpackEvent: Create instance " << name << endl;

    addEventElement(new TnfisTrigAccClk("Trig",0));
    addEventElement(new TnfisTrigAccClk("Acc",1));
    addEventElement(new TnfisEvtCounter("EvtCounter",2));
//    addEventElement(new TnfisTrigAccClk("Clock",2));
    addEventElement(new TnfisTimer("Timer",3));
    addEventElement(new TnfisTrigAccClk("Veto",4));
//    addEventElement(new TnfisAbsorber("Absorber",4));
    addEventElement(new TnfisTrigAccClk("notVeto",5));
    addEventElement(new TnfisUnpFC("HZDR_FC",6));
//    addEventElement(new TnfisVetoLength("VetoLength",6));

}

TnfisUnpackEvent::~TnfisUnpackEvent()
{
  cout << "**** TnfisEventUnp: Delete instance " << endl;
}

void  TnfisUnpackEvent::Clear(Option_t *t)
{

}


//----------------------------END OF GO4 SOURCE FILE ---------------------
