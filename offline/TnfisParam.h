#ifndef SPAR_H
#define SPAR_H


#include "nfisGlobals.h"

#include "TGo4Parameter.h"
#include "TTimeStamp.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//                            Global Parameters                               //
////////////////////////////////////////////////////////////////////////////////

class TnfisParamGlobal : public TGo4Parameter {
   public:
      TnfisParamGlobal();
      TnfisParamGlobal(const char* name);
      virtual ~TnfisParamGlobal();
      virtual Int_t  PrintParameter(Text_t * n, Int_t);
      Int_t  PrintParameter();
      virtual Bool_t UpdateFrom(TGo4Parameter *);

      TTimeStamp *ExpStart, *ExpEnd;

      Bool_t    fFill;          // enable filling histograms
      Bool_t    fOutput;        // enable filling Tnfis event

      UInt_t    s_year;         //year of exp. start
      UInt_t    s_month;        //month of exp. start
      UInt_t    s_day;          //day of exp. start
      UInt_t    s_hour;         //hour of exp. start
      UInt_t    s_min;          //minute of exp. start
      UInt_t    s_sec;          //second of exp. start

      UInt_t    e_year;         //year of exp. end
      UInt_t    e_month;        //month of exp. end
      UInt_t    e_day;          //day of exp. end
      UInt_t    e_hour;         //hour of exp. end
      UInt_t    e_min;          //minute of exp. end
      UInt_t    e_sec;          //second of exp. end

      Float_t   ReadoutPeriod;  //Scaler readout period

      Bool_t    fShadowBar;     //analysis type. True: Underground measurement with shadow bar.
                                //               False: Neutron induced fission measurement.

   ClassDef(TnfisParamGlobal,1)
};


////////////////////////////////////////////////////////////////////////////////
//                            QDC Parameter                                   //
////////////////////////////////////////////////////////////////////////////////

class TnfisParamQDC : public TGo4Parameter {
   public:
      TnfisParamQDC();
      TnfisParamQDC(const char* name);
      virtual ~TnfisParamQDC();
      virtual Int_t  PrintParameter(Text_t * n, Int_t);
      Int_t  PrintParameter();
      virtual Bool_t UpdateFrom(TGo4Parameter *);

      Float_t  fOffsetH[NumHZDRFC]; // QDC high gain offset
      Float_t  fFactorH[NumHZDRFC]; // QDC high gain factor
      Float_t  fOffsetL[NumHZDRFC]; // QDC low gain factor
      Float_t  fFactorL[NumHZDRFC]; // QDC low gain factor

   ClassDef(TnfisParamQDC,1)
};

////////////////////////////////////////////////////////////////////////////////
//                            ToF Parameter                                   //
////////////////////////////////////////////////////////////////////////////////

class TnfisParamToF : public TGo4Parameter {
   public:
      TnfisParamToF();
      TnfisParamToF(const char* name);
      virtual ~TnfisParamToF();
      virtual Int_t  PrintParameter(Text_t * n, Int_t);
      Int_t  PrintParameter();
      virtual Bool_t UpdateFrom(TGo4Parameter *);

      Float_t   TDC_Dispersion;             //TDC Dispersion Time [ch/ps]

      //HZDR UFC
      Float_t   UFC_DepDist;                // Distance between U(235) deposits in the FC [cm]
      Float_t   UFC_FlightPath;             // neutron flight path to fission layer 1 of HZDR U(235) FC [cm]
      Int_t     UFC_Offset[NumHZDRFC];      // gamma-flash offset HZDR U(235) FC [ch]

   ClassDef(TnfisParamToF,1)
};


#endif //SPAR_H

//----------------------------END OF GO4 SOURCE FILE ---------------------
//*/
