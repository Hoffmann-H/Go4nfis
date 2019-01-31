#include "TnfisParam.h"
#include "Riostream.h"

////////////////////////////////////////////////////////////////////////////////
//                            Global Parameters                               //
////////////////////////////////////////////////////////////////////////////////

//******************************************************************************
TnfisParamGlobal::TnfisParamGlobal() : TGo4Parameter("Parameter")
{
    fFill=kTRUE;
    fOutput=kFALSE;
}

TnfisParamGlobal::TnfisParamGlobal(const char* name) : TGo4Parameter(name)
{
    fFill=kTRUE;
    fOutput=kFALSE;

    s_year    = 2014;   //year of exp. start
    s_month   = 05;     //month of exp. start
    s_day     = 22;     //day of exp. start
    s_hour    = 12;     //hour of exp. start
    s_min     = 0;      //minute of exp. start
    s_sec     = 0;      //second of exp. start

    ExpStart = new TTimeStamp(s_year, s_month, s_day,
                              s_hour, s_min, s_sec, 0, false);

//    e_year    = 2014;   //year of exp. end
//    e_month   = 05;     //month of exp. end
//    e_day     = 28;     //day of exp. end
//    e_hour    = 11;     //hour of exp. end
//    e_min     = 0;      //minute of exp. end
//    e_sec     = 0;      //second of exp. end

//    e_year    = 2014;   //year of exp. end
//    e_month   = 05;     //month of exp. end
//    e_day     = 28;     //day of exp. end
//    e_hour    = 22;     //hour of exp. end
//    e_min     = 0;      //minute of exp. end
//    e_sec     = 0;      //second of exp. end

//    s_year    = 2014;   //year of exp. start
//    s_month   = 05;     //month of exp. start
//    s_day     = 28;     //day of exp. start
//    s_hour    = 7;     //hour of exp. start
//    s_min     = 0;      //minute of exp. start
//    s_sec     = 0;      //second of exp. start

//    ExpStart = new TTimeStamp(s_year, s_month, s_day,
//                              s_hour, s_min, s_sec, 0, false);

    e_year    = 2014;   //year of exp. end
    e_month   = 06;     //month of exp. end
    e_day     = 02;     //day of exp. end
    e_hour    = 12;     //hour of exp. end
    e_min     = 0;      //minute of exp. end
    e_sec     = 0;      //second of exp. end


    ReadoutPeriod = 15.;//scaler readout period

    fShadowBar = kFALSE;

    ExpEnd = new TTimeStamp(e_year, e_month, e_day,
                            e_hour, e_min, e_sec, 0, false);
}

TnfisParamGlobal::~TnfisParamGlobal()
{
}
//******************************************************************************

Int_t TnfisParamGlobal::PrintParameter()
{
    cout << "Parameter " << GetName() << ":";
    cout << " fill="     << fFill;
    cout << " output="   << fOutput << endl;

    cout << " ReadoutPeriod="   << ReadoutPeriod;

    cout << " Experiment was started @: " << s_year << "/" << s_month << "/" << s_day
         << " @ UTC:" << s_hour << ":" << s_min << ":" << s_sec << endl;

    cout << " Experiment will end @: " << e_year << "/" << e_month << "/" << e_day
         << " @ UTC:" << e_hour << ":" << e_min << ":" << e_sec << endl;

    return 0;
}

//______________________________________________________________________________
Int_t TnfisParamGlobal::PrintParameter(Text_t * n, Int_t)
{
    PrintParameter();

    return 0;
}

//______________________________________________________________________________
Bool_t TnfisParamGlobal::UpdateFrom(TGo4Parameter *pp)
{
    if(pp->InheritsFrom("TnfisParamGlobal"))
    {   TnfisParamGlobal * from;
        from = (TnfisParamGlobal *) pp;

        cout << "**** TnfisParamGlobal " << GetName()
             << " updated from auto save file" << endl;

        fFill   = from->fFill;
        fOutput = from->fOutput;

        ReadoutPeriod = from->ReadoutPeriod;

        s_year    = from->s_year;   //year of exp. start
        s_month   = from->s_month;  //month of exp. start
        s_day     = from->s_day;    //day of exp. start
        s_hour    = from->s_hour;   //utc hour of exp. start
        s_min     = from->s_min;    //utc minute of exp. start
        s_sec     = from->s_sec;    //utc second of exp. start

        e_year    = from->e_year;   //year of exp. end
        e_month   = from->e_month;  //month of exp. end
        e_day     = from->e_day;    //day of exp. end
        e_hour    = from->e_hour;   //utc hour of exp. end
        e_min     = from->e_min;    //utc minute of exp. end
        e_sec     = from->e_sec;    //utc second of exp. end

        PrintParameter();
    }
    else
        cout << "Wrong parameter object: " << pp->ClassName() << endl;

    return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
//                            QDC Parameter                                   //
////////////////////////////////////////////////////////////////////////////////


//******************************************************************************
TnfisParamQDC::TnfisParamQDC() : TGo4Parameter("Parameter")
{

}

TnfisParamQDC::TnfisParamQDC(const char* name) : TGo4Parameter(name)
{
    //init paramters with std. values
    //low gain      //CH:1       2       3       4       5       6       7       8
    Float_t OffsetL[] = {-38.13, -18.86, -55.40, -50.49, -26.83, -31.96, -33.86, -15.03}; // shifting pedastal of ch 1 to 0
    Float_t FactorL[] = {1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.};

    //high gain
    Float_t OffsetH[] = {-1112., -1001.,-817.,-1129.,-1024., 0.,-1177., 0.};
    Float_t FactorH[] = {1., 1.02,  1.10,  1.14,  1.05, 0.,  1.08, 0.};

    for (int i_PreAmp=0; i_PreAmp<8; i_PreAmp++)
    {   //QDC low gain
        fOffsetL[i_PreAmp] = OffsetL[i_PreAmp];
        fFactorL[i_PreAmp] = FactorL[i_PreAmp];

        //QDC high gain
        fOffsetH[i_PreAmp] = OffsetH[i_PreAmp];
        fFactorH[i_PreAmp] = FactorH[i_PreAmp];
    }

}

TnfisParamQDC::~TnfisParamQDC()
{
}
//******************************************************************************

Int_t TnfisParamQDC::PrintParameter()
{
    cout << "Parameter " << GetName() << ":";


    return 0;
}

//______________________________________________________________________________
Int_t TnfisParamQDC::PrintParameter(Text_t * n, Int_t)
{
    PrintParameter();

    return 0;
}

//______________________________________________________________________________
Bool_t TnfisParamQDC::UpdateFrom(TGo4Parameter *pp)
{
    if(pp->InheritsFrom("TnfisParamQDC"))
    {   TnfisParamQDC * from;
        from = (TnfisParamQDC *) pp;

        cout << "**** TnfisParamQDC " << GetName()
             << " updated from auto save file" << endl;

        for (int i_PreAmp=0; i_PreAmp<8; i_PreAmp++)
        {   //QDC high gain
            fOffsetH[i_PreAmp] = from->fOffsetH[i_PreAmp];
            fFactorH[i_PreAmp] = from->fFactorH[i_PreAmp];

            //QDC low gain
            fOffsetL[i_PreAmp] = from->fOffsetL[i_PreAmp];
            fFactorL[i_PreAmp] = from->fFactorL[i_PreAmp];
        }

        PrintParameter();
    }
    else
        cout << "Wrong parameter object: " << pp->ClassName() << endl;

    return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
//                            TDC Parameter                                   //
////////////////////////////////////////////////////////////////////////////////


//******************************************************************************
TnfisParamToF::TnfisParamToF() : TGo4Parameter("Parameter")
{

}

TnfisParamToF::TnfisParamToF(const char* name) : TGo4Parameter(name)
{
    //init paramters with std. values
    //n flight Path
    UFC_FlightPath  =   624.2;  //[cm]
    //distance between uranium deposits in the fission chamber
    UFC_DepDist     =   1.0;    //[cm]
    //TDC dispersion
    TDC_Dispersion  =   25.0;   //[ps/ch]

    //array containing the gamma-peak positions [ch] of the time-diff spectra
                        //CH:   1       2       3       4
    Int_t UFC_gPeak[] = {       14583,  14070,  13916,  13950,
                        //      5       6       7       8
                                13971,  13979,  14120,  14014};

    for (int i_preAmp=0; i_preAmp < NumHZDRFC; i_preAmp++ )
        UFC_Offset[i_preAmp] = UFC_gPeak[i_preAmp];


}

TnfisParamToF::~TnfisParamToF()
{
}

//******************************************************************************

Int_t TnfisParamToF::PrintParameter()
{
    cout << "Parameter " << GetName() << ":" << endl;

    cout << "TDC dispersion: " << TDC_Dispersion <<  " ps/ch" << endl;

    cout << "Distance of deposits in HZDR UFC: " << UFC_DepDist << " cm" << endl;
    cout << "Neutron Flight Path HZDR FC ch1:"<< UFC_FlightPath << " cm" << endl;
    for (int i_PreAmp=0; i_PreAmp<8; i_PreAmp++)
        cout << "Offset gamma-flash HZDR UFC PreAmp " << i_PreAmp << ": "
             << UFC_Offset[i_PreAmp] << endl;

    return 0;
}

//______________________________________________________________________________
Int_t TnfisParamToF::PrintParameter(Text_t * n, Int_t)
{
  PrintParameter();

  return 0;
}

//______________________________________________________________________________
Bool_t TnfisParamToF::UpdateFrom(TGo4Parameter *pp)
{
    if(pp->InheritsFrom("TnfisParamToF"))
    {   TnfisParamToF * from;
        from = (TnfisParamToF *) pp;

        cout << "**** TnfisParamToF " << GetName()
             << " updated from auto save file" << endl;

        TDC_Dispersion  = from->TDC_Dispersion;
        UFC_DepDist     = from->UFC_DepDist;
        UFC_FlightPath  = from->UFC_FlightPath;

        for (int i_PreAmp=0; i_PreAmp<NumHZDRFC; i_PreAmp++)
            UFC_Offset[i_PreAmp]  = from->UFC_Offset[i_PreAmp];

        PrintParameter();
    }
    else
        cout << "Wrong parameter object: " << pp->ClassName() << endl;

    return kTRUE;
}

//----------------------------END OF GO4 SOURCE FILE ---------------------
//*/
