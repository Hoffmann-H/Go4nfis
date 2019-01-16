//----------------------start of nfisHistograms.cxx-----------------------------

#include "nfisHistograms.h"


////////////////////////////////////////////////////////////////////////////////

nfisHistograms::nfisHistograms()
{   //init parameter
    pParGlobal= (TnfisParamGlobal *) MakeParameter("GlobalPar", "TnfisParamGlobal");
}

nfisHistograms::~nfisHistograms()
{
}

////////////////////////////////////////////////////////////////////////////////

//define all nfis histogrmas
void nfisHistograms::DefineHistograms(const char* step_name, const char* det)
{
    //define histograms for the unpacking step
    if (strcmp(step_name, "Raw")==0)
    {   DefineRawScaler(step_name,      "Scaler");
        DefineRawTDC(step_name,         "TDC");
        DefineRawQDC(step_name,         "QDC");
    }

    //define histograms for the fission chamber analysis step
    if (strcmp(step_name, "Analysis")==0)
    {   DefineAnaCommon(step_name,  "Common");
        if (strcmp(det, "FC")==0)
        {   DefineAnaTrigAcc(step_name, "Acc&Trig");
            DefineAnaHZDR(step_name,    "PuFC");
        }
    }

    /*//define histograms for the calibration step
    if (strcmp(step_name, "Calibration")==0)
    {   DefineCalHZDR(step_name,        "nELBE_UFC");
    }

    //define histograms for the physics step
    if (strcmp(step_name, "Physics")==0)
    {
    }
    //*/

}


////////////////////////////////////////////////////////////////////////////////
////                    Histograms for the raw step                         ////
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________

//common histograms
void  nfisHistograms::DefineRawScaler(const char *step_name,
                                      const char *dir_name,
                                      unsigned int NumHist)
{//method to define common histograms for the unpacking step
char    path_name[256] = "",
        obj_name[256]  = "",
        obj_title[256] = "";
//load a vector containing the descriptions of each scaler channel
    vector<string> *ChannelNames = ScalerChannelNames();

//define directory structure of histograms
if (strcmp(dir_name, "")==0)
    sprintf(path_name,"%s", step_name);
else
    sprintf(path_name,"%s/%s", step_name, dir_name);


//determine start and stop time of experiment
UInt_t t_start = pParGlobal->ExpStart->GetSec();
UInt_t t_end   = pParGlobal->ExpEnd->GetSec();
//calculate the beam time in seconds
UInt_t tdiff = t_end-t_start;
//calculate number of bins
UInt_t tbins = tdiff/pParGlobal->ReadoutPeriod;

cout << "tdiff=" << tdiff << endl;
cout << "tbins=" << tbins << endl;

for (int i_ch=0; i_ch<NumHist; i_ch++)
{   ////Rate
    sprintf(obj_name,"%s/Rates/H1%sRate_%i", path_name, step_name, i_ch+1);
    sprintf(obj_title,"Scaler rate, %s", ChannelNames->at(i_ch).c_str());

    pH1RawRate[i_ch]   = (TH1D *) MakeTH1('D',obj_name, obj_title,
//                                          32500,0.,32500.,
                                            tbins,t_start, t_end,
                                            "", "scaler rate / 1/s");
    ////Veto
    if (i_ch<NumVeto){
        sprintf(obj_name,"%s/Veto/H1%sVeto_%i", path_name, step_name, i_ch+1);
        sprintf(obj_title,"Veto length, channel %i", i_ch+1);

        pH1RawVeto[i_ch]   = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                              262144,0.,6.55360e6,
                                              "#font[12]{t}_{veto} / ns",
                                              "counts");
    }

    ////Veto Length
    if (i_ch==0){
        sprintf(obj_name,"%s/Veto/H1%sVetoTot", path_name, step_name);
        sprintf(obj_title,"Veto length, channel 4");

        pH1RawVetoTot   = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            6.55360e6,0.,6.55360e6,
                                            "#font[12]{t}_{veto} / ns",
                                           "veto channel");
    }
}

    ////Trigger time
    sprintf(obj_name,"%s/H1%sTRGTime", path_name, step_name);
    pH1RawTRGTime      = (TH1I *) MakeTH1('I',obj_name, "Trigger time tag",
                                        50000,0,50000,"Channels", "counts");

    ////Rates 2D
    sprintf(obj_name,"%s/Rates/H2%sRate", path_name, step_name);
    pH2RawRate         = (TH2D *) MakeTH2('D',obj_name,"Scaler rate",
                                        32500,0.,32500,
                                        2*NumScaler,0.,NumScaler,
                                        "#font[12]{t} / readout intervall",
                                        "scaler channel");
    ////Veto 2D
    sprintf(obj_name,"%s/Veto/H2%sVeto", path_name, step_name);
    pH2RawVeto         = (TH2I *) MakeTH2('I',obj_name,"Veto length",
                                        262144,0.,6.55360e6,
                                        10,0,5,
                                       "#font[12]{t}_{veto} / ns",
                                       "Channels");
//*/
}

//______________________________________________________________________________

//TDC histograms
void nfisHistograms::DefineRawTDC(const char* step_name, const char* dir_name,
                                  unsigned int NumHist)
{//method to define TDC histograms for the unpacking step
    char    path_name[256] = "",
            obj_name[256] = "",
            obj_title[256] = "";

    //define directory structure of histograms
    if (strcmp(dir_name, "")==0)
        sprintf(path_name,"%s", step_name);
    else
        sprintf(path_name,"%s/%s", step_name, dir_name);

    ////TDC
    for (int i_ch=0; i_ch<NumHist; i_ch++)
    {   sprintf(obj_name,"%s/H1%sTDC_%i", path_name, step_name, i_ch+1);
        sprintf(obj_title,"TDC raw data, channel %i", i_ch+1);
        pH1RawTDC[i_ch] = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                        53500,0.,535000.,
                                        "#font[12]{t} / ch", "counts");
    }

    ////TDC 2D (channel independend)
    sprintf(obj_name,"%s/H2%sTDC", path_name, step_name);
    pH2RawTDC        = (TH2I *) MakeTH2('I',obj_name,
                                       "TDC raw data",
                                        53500,0.,535000.,
                                        2*NumTDC,0,NumTDC,"#font[12]{t} / ch",
                                       "TDC channel");
 }

//______________________________________________________________________________

//QDC histograms
void nfisHistograms::DefineRawQDC(const char* step_name, const char* dir_name,
                                  unsigned int NumHist)
{   //method to define QDC histograms for the unpacking step
    char    path_name[256] = "",
            obj_name[256] = "",
            obj_title[256] = "";

    //define directory structure of histograms
    if (strcmp(dir_name, "")==0)
        sprintf(path_name,"%s", step_name);
    else
        sprintf(path_name,"%s/%s", step_name, dir_name);

    //determine start and stop time of experiment
    UInt_t t_start = pParGlobal->ExpStart->GetSec();
    UInt_t t_end   = pParGlobal->ExpEnd->GetSec();
    //calculate the beam time in seconds
    UInt_t tdiff = t_end-t_start;
    //calculate number of bins
    UInt_t tbins = tdiff/pParGlobal->ReadoutPeriod;

    //create channel dependend historgrams
    for (int i_ch=0; i_ch<NumHist; i_ch++)
    {
        ////QDC (high gain)
        sprintf(obj_name,"%s/high/H1%sQDCh_%i",
                path_name, step_name, i_ch+1);
        sprintf(obj_title,"QDC (high gain), channel %i", i_ch+1);


        pH1RawQDCh[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                            "counts");
        ////QDC (low gain)
        sprintf(obj_name,"%s/low/H1%sQDCl_%i",
                path_name, step_name, i_ch+1);
        sprintf(obj_title,"QDC (low gain), channel %i", i_ch+1);

        pH1RawQDCl[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");

        ////QDC (low gain) vs Time
        sprintf(obj_name,"%s/low/QDCvsTime/H2%sQDCvsTime_%i", path_name,
                step_name, i_ch);
        sprintf(obj_title,"QDCS(low) vs Time channel %i", i_ch+1);
        pH2QDCTime[i_ch]   = (TH2I *) MakeTH2('I',obj_name,obj_title,
                                            350,0.,350.,
                                            0.25*tbins,t_start, t_end,
                                           "#font[12]{Q} / ch",
                                           "");
//*/
    }

    ////QDC (high gain) 2D
    sprintf(obj_name,"%s/high/H2%sQDCh", path_name, step_name);
    pH2RawQDCh     = (TH2I *) MakeTH2('I',obj_name,
                                   "QDC (high)",
                                    NumQDCReg,0,NumQDCReg,
                                    2*NumQDC,0,NumQDC,"#font[12]{Q} / ch",
                                   "QDC channel");
    ////QDC (low gain) 2D
    sprintf(obj_name,"%s/low/H2%sQDCl", path_name, step_name);
    pH2RawQDCl     = (TH2I *) MakeTH2('I',obj_name,
                                   "QDC (low gain)",
                                    NumQDCReg,0,NumQDCReg,
                                    2*NumQDC,0,NumQDC,"#font[12]{Q} / ch",
                                   "QDC channel");
//*/
}

////////////////////////////////////////////////////////////////////////////////
////          Histograms for the analysis steps                             ////
////////////////////////////////////////////////////////////////////////////////

void  nfisHistograms::DefineAnaCommon(const char* step_name,
                                      const char* dir_name)
{   //method to define common histograms for the analysis step
    char    path_name[256] = "",
            obj_name[256]  = "";

    //define directory structure of histograms
    sprintf(path_name,"%s", step_name);

    ////Hit 2D
    sprintf(obj_name,"%s/H2AnaHit", path_name);
    pH2AnaHit = (TH2I *) MakeTH2(
                    'I',obj_name, "Number of hits per event",
                    MAXHITS,0,MAXHITS,
                    2*(NumHZDRFC),0,NumHZDRFC,
                    "hits per event", "hit channel");

    ////TimeDiff to acc 2D
    sprintf(obj_name,"%s/H2AnaDt", path_name);
    pH2AnaDt = (TH2I *) MakeTH2(
                    'I',obj_name, "Time Diff to Acc (coarse bin)",
                    107000,-535000.,535000.,
                    2*NumTDC,0,NumTDC,
                   "#font[12]{t} / ch","TDC channel");

    ////TimeDiff to acc 2D (pulse height gated)
    sprintf(obj_name,"%s/H2AnaDtG", path_name, step_name);
    pH2AnaDt_g  = (TH2I *) MakeTH2(
                    'I',obj_name, "Time Diff to Acc (coarse bin, qdc gate on FF)",
                    107000,-535000.,535000.,
                    2*NumTDC,0,NumTDC,
                   "#font[12]{t} / ch", "TDC channel");
    //*/

}

//______________________________________________________________________________
void  nfisHistograms::DefineAnaTrigAcc(const char* step_name,
                                       const char* dir_name)
{   //method to define common histograms for the analysis step
    char    path_name[256] = "",
            obj_name[256]  = "",
            obj_title[256] = "";

    //define directory structure of histograms
    if (strcmp(dir_name, "")==0)
        sprintf(path_name,"%s", step_name);
    else
        sprintf(path_name,"%s/%s", step_name, dir_name);

    ////AccHits
        sprintf(obj_name,"%s/Hits/H1AnaAccHit", path_name);
        sprintf(obj_title,"Number of hits per event");

        pH1AnaHitAcc    = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            32,0,32,"hits per event",
                                           "counts");

    ////TrigHits
        sprintf(obj_name,"%s/Hits/H1AnaTrigHit", path_name);
        sprintf(obj_title,"Number of hits per event");

        pH1AnaHitTrig   = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            32,0,32,"hits per event",
                                           "counts");
//*/
}


//______________________________________________________________________________
void nfisHistograms::DefineAnaHZDR(const char* step_name, const char* dir_name,
                                   unsigned int NumHist)
{   //method to define analysis histograms for the HZDR FC
    char    path_name[256] = "",
            obj_name[256] = "",
            obj_title[256] = "";

    //define directory structure of histograms
    if (strcmp(dir_name, "")==0)
        sprintf(path_name,"%s", step_name);
    else
        sprintf(path_name,"%s/%s", step_name, dir_name);

    for (int i_ch=0; i_ch<NumHist; i_ch++)
    {
        ////Hits
        sprintf(obj_name,"%s/Hits/H1AnaHZDRHit_%i", path_name, i_ch+1);
        sprintf(obj_title,"Number of hits per event Channel %i", i_ch+1);


        pH1AnaHitHZDR[i_ch] = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                                32,0,32,"hits per event",
                                                "counts");
        //*/
        ////TimeDiff to acc
        sprintf(obj_name,"%s/TimeDiff/H1AnaHZDRDt_%i", path_name, i_ch+1);
        sprintf(obj_title,"Time Diff to Acc (coarse bin), channel %i", i_ch+1);
        pH1AnaDtHZDR[i_ch]  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            107000,-535000.,535000.,
                                           "#font[12]{t} / ch", "counts");

        ////TimeDiff to acc (pulse heigt gated)
        sprintf(obj_name,"%s/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", path_name, i_ch+1);
        sprintf(obj_title,"Time Diff to Acc (coarse bin, pulse height gate),"
                          " channel %i", i_ch+1);
        pH1AnaDtHZDR_g[i_ch]  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            107000,-535000.,535000.,
                                           "#font[12]{t} / ch", "counts");

        ////QDC (low gain)
        sprintf(obj_name,"%s/QDC/low/H1AnaQDCl_%i", path_name, i_ch+1);
        sprintf(obj_title,"QDC (low gain), channel %i", i_ch+1);

        pH1AnaQDCl[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");

        ////QDC (low gain, ToF gated, neutron-induced fissions)
        sprintf(obj_name,"%s/QDC/low/ToF_gated/NIF/H1AnaQDCl_NIF_%i", path_name, i_ch+1);
        sprintf(obj_title,"QDC (low gain, ToF gated, n-induced fission), channel %i", i_ch+1);

        pH1AnaQDCl_NIF[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");

        ////QDC (low gain, ToF gated, spontaneous fission)
        sprintf(obj_name,"%s/QDC/low/ToF_gated/SF/H1AnaQDCl_SF_%i", path_name, i_ch+1);
        sprintf(obj_title,"QDC (low gain, ToF gated, spontaneous), channel %i", i_ch+1);

        pH1AnaQDCl_SF[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");

        ////QDC (low gain, ungated, self triggered )
        sprintf(obj_name,"%s/QDC/low/trig/H1AnaQDCl_trig_%i", path_name, i_ch+1);
        sprintf(obj_title,"QDC (low gain, ungated, self triggered), channel %i", i_ch+1);

        pH1AnaQDCl_trig[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");
        //*/
    }

 }


////////////////////////////////////////////////////////////////////////////////
////          Histograms for the calibration steps                          ////
////////////////////////////////////////////////////////////////////////////////
/*
void nfisHistograms::DefineCalHZDR(const char* step_name, const char* dir_name,
                                   unsigned int NumHist)
{   //method to define analysis histograms for the HZDR FC
    char    path_name[256] = "",
            obj_name[256] = "",
            obj_title[256] = "";

    //define directory structure of histograms
    if (strcmp(dir_name, "")==0)
        sprintf(path_name,"%s", step_name);
    else
        sprintf(path_name,"%s/%s", step_name, dir_name);

    for (int i_ch=0; i_ch<NumHist; i_ch++)
    {
        ////QDC (low gain)
        sprintf(obj_name,"%s/QDC/low/H1CalQDCl_%i", path_name, i_ch+1);
        sprintf(obj_title,"QDC (low gain), channel %i", i_ch+1);

        pH1CalQDCl[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");

        ////QDC (high gain)
        sprintf(obj_name,"%s/QDC/high/H1CalQDCh_%i", path_name, i_ch+1);
        sprintf(obj_title,"QDC (high gain), channel %i", i_ch+1);

        pH1CalQDCh[i_ch]   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                            NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                           "counts");

        ////Time-of-Flight
        sprintf(obj_name,"%s/ToF/H1CalHZDRDt_%i", path_name, i_ch+1);
        sprintf(obj_title,"Time-of-Flight (coarse bin), channel %i", i_ch+1);
        pH1CalDtHZDR[i_ch]  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            107000,-13375.,13375.,
                                           "#font[12]{t} / ns", "counts");
        ////Time-of_Flight (pule heigt gated)
        sprintf(obj_name,"%s/ToF/PH-Gated/H1CalHZDRDtG_%i", path_name, i_ch+1);
        sprintf(obj_title,"Time-of-Flight (coarse bin, pulse height gate),"
                          " channel %i", i_ch+1);
        pH1CalDtHZDR_g[i_ch]  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            107000,-13375.,13375.,
                                           "#font[12]{t} / ns", "counts");

        ////kinetic energy
        sprintf(obj_name,"%s/E_kin/H1CalEkin_%i", path_name, i_ch+1);
        sprintf(obj_title,"Neutron kinetic energy,"
                          " channel %i", i_ch+1);
        pH1CalEkin[i_ch]  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                            50000,0.01,50.,
                                           "#font[12]{E} / MeV", "counts");
    }

    ////ToF, Sum
    sprintf(obj_name,"%s/ToF/H1CalHZDRDtSum", path_name);
    sprintf(obj_title,"Time-of-Flight (coarse bin), all channel");
    pH1CalDtHZDRSum  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                        107000,-13375.,13375.,
                                       "#font[12]{t} / ns", "counts");

    ////ToF, Sum (QDC-gated)
    sprintf(obj_name,"%s/ToF/PH-Gated/H1CalHZDRDtSum_g", path_name);
    sprintf(obj_title,"Time-of-Flight (coarse bin), all channel (QDC gated)");
    pH1CalDtHZDRSum_g  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                        107000,-13375.,13375.,
                                       "#font[12]{t} / ns", "counts");

    ////kinetic energy, Sum
    sprintf(obj_name,"%s/H1CalEkin", path_name);
    sprintf(obj_title,"Neutron kinetic energy");
    pH1CalEkinSum  = (TH1I *) MakeTH1('I',obj_name, obj_title,
                                        50000,0.01,50.,
                                       "#font[12]{E} / MeV", "counts");

    ////QDC (low gain) Sum
    sprintf(obj_name,"%s/QDC/low/H1CalQDClSum", path_name);
    sprintf(obj_title,"QDC (low gain), all channel");

    pH1CalQDClSum   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                        NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                       "counts");

    ////QDC (high gain) Sum
    sprintf(obj_name,"%s/QDC/high/H1CalQDChSum", path_name);
    sprintf(obj_title,"QDC (high gain), all channels");

    pH1CalQDChSum   = (TH1I *) MakeTH1('I',obj_name,obj_title,
                                        NumQDCReg,0,NumQDCReg,"#font[12]{Q} / ch",
                                       "counts");

    ////QDC (low gain) 2D
    sprintf(obj_name,"%s/QDC/low/H2CalQDCl", path_name);
    pH2CalQDCl     = (TH2I *) MakeTH2('I',obj_name,
                                   "QDC (low gain)",
                                    NumQDCReg,0,NumQDCReg,
                                    2*NumQDC,0,NumQDC,"#font[12]{Q} / ch",
                                   "QDC channel");


    ////QDC (high gain) 2D
    sprintf(obj_name,"%s/QDC/high/H2CalQDCh", path_name);
    pH2CalQDCh     = (TH2I *) MakeTH2('I',obj_name,
                                   "QDC (high)",
                                    NumQDCReg,0,NumQDCReg,
                                    2*NumQDC,0,NumQDC,"#font[12]{Q} / ch",
                                   "QDC channel");

    ////ToF vs QDC (low gain) 2D
    sprintf(obj_name,"%s/H2CalToFvsQDCl", path_name);
    pH2CalToFvsQDC     = (TH2I *) MakeTH2('I',obj_name,
                                   "ToF vs. QDC (low)",
                                   250,20.,45.,
                                   256,0., NumQDCReg,
                                   "#font[12]{t} / ns",
                                   "QDC channel");

 }//*/

//______________________________________________________________________________
vector<string> *nfisHistograms::ScalerChannelNames()
{   //creates a QList containing the names of the scaler channels

string strarray[] = {
    "FC ch 1"   ,"FC ch 2"      ,"FC ch 3"      ,"FC ch 4"      ,"FC ch 5",     //5
    "FC ch 6"   ,"FC ch 7"      ,"FC ch 8"      ,"LaBr3 1"      ,"LaBr3 2",     //10
    "LaBr3 3"   ,"LaBr3 4"      ,"LaBr3 5"      ,"unused"       ,"unused",      //15
    "unused"    ,"HPGe 1"       ,"HPGe 2"       ,"HPGe 3"       ,"HPGe 4",      //20
    "HPGe 5"    ,"unused"       ,"unused"       ,"unused"       ,"Accelerator", //25
    "SOR"       ,"unused"       ,"unused"       ,"unused"       ,"unused",      //30
    "unused"    ,"unused"       ,"in_or_00-07"  ,"in_or_08-15"  ,"in_or_16-23", //35
    "pl_coin_0" ,"pl_coin_1"    ,"global_or"    ,"trigger_raw"  ,"trigger",     //40
    "trigger_ds","LCLK"         ,"LCLK & veto"  ,"LCLK & nveto" ,"unused",      //45
    "t_{dead}"  ,"t_{live}"     ,"t_{real}"     ,"ADC words"    ,"ADC events",  //50
    "QDC1 words","QDC1 events"  ,"QDC2 words"   ,"QDC2 events"  ,"TDC words",   //55
    "TDC events","veto words"   ,"lmd words"    ,"unused"       ,"unused",      //60
    "unused"    ,"unused"       ,"unused"       ,"target pos"};


    vector<string> *ChannelNames = new vector<string>(strarray, strarray + 64);

    return ChannelNames;
}
//*/
//------------------------end of nfisHistograms.cxx-----------------------------


//////QDC (high gain)
//if (i_ch<NumHist-1)         // histogram including one Pre-Amp
//{   sprintf(obj_name,"%s/high/%sH1QDCh_%i",
//            path_name, step_name, i_ch+1);
//    sprintf(obj_title,"QDC (high gain), channel %i", i_ch+1);
//}
//else                        // histogram including all Pre-Amps
//{   sprintf(obj_name,"%s/high/%sH1QDCh_%i",
//            path_name, step_name, i_ch+1);
//    sprintf(obj_title,"QDC (high gain)");

//}
//*/
