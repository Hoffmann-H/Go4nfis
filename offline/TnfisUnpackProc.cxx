#include "TnfisUnpackProc.h"
#include "Riostream.h"
#include "TH1.h"

//#include "TnfisAnaDataWord.h"
#include "TnfisUnpackEvent.h"


//***********************************************************
TnfisMakeUnp::TnfisMakeUnp() : TGo4EventProcessor("Proc")
{
  cout << "**** TnfisMakeUnp: Create std. instance " << endl;
  GlobalScaler = 0;        //should not be reseted
}
//***********************************************************
TnfisMakeUnp::~TnfisMakeUnp()
{
  cout << "**** TnfisMakeUnp: Delete instance " << endl;
}
//***********************************************************
// this one is used in standard factory
TnfisMakeUnp::TnfisMakeUnp(const char* name) : TGo4EventProcessor(name)
{
    cout << "**** TnfisMakeUnp: Create par instance " << name << endl;
    //init flags and parameters
    CommentFlag[0]  = false;    //print current data word
    CommentFlag[1]  = false;    //print measurement times
    CommentFlag[2]  = false;    //print time flags
    CommentFlag[3]  = false;    //print scaler values
    CommentFlag[4]  = false;    //print veto scaler values
    CommentFlag[5]  = false;    //print adc values
    CommentFlag[6]  = false;    //print qdc values
    CommentFlag[7]  = false;    //print tdc values
    CommentFlag[8]  = false;    //print trigger values
    CommentFlag[9]  = false;    //print current absorber values
    CommentFlag[10] = false;    //print absolute time
    CommentFlag[11] = false;    //print loop information

    PrintScalerFlag = 0;
    l_scaler        = 0;

    LNE_time        = 0.;

    t_live = t_real = t_dead = 0;
    //create time stamp object for absolute time
    AbsTime = new TTimeStamp();

    //load a vector containing the descriptions of each scaler channel
    pChannelNames = new vector<string>();
    pChannelNames = ScalerChannelNames();

    //init parameter
    pParGlobal= (TnfisParamGlobal *) MakeParameter("GlobalPar", "TnfisParamGlobal");

    // Create of parameters and histograms (check if restored from auto save file):
    pHist = new nfisHistograms();
    pHist->DefineHistograms("Raw");

    //pictures
    MakePictures();

}
//-----------------------------------------------------------
// event function
Bool_t TnfisMakeUnp::BuildEvent(TGo4EventElement* output)
{
    Int_t lwords;
    Int_t *pdata;
    TGo4MbsSubEvent * rawsub;

    //load nELBE data structure
    nELBE_data *u_data = new nELBE_data();

    // define and init event pointers
    TGo4MbsEvent * raw    = (TGo4MbsEvent* ) GetInputEvent();

    //get pointer to the event source
    TGo4EventSource *EvtSource = raw->GetEventSource();

    //define output event structure pointer
    TnfisUnpackEvent *pTarget = (TnfisUnpackEvent*) output;
    TnfisUnpackEvent& ev=*pTarget; // ref instead pointer for array syntax:
    //init all items of the target classes, clear their content
    InitEventStructure(pTarget);

    //draw line to pH2TDCTime histogram, if absorber changes
/*    if ((raw->GetTrigger()==14) ||
        (raw->GetTrigger()==15))
    {   for (int i_ch=0; i_ch<NumHZDRFC; i_ch++)
            DrawLine(pPicQDCvsTime[i_ch], pHist->pH2QDCTime[i_ch],
                     pTimer->GetRealTime(), "trig", raw->GetTrigger());
    }
*/
    //check, if there is a valid input event
    if (raw==0)
    {  cout << "TnfisUnpackProc: no input event !"<< endl;
       return kFALSE;
    }

    // prepare loop over sub-events
    raw->ResetIterator();

    // loop over subevents
    if (CommentFlag[11] == true)
        cout << "Starting loop over sub-events" << endl;
    while((rawsub = raw->NextSubEvent()) != 0)
    {   //clear the event counter
        ResetEventStructure();

//        time_t AbsTime = (time_t)raw->GetMbsBufferHeader()->l_time;



        if (CommentFlag[11] == true)
            cout << "Beginning another sub-event loop" << endl;

        pdata  = rawsub->GetDataField(); // pointer to first data longword
        lwords = rawsub->GetIntLen();    // number of longwords

        if ((raw->GetTrigger()==14) || (raw->GetTrigger()==15))
            cout << "Type: " << raw->GetTrigger() << endl;

        while (lwords>0)
        {   //fill nELBE data structure with MBS formated data
            if (CommentFlag[11] == true)
                cout << "Analyze word number " << lwords << endl;
            u_data->value= (*pdata++) & 0xffffffff; lwords--;

            if (CommentFlag[0] == true)
                cout <<"<DATAWORD>: " << std::hex << u_data->value << endl;

            /*** analyze dataword and put into event array ******************/
            if (CommentFlag[11] == true)
                cout << "u_data->common.geo: " << std::hex << u_data->common.geo << endl;
            //classify the unpacking after geo id
            switch(u_data->common.geo)
            {
                case 0: // measurement times
                    AnalyzeMeasurementTimes(u_data);
                break;

                case 1: // time Flag
                    AnalyzeTimeFlag(u_data);
                break;

                case 4: // Scaler
                    AnalyzeScaler(u_data, &pdata, &lwords, pChannelNames,
                                  pTimer);
                break;

                case 6: // Veto Scaler
                    AnalyzeVetoScaler(u_data, &pdata, &lwords, pTarget,
                                      pEvtCounter);
                break;

                case 8: // TDC
                    AnalyzeTDC(u_data, pTarget, pEvtCounter);
                break;

                case 9: // Trigger time tag
                    AnalyzeTriggerTimeTag(u_data);
                break;

//                case 10: // HPGe ADC

//                break;

                case 12: // QDC FC 1
                    AnalyzeQDC(u_data, pTarget, pEvtCounter);
                break;

//                case 14: // QDC FC 2

//                break;

                case 28: // Absorber change flag
                    AnalyzeAbsorberFlag(u_data);
                break;

                case 31: // Absolute time
                    AnalyzeAbsoluteTime(u_data, pEvtCounter);
                break;

            default:
                cerr << endl
                     << "<ERROR>: Unknown GEO: "
                     << std::hex << u_data->common.geo
                     << ", data " << std::hex << u_data->value
                     << endl;
            break;
            } // switch (u_data->common.geo)

        } // while lwords>0

    } // while subevents

  pTarget->SetValid(kTRUE);
//  FillHistograms(pTarget); // fill histograms from output event

  return kTRUE;
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeMeasurementTimes(nELBE_data *u_data)
{// prints the real and and life time of analyzed file(s) after finishing
    if (CommentFlag[1]==true)
    {   cout << "\tMeasurement time" << endl;
        if (u_data->time.flag == 0)
            cout    <<"<I>: Measurement was running for \t"
                    << (Float_t)u_data->time.data/10
                    <<" s."
                    << endl;
        if (u_data->time.flag == 1)
            cout    << endl
                    <<"<I>: Live Time was \t \t"
                    << (Float_t)u_data->time.data/10
                    <<" s."
                    << endl;
    }
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeTimeFlag(nELBE_data *u_data)
{   //prints the scaler time flag if
    LNE_time = (Double_t)u_data->time_flag.data; // in ms

    if (CommentFlag[2] == true)
        cout    << "\tscaler time flag\tTime Flag "
                << l_scaler << " after "
                << LNE_time/1000. <<" sec."
                << endl;
    if (PrintScalerFlag == 1)
        cout    << "<I>: Time Flag " << l_scaler
                << " after " << LNE_time/1000.
                << endl;

    l_scaler+=1;
//    if (l_scaler==MAXSCALERCHANNEL) l_scaler=1;
//    else                            l_scaler+=1;

    GlobalScaler++;
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeScaler(nELBE_data *u_data,
                                     Int_t **pdata, Int_t *lwords,
                                     vector<string> *ChannelNames,
                                    TnfisTimer *Timer)
{   //prints out veto scaler values and fills the veto histograms
    Double_t rate;

    //determine the number of data word from the first word (header)
    dat[0] = u_data->scaler.words;

    if (CommentFlag[3] == true)
        cout << "\tScaler\t" << dat[0]
             << " words to read" << endl;

    for (UInt_t i_ch = 0; i_ch < dat[0]; i_ch++)
    {   // read one more data word = Scaler data
        u_data->value = **pdata;
        // scaler value
        dat[1] = u_data->value;

        (*pdata)++; (*lwords)--;


        //calculate the rate in 1/s (because [LNE_time] = ms
        if (LNE_time!=0)
            rate = (Double_t)dat[1]*1000./(Double_t)LNE_time;

        if (PrintScalerFlag == 1)
            cout    << "<I> Scaler Channel " << i_ch << " "
                    << ChannelNames->at(i_ch)
                    << ":   " << dat[1]
                    << "counts = " << rate/1000. << "kHz"
                    << endl;

        if (dat[1] > 0)
        {   //calculate dead, live and real time
            if ((i_ch == 45) || (i_ch == 46) || (i_ch == 47))
            {   // $t_dead/_live/_real = counts / 40 MHz
                Double_t time = rate * LNE_time/1000. / 40.0e6;

                switch(i_ch)
                {
                    case 45:    //dead time
                        t_dead += time;
                        Timer->SetDeadTime(t_dead);
                    break;

                    case 46:    //live time
                        t_live += time;
                        Timer->SetLiveTime(t_live);
                    break;

                    case 47:    //real time
                        t_real += time;
                        Timer->SetRealTime(t_real);
                    break;
                }
                rate = time;

            } //end if ((i_ch == 45) || (i_ch == 46) || (i_ch == 47))

            //Fill histograms
            // eiter...
//            pHist->pH1RawRate[i_ch]->Fill(l_scaler, rate);
//            pHist->pH2RawRate->Fill(l_scaler, i_ch, rate);
            // ...or...
            pHist->pH1RawRate[i_ch]->Fill(AbsTime->GetSec(), rate);
            pHist->pH2RawRate->Fill(AbsTime->GetSec(), i_ch, rate);

        } // end if (dat[1] > 0)
    }

}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeVetoScaler(nELBE_data *u_data,
                                     Int_t **pdata, Int_t *lwords,
                                     TnfisUnpackEvent *pTarget,
                                     TnfisEvtCounter *EvtCounter)
{   /*/prints out veto scaler values and fills the veto histograms

    TnfisUnpackEvent& ev=*pTarget; // ref instead pointer for array syntax:
//    TnfisVetoLength *pVetoLength = dynamic_cast<TnfisVetoLength*>(&ev[6]);

    //determine the number of data words from the first word (header)
    dat[0] = u_data->veto.words;

    if (CommentFlag[4] == true)
        cout << "\tVeto scaler\t" << dat[0]
             << " words to read" << endl;

    for (UInt_t i_ch = 0; i_ch < dat[0]; i_ch++)
    {   // read one more data word = veto data
        u_data->value = **pdata;
        // veto time in ns
        dat[1] = u_data->value * 25;

        (*pdata)++; (*lwords)--;


        if (CommentFlag[4] == true)
            cout    << "\tVeto data word " << i_ch
                    << ": Event " << EvtCounter->GetVeto()
                    <<" Channel " << i_ch  << " Value " << dat[1] << endl;

        pHist->pH1RawVeto[i_ch]->Fill(dat[1]);
        pHist->pH2RawVeto->Fill(dat[1],i_ch+1);

        Double_t vetolength = ((Double_t)u_data->value +
                               randtriang()) * 25.;

        if (i_ch == 5) // total veto
        {   //pVetoLength->SetVetoLength(EvtCounter->GetVeto(), vetolength);
            pHist->pH1RawVetoTot->Fill(vetolength);
        }
    }
    //increase veto counter
    EvtCounter->IncreaseVeto();

    if (EvtCounter->GetVeto() >= MAXEVENTS)
    {   cerr    << "<ERROR>: Reached maximum number of VETO "
                << "events per readout event"
                << endl;
        EvtCounter->SetVeto(MAXEVENTS-1);
    }

//*/
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void TnfisMakeUnp::AnalyzeQDC(nELBE_data *u_data, TnfisUnpackEvent *pTarget,
                              TnfisEvtCounter *EvtCounter)
{   if (CommentFlag[6] == true) cout <<"\tQDC" << endl;
    switch (u_data->qdc.flag)
    {
        case 0x0 : // data word
            AnalyzeQDCDataWord(u_data, pTarget, EvtCounter);
            break;

        case 0x2 : // header
            AnalyzeQDCHeader(u_data, EvtCounter);
        break;

        default :
            cerr    << "<ERROR>: user readout: UNKNOWN QDC HEADER "
                    << u_data->qdc.flag
                    << " data word: " << u_data->value
                    << endl;
        break;
    } //switch (u_data->qdc.flag)
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeQDCDataWord(nELBE_data *u_data,
                                      TnfisUnpackEvent *pTarget,
                                      TnfisEvtCounter *EvtCounter)
{   ULong_t data[6];

    //analyzes the qdc data word
    data[0] = u_data->qdc_data.geo;  // GEO
    data[1] = u_data->qdc_data.ch;   // Channel
    data[2] = u_data->qdc_data.r;    // Range
    data[3] = u_data->qdc_data.un;   // below threshold
    data[4] = u_data->qdc_data.ov;   // overflow
    data[5] = u_data->qdc_data.data; // QDC value

    if (CommentFlag[6] == true)
    cout << "\tData\tGEO = " << data[0]
    << " Channel = " << data[1]
    << " Range = " << data[2]
    << " UN/OV = " << data[3] << "/" << data[4]
    << " Value = " << data[5]
    << endl;
    //set count to QDC value of the event counter
    count = EvtCounter->GetQDC();

    AnalyzeQDCRange(data, pTarget, EvtCounter);

}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeQDCRange(ULong_t *data, TnfisUnpackEvent *pTarget,
                                   TnfisEvtCounter *EvtCounter)
{//analyzes the qdc range dependence (high and low gain range)

    ULong_t     ch          = data[1];  //channel
    ULong_t     GainMode    = data[2];  //gain mode (high/low)
     Long_t     Value       = data[5];  //QDC value

    Bool_t      Underflow   = data[3];  //event below threshold
    Bool_t      Overflow    = data[4];  //overflow

    if (Underflow)  Value   = -1;
    if (Overflow)   Value   = NumQDCReg+1;

    switch (GainMode) // high or low gain mode
    {
        case 0: // low gain
            pHist->pH1RawQDCl[ch]->Fill(Value);
            pHist->pH2RawQDCl->Fill(Value, ch);

            pHist->pH2QDCTime[ch]->Fill(Value, AbsTime->GetSec());

            if ((ch >= 0) && (ch < NumHZDRFC))
            {   if (Underflow)
                    pHZDRPreAmp[ch]->SetQDCl(EvtCounter->GetQDC(), -10);
                else if (Overflow)
                    pHZDRPreAmp[ch]->SetQDCl(EvtCounter->GetQDC(), NumQDCReg+10);
                else
                pHZDRPreAmp[ch]->SetQDCl(EvtCounter->GetQDC(), Value);
            }
            else
                cerr    << "<ERROR>: Unknown QDC channel "
                        << ch << endl;


        break;

        case 1: // high gain
            pHist->pH1RawQDCh[ch]->Fill(Value);
            pHist->pH2RawQDCh->Fill(Value, ch);

            if ((ch >= 0) && (ch < NumHZDRFC))
            {   if (Underflow)
                    pHZDRPreAmp[ch]->SetQDCh(EvtCounter->GetQDC(), -10);
                else if (Overflow)
                    pHZDRPreAmp[ch]->SetQDCh(EvtCounter->GetQDC(), NumQDCReg+10);
                else
                pHZDRPreAmp[ch]->SetQDCh(EvtCounter->GetQDC(), Value);
            }
            else
                cerr << "<ERROR>: Unknown QDC channel "
                     << ch << endl;
        break;

        default: break;
    }

}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeQDCHeader(nELBE_data *u_data,
                                    TnfisEvtCounter *EvtCounter)
{//analyzes the qdc header
    ULong_t data[3];

    data[0] = u_data->qdc_header.geo; // GEO
    data[1] = u_data->qdc_header.crt; // Crate
    data[2] = u_data->qdc_header.cnt; // Cnt

    if (CommentFlag[6] == true)
    {   cout << "\tHeader\tGEO = " << data[0]
        << " Crate = " << data[1]
        << " Cnt = " << data[2]
        << endl;
        cout <<"\tEventnumber = " << EvtCounter->GetQDC()
        << endl;
    }

    //increase QDC event counter
    EvtCounter->IncreaseQDC();

    if (EvtCounter->GetQDC() >= MAXEVENTS)
    {   cerr    << "<ERROR>: Reached maximum number of QDC "
                << "events per readout event" << endl;
        EvtCounter->SetQDC(MAXEVENTS-1);
    }

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void TnfisMakeUnp::AnalyzeTDC(nELBE_data *u_data, TnfisUnpackEvent *pTarget,
                              TnfisEvtCounter *EvtCounter)
{   if (CommentFlag[7] == true) cout << "\tTDC" << endl;

    switch(u_data->tdc.flag)
    {
        case 0 : //global trailer
            AnalyzeTDCGlobalTrailer(EvtCounter);
        break;

        case 1 : // TDC measurement
            if (CommentFlag[7] == true)
                cout << "\tData\t=>" << endl;

            ULong_t data[2];
            data[0] = u_data->tdc_data.data;    //measurement
            data[1] = u_data->tdc_data.ch;      //channel

            if (CommentFlag[7] == true)
                cout    <<"GEO = " << u_data->common.geo
                        <<" Channel = " << data[1]
                        <<" Value = " << data[0]
                        << endl;

            channelok = 1;

            switch (data[1])
            {
                case 0: //HZDR FC 1
                case 1: //HZDR FC 2
                case 2: //HZDR FC 3
                case 3: //HZDR FC 4
                case 4: //HZDR FC 5
                case 5: //HZDR FC 6
                case 6: //HZDR FC 7
                case 7: //HZDR FC 8
                    if (CommentFlag[7] == true) cout << "HZDR FC" << endl;

                    count = pHZDRPreAmp[data[1]]->GetHitCounter(EvtCounter->GetTDC());

                    if (count < MAXHITS-1)
                        // to avoid array overflow all hits beyond number 31
                        // will be put into the same memory
                        pHZDRPreAmp[data[1]]->SetHit(EvtCounter->GetTDC(), count, data[0]);
                    else
                        pHZDRPreAmp[data[1]]->SetHit(EvtCounter->GetTDC(), MAXHITS-1, data[0]);

                    pHZDRPreAmp[data[1]]->IncreaseHitCounter(EvtCounter->GetTDC());
                break;

                case  8: // LaBr3 1
                case  9: // LaBr3 2
                case 10: // LaBr3 3
                case 11: // LaBr3 4
                case 12: // LaBr3 5
                break;


                case 24: //accelerator
                {   if (CommentFlag[7] == true) cout << "Accelerator" << endl;

                    AnalyzeTDCTrigAccClkVeto(data, EvtCounter, pAcc);
                }
                break;

                case 16: // HPGe 1
                case 17: // HPGe 2
                case 18: // HPGe 3
                case 19: // HPGe 4
                case 20: // HPGe 5
                break;

                case 28: //TDC trigger
                {   if (CommentFlag[7] == true) cout << "TDC trigger" << endl;

                }
                break;

                case 29: //trigger
                {   if (CommentFlag[7] == true) cout << "Trigger" << endl;

                    AnalyzeTDCTrigAccClkVeto(data, EvtCounter, pTrig);
                }
                break;

                case 30: //veto
                {   if (CommentFlag[7] == true) cout << "Veto" << endl;

                    AnalyzeTDCTrigAccClkVeto(data, EvtCounter, pVeto);
                }
                break;

                case 31: //not veto
                {   if (CommentFlag[7] == true) cout << "Not Veto" << endl;

                    AnalyzeTDCTrigAccClkVeto(data, EvtCounter, pNotVeto);
                }
                break;


                default:
                {
                    channelok = 0;
                    cerr << "<ERROR> Unknown TDC 0 channel " << data[1]
                         << " had data: " << u_data->value
                         << "!" << endl;
                }
                break;
            } // switch (data[1])

            if (channelok == 1)
            {   pHist->pH1RawTDC[data[1]]->Fill(data[0]);    //TDC raw data
                pHist->pH2RawTDC->Fill(data[0], data[1]);    //channel bin
            }
    break;

    default:
        cerr << "<ERROR> Unknown TDC header: geo: " << u_data->common.geo
             << " header: " << std::hex << u_data->tdc.flag
             << " data: " << std::hex << u_data->value
             << endl;
    break;
    } // switch (u_data->tdc.flag)
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeTDCGlobalTrailer(TnfisEvtCounter *EvtCounter)
{   //analyzes the global trailer of the tdc
    if (CommentFlag[7] == true)
        cout << "\tTrailer (eventnumber="
             << EvtCounter->GetTDC()
             << ")" << endl;

    //increase plastic TDC event counter
    EvtCounter->IncreaseTDC();

    if (EvtCounter->GetTDC() >= MAXEVENTS)
    {   cerr    << "<ERROR>: Reached maximum number of TDC"
                << "events per readout event"
                << endl;
        EvtCounter->SetTDC(MAXEVENTS-1);
    }

    // total number of events -> for comparison reasons
    EvtCounter->SetTotal(EvtCounter->GetTDC());
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeTDCTrigAccClkVeto(ULong_t *data,
                                            TnfisEvtCounter *EvtCounter,
                                            TnfisTrigAccClk *pTrigAccClk)
{//handles nELBE Event structure for Trigger, Accelerator, Clock and Vetos
    Long_t CounterTDC = EvtCounter->GetTDC();

    count = pTrigAccClk->GetHitCounter(CounterTDC);
    // to avoid array overflow all hits beyond number 31
    // will be put into the same memory
    if (count < MAXHITS-1)
        pTrigAccClk->SetHit(CounterTDC, count, data[0]);
    else
    {   channelok = 0;
        pTrigAccClk->SetHit(CounterTDC, MAXHITS-1, data[0]);
    }
    pTrigAccClk->IncreaseHitCounter(CounterTDC);
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeTriggerTimeTag(nELBE_data *u_data)
{//analyzes the trigger time tag data words
    dat[0] = u_data->tdc_trg_time.geo; // geo
    dat[1] = u_data->tdc_trg_time.data; // trigger time tag

    if (CommentFlag[8] == true)
        printf("\tTrigger Time Tag: 0x%07lx (last was 0x%07lx)",
               dat[1],lasttriggertime);

    if (dat[1] < lasttriggertime)
      triggertimediff = dat[1] + 0x7ffffff - lasttriggertime;
    else
      triggertimediff = dat[1] - lasttriggertime;

    if (CommentFlag[8] == true)
        printf("\tTime to last tigger = %ld ns\n",triggertimediff*800);

    lasttriggertime = dat[1];

    if (triggertimediff >= 40000)
    {   if (triggertimediff >= 1039900)
        {   cerr    <<"<ERROR>: Value out of bounds of TRGTIME: "
                    << triggertimediff << endl;
            triggertimediff = 49999;
        }
        else
            triggertimediff = (triggertimediff - 40000)/100 + 40000;
    }
    pHist->pH1RawTRGTime->Fill(triggertimediff);

}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeAbsorberFlag(nELBE_data *u_data)
{//analyszes the absorber change flag data words

    if (CommentFlag[9] == true) cout << "\tAbsorber flag" << endl;

    dat[0] = u_data->absorber.geo; // geo
    dat[1] = u_data->absorber.abs; // absorber

    if (CommentFlag[9] == true) cout << "\t" << dat[1] << endl;
    else cout << "<I> New Absorber " << dat[1] << endl;


    pAbs->SetAbsorber(dat[1]);

    //draw line to pH2TDCTime histogram, if absorber changes
    for (int i_ch=0; i_ch<NumHZDRFC; i_ch++)
        DrawLine(pPicQDCvsTime[i_ch], pHist->pH2QDCTime[i_ch],
                 pTimer->GetRealTime(), "abs");
}

//______________________________________________________________________________
void TnfisMakeUnp::AnalyzeAbsoluteTime(nELBE_data *u_data,
                                       TnfisEvtCounter *EvtCounter)
{
    switch (u_data->abstime.id)
    {
        case 0:
            AbsTime->SetSec(0); AbsTime->SetNanoSec(0);
            AbsTime->SetSec((u_data->abstime.data & 0xffff) << 16);
            if (CommentFlag[10] == true)
              cout << "Time, case 0: " << AbsTime->AsString() << endl;

        break;

        case 1:
            AbsTime->SetSec( AbsTime->GetSec() | (u_data->abstime.data & 0xffff));
            if (CommentFlag[10] == true)
              cout << "Time, case 1: " << AbsTime->AsString() << endl;

        break;

        case 2:
            AbsTime->SetNanoSec(1000*((u_data->abstime.data & 0xffff) << 16));
            if (CommentFlag[10] == true)
              cout << "Time, case 2: " << AbsTime->AsString() << endl;

        break;

        case 3:
            AbsTime->SetNanoSec( 1000*((AbsTime->GetNanoSec()/1000) | (u_data->abstime.data & 0xffff)));

            //save time stamp in output event
            EvtCounter->SetAbsTime(*AbsTime);

            if (CommentFlag[10] == true)
                cout << "Absolute time: " << AbsTime->AsString() << endl;

        break;

        default: break;
    }
}

//______________________________________________________________________________
void TnfisMakeUnp::InitEventStructure(TnfisUnpackEvent *pTarget)
{   //clears all items of the target event structure
    TnfisUnpackEvent& ev=*pTarget;
    pTrig           = dynamic_cast<TnfisTrigAccClk*>(&ev[0]);
    pAcc            = dynamic_cast<TnfisTrigAccClk*>(&ev[1]);
    pEvtCounter     = dynamic_cast<TnfisEvtCounter*>(&ev[2]);
//    pClock          = dynamic_cast<TnfisTrigAccClk*>(&ev[2]);
    pTimer          = dynamic_cast<TnfisTimer*>(&ev[3]);
    pVeto           = dynamic_cast<TnfisTrigAccClk*>(&ev[4]);
//    pAbs            = dynamic_cast<TnfisAbsorber*>(&ev[4]);
    pNotVeto        = dynamic_cast<TnfisTrigAccClk*>(&ev[5]);
    pHZDRFC         = dynamic_cast<TnfisUnpFC*>(&ev[6]);
        for (int i_preAmp = 0; i_preAmp < NumHZDRFC; i_preAmp++)
            pHZDRPreAmp[i_preAmp] = dynamic_cast<TnfisUnpPreAmp*>(
                    pHZDRFC->GetPreAmp(i_preAmp));
//    pVetoLength     = dynamic_cast<TnfisVetoLength*>(&ev[6]);

    ResetEventStructure("init");

}

//______________________________________________________________________________
void TnfisMakeUnp::ResetEventStructure(Option_t *t)
{
    //clear the event structures
    pTrig->Clear(t);
    pAcc->Clear(t);
//    pClock->Clear(t);
    pVeto->Clear(t);
    pNotVeto->Clear(t);
//    pVetoLength->Clear(t);
    pEvtCounter->Clear(t);
//    pAbs->Clear(t);
    pHZDRFC->Clear(t);
    for (int i_preAmp = 0; i_preAmp < NumHZDRFC; i_preAmp++)
        pHZDRPreAmp[i_preAmp]->Clear(t);
}

//______________________________________________________________________________
void TnfisMakeUnp::DrawLine(TGo4Picture *pPic, TH2I *pHist, Double_t time,
                             const char *kind, UInt_t value)
{//draws a horizontal line into a histogram within a picture
    //get number of bins and lower and upper edge of the histogram
    UInt_t nbins = pHist->GetNbinsX();
    Double_t xmin, xmax;
    xmin = pHist->GetXaxis()->GetBinLowEdge(1);
    xmax = pHist->GetXaxis()->GetBinLowEdge(nbins);

    //create horizontal line at time from one (x)end to the other
    TLine *pLine = new TLine(xmin, time, xmax, time);

    if (strcmp(kind,"abs")==0)
    {   //set dot line dot style
        pLine->SetLineStyle(5);
        pLine->SetLineWidth(2);
        TText *pAbsText = new TText(0.75*xmax,time+20,
                                    Form("Abs: %i", value));
            pAbsText->SetTextFont(63);
            pAbsText->SetTextSize(18);

        pPic->AddSpecialObject(pAbsText);

    }
    else if (strcmp(kind,"run")==0)
    {   //set thick line
        pLine->SetLineWidth(3);

        TText *pRunText = new TText(1.01*xmax,time,
                                    Form("%s", CurrentRun.c_str()));
            pRunText->SetTextFont(63);
            pRunText->SetTextSize(23);

        pPic->AddSpecialObject(pRunText);
    }
    else if (strcmp(kind,"trig")==0)
    {   pLine->SetLineWidth(5);

        if (value==14)
            pLine->SetLineColor(632);
        else if (value==15)
            pLine->SetLineColor(600);
    }

    pPic->AddSpecialObject(pLine);
}

//______________________________________________________________________________
string TnfisMakeUnp::ExtractRun (const string& str)
{ // extracts the run number from filename
  unsigned found = str.find_last_of("/\\");

  string file = str.substr(found+1);
  string run  = file.substr(9,4);

  return run;
}


//______________________________________________________________________________
vector<string> *TnfisMakeUnp::ScalerChannelNames()
{   //creates a QList containing the name of the scaler channels

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

//______________________________________________________________________________
Double_t TnfisMakeUnp::randtriang()
{   // method to create an triang distr to smear cal. histograms
  Double_t xx1,xx2;

  xx1 = (Double_t)rand() / (Double_t)RAND_MAX;
  xx2 = (Double_t)rand() / (Double_t)RAND_MAX;

  return(xx1-xx2);
}


//-----------------------------------------------------------
// histogramming function
Bool_t TnfisMakeUnp::FillHistograms(TnfisUnpackEvent * event)
{

  return kTRUE;
}

//______________________________________________________________________________
void TnfisMakeUnp::MakePictures()
{//set up all pictures
    char    obj_name[256] = "",
            obj_title[256] = "";


    ////QDC vs Time picture (for each channel one)
    for (int i_ch=0; i_ch<NumHZDRFC; i_ch++)
    {   sprintf(obj_name,"PicQDCvsTime_%i", i_ch);
        sprintf(obj_title,"Picture of QDC(low) vs Time channel %i", i_ch+1);
        pPicQDCvsTime[i_ch] = new TGo4Picture(obj_name, obj_title);
        pPicQDCvsTime[i_ch]->AddH1(pHist->pH2QDCTime[i_ch]);
        pPicQDCvsTime[i_ch]->SetLogScale(2);
        AddPicture(pPicQDCvsTime[i_ch],"Raw/QDCvsTime"); // add picture to frame work
    }

    Int_t xPic = 0;
    Int_t yPic = 0;
    pPicQDCvsTimeAll = new TGo4Picture("PicQDCvsTime",
                                       "Picture of QDC(low) vs Time",2,4);
    for (int i_ch=0; i_ch<NumHZDRFC; i_ch++)
    {   xPic = i_ch%4; yPic = i_ch/4;
        pPicQDCvsTimeAll->Pic(yPic,xPic)->AddH1(pHist->pH2QDCTime[i_ch]);
        pPicQDCvsTimeAll->Pic(yPic,xPic)->SetLogScale(2);
    }
    AddPicture(pPicQDCvsTimeAll, "Raw");


    ////QDC vs Time picture (combine all sub-pictures to one)
    pPicQDCl = new TGo4Picture("PicQDCl","Picture of low gain QDC values of "
                               " all HZDR FC PreAmps");
    pPicQDCh = new TGo4Picture("PicQDCh","Picture of high gain QDC values of "
                               " all HZDR FC PreAmps");

    //get qdc high and low histogram of each pre-Amp
    for (int i_preAmp=0; i_preAmp<NumHZDRFC; i_preAmp++)
    {   pPicQDCl->AddH1(pHist->pH1RawQDCl[i_preAmp]);
        pPicQDCh->AddH1(pHist->pH1RawQDCh[i_preAmp]);
    }

    //set log scale to pictures
    pPicQDCl->SetLogScale(1);   pPicQDCh->SetLogScale(1);
    //add picture to frame work
    AddPicture(pPicQDCl, "Raw");       AddPicture(pPicQDCh, "Raw");


    ////picture of count-rates
    pPicRate = new TGo4Picture("PicRate","Picture of count-rates of "
                               " all HZDR FC PreAmps");

    //add rate histograms of all FC channel
    for (int i_preAmp=0; i_preAmp<NumHZDRFC; i_preAmp++)
        pPicRate->AddH1(pHist->pH1RawRate[i_preAmp]);

    char TimeFormat[100];

    sprintf(TimeFormat, "#splitline{%%d.%%m.}{%%H:%%M}%%F1970-01-01 00:00:00s0");

    pPicRate->SetXAxisTimeDisplay(kTRUE);
    pPicRate->SetXAxisTimeFormat(TimeFormat);

    //add picture to frame work
    AddPicture(pPicRate, "Raw");

    ////picture of scaler rates
    pPicScalerRate = new TGo4Picture("PicScalerRate","Picture of count-rates of "
                               " all scaler channels");

    //add rate histogram
    pPicScalerRate->AddObject((TObject*)pHist->pH2RawRate);

//    char TimeFormat[100];

    sprintf(TimeFormat, "#splitline{%%d.%%m.}{%%H:%%M}%%F1970-01-01 00:00:00s0");

    pPicScalerRate->SetXAxisTimeDisplay(kTRUE);
    pPicScalerRate->SetXAxisTimeFormat(TimeFormat);

    //add picture to frame work
    AddPicture(pPicScalerRate, "Raw");

}

//----------------------------END OF GO4 SOURCE FILE ---------------------
