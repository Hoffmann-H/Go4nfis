#include <time.h>
#include "TnfisAnalysisFC.h"
#include "TnfisUnpackEvent.h"

//***********************************************************
TnfisFCAnalysis::TnfisFCAnalysis() : TGo4EventProcessor("Proc")
{
  cout << "**** TnfisFCAnalysis: Create instance " << endl;


}
//***********************************************************
TnfisFCAnalysis::~TnfisFCAnalysis()
{
  cout << "**** TnfisFCAnalysis: Delete instance " << endl;
}
//***********************************************************
// this one is used in standard factory
TnfisFCAnalysis::TnfisFCAnalysis(const char* name) : TGo4EventProcessor(name)
{   cout << "**** TnfisFCAnalysis: Create instance " << name << endl;

    //set level of comment output
    CommentFlag = false;

    //// init user analysis objects:
    //histograms
    pHistFC = new nfisHistograms();
    pHistFC->DefineHistograms("Analysis", "FC");

    //conditions
    MakeConditions();

    //pictures
    //MakePictures();

}
//-----------------------------------------------------------
// event function
Bool_t TnfisFCAnalysis::BuildEvent(TGo4EventElement * output)
{
    if (CommentFlag == true) cout << "<I>: Analyze fission chamber data" << endl;

    //define input and output object

    //input
    TnfisUnpackEvent  *pInput = (TnfisUnpackEvent *) GetInputEvent();
//    TnfisUnpackEvent  *pInput = (TnfisUnpackEvent *) GetOutputEvent("Unpacking");
    TnfisUnpackEvent  &UnEv    = *pInput;

    //get pointers to the input event classes
    TnfisTrigAccClk *pTrigInput = dynamic_cast<TnfisTrigAccClk*>(&UnEv[0]);
    TnfisTrigAccClk *pAccInput  = dynamic_cast<TnfisTrigAccClk*>(&UnEv[1]);
    TnfisEvtCounter *pCounter   = dynamic_cast<TnfisEvtCounter*>(&UnEv[2]);
//    TnfisTimer      *pTimer     = dynamic_cast<TnfisTimer*>(&UnEv[3]);
    TnfisUnpFC      *pHZDRFC    = dynamic_cast<TnfisUnpFC*>(&UnEv[6]);

    TnfisUnpPreAmp  *pHZDRPreAmp[NumHZDRFC];
    for (int i=0; i<NumHZDRFC; i++)
        pHZDRPreAmp[i] = dynamic_cast<TnfisUnpPreAmp*>(pHZDRFC->GetPreAmp(i));

    //output
    TnfisAnaFCEvent    *pTarget    = (TnfisAnaFCEvent*) output;

    InitEventStructure(pTarget);

    //get number of trigger events in readout event
    UInt_t NumTrigEvt = pCounter->GetTotal();
    pTarget->SetNumTrigEvents(NumTrigEvt);
    Int_t time = (Int_t) pCounter->GetAbsTime().GetSec();
//    cout << time << endl;

    //loop over all trigger events in the readout event
    for (UInt_t event=0; event < pCounter->GetTotal(); event++)
    {
    // 1st filter ==> check, if FCs have got a hit /////////////////////////////
        multiplicity_struct mult;
        memset(&mult,0,sizeof(mult));

        // loop over HZDR preAmps
        for (int i_preAmp = 0; i_preAmp < NumHZDRFC; i_preAmp++)
        {   // check, if there was at least one hit in the preAmp channel
            if (pHZDRPreAmp[i_preAmp]->GetHit(event, 0) != 0)
            {   // determine the preAmps, which were hit
                mult.fchzdr.ch[mult.fchzdr.cnt] = i_preAmp;
                mult.fchzdr.cnt++;
            }
        } //end of preAmp loop

        //get number of Trigger-Hits in trigger event
        Long_t TrigHitCounter = pTrigInput->GetHitCounter(event);

    // analyze trigger /////////////////////////////////////////////////////////
        // if no trigger was found
        if (TrigHitCounter == 0)
            cerr << "<ERROR> no trigger found in TDC" << endl;

        // if more than one trigger found
        for (int i_trig=2; i_trig<TrigHitCounter; i_trig++)
            cerr << "<ERROR> " << i_trig+1
                 << " trigger found " << pTrigInput->GetHit(event, i_trig) -
                    pTrigInput->GetHit(event, 0)
                 << endl;

        // number of trigger hits //////////////////////////////////////////////////
        if (TrigHitCounter > 0)
        {   if (CommentFlag == true)
                cout << "<I>: increment HIT y="<< 10
                     << " x= " << TrigHitCounter
                     << endl;
            pHistFC->pH1AnaHitTrig->Fill(TrigHitCounter);
            pHistFC->pH2AnaHit->Fill(TrigHitCounter, 10);
        }
//*/

    // analyze data if only one trigger found //////////////////////////////////
        if (TrigHitCounter == 1)
        {   // analyze accelerator

            //Get acc hits in readout event
            Long_t AccHitCounter = pAccInput->GetHitCounter(event);

            //check, if there was at least one trigger from accelerator
            if (AccHitCounter > 0)
            {   // number of accelerator hits
                if (CommentFlag == true)
                    cout << "<I>: increment HIT y=" << 9
                         << " x=" << AccHitCounter
                         << endl;
                pHistFC->pH1AnaHitAcc->Fill(AccHitCounter);
                pHistFC->pH2AnaHit->Fill(AccHitCounter, 9);

                // search for the first accelerator hit
                pAccInput->SetFirstHit(event, pAccInput->GetHit(event,0));

                // to avoid that the following loop will try to touch unused memory
                if (AccHitCounter > 30)
                    pAccInput->SetHitCounter(event, 30);

                // search through all hits and find first hit
                for(int i_hit=1; i_hit<AccHitCounter; i_hit++)
                {   if (pAccInput->GetHit(event, i_hit) <
                        pAccInput->GetFirstHit(event))
                        pAccInput->SetFirstHit(event,
                            pAccInput->GetHit(event, i_hit));
                }
            } //if (AccHitCounter > 0)
//*/
            UInt_t channel=0;

            // analyze HZDR fission chamber
            for (int i_preAmp=0; i_preAmp<mult.fchzdr.cnt; i_preAmp++)
            {   channel = mult.fchzdr.ch[i_preAmp];

                if (CommentFlag == true)
                  cout << "<I>: increment HIT y=" << channel
                       << " x="<< pHZDRPreAmp[channel]->GetHitCounter(event)
                        << endl;
                pHistFC->pH1AnaHitHZDR[channel]->Fill(pHZDRPreAmp[channel]->GetHitCounter(event));
                pHistFC->pH2AnaHit->Fill(pHZDRPreAmp[channel]->GetHitCounter(event), channel);

                //set hit no. 0 to first hit
                Long_t FirstHit = pHZDRPreAmp[channel]->GetHit(event,0);

                if (pHZDRPreAmp[channel]->GetHitCounter(event) > 30)
                  pHZDRPreAmp[channel]->SetHitCounter(event, 30);

                // search through all hits and find first hit
                for(int i_hit=1; i_hit<pHZDRPreAmp[channel]->GetHitCounter(event); i_hit++)
                {   if (pHZDRPreAmp[channel]->GetHit(event, i_hit) < FirstHit)
                        FirstHit = pHZDRPreAmp[channel]->GetHit(event, i_hit);
                }

                //set the earliest non-zero hit to first hit
                pHZDRPreAmp[channel]->SetFirstHit(event, FirstHit);

            } //for (int i_preAmp=0;i_preAmp<mult.fchzdr.cnt;i_preAmp++)
//*/

    // process coincidence if Reference-hit found //////////////////////////////
            // there has to be an accelerator hit
            if (pAccInput->GetFirstHit(event) != 0)
            {   // calculate time difference to reference for each channel

                //HZDR FC
                for (int i_preAmp=0; i_preAmp<mult.fchzdr.cnt; i_preAmp++)
                {   channel = mult.fchzdr.ch[i_preAmp];

                    if (pHZDRPreAmp[channel]->GetFirstHit(event) != 0)
                    {   Double_t TimeDiff = ( pHZDRPreAmp[channel]->GetFirstHit(event)-
                                            pAccInput->GetFirstHit(event) - 62000) / 40.96
                                            + nfisHistograms::tof_bias(channel); // bias according to simulation

                        Long_t QDCl = pHZDRPreAmp[channel]->GetQDCl(event);
                        Long_t QDCh = pHZDRPreAmp[channel]->GetQDCh(event);

                        pHistFC->pH1AnaQDCl_trig[channel]->Fill(QDCl);

                        //fill qdc and time diff to the output event
                        TnfisAnaPreAmp *pPreAmpOut =
                                dynamic_cast<TnfisAnaPreAmp*>(
                                        pHZDRFCOut->GetPreAmp(channel));
                            pPreAmpOut->SetTime(TimeDiff);
                            pPreAmpOut->SetQDCl(QDCl);
                            pPreAmpOut->SetQDCh(QDCh);

                        //fill time diff into corresponding histo
                        pHistFC->pH1AnaDtHZDR[channel]->Fill(TimeDiff);

                        //fill time diff and channel into 2D histo
                        pHistFC->pH2AnaDt->Fill(TimeDiff, channel);

                        //fill in time difference with qdc gate on fission frag
                        if (pConQDC[channel]->Test(QDCl))  // Test if the measured energy fits better as alpha or as fission fragment
                        {   // Fission fragments
                            if (CommentFlag)
                                cout << "ch " << channel << ": QDC pulse height " << QDCl << " within gate. Filling pH2AnaDt_g." << endl;
                            pHistFC->pH1AnaDtHZDR_g[channel]->Fill(TimeDiff);
                            pHistFC->pH2DtTime[channel]->Fill(TimeDiff, (Int_t)time);
                            //pHistFC->pH2AnaDt_g->Fill(TimeDiff, channel);
                            pHistFC->pH1AnaQDCl[channel]->Fill(QDCl);

                            //FF-QDC spectra with gate on ToF
//                            if (pConToF[channel]->Test(TimeDiff))
//                              //n-induced fission
//                                pHistFC->pH1AnaQDCl_NIF[channel]->Fill(QDCl);
//                            else
//                              //spontaneous fission
//                                pHistFC->pH1AnaQDCl_SF[channel]->Fill(QDCl);
                        }
//                        else
//                        { // pulse-height refused
//                            if (CommentFlag)
//                                cout << "ch " << channel << ": QDC pulse height " << QDCl << " refused." << endl; //*/
//                            pHistFC->pH1AnaDtHZDR_r[channel]->Fill(TimeDiff);
//                        } // if (pConQDC[channel]->Test(QDCl))
                    } //if (pPreAmp[channel]->GetFirstHit(event) != 0)
//                    else
//                        cout << "No Fission Chamber hit!" << endl;

                } // for (int i_preAmp=0; i_preAmp<mult.fchzdr.cnt; i_preAmp++)

            } // if (pAccInput->GetFirstHit(event) != 0)
//            else
//                cout << "No Accalerator hit!" << endl;

        } //if (HitCounter == 1)

    } //for (int event=0; event<MAXHITS-1; event++)


    if (CommentFlag == true)
        cout << "<I>: Fission chamber data analyzed" << endl;

    //cout << n_NIF << " " << n_SF << endl;

    return kTRUE;
}

//-----------------------------------------------------------
// histogramming function
Bool_t TnfisFCAnalysis::FillHistograms()
{

  return kTRUE;
}

//______________________________________________________________________________
void TnfisFCAnalysis::InitEventStructure(TnfisAnaFCEvent *pTarget)
{   //clears all items of the target event structure
    TnfisAnaFCEvent    &FCEv       = *pTarget;

    pHZDRFCOut         = dynamic_cast<TnfisAnaFC*>(&FCEv[1]);
    for (int i_preAmp = 0; i_preAmp < NumHZDRFC; i_preAmp++)
        pHZDRPreAmpOut[i_preAmp] = dynamic_cast<TnfisAnaPreAmp*>(
                    pHZDRFCOut->GetPreAmp(i_preAmp));

    ResetEventStructure("init");

}

//______________________________________________________________________________
void TnfisFCAnalysis::ResetEventStructure(Option_t *t)
{
    //clear the event structures
    pHZDRFCOut->Clear(t);
    for (int i_preAmp = 0; i_preAmp < NumHZDRFC; i_preAmp++)
        pHZDRPreAmpOut[i_preAmp]->Clear(t);
}

Double_t TnfisFCAnalysis::ChToNanosec(Double_t channel)
{
    return (channel - 62000) / 40.96;
}

//______________________________________________________________________________
void TnfisFCAnalysis::MakeConditions()
{
    //dummies to form right object names
    char CondNameDummy[50] ="";
    char HistNameDummy[20] ="";

    // define qdc and ToF conditions
    // ToF conditions estimated per Gaussian fit with constant background applied to H1AnaHZDRDtG[Channel]
    // unit: QDC channels / TDC channels

        // Channel              1        2        3        4        5        6        7        8
//        Double_t qdc_min[]   = {899.24,  853.668, 895.652, 849.393, 1046.41, 891.396, 906.123, 837.486}; // PuFC
//        Double_t tof_mean[]  = {67725,   67252,   67156,   67311,   67307,   67267,   67221,   67222}; // PuFC
//        Double_t tof_width[] = {115,     94,      90,      84,      90,      96,      86,      94}; // PuFC
/// Auch in nfisHistograms.cxx, ll 341-344 auswählen!
        // Channel              1        2        3        4        5        6        7        8
        Double_t qdc_min[]   = {220.396, 500.964, 219.016, 214.101, 199.19,  171.929, 169.046, 121.56}; // UFC
        Double_t tof_mean[]  = {73202,   73016,   72674,   72864,   72863,   72828,   72793,   72776}; // UFC
        Double_t tof_width[] = {134,     119,     117,     118,     118,     116,     115,     116}; // UFC
        
    Double_t qdc_max[]      = {4096, 4096, 4096, 4096, 4096, 4096, 4096, 4096};

    //make qdc conditions for each pre-Amp
    for (int i_PMT=0; i_PMT<NumHZDRFC; i_PMT++)
    {   sprintf(CondNameDummy, "Analysis/QDC/QDCl_Thr_Ch%i", i_PMT+1);
        sprintf(HistNameDummy, "H1RawQDCl_%i", i_PMT+1);

        pConQDC[i_PMT] = MakeWinCond(CondNameDummy,
                                     qdc_min[i_PMT], qdc_max[i_PMT],
                                     HistNameDummy);
        //set initial thresholds
        pConQDC[i_PMT]->ResetCounts();
        pConQDC[i_PMT]->Enable();
    }

    //define time window for underground estimation
    Double_t Dt_min = 63600;
    Double_t Dt_max = 78500;

    //make qdc conditions for each pre-Amp
    for (int i_PMT=0; i_PMT<NumHZDRFC; i_PMT++)
    {   sprintf(CondNameDummy, "Analysis/ToF/ToF_Thr_Ch%i", i_PMT+1);
        sprintf(HistNameDummy, "H1AnaHZDRDtG_%i", i_PMT+1);

        pConToF[i_PMT] = MakeWinCond(CondNameDummy,
                                     ChToNanosec(tof_mean[i_PMT] - 3 * tof_width[i_PMT]), ChToNanosec(tof_mean[i_PMT] + 3 * tof_width[i_PMT]),
                                     HistNameDummy);

        sprintf(CondNameDummy, "Analysis/ToF/ToF_Dt_Ch%i", i_PMT+1);
        pConDt[i_PMT] = MakeWinCond(CondNameDummy, ChToNanosec(Dt_min), ChToNanosec(Dt_max), HistNameDummy);

        //set initial thresholds
        pConToF[i_PMT]->ResetCounts();
        pConToF[i_PMT]->Enable();

        pConDt[i_PMT]->ResetCounts();
        pConDt[i_PMT]->Enable();
    }


}

//______________________________________________________________________________
void TnfisFCAnalysis::MakePictures()
{
    //pictures
    pPicToF   = new TGo4Picture("AnaPicToF","Picture of uncalibrated time-of-flight"
                                "spectra of all HZDR FC PreAmps and H19 FC");
//    pPicToF_g = new TGo4Picture("AnaPicToF_g","Picture of uncalibrated time-of-flight"
//                                "spectra of all HZDR FC PreAmps and H19 FC with"
//                                "QDC gate");

    //get time-of-flight single spectra and add them to the picture
    for (int i_preAmp=0; i_preAmp<NumHZDRFC; i_preAmp++)
    {   pPicToF->AddH1(pHistFC->pH1AnaDtHZDR[i_preAmp]);
//        pPicToF_g->AddH1(pHistFC->pH1AnaDtHZDR_g[i_preAmp]);
    }

    //set log scale to pictures
    pPicToF->SetLogScale(1); //pPicToF_g->SetLogScale(1);
    //add picture to frame work
    AddPicture(pPicToF, "Analysis"); //AddPicture(pPicToF_g, "Analysis");

}

////////////////////////////////////////////////////////////////////////////////
//                      Fission Chamber Event Class                           //
////////////////////////////////////////////////////////////////////////////////

TnfisAnaFCEvent::TnfisAnaFCEvent() : TGo4CompositeEvent()
{ // void constructor neccessary for the TTree lib
}

TnfisAnaFCEvent::TnfisAnaFCEvent(const char* name) : TGo4CompositeEvent()
{
  cout << "**** TnfisAnaFCEvent: Create instance " << name << endl;

  Clear();

    addEventElement(new TnfisAnaFC("HZDR_FC",1));

}

TnfisAnaFCEvent::~TnfisAnaFCEvent()
{
  cout << "**** TnfisAnaFCEvent: Delete instance " << endl;
}

//______________________________________________________________________________
void TnfisAnaFCEvent::Clear(Option_t *t)
{
  // all members should be cleared.
  NumTrigEvents = 0;
}

////////////////////////////////////////////////////////////////////////////////

TnfisAnaFC::TnfisAnaFC():TGo4CompositeEvent()
{ // void constructor neccessary for the TTree lib
}

TnfisAnaFC::TnfisAnaFC(const char* name, Short_t id) :
    TGo4CompositeEvent(Form("%s_Common", name),name,id)
{
  cout << "**** TnifsAnaFC: Create instance " << name << endl;

  TString chname;

  //loop over all PreAmps in FC
  for(Int_t i = 0; i < NumHZDRFC; i++)
  {   chname.Form("%s_PreAmpCh%i",name, i);
      PreAmp[i] = new TnfisAnaPreAmp(chname.Data(), i);
      addEventElement(PreAmp[i]);
  }
  Clear();
}
//______________________________________________________________________________
TnfisAnaFC::~TnfisAnaFC()
{
  cout << "**** TnifsAnaFC: Delete instance " << endl;
}

//______________________________________________________________________________
void  TnfisAnaFC::Clear(Option_t *t)
{
  // all members should be cleared.
}

////////////////////////////////////////////////////////////////////////////////

TnfisAnaPreAmp::TnfisAnaPreAmp():TGo4EventElement()
{ // void constructor neccessary for the TTree lib
}

TnfisAnaPreAmp::TnfisAnaPreAmp(const char* name, Short_t id) :
    TGo4EventElement(name, name, id)
{
    cout << "\t**** TnfisAnaPreAmp: Create instance " << name << endl;
    //init arrays with 0
    pVecQDCh = new vector<Long_t>();    //pVecQDCh->resize(MAXHITS);
    pVecQDCl = new vector<Long_t>();    //pVecQDCl->resize(MAXHITS);
    pVecTime = new vector<Long_t>();    //pVecTime->resize(MAXHITS);

    Clear();
}

//______________________________________________________________________________
TnfisAnaPreAmp::~TnfisAnaPreAmp()
{
    cout << "\t**** TnfisAnaPreAmp:" << this->GetName()
         << " Delete instance " << endl;
//    Clear();
}

//______________________________________________________________________________
void  TnfisAnaPreAmp::Clear(Option_t *t)
{
    // all members should be cleared.
//      memset(ADC,0, sizeof(ADC));
//      memset(QDCh,0, sizeof(QDCh));
//      memset(QDCl,0, sizeof(QDCl));
//      memset(Time,TDCCUTOFF, sizeof(Time));

    pVecQDCh->clear();
    pVecQDCl->clear();
    pVecTime->clear();

}

////////////////////////////////////////////////////////////////////////////////


//----------------------------END OF GO4 SOURCE FILE ---------------------
//*/
