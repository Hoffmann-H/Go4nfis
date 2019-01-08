// This macro setup.C is called in TnfisAnalysis
void setup()
{
  TGo4AnalysisStep        * step;
  TGo4StepFactory         * fact;
  TGo4FileStoreParameter  * store;
  TGo4FileSourceParameter * input;

  Text_t autosave[]     ="nfisAutoSave";    // name of autosave file (.root)
  Text_t UnpOut[]       ="RawOutput";       // name of output file step 1  (.root)
  Text_t AnaFCOut[]     ="AnaFCOut";        // name of output file step 2  (.root)
//  Text_t CalFCOut[]     ="CalFCOut";        // name of output file step 2.1(.root)

//  Text_t node[]="rio42";                    // name of MBS or event server node
  Text_t node[]="File-Lists/nfis2014_NIF.lml";// name of MBS list mode file list
  Int_t  port=0;                            // depends on the remote event server

  go4->SetAutoSaveFile(autosave,kTRUE);   // optional, overwrite on
  go4->SetAutoSaveInterval(0);      // after n seconds , 0 = at end of event loop
  go4->SetAutoSave(kTRUE);          // optional
//  go4->SetAutoSaveFileChange(kTRUE);// create new autosave file for each file

  TGo4StepFactory*  facUnp      = new TGo4StepFactory("nfisUnp");
  TGo4AnalysisStep* stepUnp     = new TGo4AnalysisStep("Unpacking",facUnp,0,0,0);
  go4->AddAnalysisStep(stepUnp);

  TGo4StepFactory*  facAnaFC    = new TGo4StepFactory("nfisAnaFC");
  TGo4AnalysisStep* stepAnaFC   = new TGo4AnalysisStep("Anal.FC",facAnaFC,0,0,0);
  go4->AddAnalysisStep(stepAnaFC);

//  TGo4StepFactory*  facCalFC    = new TGo4StepFactory("nfisCalFC");
//  TGo4AnalysisStep* stepCalFC   = new TGo4AnalysisStep("Cali.FC",facCalFC,0,0,0);
//  go4->AddAnalysisStep(stepCalFC);

  // some processors do not use output of previous as input.
  // Therefore we must switch off step checking.
  // We must be sure that all processor get their input event properly
  // In your case means that:
  // MakeUnp is needed by FCAnalysis and PlAnalysis
  // In TnfisAnalysis::UserPreLoop() this is checked.

  go4->SetStepChecking(kFALSE);

// Analysis has the step, which has the factory, 
// which handles the event processor and event objects

  //=========================================================
  // Set up step 1 Unpack
  //=========================================================
  step = go4->GetAnalysisStep("Unpacking");
  fact=(TGo4StepFactory *)step->GetStepFactory();
  fact->DefEventProcessor("nfisMakeUnp","TnfisMakeUnp");
  fact->DefOutputEvent("nfisUnpackEvent","TnfisUnpackEvent");


  // activate one of the following MBS event sources:
  //--------------------------------------------------
  // TGo4MbsRandomParameter * source = new TGo4MbsRandomParameter("Random");
  //--------------------------------------------------
  // TGo4MbsEventServerParameter * source = new TGo4MbsEventServerParameter(node);
  //--------------------------------------------------
  // TGo4MbsStreamParameter * source = new TGo4MbsStreamParameter(node);
  //--------------------------------------------------
//   TGo4MbsTransportParameter * source = new TGo4MbsTransportParameter(node);
  //--------------------------------------------------
   TGo4MbsFileParameter * source = new TGo4MbsFileParameter(node);
  //--------------------------------------------------
  // TGo4RevServParameter * source=new TGo4RevServParameter(node);
  // source->SetPort(port);
  //--------------------------------------------------

  step->SetEventSource(source);  // register event source
  step->SetSourceEnabled(kTRUE);
  step->SetProcessEnabled(kTRUE);
  step->SetErrorStopEnabled(kTRUE);
  store = new TGo4FileStoreParameter(UnpOut);
  store->SetOverwriteMode(kTRUE);
  step->SetEventStore(store);
  step->SetStoreEnabled(kFALSE);  // en-disable output


  //=========================================================
  // Set up step 2 analysze fission chamber data
  //=========================================================
  step = go4->GetAnalysisStep("Anal.FC");
  fact=(TGo4StepFactory *)step->GetStepFactory();
  fact->DefInputEvent("nfisUnpackEvent","TnfisUnpackEvent");
  fact->DefEventProcessor("nfisFCAnalysis","TnfisFCAnalysis");
  fact->DefOutputEvent("nfisAnaFCEvent","TnfisAnaFCEvent");

  step->SetProcessEnabled(kTRUE);
  step->SetErrorStopEnabled(kTRUE);

  input = new TGo4FileSourceParameter(UnpOut);
  step->SetEventSource(input);
  step->SetSourceEnabled(kFALSE);

  store = new TGo4FileStoreParameter(AnaFCOut);
  store->SetOverwriteMode(kTRUE);
  step->SetEventStore(store);
  step->SetStoreEnabled(kFALSE);  // en-disable output

  /*/=========================================================
  // Set up step 2.1 calibrate fission chamber data
  //=========================================================
  step = go4->GetAnalysisStep("Cali.FC");
  fact=(TGo4StepFactory *)step->GetStepFactory();
  fact->DefInputEvent("nfisAnaFCEvent","TnfisAnaFCEvent");
  fact->DefEventProcessor("nfisFCCalibrate","TnfisFCCalibrate");
  fact->DefOutputEvent("nfisCaliFCEvent","TnfisCalFCEvent");

  step->SetProcessEnabled(kTRUE);
  step->SetErrorStopEnabled(kTRUE);

  input = new TGo4FileSourceParameter(AnaFCOut);
  step->SetEventSource(input);
  step->SetSourceEnabled(kFALSE);

  store = new TGo4FileStoreParameter(CalFCOut);
  store->SetOverwriteMode(kTRUE);
  step->SetEventStore(store);
  step->SetStoreEnabled(kFALSE);  // en-disable output
//*/

  printf("Setup done!\n");
}
