
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum f�r Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#include "TnfisAnalysis.h"

#include <stdlib.h>
#include "Riostream.h"

#include "Go4EventServer.h"
#include "TGo4StepFactory.h"
#include "TGo4AnalysisStep.h"
#include "TGo4Version.h"

//***********************************************************
TnfisAnalysis::TnfisAnalysis(){}

//***********************************************************
// this constructor is called by go4analysis executable
TnfisAnalysis::TnfisAnalysis(int argc, char** argv) : TGo4Analysis(argc, argv),
      fPreLoopDone(kFALSE)
{
   cout << "**** Create TnfisAnalysis " << argv[0] << endl;

   if (!TGo4Version::CheckVersion(__GO4BUILDVERSION__)) {
      cout << "****  Go4 version mismatch" << endl;
      exit(-1);
   }

// More settings are done in macro setup.C

   gROOT->ProcessLine(".x setup.C()");
   if(argc > 1) {
     cout << "**** TnfisAnalysis Argument "<< argv[1]<<endl;
   }
   Print(); // print setup

   // Define custom passwords for analysis server
   DefineServerPasswords("nfisadmin", "nfiscontrol", "nfisview");

}

//***********************************************************
TnfisAnalysis::~TnfisAnalysis()
{
   cout << "**** TnfisAnalysis: Delete instance" << endl;
}

//-----------------------------------------------------------
Int_t TnfisAnalysis::UserPreLoop()
{
   // called after start analysis before first event
   cout << "**** TnfisAnalysis: PreLoop" << endl;
   Print(); // print setup


   // At this point all step objects have beeen created.
   // We check here that the necessary steps are enabled.
   // If not, we throw exception to not to start the event loop.

   Bool_t OK=kTRUE;
   TGo4AnalysisStep * step;

    if(fPreLoopDone)
        cout << "**** TnfisAnalysis: UserPreLoop already done!" << endl;
    else
    {
        cout << "**** TnfisAnalysis: UserPreLoop, checking step logic" << endl;

//        // Both Unpack steps need input from Unpacking:
//        step = GetAnalysisStep("Anal.FC");
//        if(step->IsProcessEnabled())
//        {
//            step=GetAnalysisStep("Unpacking");
//            if(!step->IsProcessEnabled())
//            {
//            cout << "**** !! TnfisAnalysis: Step Anal. FC needs step Unpacking"
//                 << endl;
//            OK=kFALSE;
//            }
//        }


//        step = GetAnalysisStep("Cali.FC");
//        if(step->IsProcessEnabled())
//        {
//            step=GetAnalysisStep("Anal.FC");
//            if(!step->IsProcessEnabled())
//            {
//                cout << "**** !! TnfisAnalysis: Step Cali.FC needs step"
//                     <<" Anal.FC" << endl;
//                OK=kFALSE;
//            }
//        }

    }

   Go4EventCounter = 0;

   fPreLoopDone=kTRUE;
   if(!OK) throw TGo4UserException(3,"TnfisAnalysis: Error in step setup!");

   return 0;
}
//-----------------------------------------------------------
Int_t TnfisAnalysis::UserPostLoop()
{
   // called after close analysis after last  event
   // all this is optional:
   cout << "**** TnfisAnalysis: PostLoop" << endl;

   cout << Go4EventCounter << " events processed" << endl;

   return 0;
}

//-----------------------------------------------------------
Int_t TnfisAnalysis::UserEventFunc()
{
    // all this is optional:
    // This function is called once for each event after all steps.
    Go4EventCounter++;

//    if (Go4EventCounter >= 1000)
//    {   cout << "**********************************************" << endl << endl;
//        cout << "Event \t" << Go4EventCounter << " processed" << endl << endl;
//        cout << "**********************************************" << endl << endl;
//    }

    return 0;
}
//*/
