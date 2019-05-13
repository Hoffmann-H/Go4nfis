//std. c++ includes
#include <iostream>
#include <stdlib.h>

//root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

//my includes
#include "Xsection.h"
#include "PuFC.h"
#include "UFC.h"
#include "Plot.h"

using namespace std;

int main(int argc, char **argv)
{

    TApplication theApp("App", &argc, argv);

    gROOT->Reset();

//    Xsection *pXs = new Xsection();

    Plot *p = new Plot("PuFC", "Open");
//    p->Source_E();

//    Hist *h = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root", "NIF");
//    h->SetDraw(p);
//    h->DoAnalyzeQDC();

//    PuFC *fc = new PuFC();
//    fc->SetDraw(p);
//    fc->ScatCorrSim();
//    fc->Corrections();
//    fc->CompareShadowCone("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_SB_5E6.root");

//    UFC *ufc = new UFC();
//    ufc->SetDraw(p);
//    ufc->ScatCorrSim();
//    ufc->Corrections();
//    ufc->CompareShadowCone("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/UFC_SB_1E6.root");

//    Sim *sim = new Sim("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_3E6.root", "PuFC", "Open", 0, 0);


    AnaSim *ana = new AnaSim("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_3E6.root", "PuFC", "Open", p);
//    ana->nProjectiles();
    ana->nFissions();
//    ana->Corrections();

    theApp.Run();

    return 0;
}
