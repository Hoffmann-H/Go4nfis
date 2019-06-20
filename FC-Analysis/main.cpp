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

//    Plot *p = new Plot("UFC", "Open", 0);
//    p->Source_E();

//    Hist *h = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root", "NIF", "UFC");
//    h->SetDraw(p);
//    h->DoAnalyzeQDC();

//    PuFC *fc = new PuFC(0);
    UFC *ufc = new UFC(1);
    ufc->Calibrate();

//    Xsection *xs = new Xsection();
//    xs->RelativeCS();

//    fc->CompareTransmission();

//    cout << "t_live(FG) = " << ufc->pHFG->t_live << endl;
//    cout << "t_live(BG) = " << ufc->pHBG->t_live << endl;
//    fc->AnalyzeDt();
//    fc->CrossSection();
//    ufc->CrossSection();
//    fc->Corrections();
//    ufc->Corrections();
//    fc->CompareShadowCone();
//    ufc->CompareTransmission();

//    Sim *sim = new Sim("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_5E7.root", "PuFC", "Open", 1, 1);
//    sim->Calculate();

//    AnaSim *Pu = new AnaSim("PuFC", 1, 0);
//    AnaSim *U = new AnaSim("UFC", 1, p);
//    Pu->Corrections();
//    U->Corrections();

//    p->SimF(Pu->T, Pu->DT, Pu->S, Pu->DS, Pu->F, Pu->DF, U->T, U->DT, U->S, U->DS, U->F, U->DF);

    theApp.Run();

    return 0;
}
