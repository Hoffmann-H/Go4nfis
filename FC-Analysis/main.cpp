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

//    Hist *h = new Hist("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root", "NIF");
//    h->SetDraw(p);
//    h->DoAnalyzeQDC();

//    PuFC *fc = new PuFC(1);
//    UFC *fc = new UFC(1);

    Xsection *xs = new Xsection();
    xs->RelativeCS();

//    fc->CompareTransmission();

//    cout << "t_live(FG) = " << ufc->pHFG->t_live << endl;
//    cout << "t_live(BG) = " << ufc->pHBG->t_live << endl;
//    ufc->AnalyzeDt();
//    ufc->CrossSection();
//    fc->CompareShadowCone();
//    ufc->CompareTransmission();

//    Sim *sim = new Sim("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_5E7.root", "PuFC", "Open", 1, 1);
//    sim->Calculate();

//    Plot *p = new Plot("PuFC", "Open", 0);
//    p->Source_E();

//    AnaSim *Pu = new AnaSim("PuFC", 1, 0);
//    AnaSim *U = new AnaSim("UFC", 1, p);
//    Pu->Corrections();
//    U->Corrections();

//    p->SimF(Pu->T, Pu->DT, Pu->S, Pu->DS, Pu->F, Pu->DF, U->T, U->DT, U->S, U->DS, U->F, U->DF);

    theApp.Run();

    return 0;
}
