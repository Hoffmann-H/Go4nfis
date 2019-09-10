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
#include "Run.h"
#include "Plot.h"
#include "ToF.h"
#include "root/QDCmin.C"
#include "root/PeakWidth.C"
#include "root/NumberOfAtoms.C"

using namespace std;

//void ProvideNotebookData()
//{
//    // Fit QDC minima
//    FindPuMinima();
//    FindUMinima();

//    // Create graphs: Peak content vs Peak width, right background vs Peak width.
//    PeakWidth();

//    // Calculate eff. number of Pu atoms from spontaneous fission
//    NumberOfPuAtoms();
//    // Calculate eff. number of U atoms from H19 calibration and Carlson
//    NumberOfUAtoms();

//    // Create an instance of ToF and calculate background
//    PuFC *fc = new PuFC(1, 0);
//    UFC *ufc = new UFC(1, 0);
//}

int main(int argc, char **argv)
{

    TApplication theApp("App", &argc, argv);

    gROOT->Reset();

//    ProvideNotebookData();

//    Plot *p = new Plot("UFC", "Open", 0);
//    p->Source_E();

//    PuFC *fc = new PuFC(0, 0);
//    fc->GetNatoms();
//    fc->DrawStability();
//    fc->CrossSection();
//    fc->Corrections();

    UFC *ufc = new UFC(0, 0);
//    ufc->IsoVec();
    ufc->CrossSection();
    ufc->Corrections();

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

//    Sim *sim = new Sim("/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results/4_ene/PuFC_Open_5E7.root", "PuFC", "Open", 1, 0);
//    sim->Calculate();
//    sim->SimToF();

//    AnaSim *Pu = new AnaSim("PuFC", 1, 0);
//    AnaSim *U = new AnaSim("UFC", 1, p);
//    Pu->Corrections();
//    U->Corrections();

//    p->SimF(Pu->T, Pu->DT, Pu->S, Pu->DS, Pu->F, Pu->DF, U->T, U->DT, U->S, U->DS, U->F, U->DF);

    theApp.Run();

    return 0;
}
