//std. c++ includes
#include <iostream>
#include <stdlib.h>

//root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

//my includes
//#include "Xsection.h"
//#include "PuFC.h"
//#include "UFC.h"
//#include "Run.h"
//#include "Plot.h"
#include "root/Target.C"
#include "root/AnaSim.C"
#include "root/FC.C"
#include "root/QDCmin.C"
#include "root/PeakWidth.C"
#include "root/NeutronField.C"
#include "root/ToF.C"
#include "root/NumberOfAtoms.C"
#include "root/CrossSection.C"
#include "root/Carlson.C"
#include "root/Correction.C"
#include "root/Runs.C"

using namespace std;


int main(int argc, char **argv)
{

    TApplication theApp("App", &argc, argv);

    gROOT->Reset();

    // Determine QDC threshold
    QDCmin();

    // Create plots: Peak content over peak width
    PeakWidth();

    // Get source spectrum and Ti(T) scattering correction
    Target();

    // Convert simulation data
    MCNPtoROOT();

    // Calculate neutron scattering correction factors
    AnaSim();

    // Count (n,f) events in time difference spectrums
    ToF();

    // Calculate neutron field intensity
    NeutronField();

    // Calculate UFC inefficiency correction
    Carlson();

    // Get UFC@nELBE and PuFC SF calibration
    NumberOfAtoms();

    // Calculate cross section
    CrossSection();

    // Apply corrections to CS measurement,
    // calibrate areal mass densities
    Correction();

//    Runs();

    theApp.Run();

    return 0;
}
