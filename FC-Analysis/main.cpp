/// Analysis of the Neutron-induced fission experiment done in 2014 at PTB.
/// Root ToF histograms from the Go4 analysis are used to determine 
/// (n,f) cross sections at 15MeV for Pu-242 and U-235. 
///
/// Main analysis scripts: ToF, NumberOfAtoms, NeutronField, CrossSection
/// Corrections: Target, MCNPtoROOT, AnaSim, Carlson, Correction
/// Necessary tools: FC, Runs, SaveToFile
/// optional checks: QDCmin, Deposit, nELBEsim, drawTL, PeakWidth, 
///                  ShadowCone, Stability, VglSim
/// 
/// First run the Go4 analysis. Save the results as:
///     UFC_NIF.root, UFC_SB.root, NIF.root, SB.root, SF.root
/// In SaveToFile.C, provide the path where to find the histogrammed data.
/// Run the QDCmin script.
/// Use the QDC minima for another Go4 analysis run. 
/// Run the PeakWidth script and 
/// set ToF limits for peak and background in the FC script. 
/// Run ToF, NumberOfAtoms, NeutronField
/// 
/// Author: Hans Hoffmann

//std. c++ includes
#include <iostream>
#include <stdlib.h>

//root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

//my includes
#include "root/Target.C"
#include "root/MCNPtoROOT.C"
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
    /// Use results to re-run the Go4 analysis!

    // Create plots: Peak content over peak width
    PeakWidth();
    /// Set  limits in FC.C!

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

    // Apply corrections to CS measurement,
    // calibrate areal mass densities
    Correction();

//    CrossSection();
//    Runs();

    theApp.Run();

    return 0;
}
