#ifndef CORRECTION_H
#define CORRECTION_H
#include "Runs.C"
#include <fstream>
#include "NeutronField.C"
#include "SaveToFile.C"

void IsoVec(string FC, Double_t &I, Double_t &DI)
{
    if (strcmp(FC.c_str(), "UFC"))
    { // PuFC
        I = 1;
        DI = 0;
    } else {
        // TODO: Calculate Isotope vector
        I = 1.0398;
        DI = 0;
    }
    cout << FC << " Isotope correction: " << I << " +- " << DI << endl;
}

void DeltaIso()
{
    Int_t nPu = 6;
    Double_t isoPu[] = {0.00002, 0.00005, 0.00022, 0.00002, 0.99967, 0.00002};
    Double_t DisoPu[] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001};
    Double_t evalPu[] = {2.74344, 2.415756, 2.3578, 2.09, 2.083152, 1.829932};
    Int_t nU = 5;
    Double_t isoU[] = {0.0000011, 0.00459, 0.904, 0.00401, 0.0912};
    Double_t DisoU[] = {0.0000005, 0.00005, 0.005, 0.00005, 0.0006};
    Double_t evalU[] = {2.09894, 2.272568, 2.1256, 1.750379, 1.2415};

    Double_t DPu1 = sqrt( pow(DisoPu[0] * evalPu[0], 2) +
                          pow(DisoPu[1] * evalPu[1], 2) +
                          pow(DisoPu[2] * evalPu[2], 2) +
                          pow(DisoPu[3] * evalPu[3], 2) +
                          pow(DisoPu[4] * evalPu[4], 2) +
                          pow(DisoPu[5] * evalPu[5], 2)
                        ) / (
                          isoPu[0] * evalPu[0] +
                          isoPu[1] * evalPu[1] +
                          isoPu[2] * evalPu[2] +
                          isoPu[3] * evalPu[3] +
                          isoPu[4] * evalPu[4] +
                          isoPu[5] * evalPu[5] );
    Double_t Pu = evalPu[4] / (
                          isoPu[0] * evalPu[0] +
                          isoPu[1] * evalPu[1] +
                          isoPu[2] * evalPu[2] +
                          isoPu[3] * evalPu[3] +
                          isoPu[4] * evalPu[4] +
                          isoPu[5] * evalPu[5] );
    cout << Pu << endl << DPu1 << endl;
    Double_t DU1 = sqrt( pow(DisoU[0] * evalU[0], 2) +
                          pow(DisoU[1] * evalU[1], 2) +
                          pow(DisoU[2] * evalU[2], 2) +
                          pow(DisoU[3] * evalU[3], 2) +
                          pow(DisoU[4] * evalU[4], 2)
                        ) / (
                          isoU[0] * evalU[0] +
                          isoU[1] * evalU[1] +
                          isoU[2] * evalU[2] +
                          isoU[3] * evalU[3] +
                          isoU[4] * evalU[4] );
    cout << endl << DU1 << endl;
}

void ScatteringUFC()
{ // Experimental scattering correction
    char name[64] = "";
    Double_t MonitorFG = 72661227.5811;
    Double_t DMonitorFG = sqrt(MonitorFG);
    Double_t MonitorBG = 28042758.0879;
    Double_t DMonitorBG = sqrt(MonitorBG);
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *geFissionFG = (TGraphErrors*)fAna->Get("UFC/ToF/Signal/UFC_NIF/InducedFission");
    if (!geFissionFG)
        cout << "Could not open " << "UFC/ToF/Signal/UFC_NIF/InducedFission" << endl;
    TGraphErrors *geFissionBG = (TGraphErrors*)fAna->Get("UFC/ToF/Signal/UFC_SB/InducedFission");
    if (!geFissionBG)
        cout << "Could not open " << "UFC/ToF/Signal/UFC_SB/InducedFission" << endl;
    TGraphErrors *geExpS = new TGraphErrors(8);

    TGraphErrors *geSimSB = (TGraphErrors*)fAna->Get("Simulation/Geant4/UFC_SB/Correction/Total");
    if (!geSimSB)
        cout << "Could not open " << "Simulation/Geant4/UFC_SB/Correction/Total" << endl;
    TGraphErrors *geSimOpen = (TGraphErrors*)fAna->Get("Simulation/Geant4/UFC_Open/Correction/Total");
    if (!geSimOpen)
        cout << "Could not open " << "Simulation/Geant4/UFC_Open/Correction/Total" << endl;
    TGraphErrors *geSimS = new TGraphErrors(8);
    for (Int_t i = 0; i < 8; i++)
    {
        // Experimental scattering
		Double_t x, C_FG, DC_FG, C_BG, DC_BG;
        geFissionFG->GetPoint(i, x, C_FG);
        DC_FG = geFissionFG->GetErrorY(i);
        geFissionBG->GetPoint(i, x, C_BG);
        DC_BG = geFissionBG->GetErrorY(i);
        Double_t S = 1 - C_BG / MonitorFG / C_FG * MonitorFG;
        Double_t DS = (1 - S) * sqrt( pow(DC_BG / C_BG, 2) +
									  pow(DMonitorFG / MonitorFG, 2) +
									  pow(DC_FG / C_FG, 2) +
									  pow(DMonitorBG / MonitorBG, 2) );
		geExpS->SetPoint(i, i+1, S);
		geExpS->SetPointError(i, 0, DS);

        // Simulated Shadow Bar correction
//        Double_t effNtotSB, DeffNtotSB, effNtotOpen, DeffNtotOpen, SimS, DSimS;
//        geSimSB->GetPoint(i, x, effNtotSB);
//        DeffNtotSB = geSimSB->GetErrorY(i);
//        geSimOpen->GetPoint(i, x, effNtotOpen);
//        DeffNtotOpen = geSimOpen->GetErrorY(i);
//        SimS = 1 - effNtotSB / effNtotOpen;
//        DSimS = effNtotSB / effNtotOpen * sqrt( pow(DeffNtotSB / effNtotSB, 2) +
//                                                pow(DeffNtotOpen / effNtotOpen, 2) );
//        geSimS->SetPoint(i, i+1, SimS);
//        geSimS->SetPointError(i, 0, DSimS);
    }
    geExpS->SetName("UFC_Exp_S");
    Save(fAna, "UFC/Correction", geExpS);
//    geSimS->SetName("UFC_Geant4_S");
//    Save(fAna, "UFC/Correction", geSimS);
    fAna->Save();
    fAna->Close();
}

/*void TransmissionPuFC()
{
	char name[64] = "";
	TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
	
	for (Int_t j = 0; j < nRunsFG; j++)
	{
		string Run = GetRunName("PuFC_FG", j);
		sprintf(name, "PuFC/ToF/Signal/%s/", Run.c_str());

	}
}//*/

void ApplyCorrectionsU(string Simulation = "Geant4", Int_t Gate = 15)
{
    Double_t evalCS = 2.0730695;
    Double_t eff = 0.946152;
//    Double_t Deposit_area[] = {42.7736, 41.5982, 43.3861, 43.4792, 43.2923, 43.1283, 42.7877, 43.6696}; // cm^2
    Double_t Deposit_area[] = {43, 43, 43, 43, 43, 43, 43, 43};
    Double_t DeltaDepositArea = 0.7 / sqrt(8); // cm^2
    Double_t MolarMass = 235.3175644086;
    Double_t u = 1.660539E-24; // [g]

    char name[64] = "";
    TFile* fAna = TFile::Open(results_file, "UPDATE");

    // Get fission signal rate
    sprintf(name, "UFC/ToF/Signal/UFC_NIF/FissionRate");
    TGraphErrors *geFissionRate = (TGraphErrors*)fAna->Get(name);
    if (!geFissionRate) cout << "Could not get " << name << endl;

    // Get effective areal mass density calibrated at nELBE
    TGraphErrors *geELBEmA = (TGraphErrors*)fAna->Get("UFC/nAtoms/UFC_emA");
    if (!geELBEmA) cout << "Could not open " << "UFC/nAtoms/UFC_emA" << endl;

    // Get neutron density rate
    sprintf(name, "UFC/NeutronField/UFC_NIF/UFC_NIF_Density");
    TGraphErrors *geFlux = (TGraphErrors*)fAna->Get(name);
    if (!geFlux) cout << "Could not get " << name << endl;

    // Get TARGET ToF-gate correction factor
    sprintf(name, "UFC/Correction/UFC_Target_Gate"); // for vacuum simulation with total source spectrum
//    sprintf(name, "UFC/Correction/UFC_Target"); // for vacuum simulation with direct source spectrum only
    TGraphErrors *gT = (TGraphErrors*)fAna->Get(name); if (!gT) cout << "Could not get " << name << endl;
//    sprintf(name, "UFC/Correction/UFC_%s_Gate_%ins", Simulation.c_str(), Gate);
//    TGraphErrors *gG = (TGraphErrors*)fAna->Get(name); if (!gG) cout << "Could not get " << name << endl;

    // Get simulated correction for neutron scattering
    sprintf(name, "UFC/Correction/UFC_%s_real_C", Simulation.c_str());
    TGraphErrors *gC = (TGraphErrors*)fAna->Get(name); if (!gC) cout << "Could not get " << name << endl;

    // Get Inefficiency correction ("Carlson-Faktor")
    TGraphErrors *gCarlson = (TGraphErrors*)fAna->Get("Carlson/Carlson_Correction");
    if (!gCarlson) cout << "Could not get " << "Carlson/Carlson_Correction" << endl;

    // Create corrected cross section graph
    TGraphErrors *gCorr = new TGraphErrors(8);
    sprintf(name, "UFC_CS_%s_corrected", Simulation.c_str());
    gCorr->SetName(name);

    // Create atom number calibration graph
    TGraphErrors *gN = new TGraphErrors(8);
    sprintf(name, "UFC_eN_cal_%s", Simulation.c_str());
    gN->SetName(name);
    gN->SetTitle("U Calibration; Deposit; #varepsilon_{#it{n}ELBE}#it{N}_{U}");

    // Create areal mass density calibration graph
    TGraphErrors *gmA = new TGraphErrors(8);
    sprintf(name, "UFC_emA_cal_%s", Simulation.c_str());
    gmA->SetName(name);
    gmA->SetTitle("U areal mass density; Deposit; #varepsilon_{#it{n}ELBE}#it{m}_{A} / mg cm^{-2}");

    // Raw cross section parameters
    Double_t FissionRate[8], DFissionRate[8];
    Double_t nELBEmA[8], DnELBEmA[8];
    Double_t Flux[8], DFlux[8];

    // Correction factors
    Double_t I, DI;
    IsoVec("UFC", I, DI);
    Double_t T[8], DT[8];
    Double_t C[8], DC[8];
    Double_t Carlson[8], DCarlson[8];

    // Calibration and cross section
    Double_t N[8], DN[8];
    Double_t mA[8], DmA[8];
    Double_t Sigma[8], DSigma[8];

    Double_t sumSigma = 0;
    Double_t D2sumSigma = 0;
    Double_t sumAtoms = 0;
    ofstream output("../results/Correction_UFC.txt");
    cout << "ch   eff.m_A(nELBE)/(mg*cm^-2)   N(cal)   Sigma/b" << endl;
    output << "# deposit cross-section" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get data
        Double_t x;
        geFissionRate->GetPoint(i, x, FissionRate[i]);
        DFissionRate[i] = geFissionRate->GetErrorY(i);
        geELBEmA->GetPoint(i, x, nELBEmA[i]);
        DnELBEmA[i] = geELBEmA->GetErrorY(i);
        geFlux->GetPoint(i, x, Flux[i]);
        DFlux[i] = geFlux->GetErrorY(i);
        gT->GetPoint(i, x, T[i]);
        DT[i] = gT->GetErrorY(i);
        gC->GetPoint(i, x, C[i]);
        DC[i] = gC->GetErrorY(i);
        gCarlson->GetPoint(i, x, Carlson[i]);
        DCarlson[i] = gCarlson->GetErrorY(i);

        /// Calculate cross section using nELBE calibration ///////////////////////////////////////
        Sigma[i] = I * T[i] * C[i] * Carlson[i] * FissionRate[i] * MolarMass * u / Deposit_area[i] / nELBEmA[i] / Flux[i] * 1.E22 * 1000;
        ///////////////////////////////////////////////////////////////////////////////////////////

        DSigma[i] = Sigma[i] * sqrt( //pow(DI / I, 2) + // Isotope: systematic
                                     pow(DC[i] / C[i], 2) +
                                     pow(DFissionRate[i] / FissionRate[i], 2) +
                                     pow(DnELBEmA[i] / nELBEmA[i], 2) +
                                     pow(DFlux[i] / Flux[i], 2) );
        output << i+1 << " " << Sigma[i] << " " << DSigma[i] << endl;
        gCorr->SetPoint(i, i+1, Sigma[i]);
        gCorr->SetPointError(i, 0, DSigma[i]);

        /// Calibrate atom number using evaluated cross section ///////////////////////////////////
        N[i] = I * T[i] * C[i] * Carlson[i] * FissionRate[i] / evalCS / Flux[i] * 1.E22;
        ///////////////////////////////////////////////////////////////////////////////////////////

        DN[i] = N[i] * sqrt( pow(DC[i] / C[i], 2) +
                             pow(DFissionRate[i] / FissionRate[i], 2) +
                             pow(DFlux[i] / Flux[i], 2) );
        gN->SetPoint(i, i+1, N[i]);
        gN->SetPointError(i, 0, DN[i]);
        cout << " " << i+1 << " " << nELBEmA[i] << "+-" << DnELBEmA[i] << "  " << N[i] << "+-" << DN[i] << "  " << Sigma[i] << "+-" << DSigma[i] << endl;

        /// Calibrate areal mass density using evaluated cross section ////////////////////////////
        mA[i] = I * T[i] * C[i] * Carlson[i] * FissionRate[i] * MolarMass * u / Deposit_area[i] / evalCS / Flux[i] * 1.E22 * 1000; // [mm^2 -> b, g -> mg]
        ///////////////////////////////////////////////////////////////////////////////////////////

        DmA[i] = mA[i] * sqrt( pow(DC[i] / C[i], 2) +
                               pow(DFissionRate[i] / FissionRate[i], 2) +
                               pow(DFlux[i] / Flux[i], 2) +
                               pow(DeltaDepositArea / Deposit_area[i], 2) );
        gmA->SetPoint(i, i+1, mA[i]);
        gmA->SetPointError(i, 0, DmA[i]);

        // Make average
        sumSigma += Sigma[i] * N[i];
        D2sumSigma += pow(DSigma[i] * N[i], 2);
        sumAtoms += N[i];
    }
    Save(fAna, "UFC/Correction", gCorr);
    Save(fAna, "UFC/nAtoms", gN);
    Save(fAna, "UFC/nAtoms", gmA);
    fAna->Save();
    fAna->Close();
    Double_t avSigma = sumSigma / sumAtoms;
    Double_t DavSigma = sqrt(D2sumSigma) / sumAtoms;
    cout << "Weighted average: " << avSigma << " $\\pm$ " << DavSigma << endl;
    Double_t Variance = 0;
    for (uint i = 0; i < 8; i++)
        Variance += pow(Sigma[i] - avSigma, 2) / 8.;
    cout << "Standard deviation: " << sqrt(Variance) << endl;
    output << "# average " << avSigma << " " << DavSigma;
    output.close();
}

void ApplyCorrectionsPu(string Simulation = "Geant4")
{
    Double_t evalCS = 2.1026472;
    Double_t eff = 0.988;
    Double_t Deposit_area[] = {43, 43, 43, 43, 43, 43, 43, 43};
    Double_t DeltaDepositArea = 0.7 / sqrt(8); // cm^2
    Double_t MolarMass = 242.058089519385;
    Double_t u = 1.660539E-24; // [g]

    char name[64] = "";
    TFile* fAna = TFile::Open(results_file, "UPDATE");

    // Get fission signal rate
    sprintf(name, "PuFC/ToF/Signal/NIF/FissionRate");
    TGraphErrors *geFissionRate = (TGraphErrors*)fAna->Get(name);
    if (!geFissionRate) cout << "Could not get " << name << endl;

    // Get effective atom number calibrated with SF rate
    TGraphErrors *geNsponFis = (TGraphErrors*)fAna->Get("PuFC/nAtoms/PuFC_eN");
    if (!geNsponFis) cout << "Could not open " << "PuFC/nAtoms/PuFC_eN" << endl;

    // Get neutron density rate
    sprintf(name, "PuFC/NeutronField/NIF/NIF_Density");
    TGraphErrors *geFlux = (TGraphErrors*)fAna->Get(name);
    if (!geFlux) cout << "Could not get " << name << endl;

    // Get TARGET ToF-gate correction factor
//    sprintf(name, "PuFC/Correction/PuFC_Target_Gate"); // for vacuum simulation with total source spectrum
    sprintf(name, "PuFC/Correction/PuFC_Target"); // for vacuum simulation with direct source spectrum only
    TGraphErrors *gT = (TGraphErrors*)fAna->Get(name); if (!gT) cout << "Could not get " << name << endl;
//    sprintf(name, "PuFC/Correction/PuFC_%s_Gate_%ins", Simulation.c_str(), Right("PuFC"));
//    TGraphErrors *gG = (TGraphErrors*)fAna->Get(name); if (!gG) cout << "Could not get " << name << endl;

    // Get simulated correction for neutron scattering
    sprintf(name, "PuFC/Correction/PuFC_%s_real_C", Simulation.c_str());
    TGraphErrors *gC = (TGraphErrors*)fAna->Get(name); if (!gC) cout << "Could not open " << name << endl;

    /// note: For the PuFC, SF and (n,f) detection efficiencies are assumed equal.
    /// Carlson correction k_epsilon = 1

    // Create corrected cross section graph
    TGraphErrors *gCorr = new TGraphErrors(8);
    sprintf(name, "PuFC_CS_%s_corrected", Simulation.c_str());
    gCorr->SetName(name);

    // Create atom number calibration graph
    TGraphErrors *gN = new TGraphErrors(8);
    sprintf(name, "PuFC_eN_cal_%s", Simulation.c_str());
    gN->SetName(name);
    gN->SetTitle("Pu Calibration; Deposit; #varepsilon#it{N}_{Pu}");

    // Create areal mass density calibration graph
    TGraphErrors *gmA = new TGraphErrors(8);
    sprintf(name, "PuFC_emA_cal_%s", Simulation.c_str());
    gmA->SetName(name);
    gmA->SetTitle("Pu areal mass density; Deposit; #varepsilon#it{m}_{A} / mg cm^{-2}");

    // Raw cross section parameters
    Double_t FissionRate[8], DFissionRate[8];
    Double_t NsponFis[8], DNsponFis[8];
    Double_t Flux[8], DFlux[8];

    // Correction factors
    Double_t I, DI;
    IsoVec("PuFC", I, DI);
    Double_t T[8], DT[8];
    Double_t C[8], DC[8];

    // Calibration and cross section
    Double_t N[8], DN[8];
    Double_t mA[8], DmA[8];
    Double_t Sigma[8], DSigma[8];

    Double_t sumSigma = 0;
    Double_t D2sumSigma = 0;
    Double_t sumAtoms = 0;
    ofstream output("../results/Correction_PuFC.txt");
    cout << "ch   N(SF)   N(cal)   Sigma/b" << endl;
    output << "# deposit corrected" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get data
        Double_t x;
        geFissionRate->GetPoint(i, x, FissionRate[i]);
        DFissionRate[i] = geFissionRate->GetErrorY(i);
        geNsponFis->GetPoint(i, x, NsponFis[i]);
        DNsponFis[i] = geNsponFis->GetErrorY(i);
        geFlux->GetPoint(i, x, Flux[i]);
        DFlux[i] = geFlux->GetErrorY(i);
        gT->GetPoint(i, x, T[i]);
        DT[i] = gT->GetErrorY(i);
        gC->GetPoint(i, x, C[i]);
        DC[i] = gC->GetErrorY(i);

        /// Calculate cross section using spontaneous fission rate ////////////////////////////////
        Sigma[i] = T[i] * C[i] * FissionRate[i] / NsponFis[i] / Flux[i] * 1.E22;
        ///////////////////////////////////////////////////////////////////////////////////////////

        DSigma[i] = Sigma[i] * sqrt( //pow(DI / I, 2) + // Isotope: systematic
                                     pow(DC[i] / C[i], 2) +
                                     pow(DFissionRate[i] / FissionRate[i], 2) +
                                     pow(DNsponFis[i] / NsponFis[i], 2) +
                                     pow(DFlux[i] / Flux[i], 2) );
        output << i+1 << " " << C[i] << " " << DC[i] << " " << Sigma[i] << " " << DSigma[i] << endl;
        gCorr->SetPoint(i, i+1, Sigma[i]);
        gCorr->SetPointError(i, 0, DSigma[i]);

        /// Calibrate atom number using evaluated cross section ///////////////////////////////////
        N[i] = I * T[i] * C[i] * FissionRate[i] / evalCS / Flux[i] * 1.E22;
        ///////////////////////////////////////////////////////////////////////////////////////////

        DN[i] = N[i] * sqrt( pow(DC[i] / C[i], 2) +
                             pow(DFissionRate[i] / FissionRate[i], 2) +
                             pow(DFlux[i] / Flux[i], 2) );
        gN->SetPoint(i, i+1, N[i] * eff);
        gN->SetPointError(i, 0, DN[i] * eff);
        cout << " " << i+1 << "  " << NsponFis[i] << "+-" << DNsponFis[i] << "  " << N[i] << "+-" << DN[i] << "  " << Sigma[i] << "+-" << DSigma[i] << endl;

        /// Calibrate areal mass density using evaluated cross section ////////////////////////////
        mA[i] = I * T[i] * C[i] * FissionRate[i] * MolarMass * u / Deposit_area[i] / evalCS / Flux[i] * 1.E22 * 1000; // [mm^2 -> b, g -> mg]
        ///////////////////////////////////////////////////////////////////////////////////////////

        DmA[i] = mA[i] * sqrt( pow(DC[i] / C[i], 2) +
                               pow(DFissionRate[i] / FissionRate[i], 2) +
                               pow(DFlux[i] / Flux[i], 2) +
                               pow(DeltaDepositArea / Deposit_area[i], 2) );
        gmA->SetPoint(i, i+1, mA[i]);
        gmA->SetPointError(i, 0, DmA[i]);

        // Make average
        sumSigma += Sigma[i] * N[i];
        D2sumSigma += pow(DSigma[i] * N[i], 2);
        sumAtoms += N[i];
    }
    Save(fAna, "PuFC/Correction", gCorr);
    Save(fAna, "PuFC/nAtoms", gN);
    Save(fAna, "PuFC/nAtoms", gmA);
    fAna->Save();
    fAna->Close();
    Double_t avSigma = sumSigma / sumAtoms;
    Double_t DavSigma = sqrt(D2sumSigma) / sumAtoms;
    cout << "Weighted average: " << avSigma << " $\\pm$ " << DavSigma << endl;
    output << "# average " << avSigma << " " << DavSigma;
    output.close();
}

void Correction()
{
//    ScatteringUFC();
    ApplyCorrectionsU("Geant4");
//    ApplyCorrectionsU("MCNP");

    ApplyCorrectionsPu("Geant4");
//    ApplyCorrectionsPu("MCNP");
}

#endif


