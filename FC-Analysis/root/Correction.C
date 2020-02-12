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
    char name[64] = "";
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gRaw = (TGraphErrors*)fAna->Get("UFC/CrossSection/UFC_NIF_CS_raw");
    if (!gRaw) cout << "Could not get " << "UFC/CrossSection/UFC_NIF_CS_raw" << endl;
    sprintf(name, "UFC/Correction/UFC_Target_Gate");
    TGraphErrors *gT = (TGraphErrors*)fAna->Get(name); if (!gT) cout << "Could not get " << name << endl;
//    sprintf(name, "UFC/Correction/UFC_%s_Gate_%ins", Simulation.c_str(), Gate);
//    TGraphErrors *gG = (TGraphErrors*)fAna->Get(name); if (!gG) cout << "Could not get " << name << endl;
    sprintf(name, "UFC/Correction/%s_UFC_real_C", Simulation.c_str());
    TGraphErrors *gC = (TGraphErrors*)fAna->Get(name); if (!gC) cout << "Could not get " << name << endl;
    TGraphErrors *gCarlson = (TGraphErrors*)fAna->Get("Carlson/Carlson_Correction");
    if (!gCarlson) cout << "Could not get " << "Carlson/Carlson_Correction" << endl;
    TGraphErrors *gAtoms = (TGraphErrors*)fAna->Get("UFC/nAtoms/UFC_effN");
    if (!gAtoms) cout << "Could not open " << "UFC/nAtoms/UFC_effN" << endl;

    // Get fission signal rate
    sprintf(name, "UFC/ToF/Signal/UFC_NIF/FissionRate");
    TGraphErrors *geFissionRate = (TGraphErrors*)fAna->Get(name);
    if (!geFissionRate) cout << "Could not get " << name << endl;

    // Get neutron density rate
    sprintf(name, "UFC/NeutronField/UFC_NIF/UFC_NIF_Density");
    TGraphErrors *geFlux = (TGraphErrors*)fAna->Get(name);
    if (!geFlux) cout << "Could not get " << name << endl;

    TGraphErrors *gCorr = new TGraphErrors(8);
    sprintf(name, "UFC_CS_%s_corrected", Simulation.c_str());
    gCorr->SetName(name);

    TGraphErrors *gN = new TGraphErrors(8);
    sprintf(name, "UFC_%s_calN", Simulation.c_str());
    gN->SetName(name);
    gN->SetTitle("U Calibration; Deposit; #varepsilon#font[12]{N}");

    Double_t I, DI;
    IsoVec("UFC", I, DI);
    Double_t Raw[8], DRaw[8];
    Double_t T[8], DT[8];
    Double_t C[8], DC[8];
    Double_t Carlson[8], DCarlson[8];
    Double_t Atoms[8], DAtoms[8];
    Double_t N[8], DN[8];
    Double_t FissionRate[8], DFissionRate[8], Flux[8], DFlux[8];
    Double_t Sigma[8], DSigma[8];
    Double_t sumSigma = 0;
    Double_t D2sumSigma = 0;
    Double_t sumAtoms = 0;
    ofstream output("../results/Correction_UFC.txt");
    cout << "ch   cal.N   Sigma" << endl;
    output << "# deposit uncorrected corrected" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get data
        Double_t x;
        gRaw->GetPoint(i, x, Raw[i]);
        DRaw[i] = gRaw->GetErrorY(i);
        gC->GetPoint(i, x, C[i]);
        DC[i] = gC->GetErrorY(i);
        gT->GetPoint(i, x, T[i]);
        DT[i] = gT->GetErrorY(i);
        gCarlson->GetPoint(i, x, Carlson[i]);
        DCarlson[i] = gCarlson->GetErrorY(i);
        gAtoms->GetPoint(i, x, Atoms[i]);
        DAtoms[i] = gAtoms->GetErrorY(i);
        geFissionRate->GetPoint(i, x, FissionRate[i]);
        DFissionRate[i] = geFissionRate->GetErrorY(i);
        geFlux->GetPoint(i, x, Flux[i]);
        DFlux[i] = geFlux->GetErrorY(i);

        // Apply corrections
        Sigma[i] = I * T[i] * C[i] * Carlson[i] * Raw[i];
        DSigma[i] = sqrt( /*pow(DI * T[i] * C[i] * Raw[i], 2) +*/ // Isotope: systematic
                          pow(I * T[i] * DC[i] * Raw[i], 2) +
                          pow(I * T[i] * C[i] * DRaw[i], 2) );
        output << i+1 << " " << Raw[i] << " " << DRaw[i] << " " << Sigma[i] << " " << DSigma[i] << endl;
        gCorr->SetPoint(i, i+1, Sigma[i]);
        gCorr->SetPointError(i, 0, DSigma[i]);

        // calibrate atom number
        N[i] = I * T[i] * C[i] / eff * FissionRate[i] / evalCS / Flux[i] * 1.E22;
        DN[i] = N[i] * sqrt( pow(DC[i] / C[i], 2) +
                             pow(DFissionRate[i] / FissionRate[i], 2) +
                             pow(DFlux[i] / Flux[i], 2) );
        gN->SetPoint(i, i+1, N[i] * eff);
        gN->SetPointError(i, 0, DN[i] * eff);
        cout << " " << i+1 << " " << Atoms[i] << "+-" << DAtoms[i] << "  " << N[i] << "+-" << DN[i] << "  " << Sigma[i] << "+-" << DSigma[i] << endl;

        // Make average
        sumSigma += Sigma[i] * Atoms[i];
        D2sumSigma += pow(DSigma[i] * Atoms[i], 2);
        sumAtoms += Atoms[i];
    }
    Save(fAna, "UFC/Correction", gCorr);
    Save(fAna, "UFC/nAtoms", gN);
    fAna->Save();
    fAna->Close();
    Double_t avSigma = sumSigma / sumAtoms;
    Double_t DavSigma = sqrt(D2sumSigma) / sumAtoms;
    cout << "Weighted average: " << avSigma << " $\\pm$ " << DavSigma << endl;
    output << "# average " << avSigma << " " << DavSigma;
    output.close();
}

void ApplyCorrectionsPu(string Simulation = "Geant4", Int_t Gate = 15)
{
    Double_t evalCS = 2.1026472;
    Double_t eff = 0.988;
    char name[64] = "";
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gRaw = (TGraphErrors*)fAna->Get("PuFC/CrossSection/NIF_CS_raw");
    if (!gRaw) cout << "Could not open " << "PuFC/CrossSection/NIF_CS_raw" << endl;
    sprintf(name, "PuFC/Correction/PuFC_Target_Gate");
    TGraphErrors *gT = (TGraphErrors*)fAna->Get(name); if (!gT) cout << "Could not get " << name << endl;
    sprintf(name, "PuFC/Correction/%s_PuFC_real_C", Simulation.c_str());
    TGraphErrors *gC = (TGraphErrors*)fAna->Get(name); if (!gC) cout << "Could not open " << name << endl;
//    sprintf(name, "PuFC/Correction/PuFC_%s_Gate_%ins", Simulation.c_str(), Gate);
//    TGraphErrors *gG = (TGraphErrors*)fAna->Get(name); if (!gG) cout << "Could not get " << name << endl;
    TGraphErrors *gAtoms = (TGraphErrors*)fAna->Get("PuFC/nAtoms/PuFC_effN");
    if (!gAtoms) cout << "Could not open " << "PuFC/nAtoms/PuFC_effN" << endl;

    // Get fission signal rate
    sprintf(name, "PuFC/ToF/Signal/NIF/FissionRate");
    TGraphErrors *geFissionRate = (TGraphErrors*)fAna->Get(name);
    if (!geFissionRate) cout << "Could not get " << name << endl;

    // Get neutron density rate
    sprintf(name, "PuFC/NeutronField/NIF/NIF_Density");
    TGraphErrors *geFlux = (TGraphErrors*)fAna->Get(name);
    if (!geFlux) cout << "Could not get " << name << endl;

    TGraphErrors *gCorr = new TGraphErrors(8);
    sprintf(name, "PuFC_CS_%s_corrected", Simulation.c_str());
    gCorr->SetName(name);

    TGraphErrors *gN = new TGraphErrors(8);
    sprintf(name, "PuFC_%s_calN", Simulation.c_str());
    gN->SetName(name);
    gN->SetTitle("Pu Calibration; Deposit; #varepsilon#font[12]{N}");

    Double_t I, DI;
    IsoVec("PuFC", I, DI);
    Double_t Raw[8], DRaw[8];
    Double_t T[8], DT[8];
    Double_t C[8], DC[8];
    Double_t Atoms[8], DAtoms[8];
    Double_t N[8], DN[8];
    Double_t FissionRate[8], DFissionRate[8], Flux[8], DFlux[8];
    Double_t Sigma[8], DSigma[8];
    Double_t sumSigma = 0;
    Double_t D2sumSigma = 0;
    Double_t sumAtoms = 0;
    ofstream output("../results/Correction_PuFC.txt");
    cout << "ch   effN   F   Sigma" << endl;
    output << "# deposit uncorrected factor corrected" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get data
        Double_t x;
        gRaw->GetPoint(i, x, Raw[i]);
        DRaw[i] = gRaw->GetErrorY(i);
        gT->GetPoint(i, x, T[i]);
        DT[i] = gT->GetErrorY(i);
        gC->GetPoint(i, x, C[i]);
        DC[i] = gC->GetErrorY(i);
        gAtoms->GetPoint(i, x, Atoms[i]);
        DAtoms[i] = gAtoms->GetErrorY(i);
        geFissionRate->GetPoint(i, x, FissionRate[i]);
        DFissionRate[i] = geFissionRate->GetErrorY(i);
        geFlux->GetPoint(i, x, Flux[i]);
        DFlux[i] = geFlux->GetErrorY(i);

        // Apply corrections
        Sigma[i] = T[i] * C[i] * Raw[i];
        DSigma[i] = sqrt( pow(T[i] * DC[i] * Raw[i], 2) +
                          pow(T[i] * C[i] * DRaw[i], 2) );
        output << i+1 << " " << Raw[i] << " " << DRaw[i] << " " << C[i] << " " << DC[i] << " " << Sigma[i] << " " << DSigma[i] << endl;
        gCorr->SetPoint(i, i+1, Sigma[i]);
        gCorr->SetPointError(i, 0, DSigma[i]);

        // calibrate atom number
        N[i] = I * T[i] * C[i] / eff * FissionRate[i] / evalCS / Flux[i] * 1.E22;
        DN[i] = N[i] * sqrt( pow(DC[i] / C[i], 2) +
                             pow(DFissionRate[i] / FissionRate[i], 2) +
                             pow(DFlux[i] / Flux[i], 2) );
        gN->SetPoint(i, i+1, N[i] * eff);
        gN->SetPointError(i, 0, DN[i] * eff);
        cout << " " << i+1 << " " << Atoms[i] << "+-" << DAtoms[i] << "  " << N[i] << "+-" << DN[i] << "  " << Sigma[i] << "+-" << DSigma[i] << endl;

        // Make average
//        if (i != 4 && i != 7) {
            sumSigma += Sigma[i] * Atoms[i];
            D2sumSigma += pow(DSigma[i] * Atoms[i], 2);
            sumAtoms += Atoms[i];
//        }
    }
    Save(fAna, "PuFC/Correction", gCorr);
    Save(fAna, "PuFC/nAtoms", gN);
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
    ApplyCorrectionsU("Geant4", RIGHT);
    ApplyCorrectionsU("MCNP", RIGHT);

    ApplyCorrectionsPu("Geant4", RIGHT);
    ApplyCorrectionsPu("MCNP", RIGHT);
}

#endif


