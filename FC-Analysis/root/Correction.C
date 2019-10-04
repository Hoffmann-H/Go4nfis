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
        Double_t effNtotSB, DeffNtotSB, effNtotOpen, DeffNtotOpen, SimS, DSimS;
        geSimSB->GetPoint(i, x, effNtotSB);
        DeffNtotSB = geSimSB->GetErrorY(i);
        geSimOpen->GetPoint(i, x, effNtotOpen);
        DeffNtotOpen = geSimOpen->GetErrorY(i);
        SimS = 1 - effNtotSB / effNtotOpen;
        DSimS = effNtotSB / effNtotOpen * sqrt( pow(DeffNtotSB / effNtotSB, 2) +
                                                pow(DeffNtotOpen / effNtotOpen, 2) );
        geSimS->SetPoint(i, i+1, SimS);
        geSimS->SetPointError(i, 0, DSimS);
    }
    geExpS->SetName("ExpS");
    Save(fAna, "UFC/Correction", geExpS);
    geSimS->SetName("SimS");
    Save(fAna, "UFC/Correction", geSimS);
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

void ApplyCorrectionsU()
{
    char name[64] = "";
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gRaw = (TGraphErrors*)fAna->Get("UFC/CrossSection/UFC_NIF_CS_raw");
    if (!gRaw)
        cout << "Could not open " << "UFC/CrossSection/UFC_NIF_CS_raw" << endl;
//    TGraphErrors *gS = (TGraphErrors*)fAna->Get("UFC/Correction/ExpS");
    TGraphErrors *gS = (TGraphErrors*)fAna->Get("Simulation/Geant4/UFC_Open/Correction/S");
    if (!gS)
        cout << "Could not open " << "Simulation/Geant4/UFC_Open/Correction/S" << endl;
    TGraphErrors *gT = (TGraphErrors*)fAna->Get("Simulation/Geant4/UFC_Open/Correction/T");
    if (!gT)
        cout << "Could not open " << "Simulation/Geant4/UFC_Open/Correction/T" << endl;
    TGraphErrors *gAtoms = (TGraphErrors*)fAna->Get("UFC/nAtoms/UFC_effN");
    if (!gAtoms)
        cout << "Could not open " << "UFC/nAtoms/UFC_effN" << endl;
    TGraphErrors *gCorr = new TGraphErrors(8);
    Double_t I, DI;
    IsoVec("UFC", I, DI);
    Double_t Raw[8], DRaw[8];
    Double_t S[8], DS[8];
    Double_t T[8], DT[8];
    Double_t F[8], DF[8];
    Double_t Atoms[8];
    Double_t Sigma[8], DSigma[8];
    Double_t sumSigma = 0;
    Double_t D2sumSigma = 0;
    Double_t sumAtoms = 0;
    ofstream output("../results/Correction_UFC.txt");
    cout << "ch   effN   S   T   Sigma" << endl;
    output << "# deposit uncorrected factor corrected" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get data
        Double_t x;
        gRaw->GetPoint(i, x, Raw[i]);
        DRaw[i] = gRaw->GetErrorY(i);
        gS->GetPoint(i, x, S[i]);
        DS[i] = gS->GetErrorY(i);
        gT->GetPoint(i, x, T[i]);
        DT[i] = gT->GetErrorY(i);
        gAtoms->GetPoint(i, x, Atoms[i]);

        // Apply corrections
        F[i] = I * S[i] / T[i];
        DF[i] = F[i] * sqrt(pow(DI / I, 2) +
                            pow(DS[i] / S[i], 2) +
                            pow(DT[i] / T[i], 2));
        Sigma[i] = F[i] * Raw[i];
        DSigma[i] = sqrt( pow(DF[i] * Raw[i], 2) +
                          pow(DRaw[i] * F[i], 2) );
        cout << " " << i+1 << "  " << Atoms[i]  << "  " << S[i] << "+-" << DS[i] << "  " << T[i] << "+-" << DT[i] << "  " << Sigma[i] << "+-" << DSigma[i] << endl;
        output << i+1 << " " << Raw[i] << " " << DRaw[i] << " " << F[i] << " " << DF[i] << " " << Sigma[i] << " " << DSigma[i] << endl;
        gCorr->SetPoint(i, i+1, Sigma[i]);
        gCorr->SetPointError(i, 0, DSigma[i]);

        // Make average
        sumSigma += Sigma[i] * Atoms[i];
        D2sumSigma += pow(DSigma[i] * Atoms[i], 2);
        sumAtoms += Atoms[i];
    }
    gCorr->SetName("UFC_CS_corrected");
    Save(fAna, "UFC/Correction", gCorr);
    Double_t avSigma = sumSigma / sumAtoms;
    Double_t DavSigma = sqrt(D2sumSigma) / sumAtoms;
    cout << "Weighted average: " << avSigma << " $\\pm$ " << DavSigma << endl;
    output << "# average " << avSigma << " " << DavSigma;
    output.close();
}

void ApplyCorrectionsPu()
{
    char name[64] = "";
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gRaw = (TGraphErrors*)fAna->Get("PuFC/CrossSection/NIF_CS_raw");
    if (!gRaw)
        cout << "Could not open " << "PuFC/CrossSection/NIF_CS_raw" << endl;
//    TGraphErrors *gS = (TGraphErrors*)fAna->Get("Simulation/MCNP/Correction/S");
//    TGraphErrors *gT = (TGraphErrors*)fAna->Get("Simulation/MCNP/Correction/T");
    TGraphErrors *gF = (TGraphErrors*)fAna->Get("Simulation/MCNP/Correction/F");
    if (!gF)
        cout << "Could not open " << "Simulation/MCNP/Correction/F" << endl;
    TGraphErrors *gAtoms = (TGraphErrors*)fAna->Get("PuFC/nAtoms/PuFC_effN");
    if (!gF)
        cout << "Could not open " << "PuFC/nAtoms/PuFC_effN" << endl;
    TGraphErrors *gCorr = new TGraphErrors(8);
    Double_t I, DI;
    IsoVec("PuFC", I, DI);
    Double_t Raw[8], DRaw[8];
//    Double_t S[8], DS[8];
//    Double_t T[8], DT[8];
    Double_t F[8], DF[8];
    Double_t Atoms[8];
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
        gF->GetPoint(i, x, F[i]);
        DF[i] = gF->GetErrorY(i);
        gAtoms->GetPoint(i, x, Atoms[i]);

        // Apply corrections
        Sigma[i] = F[i] * Raw[i];
        DSigma[i] = sqrt( pow(DF[i] * Raw[i], 2) +
                          pow(DRaw[i] * F[i], 2) );
        cout << " " << i+1 << "  " << Atoms[i]  << "  " << F[i] << "+-" << DF[i] << "  " << Sigma[i] << "+-" << DSigma[i] << endl;
        output << i+1 << " " << Raw[i] << " " << DRaw[i] << " " << F[i] << " " << DF[i] << " " << Sigma[i] << " " << DSigma[i] << endl;
        gCorr->SetPoint(i, i+1, Sigma[i]);
        gCorr->SetPointError(i, 0, DSigma[i]);

        // Make average
        sumSigma += Sigma[i] * Atoms[i];
        D2sumSigma += pow(DSigma[i] * Atoms[i], 2);
        sumAtoms += Atoms[i];
    }
    gCorr->SetName("PuFC_CS_corrected");
    Save(fAna, "PuFC/Correction", gCorr);
    Double_t avSigma = sumSigma / sumAtoms;
    Double_t DavSigma = sqrt(D2sumSigma) / sumAtoms;
    cout << "Weighted average: " << avSigma << " $\\pm$ " << DavSigma << endl;
    output << "# average " << avSigma << " " << DavSigma;
    output.close();
}

void Correction()
{
//    ScatteringUFC();
    ApplyCorrectionsU();
    ApplyCorrectionsPu();
}

#endif
