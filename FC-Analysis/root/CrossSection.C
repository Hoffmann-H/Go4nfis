#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H
#include "Runs.C"

void CrossSection(string Run)
{
    char name[64] = "";
    // choose fission chamber
    string FC = (Run[0] == 'U') ? "UFC" : "PuFC";
    cout << "Cross section for run " << Run << ", auto-assigned to " << FC << endl;
    TFile *fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

    // Get fission signal rate
    sprintf(name, "%s/ToF/Signal/%s/FissionRate", FC.c_str(), Run.c_str());
    TGraphErrors *geFissionRate = (TGraphErrors*)fAna->Get(name);
    if (geFissionRate == 0)
        cout << "Could not open " << name << endl;

    // Get number of target atoms
    sprintf(name, "%s/nAtoms/%s_effN", FC.c_str(), FC.c_str());
    TGraphErrors *geAtoms = (TGraphErrors*)fAna->Get(name);
    if (geAtoms == 0)
        cout << "Could not open " << name << endl;

    // Get neutron density rate
    sprintf(name, "%s/NeutronField/%s/%s_Density", FC.c_str(), Run.c_str(), Run.c_str());
    TGraphErrors *geFlux = (TGraphErrors*)fAna->Get(name);
    if (geFlux == 0)
        cout << "Could not open " << name << endl;

    // Create cross section graph
    TGraphErrors *geCrossSection = new TGraphErrors(8);
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t x, FissionRate, Atoms, Flux, DFissionRate, DAtoms, DFlux, CrossSection, DCrossSection;
        // Extract values from graphs
        geFissionRate->GetPoint(i, x, FissionRate);
        DFissionRate = geFissionRate->GetErrorY(i);
        geAtoms->GetPoint(i, x, Atoms);
        DAtoms = geAtoms->GetErrorY(i);
        geFlux->GetPoint(i, x, Flux);
        DFlux = geFlux->GetErrorY(i);
        // Calculate cross section
        CrossSection = FissionRate / (Atoms * Flux);
        DCrossSection = CrossSection * sqrt(
                    pow(DFissionRate / FissionRate, 2) +
                    pow(DAtoms / Atoms, 2) +
                    pow(DFlux / Flux, 2) );
        // Write cross section to graph
        geCrossSection->SetPoint(i, i+1, CrossSection);
        geCrossSection->SetPointError(i, 0, DCrossSection);
    }
    Save(fAna, FC+"/CrossSection", geCrossSection, Run+"_CS_raw");
}

void CrossSection()
{
    Int_t nRuns = GetnRuns();
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(j);
        CrossSection(Run);
    }
}

#endif
