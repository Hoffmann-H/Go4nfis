#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H
#include "Runs.C"
using namespace std;

string CrossSectionRun(string Run)
{
    TFile *fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    char name[64] = "";

    // choose fission chamber
    string FC = (Run[0] == 'U') ? "UFC" : "PuFC";
//    cout << "Cross section for run " << Run << ", auto-assigned to " << FC << endl;

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

    // Initialize returned string
    std::stringstream line;
    line << Run;

    // loop over deposits
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t x, FissionRate, Atoms, Flux, CrossSection;
        Double_t DFissionRate, DAtoms, DFlux, DCrossSection;
        // Extract values from graphs
        geFissionRate->GetPoint(i, x, FissionRate);
        DFissionRate = geFissionRate->GetErrorY(i);
        geAtoms->GetPoint(i, x, Atoms);
        DAtoms = geAtoms->GetErrorY(i);
        geFlux->GetPoint(i, x, Flux);
        DFlux = geFlux->GetErrorY(i);
        // Calculate cross section
        CrossSection = FissionRate / (Atoms * Flux) * 1.E22; // mm^2 -> barn
        DCrossSection = CrossSection * sqrt(
                    pow(DFissionRate / FissionRate, 2) +
                    pow(DAtoms / Atoms, 2) +
                    pow(DFlux / Flux, 2) );
        // Write cross section to graph
        geCrossSection->SetPoint(i, i+1, CrossSection);
        geCrossSection->SetPointError(i, 0, DCrossSection);

        line << " " << CrossSection << " " << DCrossSection;
    }
    sprintf(name, "%s_CS_raw", Run.c_str());
    geCrossSection->SetName(name);
    Save(fAna, FC+"/CrossSection", geCrossSection);
//    cout << line.str() << endl;
    fAna->Save();
    fAna->Close();
    return line.str();
}

void DoCrossSection(string FC)
{
    string FileName = "CrossSection_"+FC+".txt";
    ofstream output("../results/"+FileName);
    output << "# nr run {cs[b] unc[b]}" << endl;
    Int_t nRuns = GetnRuns(FC);
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(FC, j);
        output << j+1 << " " << CrossSectionRun(Run) << endl;
    }
    output.close();
}

void DoCrossSection()
{
    string FileName = "CrossSection.txt";
    ofstream output("../results/"+FileName);
    output << "# nr run {cs[b] unc[b]}" << endl;
    Int_t nRuns = GetnRuns();
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(j);
        output << j+1 << " " << CrossSectionRun(Run) << endl;
    }
    output.close();
}

void CrossSection()
{
    DoCrossSection("UFC");
    DoCrossSection("UFC_FG");
    DoCrossSection("UFC_BG");
    CrossSectionRun("UFC_NIF");
    CrossSectionRun("UFC_SB");

    DoCrossSection("PuFC");
    DoCrossSection("PuFC_FG");
    DoCrossSection("PuFC_BG");
    CrossSectionRun("NIF");
    CrossSectionRun("SB");
}

#endif
