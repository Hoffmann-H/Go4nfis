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
    if (!geFissionRate) cout << "Could not get " << name << endl;

    // Get number of target atoms
    sprintf(name, "%s/nAtoms/%s_effN", FC.c_str(), FC.c_str());
    TGraphErrors *geAtoms = (TGraphErrors*)fAna->Get(name);
    if (!geAtoms) cout << "Could not get " << name << endl;

    // Get neutron density rate
    sprintf(name, "%s/NeutronField/%s/%s_Density", FC.c_str(), Run.c_str(), Run.c_str());
    TGraphErrors *geFlux = (TGraphErrors*)fAna->Get(name);
    if (!geFlux) cout << "Could not get " << name << endl;

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

//        cout << DCrossSection / CrossSection << " = " << DFissionRate / FissionRate << " + " << DAtoms / Atoms << " + " << DFlux / Flux << endl;
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

void PrintAnalysis(string StrAnaFile)
{
    cout << "----------------------------------------" << endl
         << "Neutron-induced fission at 15MeV" << endl
         << "Experiment 05.2014 at PTB Braunschweig" << endl
         << "----------------------------------------" << endl;

    char name[128] = "";
    TFile *fAna = TFile::Open(StrAnaFile.c_str(), "READ");
    string FCs[] = {"UFC", "PuFC"};
    Double_t EvalCS[] = {2.073, 2.103};
    string FGRuns[] = {"UFC_NIF", "NIF"};
    Double_t k_I[] = {1.0398, 1.0};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        cout << endl << "----------------------------------------" << endl
             << "" << FC << endl
             << "Evaluated cross section (ENDF/B-VIII.0): " << EvalCS[i_FC] << endl
             << "----------------------------------------" << endl << endl;
        cout << "\tConstant correction factors" << endl
             << "\tIsotope vector: k_I = " << k_I[i_FC] << endl;

        cout << endl << "\tChannel-dependant data and corrections" << endl
             << "\t----------------------------------------" << endl
             << "\tChannel \t(n,f) rate  \tPhi_n \teN(nELBE)  \tk_T  \tk_S \tk_eps \tmeas.CS \tcalib.eN \tdCS \t dN" << endl
             << "\t----------------------------------------" << endl;

        sprintf(name, "%s/ToF/Signal/%s/FissionRate", FC.c_str(), FGRuns[i_FC].c_str());
        TGraphErrors *pGrFissionRate = (TGraphErrors*)fAna->Get(name); if (!pGrFissionRate) cout << "Could not get " << name << endl;
        sprintf(name, "%s/NeutronField/%s/%s_Density", FC.c_str(), FGRuns[i_FC].c_str(), FGRuns[i_FC].c_str());
        TGraphErrors *pGrFluxDensity = (TGraphErrors*)fAna->Get(name); if (!pGrFluxDensity) cout << "Could not get " << name << endl;
        sprintf(name, "%s/nAtoms/%s_eN", FC.c_str(), FC.c_str());
        TGraphErrors *pGrEffN = (TGraphErrors*)fAna->Get(name); if (!pGrEffN) cout << "Could not get " << name << endl;
        sprintf(name, "%s/Correction/%s_Target_Gate", FC.c_str(), FC.c_str());
        TGraphErrors *pGrTarget = (TGraphErrors*)fAna->Get(name); if (!pGrTarget) cout << "Could not get " << name << endl;
        sprintf(name, "%s/Correction/%s_Geant4_real_C", FC.c_str(), FC.c_str());
        TGraphErrors *pGrG4 = (TGraphErrors*)fAna->Get(name); if (!pGrG4) cout << "Could not get " << name << endl;
        sprintf(name, "Carlson/Carlson_Correction");
        TGraphErrors *pGrCarlson = (TGraphErrors*)fAna->Get(name); if (!pGrCarlson) cout << "Could not get " << name << endl;
        sprintf(name, "%s/Correction/%s_CS_Geant4_corrected", FC.c_str(), FC.c_str());
        TGraphErrors *pGrCS = (TGraphErrors*)fAna->Get(name); if (!pGrCS) cout << "Could not get " << name << endl;
        sprintf(name, "%s/nAtoms/%s_eN_cal_Geant4", FC.c_str(), FC.c_str());
        TGraphErrors *pGrEffNcal = (TGraphErrors*)fAna->Get(name); if (!pGrEffNcal) cout << "Could not get " << name << endl;

        Double_t x, FissionRate, FluxDensity, EffN, Target, G4, Carlson, CS, EffNcal;
        for (Int_t i = 0; i < 8; i++)
        {
            pGrFissionRate->GetPoint(i, x, FissionRate);
            pGrFluxDensity->GetPoint(i, x, FluxDensity);
            pGrEffN->GetPoint(i, x, EffN);
            pGrTarget->GetPoint(i, x, Target);
            pGrG4->GetPoint(i, x, G4);
            if (i_FC == 0) pGrCarlson->GetPoint(i, x, Carlson); else Carlson = 1.0;
            pGrCS->GetPoint(i, x, CS);
            pGrEffNcal->GetPoint(i, x, EffNcal);

            cout << "\t   " << i+1 << "\t   " << FissionRate << "\t   " << FluxDensity << "\t   " << EffN << " \t   " << Target << " \t   " << G4 << " \t   " << Carlson << "\t "
                 << CS << "\t " << EffNcal << "\t " << 1.E22 * FissionRate / FluxDensity / EffN * k_I[i_FC] * Target * G4 * Carlson << " \t" << 1.E22 * FissionRate / FluxDensity / EvalCS[i_FC] * k_I[i_FC] * Target * G4 * Carlson << endl;
                    //CS / EvalCS[i_FC] << "\t " << EffNcal / EffN << endl;
        }
        cout << "\t----------------------------------------" << endl;
    }
}

#endif
