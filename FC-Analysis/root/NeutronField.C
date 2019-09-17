#ifndef NEUTRON_FIELD_H
#define NEUTRON_FIELD_H
#include "Runs.C"
#include <fstream>

string NeutronFieldRun(string Run)
{ // Find 1 run's neutron flux. Return output string
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Double_t Yield = 2.158821152E4;
    Double_t DYield = 543.3;
    char name[64] = "";

    // choose fission chamber
    string FC = (Run[0] == 'U') ? "UFC" : "PuFC";
    // Get monitor data
    Double_t Monitor, DeltaRel, tReal;
    if (!GetRun(Run, Monitor, DeltaRel, tReal))
        cout << "Could not find run " << Run << endl;

    // Create graphs
    TGraphErrors *geFluence = new TGraphErrors(8);
    TGraphErrors *geFlux = new TGraphErrors(8);
    TGraphErrors *geDensity = new TGraphErrors(8);

    // Initialize returned string
    std::stringstream line;
    line << Run << " " << Monitor / tReal << " " << DeltaRel * Monitor / tReal;

    // lop over deposits
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t w = SolidAngle(i, FC);
        Double_t area = DepositArea(i, FC);
        Double_t nFluence = Yield * Monitor * w;
        Double_t DnFluence = sqrt( pow(DYield * Monitor * w, 2) +
                                   pow(Yield * DeltaRel * Monitor * w, 2) +
                                   pow(Yield * Monitor * 0, 2) );
        Double_t nFlux = nFluence / area;
        Double_t DnFlux = sqrt( pow(DnFluence / area, 2) +
                                pow(0 * nFluence / (area*area), 2) );
        Double_t nDensity = nFlux / tReal;
        Double_t DnDensity = DnFlux / tReal;

        geFluence->SetPoint(i, i+1, nFluence);
        geFluence->SetPointError(i, 0, DnFluence);
        geFlux->SetPoint(i, i+1, nFlux);
        geFlux->SetPointError(i, 0, DnFlux);
        geDensity->SetPoint(i, i+1, nDensity);
        geDensity->SetPointError(i, 0, DnDensity);

        line << " " << nDensity << " " << DnDensity;
    }
    // Save
    Save(fAna, FC+"/NeutronField/"+Run, geFluence, Run+"_Fluence");
    Save(fAna, FC+"/NeutronField/"+Run, geFlux, Run+"_Flux");
    Save(fAna, FC+"/NeutronField/"+Run, geDensity, Run+"_Density");
    fAna->Save();
    fAna->Close();
    return line.str();
}

void NeutronField(string FC)
{
    string FileName = "NeutronField_"+FC+".txt";
    ofstream output("../results/"+FileName);
    output << "# nr run monitor[1/s] unc[1/s] {flux[1/(s*mm^2)] unc[1/(s*mm^2)]}" << endl;
    Int_t nRuns = GetnRuns(FC);
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(FC, j);
        // Analyze monitor data
        cout << Run << endl;
        output << j+1 << " " << NeutronFieldRun(Run) << endl;
    }
}

void NeutronField()
{
    string FileName = "NeutronField.txt";
    ofstream output("../results/"+FileName);
    output << "# nr run monitor[1/s] unc[1/s] {flux[1/(s*mm^2)] unc[1/(s*mm^2)]}" << endl;
    Int_t nRuns = GetnRuns();
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(j);
        // Analyze monitor data
        output << j+1 << " " << NeutronFieldRun(Run) << endl;
    }
}
#endif
