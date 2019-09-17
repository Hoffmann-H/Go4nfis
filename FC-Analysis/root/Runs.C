#ifndef RUNS_H
#define RUNS_H
#include "FC.C"
#include "SaveToFile.C"
#include <fstream>

/// Manage runs here
Int_t GetnRuns(string FC)
{
    if (strcmp(FC.c_str(), "PuFC"))
        return 7; // UFC
    else
        return 7; // PuFC
}
string GetRunName(string FC, Int_t j)
{
    if (strcmp(FC.c_str(), "PuFC"))
    { // UFC
        string RunNames[] = {"UFC_FG_MS20_2",
                             "UFC_FG_MS20_3",
                             "UFC_FG_MS20_4",
                             "UFC_FG_MS21_2",
                             "UFC_FG_MS21_3",
                             "UFC_BG_MS20_5",
                             "UFC_BG_MS21_4"};
        return RunNames[j];
    }
    else
    { // PuFC
        string RunNames[] = {"PuFC_FG_MS4",
                             "PuFC_FG_MS5",
                             "PuFC_FG_MS6",
                             "PuFC_FG_MS7",
                             "PuFC_BG_MS9",
                             "PuFC_BG_MS10",
                             "PuFC_BG_MS11"};
        return RunNames[j];
    }
}
Int_t GetnRuns()
{
    return 18;
}
string GetRunName(Int_t j)
{
    string RunNames[] = {
                        "UFC_FG_MS20_2",
                        "UFC_FG_MS20_3",
                        "UFC_FG_MS20_4",
                        "UFC_FG_MS21_2",
                        "UFC_FG_MS21_3",
                        "UFC_BG_MS20_5",
                        "UFC_BG_MS21_4",
                        "UFC_NIF",
                        "UFC_SB",
                        "PuFC_FG_MS4",
                        "PuFC_FG_MS5",
                        "PuFC_FG_MS6",
                        "PuFC_FG_MS7",
                        "PuFC_BG_MS9",
                        "PuFC_BG_MS10",
                        "PuFC_BG_MS11",
                        "NIF",
                        "SB"};
    return RunNames[j];
}

Bool_t GetRun(string RunName, Double_t &monitor, Double_t &delta_rel, Double_t &t_real)
{
    // Open txt tabular, break if not successful
    std::ifstream input("/home/hoffma93/Programme/Go4nfis/FC-Analysis/data/runs.txt");
    if (!input.is_open())
    {
        cout << "Could not open " << "/home/hoffma93/Programme/Go4nfis/FC-Analysis/data/runs.txt" << endl;
        return 0;
    }
    string Run = "", lastRun;
    Double_t Monitor, DeltaRel, tReal;
    do {
        lastRun = Run; // keep last run name
        // Read line
        input >> Run >> Monitor >> DeltaRel >> tReal;
        // Break if last line was reached
        if (!strcmp(Run.c_str(), lastRun.c_str()))
            return 0;
        // repeat until names match
    } while (strcmp(Run.c_str(), RunName.c_str()));
    // Save numbers to reference
    monitor = Monitor;
    delta_rel = DeltaRel;
    t_real = tReal;
    return 1;
}

void GetVariable(TFile *fAna, string var_name, string run_name, Int_t ch, Double_t &var, Double_t &d_var)
{
    char name[128] = "";
    string FC = (run_name[0] == 'U') ? "UFC" : "PuFC";
    sprintf(name, "%s/NeutronField/%s/%s_%s", FC.c_str(), run_name.c_str(), run_name.c_str(), var_name.c_str());
    TGraphErrors *geVar = (TGraphErrors*)fAna->Get(name);
    if (geVar == 0)
        cout << "Could not get TGraphErrors " << name << endl;
    Double_t x, y;
    geVar->GetPoint(ch, x, y);
    var = y;
    d_var = geVar->GetErrorY(ch);
    return;
}

void Runs(string FC)
{
    char name[128] = "";
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Int_t nRuns = GetnRuns(FC);

    TGraphErrors *geMonitor = new TGraphErrors(nRuns);
    TGraphErrors *geMonitorRate = new TGraphErrors(nRuns);
    TGraphErrors *geFluence[8];
    TGraphErrors *geFlux[8];
    TGraphErrors *geDensity[8];


    for (Int_t i = 0; i < 8; i++)
    {
        geFluence[i] = new TGraphErrors(nRuns);
        geFlux[i] = new TGraphErrors(nRuns);
        geDensity[i] = new TGraphErrors(nRuns);
    }

    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(FC, j);
        // Get monitor counts
        Double_t Monitor, DeltaRel, tReal;
        if (!GetRun(Run, Monitor, DeltaRel, tReal))
            cout << "Could not find run " << Run << endl;
        geMonitor->SetPoint(j, j+1, Monitor);
        geMonitor->SetPointError(j, 0, DeltaRel * Monitor);
        geMonitorRate->SetPoint(j, j+1, Monitor / tReal);
        geMonitorRate->SetPointError(j, 0, DeltaRel * Monitor / tReal);
        for (Int_t i = 0; i < 8; i++)
        {
            Double_t Var, DVar;
            // Get neutron fluence
            GetVariable(fAna, "Fluence", Run, i, Var, DVar);
            geFluence[i]->SetPoint(j, j+1, Var);
            geFluence[i]->SetPointError(j, 0, DVar);
            // Get neutron flux
            GetVariable(fAna, "Flux", Run, i, Var, DVar);
            geFlux[i]->SetPoint(j, j+1, Var);
            geFlux[i]->SetPointError(j, 0, DVar);
            // Get neutron flux density
            GetVariable(fAna, "Density", Run, i, Var, DVar);
            geDensity[i]->SetPoint(j, j+1, Var);
            geDensity[i]->SetPointError(j, 0, DVar);

        }

    }
    Save(fAna, FC+"/NeutronField", geMonitor, FC+"_Monitor");
    Save(fAna, FC+"/NeutronField", geMonitorRate, FC+"_MonitorRate");
    for (Int_t i = 0; i < 8; i++)
    {
        Save(fAna, FC+"/NeutronField/Fluence", geFluence[i], FC+"_Fluence_"+to_string(i+1));
        Save(fAna, FC+"/NeutronField/Flux", geFlux[i], FC+"_Flux_"+to_string(i+1));
        Save(fAna, FC+"/NeutronField/Density", geDensity[i], FC+"_Density_"+to_string(i+1));
    }
    fAna->Save();
    fAna->Close();
}//*/

void Runs()
{
    Runs("UFC");
    Runs("PuFC");
}
#endif
