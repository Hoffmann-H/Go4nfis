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

void GetVariable(TFile *fAna, string var_name, string run_name, Int_t ch, Double_t &var, Double_t &d_var)
{ // unused
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

string Format(Int_t skip)
{
    std::stringstream s;
    s << "%lg %*s";
    for (Int_t i = 0; i < skip; i++)
        s << " %*lg %*lg";
    s << " %lg %lg";
    return s.str();
}

void Runs(string FC)
{
    char name[128] = "";
    TFile* fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Int_t nRuns = GetnRuns(FC);

    // treat neutron field data
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/NeutronField_%s.txt", FC.c_str());
    TGraphErrors *geMonitorRate = new TGraphErrors(name, Format(0).c_str());
    if (!geMonitorRate)
        cout << "Could not open " << name << endl;
    Save(fAna, FC+"/NeutronField", geMonitorRate, FC+"_MonitorRate");

    // Neutron flux
    TGraphErrors *geFlux[8];
    for (Int_t i = 0; i < 8; i++)
    {
        geFlux[i] = new TGraphErrors(name, Format(i+1).c_str());
        Save(fAna, FC+"/NeutronField", geFlux[i], FC+"_NeutronFlux_"+to_string(i+1));
    }

    // Uncorrected cross section
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/CrossSection_%s.txt", FC.c_str());
    TGraphErrors *geCrossSection[8];
    for (Int_t i = 0; i < 8; i++)
    {
        geCrossSection[i] = new TGraphErrors(name, Format(i).c_str());
        Save(fAna, FC+"/CrossSection", geCrossSection[i], FC+"_CrossSection_"+to_string(i+1));
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
