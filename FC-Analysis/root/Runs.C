#ifndef RUNS_H
#define RUNS_H
#include "FC.C"
#include "SaveToFile.C"
#include <fstream>
#include <sstream>
using namespace std;

/// Manage runs here
Int_t GetnRuns(string FC)
{
    if (!strcmp(FC.c_str(), "UFC"))
        return 7;
    else if (!strcmp(FC.c_str(), "PuFC"))
        return 7;
    else if (!strcmp(FC.c_str(), "UFC_FG"))
        return 5;
    else if (!strcmp(FC.c_str(), "PuFC_FG"))
        return 4;
    else if (!strcmp(FC.c_str(), "UFC_BG"))
        return 2;
    else if (!strcmp(FC.c_str(), "PuFC_BG"))
        return 2;
    return 0;
}
string GetRunName(string FC, Int_t j)
{
    if (!strcmp(FC.c_str(), "UFC"))
    { // UFC
        string RunNames[] = {"UFC_FG_MS20_2",
                             "UFC_FG_MS20_3",
                             "UFC_FG_MS20_4",
                             "UFC_FG_MS21_2",
                             "UFC_FG_MS21_3",
                             "UFC_BG_MS20_5",
                             "UFC_BG_MS21_4"};
        return RunNames[j];
    } else if (!strcmp(FC.c_str(), "PuFC"))
    { // PuFC
        string RunNames[] = {"PuFC_FG_MS4",
                             "PuFC_FG_MS5",
                             "PuFC_FG_MS6",
                             "PuFC_FG_MS7",
                             "PuFC_BG_MS9",
                             "PuFC_BG_MS10",
                             "PuFC_BG_MS11"};
        return RunNames[j];
    } else if (!strcmp(FC.c_str(), "UFC_FG"))
    { // PuFC
        string RunNames[] = {"UFC_FG_MS20_2",
                             "UFC_FG_MS20_3",
                             "UFC_FG_MS20_4",
                             "UFC_FG_MS21_2",
                             "UFC_FG_MS21_3"};
        return RunNames[j];
    } else if (!strcmp(FC.c_str(), "PuFC_FG"))
    { // PuFC
        string RunNames[] = {"PuFC_FG_MS4",
                             "PuFC_FG_MS5",
                             "PuFC_FG_MS6",
                             "PuFC_FG_MS7"};
        return RunNames[j];
    } else if (!strcmp(FC.c_str(), "UFC_BG"))
    { // PuFC
        string RunNames[] = {"UFC_BG_MS20_5",
                             "UFC_BG_MS21_4"};
        return RunNames[j];
    } else if (!strcmp(FC.c_str(), "PuFC_BG"))
    { // PuFC
        string RunNames[] = {"PuFC_BG_MS9",
                             "PuFC_BG_MS10",
                             "PuFC_BG_MS11"};
        return RunNames[j];
    }
    return "";
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

TGraphErrors* Divide(TGraphErrors *g1, TGraphErrors *g2)
{
    Int_t N = g1->GetN();
    TGraphErrors *g3 = new TGraphErrors(N);
    cout << g1->GetName() << " / " << g2->GetName() << endl;
    for (Int_t j = 0; j < N; j++)
    {
        Double_t x, y1, y2, yerr1, yerr2;
        g1->GetPoint(j, x, y1);
        yerr1 = g1->GetErrorY(j);
        g2->GetPoint(j, x, y2);
        yerr2 = g2->GetErrorY(j);
        Double_t y3 = y1 / y2;
        Double_t yerr3 = y3 * sqrt(pow(yerr1 / y1, 2) + pow(yerr2 / y2, 2));
        cout << " " << y1 << " / " << y2 << " = " << y3 << ", " << yerr1 << " + " << yerr2 << " -> " << yerr3 << endl;
        g3->SetPoint(j, x, y3);
        g3->SetPointError(j, 0, yerr3);
    }
    return g3;
}

void Runs(string FC, string Setup = "")
{
    char name[128] = "";
    TFile* fAna = TFile::Open(results_file, "UPDATE");
    Int_t nRuns = GetnRuns(FC+Setup);

    // treat neutron field data
    sprintf(name, "../results/NeutronField_%s%s.txt", FC.c_str(), Setup.c_str());
    TGraphErrors *geMonitor = new TGraphErrors(name, Format(0).c_str());
    if (!geMonitor)
        cout << "Could not open " << name << endl;
    sprintf(name, "%s%s_Monitor", FC.c_str(), Setup.c_str());
    geMonitor->SetName(name);
    Save(fAna, FC+"/NeutronField", geMonitor);
    sprintf(name, "../results/NeutronField_%s%s.txt", FC.c_str(), Setup.c_str());
    TGraphErrors *geMonitorRate = new TGraphErrors(name, Format(1).c_str());
    sprintf(name, "%s%s_MonitorRate", FC.c_str(), Setup.c_str());
    geMonitorRate->SetName(name);
    Save(fAna, FC+"/NeutronField", geMonitorRate);

    // Neutron flux
    TGraphErrors *geFlux[8];
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "../results/NeutronField_%s%s.txt", FC.c_str(), Setup.c_str());
        geFlux[i] = new TGraphErrors(name, Format(i+2).c_str());
//        cout << "Saving " << FC << "/CrossSection/" << FC << Setup << "_NeutronFlux_" << i+1 << endl;
        Save(fAna, FC+"/NeutronField", geFlux[i], FC+Setup+"_NeutronFlux_"+to_string(i+1));
    }

    // Ratio to monitor
    TGraphErrors *geFission[8];
    TGraphErrors *geFissionRate[8];
    TGraphErrors *geRatio[8];
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "../results/InducedFission_%s%s.txt", FC.c_str(), Setup.c_str());
        geFission[i] = new TGraphErrors(name, Format(2*i).c_str());
        sprintf(name, "%s%s_InducedFission_%i", FC.c_str(), Setup.c_str(), i+1);
        geFission[i]->SetName(name);
//        cout << "Saving " << FC << "/ToF/Signal/" << name << endl;
        Save(fAna, FC+"/ToF/Signal", geFission[i]);
        
        sprintf(name, "../results/InducedFission_%s%s.txt", FC.c_str(), Setup.c_str());
        geFissionRate[i] = new TGraphErrors(name, Format(2*i+1).c_str());
        sprintf(name, "%s%s_FissionRate_%i", FC.c_str(), Setup.c_str(), i+1);
        geFissionRate[i]->SetName(name);
//        cout << "Saving " << FC << "/ToF/Signal/" << name << endl;
        Save(fAna, FC+"/ToF/Signal", geFissionRate[i]);
        
        geRatio[i] = Divide(geFission[i], geMonitor);
        sprintf(name, "%s%s_MonitorRatio_%i", FC.c_str(), Setup.c_str(), i+1);
        geRatio[i]->SetName(name);
//        cout << "Saving " << FC << "/CrossSection/" << name << endl;
        Save(fAna, FC+"/CrossSection", geRatio[i]);
    }

    // Uncorrected cross section
    TGraphErrors *geCrossSection[8];
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "../results/CrossSection_%s%s.txt", FC.c_str(), Setup.c_str());
        geCrossSection[i] = new TGraphErrors(name, Format(i).c_str());
        if (!geCrossSection[i])
            cout << "Could not open " << name << endl;
        sprintf(name, "%s%s_CS_raw_%i", FC.c_str(), Setup.c_str(), i+1);
        geCrossSection[i]->SetName(name);
//        cout << "Saving " << FC << "/CrossSection/" << name << endl;
        Save(fAna, FC+"/CrossSection", geCrossSection[i]);
    }

    fAna->Save();
    fAna->Close();
}//*/

void Runs()
{
    Runs("UFC");
    Runs("PuFC");
    Runs("UFC", "_FG");
    Runs("UFC", "_BG");
    Runs("PuFC", "_FG");
    Runs("PuFC", "_BG");
}
#endif
