#ifndef TRANSMISSION_H
#define TRANSMISSION_H

#include <fstream>
#include <sstream>
#include "TGraphErrors.h"
#include "SaveToFile.C"

TGraphErrors *GetCrossSection(string File)
{
//    cout << "Reading " << File << endl;

    // open cross section data as filestream
    ifstream ifs(File);
    if (!ifs) {
//        cout << "Error while opening " << File << endl;
        return 0;
    }
    Double_t A, B; Int_t N;
    ifs >> A >> B >> N;
//    cout << "First line of file: N = " << N << endl;

    // Prepare graph
    TGraphErrors *ge = new TGraphErrors();
    for (Int_t i = 0; i < N; i++)
    {
        ifs >> A >> B;
        ge->SetPoint(i, A, B);
//        cout << "Point " << i << ":  " << A << ", " << B << endl;
    }
    return ge;
}

TGraphErrors *AddGraphs(TGraphErrors *ge1, TGraphErrors *ge2, Double_t c1 = 1.0, Double_t c2 = 1.0)
{/// Add graphs' y values, use x values of graph with most points
    Double_t x, y;

    // handle nullpointers
    if (!ge1 && !ge2) return 0;
    if (!ge2) {
        for (uint i = 0; i < ge1->GetN(); i++)
        {
            ge1->GetPoint(i, x, y);
            ge1->SetPoint(i, x, c1 * y);
        }
        return ge1;
    }
    if (!ge1) {
        for (uint i = 0; i < ge1->GetN(); i++)
        {
            ge2->GetPoint(i, x, y);
            ge2->SetPoint(i, x, c2 * y);
        }
        return ge2;
    }

    // N: maximum number of points
    Int_t N = ge1->GetN() > ge2->GetN() ? ge1->GetN() : ge2->GetN();

    TGraphErrors *geSum = new TGraphErrors(N);

    if (ge1->GetN() > ge2->GetN()){
        // if ge1 has more points than ge2, use ge1 x values
        for (uint i = 0; i < N; i++)
        {
            ge1->GetPoint(i, x, y);
            y = c1 * y + c2 * ge2->Eval(x);
            geSum->SetPoint(i, x, y);
        }
    } else {
        // if ge2 has more points than ge1, use ge2 x values
        N = ge2->GetN();
        for (uint i = 0; i < N; i++)
        {
            ge2->GetPoint(i, x, y);
            y = c1 * ge1->Eval(x) + c2 * y;
            geSum->SetPoint(i, x, y);
        }
    }

    cout << ge1->Eval(1.E6) << " * " << c1 << " + " << ge2->Eval(1.E6) << " * " << c2 << " = " << geSum->Eval(1.E6) << endl;
    return geSum;
}

TGraphErrors *GetTotalCrossSection(string EvalDir, string Isotope)
{ // isotope name format: 94_242_Plutonium
    string Name;
//    string EvalDir = "/home/hoffma93/Programme/Geant4-Work/ENDF-VIII.0/";
//    cout << Isotope << " total cross section from " << EvalDir << endl;
    TGraphErrors *ge[] = {0, 0, 0, 0, 0};

    Name = EvalDir + "Inelastic/CrossSection/" + Isotope;
    ge[0] = GetCrossSection(Name);
    if (!ge[0]) {
//        cout << "Inelastic cross section for " << Isotope << " not found. Exit" << endl;
        return 0;
    }

    Name = EvalDir + "Elastic/CrossSection/" + Isotope;
    ge[1] = GetCrossSection(Name);
    if (!ge[1]) {
//        cout << "Elastic cross section for " << Isotope << " not found" << endl;
    }
    Name = EvalDir + "IsotopeProduction/CrossSection/" + Isotope;
    ge[2] = GetCrossSection(Name);
    if (!ge[2]) {
//        cout << "IsotopeProduction cross section for " << Isotope << " not found" << endl;
    }
    Name = EvalDir + "Fission/CrossSection/" + Isotope;
    ge[3] = GetCrossSection(Name);
    if (!ge[3]) {
//        cout << "Fission cross section for " << Isotope << " not found" << endl;
    }
    Name = EvalDir + "Capture/CrossSection/" + Isotope;
    ge[4] = GetCrossSection(Name);
    if (!ge[4]) {
//        cout << "Capture cross section for " << Isotope << " not found" << endl;
    }

    // find graph with most points
    Int_t N = 0;
    uint j = -1;
    for (uint graph_index = 0; graph_index < 5; graph_index++) {
        if (ge[graph_index] != 0) if (ge[graph_index]->GetN() > N) {
            j = graph_index;
            N = ge[graph_index]->GetN();
        }
    }
    if (j == -1) {
        cout << "Cross section data for " << Isotope << " not found in " << EvalDir << endl;
        return 0;
    }

    // Add graphs' y values, use x values of graph with most points
    TGraphErrors *geTot = new TGraphErrors(N);
    geTot->SetTitle((Isotope + " (n,tot)").c_str());
    Double_t x, y;
    for (uint i = 0; i < N; i++)
    {
        ge[j]->GetPoint(i, x, y);
        y = 0;
        for (uint i = 0; i < 5; i++) {
            if (ge[i])
                y += ge[i]->Eval(x);
        }
        geTot->SetPoint(i, x, y);
    }

    // DEBUG
//    cout << Isotope << "   " << EvalDir << endl
//         << "   Inelastic: " << ge[0]->Eval(1.5E7) << endl
//         << "   Elastic:   " << ge[1]->Eval(1.5E7) << endl
//         << "   Isotope:   " << ge[2]->Eval(1.5E7) << endl
//         << "   Fission:   " << ge[3]->Eval(1.5E7) << endl
//         << "   Capture:   " << ge[4]->Eval(1.5E7) << endl
//         << "   Total:     " << geTot->Eval(1.5E7) << endl;

    return geTot;
}

void CompareTransmission(string StrMatDim, string Evaluation1 = "/home/hoffma93/Programme/Geant4-Work/ENDF-VIII.0/", string Evaluation2 = "/home/hoffma93/Programme/Geant4-Work/ENDF-VII.1/")
{
    ifstream ifsMatDim(StrMatDim);
    string line = "";
    string StrIsotope = "";
    Double_t ArealDensity = 0;

    TGraphErrors *geMacroXS_1 = 0;
    TGraphErrors *geMacroXS_2 = 0;

    // Linewise read isotopes and areal number densities
    while(1) {
        getline(ifsMatDim, line);
        if (!ifsMatDim) break; // stop at end of file or in case of an error
        if (line == "") continue; // ignore empty lines
        if (line.substr(0, 1) == "#") continue; // ignore lines marked with a #

        // convert line string into isotope name and areal number density
        std::stringstream ss(line);
        ss >> StrIsotope >> ArealDensity;
        cout << StrIsotope << " " << ArealDensity << endl;

        // Get total cross section and add it to macroscopic cross sections
        TGraphErrors *geXS_1 = GetTotalCrossSection(Evaluation1, StrIsotope);
        geMacroXS_1 = AddGraphs(geMacroXS_1, geXS_1, 1.0, ArealDensity);
        TGraphErrors *geXS_2 = GetTotalCrossSection(Evaluation2, StrIsotope);
        geMacroXS_2 = AddGraphs(geMacroXS_2, geXS_2, 1.0, ArealDensity);
//        cout << geMacroXS_1->GetN() << " / " << geMacroXS_2->GetN() << endl;
    }

    Double_t x, macroXS, T; // energy, macroscopic cross section, transmission

    TGraphErrors *geTransmission1 = new TGraphErrors(geMacroXS_1->GetN());
    for (uint i = 0; i < geMacroXS_1->GetN(); i++)
    {
        geMacroXS_1->GetPoint(i,  x, macroXS);
        T = exp( - macroXS );
        geTransmission1->SetPoint(i, x, T);
    }
    TGraphErrors *geTransmission2 = new TGraphErrors(geMacroXS_2->GetN());
    for (uint i = 0; i < geMacroXS_2->GetN(); i++)
    {
        geMacroXS_2->GetPoint(i,  x, macroXS);
        T = exp( - macroXS );
        geTransmission2->SetPoint(i, x, T);
    }

    TCanvas *c2 = new TCanvas("c2");
    geTransmission1->Draw("APL");
    geTransmission2->SetLineColor(kRed);
    geTransmission2->Draw("sameLP");
    c2->Update();
    c2->Draw();

}

void Transmission(string Isotope = "94_242_Plutonium")
{
    cout << "Transmission..." << endl;
    string Evaluation1 = "/home/hoffma93/Programme/Geant4-Work/ENDF-VIII.0/";
    string Evaluation2 = "/home/hoffma93/Programme/Geant4-Work/ENDF-VII.1/";

    TGraphErrors *pGr = GetTotalCrossSection(Evaluation1, Isotope);
    TGraphErrors *g = GetTotalCrossSection(Evaluation2, Isotope);

    TCanvas *c1 = new TCanvas("c1");
//    pGr->SetMarkerStyle(20);
//    pGr->SetMarkerSize(2.0);
    pGr->Draw("AL");
//    g->SetMarkerStyle(20);
//    g->SetMarkerSize(2.0);
    g->SetLineColor(kRed);
    g->Draw("sameL");
    c1->Update();
    c1->Draw();
}



#endif
