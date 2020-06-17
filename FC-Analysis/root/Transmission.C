#ifndef TRANSMISSION_H
#define TRANSMISSION_H

#include <fstream>
#include "TGraphErrors.h"
#include "SaveToFile.C"

TGraphErrors *GetCrossSection(string File)
{
    cout << "Reading " << File << endl;

    // open cross section data as filestream
    ifstream ifs(File);
    if (!ifs) {
        cout << "Error while opening " << File << endl;
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

TGraphErrors *GetTotalCrossSection(string Isotope)
{ // isotope name format: 94_242_Plutonium
    string Name;
    string EvalDir = "/home/hoffma93/Programme/Geant4-Work/ENDF-VIII.0/";
    TGraphErrors *geInelastic, *geElastic, *geIsotopeProduction, *geFission, *geCapture;

    Name = EvalDir + "Inelastic/CrossSection/" + Isotope;
    geInelastic = GetCrossSection(Name);
    if (!geInelastic) {
        cout << "Inelastic cross section for " << Isotope << " not found. Exit" << endl;
        return 0;
    }

    Name = EvalDir + "Elastic/CrossSection/" + Isotope;
    geElastic = GetCrossSection(Name);
    if (!geElastic) {
        cout << "Elastic cross section for " << Isotope << " not found" << endl;
        return 0;
    }
//    Name = EvalDir + "IsotopeProduction/CrossSection/" + Isotope;
//    geIsotopeProduction = GetCrossSection(Name);
//    if (!geIsotopeProduction) {
//        cout << "IsotopeProduction cross section for " << Isotope << " not found" << endl;
//        return 0;
//    }
    Name = EvalDir + "Fission/CrossSection/" + Isotope;
    geFission = GetCrossSection(Name);
    if (!geFission) {
        cout << "Fission cross section for " << Isotope << " not found" << endl;
        return 0;
    }
    Name = EvalDir + "Capture/CrossSection/" + Isotope;
    geCapture = GetCrossSection(Name);
    if (!geCapture) {
        cout << "Capture cross section for " << Isotope << " not found" << endl;
        return 0;
    }


    TGraphErrors *geTot = new TGraphErrors(geElastic->GetN());
    geTot->SetTitle((Isotope + " (n,tot)").c_str());
    Double_t x, y0, y1, y2, y3, y4;
    for (uint i = 0; i < geElastic->GetN(); i++)
    {
        if (geInelastic) geInelastic->GetPoint(i, x, y0);
        else y0 = 0;
        if (geElastic) y1 = geElastic->Eval(x);
        else y1 = 0;
//        if (geIsotopeProduction) y2 = geIsotopeProduction->Eval(x);
//        else y2 = 0;
        if (geFission) y3 = geFission->Eval(x);
        else y3 = 0;
        if (geCapture) y4 = geCapture->Eval(x);
        else y4 = 0;
        geTot->SetPoint(i, x, y0 + y1 + y2 + y3 + y4);
    }

    return geTot;
}

void Transmission()
{
    cout << "Transmission..." << endl;
    TGraphErrors *pGr = GetTotalCrossSection("94_242_Plutonium");

    TGraphErrors *g = new TGraphErrors("/home/hoffma93/Programme/Go4nfis/FC-Analysis/data/Pu-242(n,tot).dat");

    TCanvas *c1 = new TCanvas("c1");
    pGr->SetMarkerStyle(20);
    pGr->SetMarkerSize(2.0);
    pGr->Draw("AP");
    g->SetMarkerStyle(20);
    g->SetMarkerSize(2.0);
    g->SetMarkerColor(kRed);
    g->Draw("sameP");
    c1->Update();
    c1->Draw();
}



#endif
