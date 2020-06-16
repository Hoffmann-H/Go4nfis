#ifndef TRANSMISSION_H
#define TRANSMISSION_H

#include <fstream>
#include "TGraphErrors.h"
#include "SaveToFile.C"

TGraphErrors *GetCrossSection(string File)
{
    cout << "Reading " << File << endl;
//    TGraphErrors *ge = new TGraphErrors((EvalDir+Name).c_str(), "%*lg %*lg %lg %lg");
//    if (!ge) {
//        cout << "Error while opening " << File << endl;
//        return 0;
//    }
//    return ge;

    // open cross section data as filestream
    ifstream ifs(File);
    if (!ifs) {
        cout << "Error while opening " << File << endl;
        return 0;
    }
    Double_t A, B; Int_t N;
    ifs >> A >> B >> N;
    cout << "First line of file: N = " << N << endl;

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

    Name = EvalDir + "Elastic/CrossSection/" + Isotope;
    TGraphErrors *geElastic = GetCrossSection(Name);
    if (!geElastic) {
        cout << "Elastic cross section for " << Isotope << " not found" << endl;
        return 0;
    }

    Name = EvalDir + "Inelastic/CrossSection/" + Isotope;
    TGraphErrors *geInelastic = GetCrossSection(Name);

    TGraphErrors *geTot = new TGraphErrors(geElastic->GetN());
    Double_t x, y;
    for (uint i = 0; i < geElastic->GetN(); i++)
    {
        geInelastic->GetPoint(i, x, y);
        geTot->SetPoint(i, x, y + geElastic->Eval(x));
    }

    return geTot;
}

void Transmission()
{
    cout << "Transmission..." << endl;
    TGraphErrors *pGr = GetTotalCrossSection("94_242_Plutonium");

    TCanvas *c1 = new TCanvas("c1");
    pGr->SetMarkerStyle(20);
    pGr->SetMarkerSize(2.0);
    pGr->Draw("APL");
    c1->Update();
    c1->Draw();
}



#endif
