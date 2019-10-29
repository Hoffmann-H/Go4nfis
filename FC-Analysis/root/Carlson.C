//#include "SaveToFile.C"
#include "NumberOfAtoms.C"

TGraph* Spectrum(string fac = "nELBE")
{
    string file_name;
    if (!strcmp(fac.c_str(), "nELBE"))
        file_name = "/home/hoffma93/Programme/ROOT/Data/nELBE_E.dat";
    else
        file_name = "/home/hoffma93/Programme/ROOT/Data/Source_E.dat";
    TGraph *g = new TGraphErrors(file_name.c_str());
    if (g == 0)
        cout << "Fehler beim Öffnen von " << file_name << endl;
    char name[64] = "";
    sprintf(name, "%s_En", fac.c_str());
    g->SetName(name);
    return g;
}

TGraphErrors* Result(Int_t isotope, string entry = "eta")
{
    // Read format: E eta(235) W(235) a2(235) eta(238) W(238) a2(238) DE Deta(235) DW(235) Da2(235) Deta(238) DW(238) Da2(238)
    char name[64] = "";
    string file_name = "/home/hoffma93/Experiment/Carlson-Korrektur/results/Carlson_nELBE.dat";
    string format;
    if (isotope == 235)
    {
        if (!strcmp(entry.c_str(), "eta"))
            format = "%lg %lg %*lg %*lg %*lg %*lg %*lg %lg %lg";
        if (!strcmp(entry.c_str(), "W"))
            format = "%lg %*lg %lg %*lg %*lg %*lg %*lg %lg %*lg %lg";
        if (!strcmp(entry.c_str(), "a2"))
            format = "%lg %*lg %*lg %lg %*lg %*lg %*lg %lg %*lg %*lg %lg";
    } else if (isotope == 238) {
        if (!strcmp(entry.c_str(), "eta"))
            format = "%lg %*lg %*lg %*lg %lg %*lg %*lg %lg %*lg %*lg %*lg %lg";
        if (!strcmp(entry.c_str(), "W"))
            format = "%lg %*lg %*lg %*lg %*lg %lg %*lg %lg %*lg %*lg %*lg %*lg %lg";
        if (!strcmp(entry.c_str(), "a2"))
            format = "%lg %*lg %*lg %*lg %*lg %*lg %lg %lg %*lg %*lg %*lg %*lg %*lg %lg";
    } else
        cout << "Unknown isotope " << isotope << endl;

    TGraphErrors *ge = new TGraphErrors(file_name.c_str(), format.c_str());
    if (ge == 0)
        cout << "Fehler beim Öffnen von " << file_name << " im Format " << format << endl;
    sprintf(name, "U%i_%s", isotope, entry.c_str());
    ge->SetName(name);
    return ge;
}

TGraphErrors* IvsE(Double_t t0, Double_t Dt0, Double_t R, Double_t DR, TGraphErrors *gEta, TGraphErrors *g_a2)
{ // Inefficiency dependent on neutron energy. Energy dependency given as TGraphErrors.

    Int_t N = g_a2->GetN();
    TGraphErrors *gI = new TGraphErrors(N);
    Double_t E, DE, Eta, DEta, a2, Da2, I, DI;
    for (Int_t j = 0;  j < N;  j++)
    {
        gEta->GetPoint(j, E, Eta);
        DEta = gEta->GetErrorY(j);
        g_a2->GetPoint(j, E, a2);
        Da2 = g_a2->GetErrorY(j);
        DE = g_a2->GetErrorX(j);
        I =            (t0 / R / 2 + Eta) * (1 - a2 / 2);
        DI = sqrt( pow(Dt0 / R / 2        * (1 - a2 / 2), 2) +
                   pow( t0 * DR/R/R / 2   * (1 - a2 / 2), 2) +
                   pow(              DEta * (1 - a2 / 2), 2) +
                   pow((t0 / R / 2 + Eta) *     Da2 / 2 , 2) );
//        cout << E << "+-" << DE << " " << t0 << " " << R << " " << I << "+-" << DI << endl;
        gI->SetPoint(j, E, I);
        gI->SetPointError(j, DE, DI);
    }
    return gI;
}

TGraphErrors* GetT(Double_t eff = 0.96)
{ // Return U3O8 layer thickness of Deposits. Unit: mg/cm^2
    Double_t M_U = 235.31765;
    Double_t M_O = 15.999;
    cout << endl << "Deposit depth [mg/cm^2]" << endl;
    TGraphErrors *geEffUmA = EffUmA(); // mg/cm^2
    TGraphErrors *geT = new TGraphErrors(8);
    Double_t x, UmA, DUmA, mA, DmA;
    for (Int_t i = 0; i < 8; i++)
    {
        geEffUmA->GetPoint(i, x, UmA);
        DUmA = geEffUmA->GetErrorY(i);
        mA = (1.0 + 8.0 / 3.0 * M_O / M_U) * UmA / eff;
        DmA = DUmA * mA / UmA;
        cout << i+1 << " " << mA << "+-" << DmA << endl;
        geT->SetPoint(i, i+1, mA);
        geT->SetPointError(i, 0, DmA);
    }
    geT->SetName("U3O8mA");
    geT->SetTitle("U3O8 areal density");
    return geT;
}

TH1F* GetR(Double_t &R, Double_t &DR, Int_t isotope = 235, string spectrum = "nELBE")
{ // Return fission fragment range distribution. Unit: mg/cm^2
    cout << endl << "Getting fission fragment range, U" << isotope << ", " << spectrum << endl;
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Experiment/Carlson-Korrektur/Geant4/G4_U%i_nf_%s_1e6.root", isotope, spectrum.c_str());
    TFile *f = TFile::Open(name);
    if (!f)
        cout << "Could not open " << name << endl;
    TH1F *hR = (TH1F*)f->Get("FF_Range");
    if (!hR)
        cout << "Could not get " << "TH1F FF_Range" << endl;
    sprintf(name, "U%i_%s_FF_Range", isotope, spectrum.c_str());
    hR->SetName(name);
    R = hR->GetMean(1);
    DR = hR->GetStdDev(1);
    cout << "<R> = " << R << ", stdev = " << DR << endl;
    return hR;
}

TGraphErrors* AddGraphs(TGraphErrors *g1, Double_t w1, Double_t Dw1, TGraphErrors *g2, Double_t w2, Double_t Dw2)
{
    Int_t N = g1->GetN();
    if (g2->GetN() != N)
        cout << "Can't add " << g1->GetName() << " and " << g2->GetName() << ": Different nummber of points" << endl;
    TGraphErrors *gS = new TGraphErrors(N);
    Double_t x1, x2, y1, y2, yS, Dx, Dy1, Dy2, DyS;
    for (Int_t j = 0; j < N; j++)
    {
        g1->GetPoint(j, x1, y1);
        g2->GetPoint(j, x2, y2);
        if (x2 != x1)
            cout << "Warning: Point " << j << ", different x values: " << x1 << ", " << x2 << endl;
        Dx = g1->GetErrorY(j);
        Dy1 = g1->GetErrorY(j);
        Dy2 = g2->GetErrorY(j);
        yS = w1 * y1 + w2 * y2;
        DyS = sqrt( pow(Dw1 * y1, 2) + pow(w1 * Dy1, 2) + pow(Dw2 * y2, 2) + pow(w2 * Dy2, 2) );
        gS->SetPoint(j, x1, yS);
        gS->SetPointError(j, Dx, DyS);
    }
    return gS;
}

Double_t WeightedIntegral(TGraph *weight, TGraph *f)
{ // assumption: weighting function has closer data points
    Double_t x1, x2, x3, w, val, sum = 0, sum_val = 0;
    Double_t N = weight->GetN();
    for (Int_t j = 1; j < N-1; j++)
    {
        weight->GetPoint(j-1, x1, w);
        weight->GetPoint(j+1, x3, w);
        weight->GetPoint(j, x2, w);

        val = f->Eval(x2);

        sum_val += w * 0.5 * (x3 - x1) * val;
        sum += w * 0.5 * (x3 - x1);
//        cout << sum_val << " " << sum << endl;
    }
    return sum_val / sum;
}

void Carlson()
{ // Calculate the efficiency correction (eff_nELBE / eff_15MeV)
    char name[64] = "";

    // Get layer thickness (8 deposits)
    Double_t approx_eff = 0.96;
    TGraphErrors *g_t0 = GetT(approx_eff);

    // Get Fission Fragment Ranges (4 scenarios)
    Double_t R[4], DR[4];
    TH1F *hR[4];
    hR[0] = GetR(R[0], DR[0], 235, "nELBE");
    hR[1] = GetR(R[1], DR[1], 235, "15MeV");
    hR[2] = GetR(R[2], DR[2], 238, "nELBE");
    hR[3] = GetR(R[3], DR[3], 238, "15MeV");

    // Get anisotropy data
    TGraphErrors *geU235eta = Result(235, "eta");
    TGraphErrors *geU235W   = Result(235, "W");
    TGraphErrors *geU235a2  = Result(235, "a2");
    TGraphErrors *geU238eta = Result(238, "eta");
    TGraphErrors *geU238W   = Result(238, "W");
    TGraphErrors *geU238a2  = Result(238, "a2");

    // Fission weighting parameters
    Double_t iso[] = {0.901989, 0.0912};
    Double_t sigma[] = {2.09, 1.2};
    Double_t Dsigma[] = {0.14, 0.54};
    Double_t w[2], Dw[2];
    w[0] = iso[0] * sigma[0] / (iso[0] * sigma[0] + iso[1] * sigma[1]);
    w[1] = iso[1] * sigma[1] / (iso[0] * sigma[0] + iso[1] * sigma[1]);
    Dw[0] = iso[0] * Dsigma[0] / (iso[0] * sigma[0] + iso[1] * sigma[1]); // uncertainty approximately
    Dw[1] = iso[1] * Dsigma[1] / (iso[0] * sigma[0] + iso[1] * sigma[1]);

    // Neutron spectra
    TGraph *g_nELBE = Spectrum("nELBE");
    TGraph *g_PTB = Spectrum("PTB");

    // Save
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Save(fAna, "Carlson/U235/nELBE", hR[0]);
    Save(fAna, "Carlson/U235/15MeV", hR[1]);
    Save(fAna, "Carlson/U238/nELBE", hR[2]);
    Save(fAna, "Carlson/U238/15MeV", hR[3]);
    Save(fAna, "Carlson/U235", geU235eta);
    Save(fAna, "Carlson/U235", geU235W);
    Save(fAna, "Carlson/U235", geU235a2);
    Save(fAna, "Carlson/U238", geU238eta);
    Save(fAna, "Carlson/U238", geU238W);
    Save(fAna, "Carlson/U238", geU238a2);
    Save(fAna, "Carlson/nELBE", g_nELBE);
    Save(fAna, "Carlson/15MeV", g_PTB);

    /// Calculate Inefficiencies, loop over deposits
    TGraphErrors *geI[6][8]; // 1st index: 235 nElbe; 235 15MeV; 238 nElbe; 238 15MeV; weighted nELBE; weighted 15MeV. 2nd index: Deposit.
    Double_t x, t0, Dt0, eff[2], Deff[2];
    TGraphErrors *gEff_nELBE = new TGraphErrors(8);
    gEff_nELBE->SetName("eff_nELBE");
    gEff_nELBE->SetTitle("Carlson efficiency nELBE");
    TGraphErrors *gEff_15MeV = new TGraphErrors(8);
    gEff_15MeV->SetName("eff_15MeV");
    gEff_15MeV->SetTitle("Carlson efficiency 15MeV");
    TGraphErrors *gC = new TGraphErrors(8);
    gC->SetName("Carlson_Correction");
    gC->SetTitle("Carlson efficiency correction factor");
    for (Int_t i = 0; i < 8; i++)
    {
        // Thickness
        g_t0->GetPoint(i, x, t0);
        Dt0 = g_t0->GetErrorY(i);

        // Inefficiency 235/238, nELBE/15MeV
        geI[0][i] = IvsE(t0, Dt0, R[0], DR[0], geU235eta, geU235a2);
        sprintf(name, "I_U235_nELBE_%i", i+1);
        geI[0][i]->SetName(name);
        geI[1][i] = IvsE(t0, Dt0, R[1], DR[1], geU235eta, geU235a2);
        sprintf(name, "I_U235_15MeV_%i", i+1);
        geI[1][i]->SetName(name);
        geI[2][i] = IvsE(t0, Dt0, R[2], DR[2], geU238eta, geU238a2);
        sprintf(name, "I_U238_nELBE_%i", i+1);
        geI[2][i]->SetName(name);
        geI[3][i] = IvsE(t0, Dt0, R[3], DR[3], geU238eta, geU238a2);
        sprintf(name, "I_U238_15MeV_%i", i+1);
        geI[3][i]->SetName(name);

        // Weighted Inefficiency
        geI[4][i] = AddGraphs(geI[0][i], w[0], Dw[0], geI[2][i], w[1], Dw[1]);
        sprintf(name, "I_nELBE_%i", i+1);
        geI[4][i]->SetName(name);
        geI[5][i] = AddGraphs(geI[1][i], w[0], Dw[0], geI[3][i], w[1], Dw[1]);
        sprintf(name, "I_15MeV_%i", i+1);
        geI[5][i]->SetName(name);

        // efficiencies for spectra
        eff[0] = 1.0 - WeightedIntegral(g_nELBE, geI[4][i]);
        gEff_nELBE->SetPoint(i, i+1, eff[0]);
        eff[1] = 1.0 - WeightedIntegral(g_PTB, geI[4][i]);
        gEff_15MeV->SetPoint(i, i+1, eff[1]);
        gC->SetPoint(i, i+1, eff[0] / eff[1]);

        // Save
        Save(fAna, "Carlson/U235/nELBE", geI[0][i]);
        Save(fAna, "Carlson/U235/15MeV", geI[1][i]);
        Save(fAna, "Carlson/U238/nELBE", geI[2][i]);
        Save(fAna, "Carlson/U238/15MeV", geI[3][i]);
        Save(fAna, "Carlson/nELBE", geI[4][i]);
        Save(fAna, "Carlson/15MeV", geI[5][i]);

    }//*/
//    new TCanvas();
//    geI[0][0]->SetLineColor(kRed);
//    geI[0][0]->Draw();
//    geI[2][0]->SetLineColor(kGreen);
//    geI[2][0]->Draw("same");
//    geI[4][0]->Draw("same");
//    Save(fAna, "Carlson/")
    Save(fAna, "Carlson/nELBE", gEff_nELBE);
    Save(fAna, "Carlson/15MeV", gEff_15MeV);
    Save(fAna, "Carlson", gC);
    fAna->Save();
    fAna->Close();
}


void f()
{
    TGraph *g = Spectrum("nELBE");
//    TFile *f = TFile::Open("/home/hoffma93/Experiment/Carlson-Korrektur/results/Carlson.root", "UPDATE");

}
