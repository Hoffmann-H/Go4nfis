#include "SaveToFile.C"
#include "FC.C"
#include <fstream>
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

void GetBinning(string file_to_read)
{
    /// Print binning used in tabular
    TGraph *gTE = new TGraph(file_to_read.c_str(), "%lg %lg %*lg %*lg");
    Int_t N = gTE->GetN();
    Double_t Tmin, Emin, T, E, Tmax, Emax;
    Int_t Tbins, Ebins;
    gTE->GetPoint(0, Tmin, Emin);
    gTE->GetPoint(N-1, Tmax, Emax);
    Int_t i = 0;
    do {
        i++;
        gTE->GetPoint(i, T, E);
    } while (T == Tmin);
    Tbins = (Tmax - Tmin) / T;
    Ebins = i - 1;
    cout << N << " entries" << endl;
    cout << "T: " << Tbins << " steps from " << Tmin << " to " << Tmax << endl;
    cout << "E: " << Ebins << " steps from " << Emin << " to " << Emax << endl;
}

TH2D* MakeEvsT(string file_to_read, Bool_t draw, string FC = "PuFC", string key = "real", Int_t ch = 0)
{
    char name[64] = "";
    Int_t N;

        // Open tabular as graph
        TGraphErrors *gTE = new TGraphErrors(file_to_read.c_str(), "%lg %lg %lg %lg");
        N = gTE->GetN();

        // Open tabular as filestream
//        std::ifstream input(file_to_read.c_str());
//        string line;
//        for (Int_t i = 0; i < 7; i++)
//            getline(input, line);
//        N = 100701;

//    cout << N << " entries" << endl;
    Double_t Tmin, Tmax, Emin, Emax;
    Int_t Tbins, Ebins;
    Tmin = 0.0, Tmax = 200.0, Emin = 0.0, Emax = 16.0;
    Tbins = 2000, Ebins = 1600;
    sprintf(name, "%s_ToFvsEkin_%s_Ch.%i", FC.c_str(), key.c_str(), ch + 1);
    TH2D *pHist = new TH2D(name, name,
                           Tbins, Tmin, Tmax, Ebins, Emin, Emax);
    sprintf(name, "%s, ToF vs Ekin, %s, ch.%i", FC.c_str(), key.c_str(), ch + 1);
    pHist->SetTitle(name);
    pHist->GetXaxis()->SetTitle("#font[12]{t} / ns");
    pHist->GetYaxis()->SetTitle("#font[12]{E} / MeV");

    Double_t T, E, C, DC;
    Int_t binT, binE;

    for (Int_t i = 0; i < N; i++)
    {
            // read from graph
            gTE->GetPoint(i, T, E);
            C = gTE->GetErrorX(i);
            DC = gTE->GetErrorY(i);

            // read from filestream
//            input >> T >> E >> C >> DC;

        // skip emty entries
        if (C == 0.0)
            continue;

        T *= 10; // unit: 1 ns
//        if (E > 2.2 && E < 2.4 && T >= 30 && T < 45)
//            cout << T << "   " << E << "   " << C << "   " << DC << "   " << binT << "   " << binE << endl;
        binT = pHist->GetXaxis()->FindBin(T - 0.001);
        binE = pHist->GetYaxis()->FindBin(E - 0.001);

        pHist->SetBinContent(binT, binE, C);
        pHist->SetBinError(binT, binE, DC * C);
    }
    if (draw)
    {
        sprintf(name, "%s_EvsT_%i", FC.c_str(), ch + 1);
        TCanvas *c1 = new TCanvas(name, name, 200, 10, 700, 500);
        pHist->SetStats(0);
        pHist->Draw("colz");
    }
//    cout << pHist->Integral() << endl;
    return pHist;
}

TH1D* TimeProjection(TH2D *pH2, Int_t ch, string FC = "PuFC", string key = "real")
{
    char name[64] = "";
    sprintf(name, "%s_ProjT_%s_%i", FC.c_str(), key.c_str(), ch + 1);
    TH1D *h = (TH1D*)pH2->ProjectionX(name);
//    h->SetName(name);
    sprintf(name, "%s, ToF, %s, ch.%i", FC.c_str(), key.c_str(), ch + 1);
    h->SetTitle(name);
    h->GetXaxis()->SetTitle("#font[12]{t} / ns");
    h->GetYaxis()->SetTitle("#font[12]{C}_{(n,f)}");
    return h;
}

TH1D *h = new TH1D();
Double_t f(Double_t *x, Double_t *par)
{ // return: histogram h, folded by gaussian (p0 ampl, p1 mean, p2 width), evaluated at x.
    Double_t s = 0;
    for (Int_t j = 1; j < h->GetNbinsX()+1; j++)
    {
        s += h->GetBinContent(j) * par[0] * exp(-pow((x[0] - h->GetBinCenter(j) - par[1]) / par[2], 2));
    }
    return s + par[3];
}

TH1D* FitPeakForm(TH1D *H, string Run = "NIF", Int_t ch = 0, string Simulation = "Geant4", string key = "real")
{ // Fit simulated ToF folding to data
    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
    char name[128] = "";

    // Experimental spectrum with Background
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *fExp = TFile::Open(name); if (!fExp) cout << "Could not open " << name << endl;
    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1D *hExp = (TH1D*)fExp->Get(name); if (!hExp) cout << "Could not get " << name << endl;
    hExp->Sumw2();
    Double_t bg = (hExp->Integral(Gate_0(ch, FC), Gate_a(ch, FC)) + hExp->Integral(Gate_b(ch, FC), Gate_3(ch, FC))) / (Gate_a(ch, FC) - Gate_0(ch, FC) + Gate_3(ch, FC) - Gate_b(ch, FC));
    // Background-subtracted
//    sprintf(name, "%s/ToF/Signal/%s/H1AnaHZDRDtG_%i", FC.c_str(), Run.c_str(), ch+1);
//    TH1D *hExp = (TH1D*)fAna->Get(name); if (!hExp) cout << "Could not get " << name << endl;
//    Double_t bg = 0;

    h = (TH1D*)H->Clone();
    Double_t Range[] = {hExp->GetBinLowEdge(Gate_0(ch, FC)), hExp->GetBinLowEdge(Gate_3(ch, FC))};
    Double_t Par[] = {(hExp->GetMaximum() - bg) / h->Integral(), 0, 3, bg};
//    cout << Par[0] << " " << Par[1] << " " << Par[2] << " " << Par[3] << endl;
    sprintf(name, "f%s_%i", Run.c_str(), ch+1);
    TF1 *fit = new TF1(name, f, Range[0], Range[1], 4);
    fit->SetParameters(Par[0], Par[1], Par[2], Par[3]);
    fit->FixParameter(3, bg);
//    fit->FixParameter(1, 0);

    /// Fit /////////////////////////////////////
    hExp->Fit(name, "LR0Q"); ////////////////////
    /////////////////////////////////////////////

//    for (Int_t p = 0; p < 4; p++)
//        cout << fit->GetParameter(p) << " +- " << fit->GetParError(p) << endl;
    cout << FC << "\t" << ch+1 << " \t" << fit->GetParameter(2) << " +- " << fit->GetParError(2) << endl;

    // Save
    Int_t N = h->GetNbinsX();
    Double_t xmin = h->GetBinLowEdge(1);
    Double_t xmax = h->GetBinLowEdge(N+1);
    sprintf(name, "%s_FitT_%s_%i", FC.c_str(), key.c_str(), ch+1);
    TH1D* hFit = new TH1D(name, name, N, xmin, xmax);
    sprintf(name, "%s, Ch. %i, Fit simulated ToF spectrum", FC.c_str(), ch+1);
    hFit->SetTitle(name);
    for (Int_t bin = 1; bin < N+1; bin++)
        hFit->SetBinContent(bin, fit->Eval(hFit->GetBinCenter(bin)));
//    new TCanvas();
//    hExp->Draw("p");
//    hFit->Draw("same");
    return hFit;
}

void MCNPtoROOT(Bool_t save, string Run = "NIF", string key = "real", string Path = "FCscat_PTB_weight/tally/FCscat_b", Long_t SimulatedN = 60000000000)
{ // Convert MCNP Track Length results to root
    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering";
    TH2D *h[8];

    TFile *fAna;
    TCanvas *c1;
    if (save)
        fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    else { // draw
        c1 = new TCanvas();
        c1->Divide(4, 2);
    }

    // Loop over files
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/%s_tally-2%i4_xyz_0.dmp", DirName.c_str(), Path.c_str(), 8 - i);

        // Create 2D histogram
        h[i] = MakeEvsT(name, 0, "PuFC", key, i);

        if (save)
        {
            Save(fAna, "Simulation/MCNP/"+FC+"_"+key+"/ToFvsEkin", h[i]);
            TH1D *hProj = TimeProjection(h[i], i, FC, key);
            Save(fAna, "Simulation/MCNP/"+FC+"_"+key+"/EffToF", hProj);
            if (strcmp(key.c_str(), "ideal")) {
                TH1D *hFit = FitPeakForm(hProj, Run, i, "MCNP", key);
                Save(fAna, "Simulation/MCNP/"+FC+"_"+key+"/FitToF", hFit);
            }
        } else { // draw
            c1->cd(i + 1);
            gPad->SetLogz(1);
            h[i]->SetStats(0);
            h[i]->GetXaxis()->SetRangeUser(0, 100);
            h[i]->GetZaxis()->SetRangeUser(0.0000000001, 0.0000001);
            h[i]->Draw("colz");
        }
    }
    if (save)
    {
//        TGraphErrors *geEmit = GetEmitMCNP(SimulatedN, FC);
//        Save(fAna, "Simulation/MCNP/Correction", geEmit, "Emit");
        fAna->Save();
        fAna->Close();
    } else {
        c1->Modified();
        c1->Update();
    }
}

void Geant4toROOT(string FileName, string Run = "NIF", string key = "real")
{
    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
    string FilePath = "/home/hoffma93/Programme/Geant4-Work/results";
    char name[128] = "";
    sprintf(name, "%s/%s", FilePath.c_str(), FileName.c_str());
    TFile *fG4 = TFile::Open(name, "READ");
    if (fG4 == 0)
        cout << "Could not open " << name << endl;
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

    // Emitted neutrons
    TH1F *pHemit = (TH1F*)fG4->Get("Source/Source_Theta");
    Double_t maxTheta = 0.08446522295019329 * 180 / TMath::Pi();
    Int_t lastBin = pHemit->FindBin(maxTheta);
    Double_t WeightLastBin = (maxTheta - pHemit->GetBinLowEdge(lastBin)) / pHemit->GetBinWidth(lastBin);
    Double_t nEmit = pHemit->Integral(0, lastBin - 1) + WeightLastBin * pHemit->GetBinContent(lastBin);
    cout << "Smeared angular distribution correction: " << nEmit / pHemit->Integral() << endl;
//    TGraphErrors *geEmit = new TGraphErrors(8);
//    for (Int_t i = 0; i < 8; i++)
//    {
//        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
//        pHemit = (TH1F*)fG4->Get(name);
//        geEmit->SetPoint(i, i+1, pHemit->Integral() / nEmit);
//        geEmit->SetPointError(i, 0, 0);
//    }
//    Save(fAna, "Simulation/Geant4/"+FC+"_"+Setup+"/Correction", geEmit, "Emit");

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH2D *pH2Tot = (TH2D*)fG4->Get(name);
        sprintf(name, "%s_ToFvsEkin_%s_Ch.%i", FC.c_str(), key.c_str(), i+1);
        pH2Tot->SetName(name);
//        cout << pH2Tot->Integral() << " ";
        if (pH2Tot->Integral() > 1)
            pH2Tot->Scale(1.0 / nEmit);
//        cout << pH2Tot->Integral() << endl;

        TH1D *hProj = TimeProjection(pH2Tot, i, FC, key);
        Save(fAna, "Simulation/Geant4/"+FC+"_"+key+"/EffToF", hProj);
        if (strcmp(key.c_str(), "ideal")) {
            TH1D *hFit = FitPeakForm(hProj, Run, i, "Geant4", key);
            Save(fAna, "Simulation/Geant4/"+FC+"_"+key+"/FitToF", hFit);
        }

        sprintf(name, "%s/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH2D *pH2Sc = (TH2D*)fG4->Get(name);
//        , FC+"_ToFvsEkin_"+key+"_Sc_Ch."+to_string(i+1)
        sprintf(name, "%s_ToFvsEkin_Sc_%s_Ch.%i", FC.c_str(), key.c_str(), i+1);
        if (pH2Sc->Integral() > 1)
            pH2Sc->Scale(1.0 / nEmit);

        Save(fAna, "Simulation/Geant4/"+FC+"_"+key+"/ToFvsEkin/Scattered", pH2Sc);
        Save(fAna, "Simulation/Geant4/"+FC+"_"+key+"/ToFvsEkin", pH2Tot);
    }
    fAna->Save();
    fAna->Close();
}

void MCNPtoROOT()
{
    Geant4toROOT("PuFC_real_c_5E7.root", "NIF", "real");
    Geant4toROOT("PuFC_ideal_c_5E7.root", "PuFC", "ideal");
    Geant4toROOT("UFC_real_c_5E7.root", "UFC_NIF", "real");
    Geant4toROOT("UFC_ideal_c_5E7.root", "UFC", "ideal");
    Geant4toROOT("UFC_SB_c_5E7.root", "UFC_SB", "SB");
    MCNPtoROOT(1, "NIF", "real", "FCscat_PTB_weight/tally/FCscat_b");
    MCNPtoROOT(1, "PuFC", "ideal", "FCscat_PTB_weight_void/tally/FCscat_d");
}

void unused_functions()
{
//string NameTag(string result_key = "2")
//{
//    if (result_key == "2")
//        return "";
//    if (result_key == "1")
//        return "Sc_";
//    if (result_key == "0")
//        return "Dir_";
//    cout << "Unknown result key " << result_key << endl;
//    return 0;
//}

//string PathTag(string result_key = "2")
//{
//    if (result_key == "2")
//        return "";
//    if (result_key == "1")
//        return "Scattered/";
//    if (result_key == "0")
//        return "Direct/";
//    cout << "Unknown result key " << result_key << endl;
//    return 0;
//}

//TGraphErrors* GetEmitMCNP(Double_t SimulatedN, string FC = "PuFC")
//{
//    Double_t SourceMaxAngle = atan(12.7/150);
//    TGraphErrors *ge = new TGraphErrors(8);
//    ge->SetTitle("Projectiles; Deposit; Neutrons");
//    for (Int_t i = 0; i < 8; i++)
//    {
//        Double_t SimulatedDistance = Distance(i, FC) - 2.0; //
//        Double_t DepositAngle = atan(DepositRadius(i, FC) / SimulatedDistance);
//        Double_t Emit = SolidAngle(DepositAngle) / SolidAngle(SourceMaxAngle);
//        Double_t DEmit = sqrt(Emit / SimulatedN);
//        ge->SetPoint(i, i+1, Emit);
//        ge->SetPointError(i, 0, DEmit);
////        cout << Emit << "+-" << DEmit << endl;
//    }
//    return ge;
//}

//void TraLenMCNPtoROOT(string result_key = "2", Bool_t save = 1, Bool_t draw = 0, string FC = "PuFC", Long_t SimulatedN = 60000000000)
//{ // Convert Track Length MCNP results to root
//    char name[128] = "";
//    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
//    TH2D *h[8];

//    TFile *fAna;
//    if (save)
//        fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

//    TCanvas *c1;
//    if (draw)
//    {
//        c1 = new TCanvas();
//        c1->Divide(4, 2);
//    }

//    // Loop over files
//    for (Int_t i = 0; i < 8; i++)
//    {
//        sprintf(name, "%s/FCscat_g_tally-2%i4_xyz_0_%s.dmp", DirName.c_str(), 8 - i, result_key.c_str());

//        // Create 2D histogram
//        h[i] = MakeEvsT(name, 0, FC, result_key, i);
//        // Track length estimator gives neutrons per started neutron and cm^2
//        Double_t DepositArea = TMath::Pi() * pow(3.7, 2);
//        h[i]->Scale(DepositArea);

//        if (save)
//        {
////            cout << "Saving " << "Simulation/MCNP/ToFvsEkin/"<<PathTag(result_key) << endl;
//            Save(fAna, "Simulation/MCNP/ToFvsEkin/"+PathTag(result_key), h[i]);
//        }
//        if (draw)
//        {
//            c1->cd(i + 1);
////            gPad->SetLogz(1);
//            h[i]->SetStats(0);
////            h[i]->GetXaxis()->SetRangeUser(0, 100);
//            h[i]->Draw("colz");
//        }
//    }
//    if (draw)
//    {
//        c1->Modified();
//        c1->Update();
//    }
//    if (save)
//    {
//        TGraphErrors *geEmit = GetEmitMCNP(SimulatedN, FC);
//        Save(fAna, "Simulation/MCNP/Correction", geEmit, "Emit");
//        fAna->Save();
//        fAna->Close();
//    }
//}

//void SimpleMCNPtoROOT(string result_key = "2", Bool_t save = 1, Bool_t draw = 0, string FC = "PuFC", Long_t SimulatedN = 6000000000)
//{ // Convert simple MCNP results to root
//    char name[128] = "";
//    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
//    TH2D *h[8];

//    TFile *fAna;
//    if (save)
//        fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");


//    TCanvas *c1;
//    if (draw)
//    {
//        c1 = new TCanvas();
//        c1->Divide(4, 2);
//    }

//    // Loop over files
//    for (Int_t i = 0; i < 8; i++)
//    {
//        sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_%s.dmp", DirName.c_str(), 8 - i, result_key.c_str());
//        string front_key, back_key;

//        // Create 2D histogram
//        if (result_key == "2")
//        { // total
//            front_key = "9";
//            back_key = "8";
//        }
//        if (result_key == "1")
//        { // scattered
//            front_key = "5";
//            back_key = "4";
//        }
//        if (result_key == "0")
//        { // direct
//            front_key = "1";
//            back_key = "0";
//        }
//        // Open 2D front
//        sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_%s.dmp", DirName.c_str(), 8-i, front_key.c_str());
//        TH2D *pH2front = MakeEvsT(name, 0, FC, result_key, i, 0);
//        // Change name to prevent overwriting
//        sprintf(name, "%s_EvsT_front_%i", FC.c_str(), i+1);
//        pH2front->SetName(name);
//        // Open 2D back
//        sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_%s.dmp", DirName.c_str(), 8-i, back_key.c_str());
//        h[i] = MakeEvsT(name, 0, FC, result_key, i, 0);
//        // Add front and back
//        h[i]->Add(pH2front);

//        if (save)
//        {
////            cout << "Saving " << "Simulation/MCNP_simple/ToFvsEkin/"<<PathTag(result_key) << endl;
//            Save(fAna, "Simulation/MCNP_simple/ToFvsEkin/"+PathTag(result_key), h[i]);
//        }
//        if (draw)
//        {
//            c1->cd(i + 1);
//            h[i]->SetStats(0);
//            h[i]->GetXaxis()->SetRangeUser(0, 100);
//            h[i]->Draw("colz");
//        }
//    }

//    if (draw)
//    {
//        c1->Modified();
//        c1->Update();
//    }
//    if (save)
//    {
//        TGraphErrors *geEmit = GetEmitMCNP(SimulatedN, FC);
//        Save(fAna, "Simulation/MCNP_simple/Correction", geEmit, "Emit");
//        fAna->Save();
//        fAna->Close();
//    }
//}//*/

//TH2D* TH2FtoTH2D(TH2F *h1, string FC, string key, Int_t ch)
//{
//    const char *name = h1->GetName();
//    const char *title = h1->GetTitle();
//    string replace_name = FC+"_"+key+"_"+to_string(ch+1);
//    h1->SetName(replace_name.c_str());
//    if (h1->GetNbinsX() != 2000)
//        cout << replace_name << " bad binning " << h1->GetNbinsX() << endl;
//    if (h1->GetNbinsY() != 1600)
//        cout << replace_name << " bad binning " << h1->GetNbinsY() << endl;
//    TH2D *h2 = new TH2D(name, title, 2000, 0, 200, 1600, 0, 16);
//    h2->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
//    h2->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
//    h2->Add(h1);
//    return h2;
//}
}
