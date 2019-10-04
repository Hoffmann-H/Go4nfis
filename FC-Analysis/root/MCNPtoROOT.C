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

string NameTag(string result_key = "2")
{
    if (result_key == "2")
        return "";
    if (result_key == "1")
        return "Sc_";
    if (result_key == "0")
        return "Dir_";
    cout << "Unknown result key " << result_key << endl;
    return 0;
}

string PathTag(string result_key = "2")
{
    if (result_key == "2")
        return "";
    if (result_key == "1")
        return "Scattered/";
    if (result_key == "0")
        return "Direct/";
    cout << "Unknown result key " << result_key << endl;
    return 0;
}

TH2D* MakeEvsT(string file_to_read, Bool_t draw, string FC = "PuFC", string key = "2", Int_t ch = 0, Bool_t track_length = 1)
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
    if (track_length)
    { // Binning for MCNP run with track length estimator
        Tmin = 0.0, Tmax = 200.0, Emin = 0.0, Emax = 16.0;
        Tbins = 2000, Ebins = 1600;
    } else { // Binning for simple MCNP run
        Tmin = 0.0, Tmax = 500.0, Emin = 0.0, Emax = 20.0;
        Tbins = 500, Ebins = 200;
    }
    sprintf(name, "%s_ToFvsEkin_%sCh.%i", FC.c_str(), NameTag(key).c_str(), ch + 1);
    TH2D *pHist = new TH2D(name, name,
                           Tbins, Tmin, Tmax, Ebins, Emin, Emax);
    sprintf(name, "%s, ToF vs Ekin, ch.%i", FC.c_str(), ch + 1);
    pHist->SetTitle(name);
    pHist->GetXaxis()->SetTitle("#font[12]{t} [ns]");
    pHist->GetYaxis()->SetTitle("#font[12]{E} [MeV]");

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

TGraphErrors* GetEmitMCNP(Double_t SimulatedN, string FC = "PuFC")
{
    Double_t SourceMaxAngle = atan(12.7/150);
    TGraphErrors *ge = new TGraphErrors(8);
    ge->SetTitle("Projectiles; Deposit; Neutrons");
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t SimulatedDistance = Distance(i, FC) - 2.0; //
        Double_t DepositAngle = atan(DepositRadius(i, FC) / SimulatedDistance);
        Double_t Emit = SolidAngle(DepositAngle) / SolidAngle(SourceMaxAngle);
        Double_t DEmit = sqrt(Emit / SimulatedN);
        ge->SetPoint(i, i+1, Emit);
        ge->SetPointError(i, 0, DEmit);
//        cout << Emit << "+-" << DEmit << endl;
    }
    return ge;
}

void TraLenMCNPtoROOT(string result_key = "2", Bool_t save = 1, Bool_t draw = 0, string FC = "PuFC", Long_t SimulatedN = 60000000000)
{ // Convert Track Length MCNP results to root
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
    TH2D *h[8];

    TFile *fAna;
    if (save)
        fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

    TCanvas *c1;
    if (draw)
    {
        c1 = new TCanvas();
        c1->Divide(4, 2);
    }

    // Loop over files
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/FCscat_g_tally-2%i4_xyz_0_%s.dmp", DirName.c_str(), 8 - i, result_key.c_str());

        // Create 2D histogram
        h[i] = MakeEvsT(name, 0, FC, result_key, i);
        // Track length estimator gives neutrons per started neutron and cm^2
        Double_t DepositArea = TMath::Pi() * pow(3.7, 2);
        h[i]->Scale(DepositArea);

        if (save)
        {
//            cout << "Saving " << "Simulation/MCNP/ToFvsEkin/"<<PathTag(result_key) << endl;
            Save(fAna, "Simulation/MCNP/ToFvsEkin/"+PathTag(result_key), h[i]);
        }
        if (draw)
        {
            c1->cd(i + 1);
//            gPad->SetLogz(1);
            h[i]->SetStats(0);
//            h[i]->GetXaxis()->SetRangeUser(0, 100);
            h[i]->Draw("colz");
        }
    }
    if (draw)
    {
        c1->Modified();
        c1->Update();
    }
    if (save)
    {
        TGraphErrors *geEmit = GetEmitMCNP(SimulatedN, FC);
        Save(fAna, "Simulation/MCNP/Correction", geEmit, "Emit");
        fAna->Save();
        fAna->Close();
    }
}

void SimpleMCNPtoROOT(string result_key = "2", Bool_t save = 1, Bool_t draw = 0, string FC = "PuFC", Long_t SimulatedN = 6000000000)
{ // Convert simple MCNP results to root
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
    TH2D *h[8];

    TFile *fAna;
    if (save)
        fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");


    TCanvas *c1;
    if (draw)
    {
        c1 = new TCanvas();
        c1->Divide(4, 2);
    }

    // Loop over files
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_%s.dmp", DirName.c_str(), 8 - i, result_key.c_str());
        string front_key, back_key;

        // Create 2D histogram
        if (result_key == "2")
        { // total
            front_key = "9";
            back_key = "8";
        }
        if (result_key == "1")
        { // scattered
            front_key = "5";
            back_key = "4";
        }
        if (result_key == "0")
        { // direct
            front_key = "1";
            back_key = "0";
        }
        // Open 2D front
        sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_%s.dmp", DirName.c_str(), 8-i, front_key.c_str());
        TH2D *pH2front = MakeEvsT(name, 0, FC, result_key, i, 0);
        // Change name to prevent overwriting
        sprintf(name, "%s_EvsT_front_%i", FC.c_str(), i+1);
        pH2front->SetName(name);
        // Open 2D back
        sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_%s.dmp", DirName.c_str(), 8-i, back_key.c_str());
        h[i] = MakeEvsT(name, 0, FC, result_key, i, 0);
        // Add front and back
        h[i]->Add(pH2front);

        if (save)
        {
//            cout << "Saving " << "Simulation/MCNP_simple/ToFvsEkin/"<<PathTag(result_key) << endl;
            Save(fAna, "Simulation/MCNP_simple/ToFvsEkin/"+PathTag(result_key), h[i]);
        }
        if (draw)
        {
            c1->cd(i + 1);
            h[i]->SetStats(0);
            h[i]->GetXaxis()->SetRangeUser(0, 100);
            h[i]->Draw("colz");
        }
    }

    if (draw)
    {
        c1->Modified();
        c1->Update();
    }
    if (save)
    {
        TGraphErrors *geEmit = GetEmitMCNP(SimulatedN, FC);
        Save(fAna, "Simulation/MCNP_simple/Correction", geEmit, "Emit");
        fAna->Save();
        fAna->Close();
    }
}//*/

TH2D* TH2FtoTH2D(TH2F *h1, string FC, string key, Int_t ch)
{
    const char *name = h1->GetName();
    const char *title = h1->GetTitle();
    string replace_name = FC+"_"+key+"_"+to_string(ch+1);
    h1->SetName(replace_name.c_str());
    if (h1->GetNbinsX() != 2000)
        cout << replace_name << " bad binning " << h1->GetNbinsX() << endl;
    if (h1->GetNbinsY() != 1600)
        cout << replace_name << " bad binning " << h1->GetNbinsY() << endl;
    TH2D *h2 = new TH2D(name, title, 2000, 0, 200, 1600, 0, 16);
    h2->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h2->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h2->Add(h1);
    return h2;
}

void Geant4toROOT(string FileName, string FC = "PuFC", string Setup = "Open")
{
    string FilePath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results";
//    string FileName = "5_ENDFVII.1/PuFC_Open_15E5.root";
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
    TGraphErrors *geEmit = new TGraphErrors(8);
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
        pHemit = (TH1F*)fG4->Get(name);
        geEmit->SetPoint(i, i+1, pHemit->Integral() / nEmit);
        geEmit->SetPointError(i, 0, 0);
    }
    Save(fAna, "Simulation/Geant4/"+FC+"_"+Setup+"/Correction", geEmit, "Emit");

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH2F *pH2tot = (TH2F*)fG4->Get(name);
        TH2D *pH2Tot = TH2FtoTH2D(pH2tot, FC, "2", i);

        if (pH2Tot->Integral() > 1)
            pH2Tot->Scale(1.0 / nEmit);

        sprintf(name, "%s/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH2F *pH2sc = (TH2F*)fG4->Get(name);
        TH2D *pH2Sc = TH2FtoTH2D(pH2sc, FC, "1", i);
        if (pH2Sc->Integral() > 1)
            pH2Sc->Scale(1.0 / nEmit);

        TH2D *pH2Dir = (TH2D*)pH2Tot->Clone();
        pH2Dir->Add(pH2Sc, -1);

        Save(fAna, "Simulation/Geant4/"+FC+"_"+Setup+"/ToFvsEkin/Direct", pH2Dir, FC+"_ToFvsEkin_Dir_Ch."+to_string(i+1));
        Save(fAna, "Simulation/Geant4/"+FC+"_"+Setup+"/ToFvsEkin/Scattered", pH2Sc, FC+"_ToFvsEkin_Sc_Ch."+to_string(i+1));
        Save(fAna, "Simulation/Geant4/"+FC+"_"+Setup+"/ToFvsEkin", pH2Tot, FC+"_ToFvsEkin_Ch."+to_string(i+1));
    }
    fAna->Save();
    fAna->Close();
}


void MCNPtoROOT()
{
    Geant4toROOT("5_ENDFVII.1/PuFC_Open_5E7_v2.root", "PuFC", "Open");
    Geant4toROOT("4_ene/UFC_Open_5E7.root", "UFC", "Open");
    Geant4toROOT("4_ene/UFC_SB_5E7.root", "UFC", "SB");
    SimpleMCNPtoROOT("0", 1, 0);
    SimpleMCNPtoROOT("1", 1, 0);
    SimpleMCNPtoROOT("2", 1, 0);
    TraLenMCNPtoROOT("0", 1, 0);
    TraLenMCNPtoROOT("1", 1, 0);
    TraLenMCNPtoROOT("2", 1, 0);
}
