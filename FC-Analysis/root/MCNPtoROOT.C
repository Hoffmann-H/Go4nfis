#include "SaveToFile.C"
#include "FC.C"
#include <fstream>

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

void TraLenMCNPtoROOT(string result_key = "2", Bool_t save = 0, Bool_t draw = 1, string FC = "PuFC")
{ // Convert Track Length MCNP results to root
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
    TH2D *h[8];

    TFile *f;
    if (save)
    {
        f = TFile::Open("/home/hoffma93/Programme/MCNP/results/MCNP.root", "UPDATE");
    }
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
            cout << "Saving " << FC<<"/ToFvsEkin/"<<PathTag(result_key) << endl;
            Save(f, FC+"/ToFvsEkin/"+PathTag(result_key), h[i]);
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
        f->Save();
        f->Close();
    }
}

void SimpleMCNPtoROOT(string result_key = "2", Bool_t save = 0, Bool_t draw = 1, string FC = "PuFC")
{ // Convert simple MCNP results to root
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
    TH2D *h[8];
    TFile *f;
    if (save)
    {
        f = TFile::Open("/home/hoffma93/Programme/MCNP/results/MCNP_simple.root", "UPDATE");
    }
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
            cout << "Saving " << FC<<"/ToFvsEkin/"<<PathTag(result_key) << endl;
            Save(f, FC+"/ToFvsEkin/"+PathTag(result_key), h[i]);
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
        f->Save();
        f->Close();
    }
}//*/


void MCNPtoROOT()
{
    TraLenMCNPtoROOT("0", 1, 0);
    TraLenMCNPtoROOT("1", 1, 0);
    TraLenMCNPtoROOT("2", 1, 0);
    SimpleMCNPtoROOT("0", 1, 0);
    SimpleMCNPtoROOT("1", 1, 0);
    SimpleMCNPtoROOT("2", 1, 0);
}
