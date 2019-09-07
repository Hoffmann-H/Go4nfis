#include "SaveToFile.C"
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

TH2D* MakeEvsT(string file_to_read, Bool_t draw, string FC = "PuFC", Int_t ch = 0, Bool_t track_length = 1)
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
    sprintf(name, "EvsT%i", ch + 1);
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

void MCNPtoROOT(string result_key = "2", Bool_t save = 0, Bool_t draw = 1, string FC = "PuFC")
{
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
    TH2D *h[8];

    TFile *f;
    TDirectory *pDir;
    if (save)
    {
        f = TFile::Open("/home/hoffma93/Programme/MCNP/results/MCNP.root", "UPDATE");
//        string path;
//        if (result_key == "2")
//            path = "PuFC/ToFvsEkin/";
//        if (result_key == "1")
//            path = "PuFC/ToFvsEkin/Scattered/";
//        if (result_key == "0")
//            path = "PuFC/ToFvsEkin/Direct/";
        pDir = Prepare(f, FC + "/ToFvsEkin/" + PathTag(result_key));
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
        h[i] = MakeEvsT(name, 0, FC, i);

        if (save)
        {
            std::stringstream s;
            if (result_key == "2")
                s << FC << "_ToFvsEkin_Ch." << i+1;
            if (result_key == "1")
                s << FC << "_ToFvsEkin_Sc_Ch." << i+1;
            if (result_key == "0")
                s << FC << "_ToFvsEkin_Dir_Ch." << i+1;
            Save(pDir, h[i], s.str());
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
{
    char name[128] = "";
    string DirName = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_PTB/tally";
    TH2D *h[8];
    TFile *f;
    TDirectory *pDir;
    if (save)
    {
        f = TFile::Open("/home/hoffma93/Programme/MCNP/results/MCNP_simple.root", "UPDATE");
        string path;
        if (result_key == "2")
            path = "PuFC/ToFvsEkin/";
        if (result_key == "1")
            path = "PuFC/ToFvsEkin/Scattered/";
        if (result_key == "0")
            path = "PuFC/ToFvsEkin/Direct/";
        pDir = Prepare(f, path);
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

        // Create 2D histogram
        if (result_key == "2")
        { // total
            sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_8.dmp", DirName.c_str(), 8 - i);
            h[i] = MakeEvsT(name, 0, FC, i, 0);
            sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_9.dmp", DirName.c_str(), 8 - i);
            TH2D *pH2front = MakeEvsT(name, 0, FC, i, 0);
            h[i]->Add(pH2front);
        }
        if (result_key == "1")
        { // scattered
            sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_4.dmp", DirName.c_str(), 8 - i);
            h[i] = MakeEvsT(name, 0, FC, i, 0);
            sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_5.dmp", DirName.c_str(), 8 - i);
            TH2D *pH2front = MakeEvsT(name, 0, FC, i, 0);
            h[i]->Add(pH2front);
        }
        if (result_key == "0")
        { // direct
            sprintf(name, "%s/FCscat_b_tally-2%i1_xyz_0_1.dmp", DirName.c_str(), 8 - i);
            h[i] = MakeEvsT(name, 0, FC, i, 0);
        }

        if (save)
        {
            std::stringstream s;
            if (result_key == "2")
                s << FC << "_ToFvsEkin_Ch." << i+1;
            if (result_key == "1")
                s << FC << "_ToFvsEkin_Sc_Ch." << i+1;
            if (result_key == "0")
                s << FC << "_ToFvsEkin_Dir_Ch." << i+1;
            Save(pDir, h[i], s.str());
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
}
