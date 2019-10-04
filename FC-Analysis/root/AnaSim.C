#include "SaveToFile.C"
#include "MCNPtoROOT.C"
#include "FC.C"

TH2D* GetSpectrum(TFile *f, string SimulationPath, Int_t ch, string key = "2", string FC = "PuFC")
{
    /// Open a histogram in a file created by MCNPtoROOT.C
    char name[64] = "";
    sprintf(name, "%s/ToFvsEkin/%s%s_ToFvsEkin_%sCh.%i", SimulationPath.c_str(), PathTag(key).c_str(), FC.c_str(), NameTag(key).c_str(), ch + 1);
    TH2D *h = (TH2D*)f->Get(name);
    if (h == 0)
        cout << "Could not open " << name << endl;
    else
        cout << "Opened " << h->GetName() << endl;
    return h;
}

TGraphErrors* GetSigma(string FC = "PuFC", string path = "/home/hoffma93/Programme/ROOT/Data")
{
    char name[64] = "";
    TGraphErrors *pSigma;
    if (!strcmp(FC.c_str(), "PuFC")) // if PuFC
    {
        sprintf(name, "%s/Pu242.dat", path.c_str());
        pSigma = new TGraphErrors(name, "%lg %lg %lg");
        if (pSigma == 0)
            cout << "Fehler beim Öffnen von " << path << endl;
        return pSigma;
    } else {
        Double_t UisoVec[] = {0.988, 0.0912}; // U-235, U-238 portion
        Double_t DUisoVec[] = {0.005, 0.0006};
        sprintf(name, "%s/U235.dat", path.c_str());
        TGraphErrors *pU235 = new TGraphErrors(name, "%lg %lg %lg");
        sprintf(name, "%s/U238.dat", path.c_str());
        TGraphErrors *pU238 = new TGraphErrors(name, "%lg %lg %lg");
        Int_t N = pU235->GetN();
        pSigma = new TGraphErrors(N);
        for (Int_t j = 0; j < N; j++)
        {
            Double_t x, y235, yerr235, y238;
            pU235->GetPoint(j, x, y235);
            yerr235 = pU235->GetErrorY(j);
            y238 = pU238->Eval(x);
            Double_t y = UisoVec[0] * y235 + UisoVec[1] * y238;
            pSigma->SetPoint(j, x, y);
            pSigma->SetPointError(j, 0, yerr235);
        }
        return pSigma;
    }
}

TH1D* GetProjection(TFile *f, Int_t ch, string key = "2", string FC = "PuFC", string SimulationPath = "Simulation/Geant4/PuFC_Open")
{
    char name[64] = "";
    sprintf(name, "%s/EffToF/%s%s_ProjT_%s%i", SimulationPath.c_str(), PathTag(key).c_str(), FC.c_str(), NameTag(key).c_str(), ch + 1);
    TH1D *h = (TH1D*)f->Get(name);
    cout << "Opened " << h->GetName() << endl;
    return h;
}

TH1D* GetProjection(TH2 *pH2, TGraphErrors *pSigma, Int_t ch, string FC = "PuFC")
{
    char name[128] = "";

    // create projection histogram
    Int_t NbinsX = pH2->GetNbinsX();
    Int_t NbinsY = pH2->GetNbinsY();
    Double_t xmin = pH2->GetXaxis()->GetBinLowEdge(1);
    Double_t xmax = pH2->GetXaxis()->GetBinLowEdge(NbinsX + 1);
    Double_t ymin = pH2->GetYaxis()->GetBinLowEdge(1);
    Double_t ymax = pH2->GetYaxis()->GetBinLowEdge(NbinsY + 1);
    sprintf(name, "%s%i", FC.c_str(), ch + 1);
    TH1D *pH1x = new TH1D(name, name, NbinsX, xmin, xmax);
    sprintf(name, "%s, Time profile, Ch.%i; #font[12]{t} / ns; Effective neutrons", FC.c_str(), ch + 1);
    pH1x->SetTitle(name);
    pH1x->Sumw2();

//    cout << "E   w   eff.N" << endl;

    // Loop over neutron energies
    for (Int_t ybin = 0; ybin < NbinsY + 2; ybin++)
    {
        Double_t E = pH2->GetYaxis()->GetBinCenter(ybin);
        Double_t w = pSigma->Eval(E) / pSigma->Eval(15);
        TH1D *row = pH2->ProjectionX("px", ybin, ybin);

//        cout << E << "  " << w << "  " << row->Integral() << endl;

        // Fill projection row by row
        pH1x->Add(row, w);
    }

    return pH1x;
}

void Projections(string key = "2", string FC = "PuFC", string SimulationPath = "Simulation/Geant4/PuFC_Open")
{
    char name[64] = "";
    TGraphErrors *pSigma = GetSigma(FC);
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

    for (Int_t i = 0; i < 8; i++)
    {
        TH2D *pH2 = GetSpectrum(fAna, SimulationPath, i, key, FC);
        TH1D *pH1 = GetProjection(pH2, pSigma, i, FC);
        sprintf(name, "%s_ProjT_%s%i", FC.c_str(), NameTag(key).c_str(), i+1);
        pH1->SetName(name);
        Save(fAna, SimulationPath+"/EffToF/"+PathTag(key), pH1);
        Double_t Int, Err;
        Int = pH1->IntegralAndError(0, -1, Err);
    }
    fAna->Save();
    fAna->Close();
}

void AnaSim(string FC, Bool_t save = 1, string SimulationPath = "Simulation/Geant4/PuFC_Open")
{
    // Create projections.
    Projections("0", FC, SimulationPath);
    Projections("1", FC, SimulationPath);
    Projections("2", FC, SimulationPath);

    Double_t S[8], DS[8], T[8], DT[8], F[8], DF[8];

    char name[64] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    sprintf(name, "%s/Correction/Emit", SimulationPath.c_str());
    TGraphErrors *gEmit = (TGraphErrors*)fAna->Get(name);

    TGraphErrors *gDirect = new TGraphErrors(8);
    gDirect->SetTitle("Direct fission; Deposit; Effective Neutrons");
    TGraphErrors *gTotal = new TGraphErrors(8);
    gTotal->SetTitle("Total fission; Deposit; Effective Neutrons");
    TGraphErrors *gS = new TGraphErrors(8);
    gS->SetTitle("Scattering correction; Deposit; S");
    TGraphErrors *gT = new TGraphErrors(8);
    gT->SetTitle("Transmission; Deposit; T");
    TGraphErrors *gF = new TGraphErrors(8);
    gF->SetTitle("Correction factor; Deposit; F");

    cout << "Ch      Emit      Direct      Total      S      T      F" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // projectiles
        Double_t x, nEmit, DnEmit;
        gEmit->GetPoint(i, x, nEmit);
        DnEmit = gEmit->GetErrorY(i);

        // effective scattered+total neutrons
        sprintf(name, "%s/EffToF/%s_ProjT_%i", SimulationPath.c_str(), FC.c_str(), i+1);
        TH1D *pTotFis = (TH1D*)fAna->Get(name); // fissions induced by backward sc neutrons over time
        if (pTotFis == 0)
            cout << "Could not open " << name << endl;
        sprintf(name, "%s/EffToF/Direct/%s_ProjT_Dir_%i", SimulationPath.c_str(), FC.c_str(), i+1);
        TH1D *pDirFis = (TH1D*)fAna->Get(name); // fissions induced by direct neutrons over time
        if (pDirFis == 0)
            cout << "Could not open " << name << endl;
        Double_t DnTotal, DnDirect;
        Double_t nTotal = pTotFis->IntegralAndError(0, -1, DnTotal);
        Double_t nDirect = pDirFis->IntegralAndError(0, -1, DnDirect);

        // Correction factors
        S[i] = nDirect / nTotal;
        DS[i] = S[i] * sqrt(pow(DnDirect / nDirect, 2) + pow(DnTotal / nTotal, 2));
        T[i] = nDirect / nEmit;
        DT[i] = T[i] * sqrt(pow(DnDirect / nDirect, 2) + pow(DnEmit / nEmit, 2));
        F[i] = nEmit / nTotal;
        DF[i] = F[i] * sqrt(pow(DnEmit / nEmit, 2) + pow(DnTotal / nTotal, 2));
        cout << " " << i+1 << "   " << nEmit << "+-" << DnEmit
                           << "   " << nDirect << "+-" << DnDirect
                           << "   " << nTotal << "+-" << DnTotal
                           << "   " << S[i] << "+-" << DS[i]
                           << "   " << T[i] << "+-" << DT[i]
                           << "   " << F[i] << "+-" << DF[i] << endl;

        // Write to graphs
        gDirect->SetPoint(i, i + 1, nDirect);
        gDirect->SetPointError(i, 0, DnDirect);
        gTotal->SetPoint(i, i + 1, nTotal);
        gTotal->SetPointError(i, 0, DnTotal);
        gS->SetPoint(i, i + 1, S[i]);
        gS->SetPointError(i, 0, DS[i]);
        gT->SetPoint(i, i + 1, T[i]);
        gT->SetPointError(i, 0, DT[i]);
        gF->SetPoint(i, i + 1, F[i]);
        gF->SetPointError(i, 0, DF[i]);
    } // end of for(Deposits)

    if (save)
    {
        Save(fAna, SimulationPath+"/Correction", gDirect, "Direct");
        Save(fAna, SimulationPath+"/Correction", gTotal, "Total");
        Save(fAna, SimulationPath+"/Correction", gS, "S");
        Save(fAna, SimulationPath+"/Correction", gT, "T");
        Save(fAna, SimulationPath+"/Correction", gF, "F");
        fAna->Save();
        fAna->Close();
    }
}

TH1D* GetPeakForm(TH1D *pH1ToF, Int_t ch = 0, string FC = "PuFC", Double_t IntExp = 1)
{
    char name[64] = "";
    char title[128] = "";

    /// simulation properties
    Double_t AccPulseLength = 7.0; // ns
    Double_t TimeResolution = 2.5; // ns
    Int_t NbinsProj = pH1ToF->GetNbinsX();
    Double_t ProjBinWidth = pH1ToF->GetXaxis()->GetBinWidth(1);
    Int_t AccPulseBins = AccPulseLength / ProjBinWidth;
    cout << endl << "Simulating ToF spectra..." << endl
         << " Acc. pulse length: " << AccPulseLength << " ns" << endl
         << " Time resolution: " << TimeResolution << " ns" << endl
         << " Original bins: " << NbinsProj << endl
         << " Original bin width: " << ProjBinWidth << endl
         << " Acc. pulse bins: " << AccPulseBins << endl;

    /// experimental properties - depending on FC
    Double_t tMaximum = PeakCenter(ch, FC);
    Double_t minToF = 0;
    Double_t maxToF = 440;
    Int_t NbinsToF = (maxToF - minToF) / ProjBinWidth; // use Tproj's bin width

    /// Prepare simulated ToF histogram
    sprintf(name, "%s_Dt_Ch.%i", FC.c_str(), ch+1);
    sprintf(title, "%s, Ch. %i, simulated ToF spectrum; #font[12]{t} [ns]; Counts", FC.c_str(), ch+1);
    TH1D *pH1Peak = new TH1D(name, title,  NbinsToF, minToF, maxToF);
    cout << " Sim. bins: " << pH1Peak->GetNbinsX() << endl
         << " Sim. ToF range: " << pH1Peak->GetBinLowEdge(1) << "-" << pH1Peak->GetBinLowEdge(pH1Peak->GetNbinsX()+1) << " ns" << endl;

    /// create folding function
    Double_t ampl = 1.0 / (Double_t)AccPulseBins / sqrt(2 * TMath::Pi()) * pH1Peak->GetBinWidth(1) / TimeResolution;
    cout << " Ampl: " << ampl << endl;
    TF1* g = new TF1("fG", "gaus", -440, 440);
    g->SetParameters(1, 0, TimeResolution);

    /// fold
    Int_t ProjMaxBin = pH1ToF->GetMaximumBin();
    Double_t ProjMax = pH1ToF->GetBinCenter(ProjMaxBin);
    //        cout << ProjMax << " -> " << tm[i] << endl;
    for (Int_t j = 1; j <= NbinsToF; j++)
    {
        Double_t t0 = pH1Peak->GetBinCenter(j); // ToF sampling point
        //            cout << t0 << endl;
        Double_t sum = 0;
        for (Int_t k = 1; k < NbinsProj - AccPulseBins; k++)
        {
            Double_t t1 = 0.5 * (pH1ToF->GetBinCenter(k) + pH1ToF->GetBinCenter(k + AccPulseBins));
            //                cout << " " << t1;
            sum += g->Eval(t1 - t0 + tMaximum - ProjMax) * pH1ToF->Integral(k, k + AccPulseBins - 1);
        }
        //            cout << j << "  " << t0 << "  " << sum << endl;
        pH1Peak->SetBinContent(j, sum);
    }
    Double_t IntSim = pH1Peak->Integral();
    pH1Peak->Scale(IntExp / IntSim);
    return pH1Peak;
}

void PeakForm(string key = "2", string FC = "PuFC", string SimulationPath = "Simulation/Geant4/PuFC_Open")
{
    char name[64] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

    for (Int_t i = 0; i < 8; i++)
    {
        TH1D *pH1ToF = GetProjection(fAna, i, key, FC, SimulationPath);
        TH1D *pH1Peak = GetPeakForm(pH1ToF, i, FC, 1.0);
        sprintf(name, "%s_FoldT_%s%i", FC.c_str(), NameTag(key).c_str(), i+1);
        pH1Peak->SetName(name);
        Save(fAna, SimulationPath+"/SimToF/"+PathTag(key), pH1Peak);
    }
    fAna->Save();
    fAna->Close();
}

void AnaSim()
{
    MCNPtoROOT();
    AnaSim("PuFC", 1, "Simulation/Geant4/PuFC_Open");
    AnaSim("UFC", 1, "Simulation/Geant4/UFC_Open");
    PeakForm("1", "UFC", "Simulation/Geant4/UFC_Open");
    PeakForm("2", "UFC", "Simulation/Geant4/UFC_Open");
    AnaSim("UFC", 1, "Simulation/Geant4/UFC_SB");
    AnaSim("PuFC", 1, "Simulation/MCNP_simple");
    AnaSim("PuFC", 1, "Simulation/MCNP");
    PeakForm("1", "PuFC", "Simulation/MCNP");
    PeakForm("2", "PuFC", "Simulation/MCNP");
}
