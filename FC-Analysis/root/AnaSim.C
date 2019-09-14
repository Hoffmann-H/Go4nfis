#include "SaveToFile.C"
#include "MCNPtoROOT.C"
#include "FC.C"

TH2D* GetSpectrum(TFile *f, Int_t ch, string key = "2", string FC = "PuFC")
{
    /// Open a histogram in a file created by MCNPtoROOT.C
    char name[64] = "";
    sprintf(name, "%s/ToFvsEkin/%s%s_ToFvsEkin_%sCh.%i", FC.c_str(), PathTag(key).c_str(), FC.c_str(), NameTag(key).c_str(), ch + 1);
    TH2D *h = (TH2D*)f->Get(name);
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
        Double_t UisoVec[] = {0.904, 0.0912}; // U-235, U-238 portion
        Double_t DUisoVec[] = {0.005, 0.0006};
        cout << "Uranium cross section not implemented!" << endl;
        // implementation needed...
        return 0;
    }
}

TH1D* GetProjection(TFile *f, Int_t ch, string key = "2", string FC = "PuFC")
{
    char name[64] = "";
    sprintf(name, "Analysis/EffToF/%s%s_ProjT_%s%i", PathTag(key).c_str(), FC.c_str(), NameTag(key).c_str(), ch + 1);
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

void MCNPprojections(string key = "2", string FC = "PuFC", string File = "/home/hoffma93/Programme/MCNP/results/MCNP.root")
{
    char name[64] = "";
    TGraphErrors *pSigma = GetSigma(FC);
    TFile *f = TFile::Open(File.c_str(), "UPDATE");
    TDirectory *pDir = Prepare(f, "Analysis/EffToF/" + PathTag(key));
    for (Int_t i = 0; i < 8; i++)
    {
        TH2D *pH2 = GetSpectrum(f, i, key, FC);
        TH1D *pH1 = GetProjection(pH2, pSigma, i);
        sprintf(name, "%s_ProjT_%s%i", FC.c_str(), NameTag(key).c_str(), i+1);
        Save(pDir, pH1, name);
        Double_t Int, Err;
        Int = pH1->IntegralAndError(0, -1, Err);
    }
    f->Save();
    f->Close();
}

void AnaMCNP(string FC = "PuFC", Bool_t save = 1, string File = "/home/hoffma93/Programme/MCNP/results/MCNP.root", Long_t SimulatedN = 60000000000)
{

    // Create projections.
    MCNPprojections("0", "PuFC", File);
    MCNPprojections("1", "PuFC", File);
    MCNPprojections("2", "PuFC", File);

    Double_t SourceMaxAngle = atan(12.7/150.0);
    cout << "Maximum simulated source angle: " << SourceMaxAngle / asin(1) * 90.0 << "°" << endl;
    Double_t DepositRadius = 37.0; // mm
    Double_t DepositDistance;
    Double_t S[8], DS[8], T[8], DT[8], F[8], DF[8];

    char name[64] = "";
    TFile *f = TFile::Open(File.c_str(), "UPDATE");

    TGraphErrors *gEmit = new TGraphErrors(8);
    gEmit->SetTitle("Projectiles; Deposit; Neutrons");
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

    cout << endl << File << endl;
    cout << "Ch      Emit      Direct      Total      S      T      F" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // projectiles
        DepositDistance = Distance(i, FC);
//        Double_t DepositMaxAngle = atan(DepositRadius / DepositDistance);
        Double_t nEmit = SolidAngle(DepositDistance, DepositRadius) / SolidAngle(150.0, 12.7);
        Double_t DnEmit = sqrt(nEmit / SimulatedN);

        // effective scattered+total neutrons
        sprintf(name, "Analysis/EffToF/%s_ProjT_%i", FC.c_str(), i+1);
        TH1D *pTotFis = (TH1D*)f->Get(name); // fissions induced by backward sc neutrons over time
        sprintf(name, "Analysis/EffToF/Direct/%s_ProjT_Dir_%i", FC.c_str(), i+1);
        TH1D *pDirFis = (TH1D*)f->Get(name); // fissions induced by direct neutrons over time
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
        gEmit->SetPoint(i, i + 1, nEmit);
        gEmit->SetPointError(i, 0, DnEmit);
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
        TDirectory *pDir = Prepare(f, "Analysis/Correction");
        Save(pDir, gEmit, "Emit");
        Save(pDir, gDirect, "Direct");
        Save(pDir, gTotal, "Total");
        Save(pDir, gS, "S");
        Save(pDir, gT, "T");
        Save(pDir, gF, "F");
        f->Save();
    }
}

void AnaGeant4(string FC = "PuFC", Bool_t save = 1, string FileName = "4_ene/PuFC_Open_5E7.root")
{ //    string FileName = "5_ENDFVII.1/PuFC_Open_15E5.root";
    /// Open G4 result file
    /// Scale TvsE hist
    /// Create time profile
    /// Calculate correction factors
    /// Save all
    char name[128] = "";
    string FilePath = "/home/hoffma93/Programme/Geant4-Work/builds/G4PuFCvsH19/results";

    // preparation...
    sprintf(name, "%s/%s", FilePath.c_str(), FileName.c_str());
    TFile *f = TFile::Open(name, "UPDATE"); // Open file

    TH2F *pH2Tot[8]; // prepare histograms
    TH2F *pH2Sc[8];
    TH1D *pH1Tot[8];
    TH1D *pH1Sc[8];
    TGraphErrors *pSigma = GetSigma(FC); // get evaluated cross section
    TH1F *pH1theta = (TH1F*)f->Get("Source/Source_Theta");
    Long_t SimulatedN = pH1theta->Integral(); // get simulation population

    // scale old results
//    SimulatedN *= (1-cos(0.08446522295019329)) / (1-cos(0.0746987));
    // Set N for an invalid file
    SimulatedN = 5.E6;

    Double_t S[8], DS[8], T[8], DT[8], F[8], DF[8]; // prepare analysis

    // create graphs
    TGraphErrors *gEmit = new TGraphErrors(8);
    gEmit->SetTitle("Projectiles");
    TGraphErrors *gDirect = new TGraphErrors(8);
    gDirect->SetTitle("Direct fission");
    TGraphErrors *gTotal = new TGraphErrors(8);
    gTotal->SetTitle("Total fission");
    TGraphErrors *gS = new TGraphErrors(8);
    gS->SetTitle("Scattering correction");
    TGraphErrors *gT = new TGraphErrors(8);
    gT->SetTitle("Transmission");
    TGraphErrors *gF = new TGraphErrors(8);
    gF->SetTitle("Correction factor");

    cout << endl << FilePath << "/" << FileName << endl;
    cout << "Ch      Emit      Direct      Total      S      T      F" << endl;
    for (Int_t i = 0; i < 8; i++)
    {
        // Open 2D T vs E histograms, give name dummies
        sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        pH2Tot[i] = (TH2F*)f->Get(name);
        sprintf(name, "TvsE_%i", i+1);
        pH2Tot[i]->SetName(name);
        sprintf(name, "%s/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        pH2Sc[i] = (TH2F*)f->Get(name);
        sprintf(name, "TvsE_Sc_%i", i+1);
        pH2Sc[i]->SetName(name);
        // scale T vs E to one simulated neutron (if not done yet)
        if (pH2Tot[i]->Integral() > 1)
            pH2Tot[i]->Scale(1.0/SimulatedN);
        if (pH2Sc[i]->Integral() > 1)
            pH2Sc[i]->Scale(1.0/SimulatedN);

        // Create projections
        pH1Tot[i] = GetProjection(pH2Tot[i], pSigma, i, FC);
        pH1Sc[i] = GetProjection(pH2Sc[i], pSigma, i, FC + "Sc");

        // projectiles
        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
        pH1theta = (TH1F*)f->Get(name);
        Double_t nEmit = pH1theta->Integral() / SimulatedN;
        Double_t DnEmit = 0.0; //sqrt(nEmit / SimulatedN);
        Double_t nDirect = pH1Tot[i]->Integral() - pH1Sc[i]->Integral();
        Double_t DnDirect = sqrt(nDirect / SimulatedN);
        Double_t nTotal = pH1Tot[i]->Integral();
        Double_t DnTotal = sqrt(nTotal / SimulatedN);

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
        gEmit->SetPoint(i, i + 1, nEmit);
        gEmit->SetPointError(i, 0, DnEmit);
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
        TDirectory *pDir;
        pDir = Prepare(f, FC + "/ToFvsEkin");
        for (Int_t i = 0; i < 8; i++)
        {
            sprintf(name, "%s_ToFvsEkin_Ch.%i", FC.c_str(), i+1);
            Save(pDir, pH2Tot[i], name);
        }
        pDir = Prepare(f, FC + "/ToFvsEkin/Scattered");
        for (Int_t i = 0; i < 8; i++)
        {
            sprintf(name, "%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), i+1);
            Save(pDir, pH2Sc[i], name);
        }
        pDir = Prepare(f, "Analysis/EffToF");
        for (Int_t i = 0; i < 8; i++)
        {
            sprintf(name, "%s_ProjT_%i", FC.c_str(), i+1);
            Save(pDir, pH1Tot[i], name);
        }
        pDir = Prepare(f, "Analysis/EffToF/Scattered");
        for (Int_t i = 0; i < 8; i++)
        {
            sprintf(name, "%s_ProjT_Sc_%i", FC.c_str(), i+1);
            Save(pDir, pH1Sc[i], name);
        }
        pDir = Prepare(f, "Analysis/Correction");
        Save(pDir, gEmit, "Emit");
        Save(pDir, gDirect, "Direct");
        Save(pDir, gTotal, "Total");
        Save(pDir, gS, "S");
        Save(pDir, gT, "T");
        Save(pDir, gF, "F");
        f->Save();
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

void PeakForm(string key = "2", string FC = "PuFC", string File = "/home/hoffma93/Programme/MCNP/results/MCNP.root")
{
    char name[64] = "";
    TFile *f = TFile::Open(File.c_str(), "UPDATE");
    TDirectory *pDir = Prepare(f, "Analysis/SimToF/" + PathTag(key));
    for (Int_t i = 0; i < 8; i++)
    {
        TH1D *pH1ToF = GetProjection(f, i, key, FC);
        TH1D *pH1Peak = GetPeakForm(pH1ToF, i, FC, 1.0);
        sprintf(name, "%s_FoldT_%i", FC.c_str(), i+1);
        Save(pDir, pH1Peak, name);
    }
    f->Save();
    f->Close();
}

void AnaSim()
{
//    MCNPtoROOT();
//    AnaMCNP("PuFC", 1, "/home/hoffma93/Programme/MCNP/results/MCNP_simple.root", 6000000000);
    AnaMCNP("PuFC", 1, "/home/hoffma93/Programme/MCNP/results/MCNP.root", 60000000000);
//    AnaGeant4("PuFC", 1, "5_ENDFVII.1/PuFC_Open_5E7.root");
}
