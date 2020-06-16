#ifndef ANASIM_H
#define ANASIM_H
#include "SaveToFile.C"
#include "MCNPtoROOT.C"
#include "FC.C"

TH2D* GetSpectrum(TFile *f, Int_t ch, string SimulationPath, string FC = "PuFC", string key = "real")
{
    /// Open a histogram in a file created by MCNPtoROOT.C
    char name[64] = "";
    sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_%s_Ch.%i", SimulationPath.c_str(), FC.c_str(), key.c_str(), ch + 1);
    TH2D *h = (TH2D*)f->Get(name);
    if (h == 0)
        cout << "Could not open " << name << endl;
//    else
//        cout << "Opened " << h->GetName() << endl;
    return h;
}

TH1D* GetProjection(TFile *f, Int_t ch, string SimulationPath = "Simulation/Geant4/PuFC_real", string FC = "PuFC", string key = "real")
{
    char name[64] = "";
    TH2D *pH2 = GetSpectrum(f, ch, SimulationPath, FC, key);
    sprintf(name, "%s_ProjT_%s_%i", FC.c_str(), key.c_str(), ch + 1);
    TH1D *h = (TH1D*)pH2->ProjectionX(name);
//    h->SetName(name);
    sprintf(name, "%s, ToF, %s, ch.%i", FC.c_str(), key.c_str(), ch + 1);
    h->SetTitle(name);
    h->GetXaxis()->SetTitle("#font[12]{t} / ns");
    h->GetYaxis()->SetTitle("#font[12]{C}_{(n,f)}");
    return h;
}

void Projections(string SimulationPath = "Simulation/Geant4/PuFC_real", string FC = "PuFC", string key = "real")
{
    char name[64] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");

    for (Int_t i = 0; i < 8; i++)
    {
        TH1D *pH1 = GetProjection(fAna, i, SimulationPath, FC, key);
        Save(fAna, SimulationPath+"/EffToF/", pH1);
        Double_t Int, Err;
        Int = pH1->IntegralAndError(0, -1, Err);
    }
    fAna->Save();
    fAna->Close();
}//*/

TH1D* GetPeakForm(TH1D *pH1ToF, Int_t ch = 0, string FC = "PuFC", Double_t IntExp = 1)
{
    char name[64] = "";
    char title[128] = "";

    /// simulation properties
    Double_t AccPulseLength = 7.0; // ns
    Double_t TimeResolution = 1.7; // ns
    Int_t NbinsProj = pH1ToF->GetNbinsX();
    Double_t ProjBinWidth = pH1ToF->GetXaxis()->GetBinWidth(1);
    Int_t AccPulseBins = AccPulseLength / ProjBinWidth;
//    cout << endl << "Simulating ToF spectra..." << endl
//         << " Acc. pulse length: " << AccPulseLength << " ns" << endl
//         << " Time resolution: " << TimeResolution << " ns" << endl
//         << " Original bins: " << NbinsProj << endl
//         << " Original bin width: " << ProjBinWidth << endl
//         << " Acc. pulse bins: " << AccPulseBins << endl;

    /// experimental properties - depending on FC
//    Double_t tMaximum = PeakCenter(ch, FC);
    Double_t minToF = 0;
    Double_t maxToF = 440;
    Int_t NbinsToF = (maxToF - minToF) / ProjBinWidth; // use Tproj's bin width

    /// Prepare simulated ToF histogram
    sprintf(name, "%s_Dt_Ch.%i", FC.c_str(), ch+1);
    sprintf(title, "%s, Ch. %i, simulated ToF spectrum; #font[12]{t} [ns]; Counts", FC.c_str(), ch+1);
    TH1D *pH1Peak = new TH1D(name, title,  NbinsToF, minToF, maxToF);
//    cout << " Sim. bins: " << pH1Peak->GetNbinsX() << endl
//         << " Sim. ToF range: " << pH1Peak->GetBinLowEdge(1) << "-" << pH1Peak->GetBinLowEdge(pH1Peak->GetNbinsX()+1) << " ns" << endl;

    /// create folding function
    Double_t ampl = 1.0 / (Double_t)AccPulseBins / sqrt(2 * TMath::Pi()) * pH1Peak->GetBinWidth(1) / TimeResolution;
//    cout << " Ampl: " << ampl << endl;
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
            Double_t t1 = 0.5 * (pH1ToF->GetBinLowEdge(k) + pH1ToF->GetBinLowEdge(k + AccPulseBins));
            //                cout << " " << t1;
            sum += g->Eval(t1 - t0) * pH1ToF->Integral(k, k + AccPulseBins - 1); // g->Eval(t1 - t0 + tMaximum - ProjMax)
        }
        //            cout << j << "  " << t0 << "  " << sum << endl;
        pH1Peak->SetBinContent(j, sum);
    }
    Double_t IntSim = pH1Peak->Integral();
    pH1Peak->Scale(IntExp / IntSim);
    return pH1Peak;
}

Double_t GatingCorrection(TH1D *pH1Real, TH1D *pH1Peak, Int_t ch, string FC = "PuFC", Double_t width = 15.0)
{
    Double_t bg = pH1Peak->GetBinContent(pH1Peak->GetNbinsX());
    Int_t tMax = (Int_t)pH1Peak->GetBinCenter(pH1Peak->GetMaximumBin());
    Int_t bl = pH1Peak->FindBin(tMax - Left(FC));
    Int_t br = pH1Peak->FindBin(tMax + Right(FC) + 1.0) - 1;
    Double_t Peak = pH1Peak->Integral() - pH1Peak->GetNbinsX() * bg;
    Double_t GatedPeak = pH1Peak->Integral(bl, br) - (br - bl + 1) * bg;

    tMax = (Int_t)pH1Real->GetBinCenter(pH1Real->GetMaximumBin());
    bl = pH1Real->FindBin(tMax - Left(FC));
    br = pH1Real->FindBin(tMax + Right(FC) + 1.0) - 1;
    Double_t Real = pH1Real->Integral();
    Double_t GatedReal = pH1Real->Integral(bl, br);
//    cout << Peak << " " << GatedPeak << " " << Real << " " << GatedReal << endl;
    return Peak / GatedPeak * GatedReal / Real;
}

void AnaSim(string Simulation, string FC = "PuFC", string key = "real")
{ // SimToF == 0: force SimToF recalculation. 1: Use SimToF if possible. 2: Use FitToF if possible. Otherwise create it.
    cout << "Analyzing " << Simulation << " " << FC << " " << key << " simulation results" << endl;
    char name[128] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    if (!fAna) cout << "Could not open " << "Analysis.root" << endl;
    sprintf(name, "%s/Correction/%s_Target_Gate", FC.c_str(), FC.c_str());
//    TGraphErrors *geT = (TGraphErrors*) fAna->Get(name); if (!geT) cout << "Could not get " << name << endl;
//    geT unused! Is that N_{FG} / N_{FG+BG} ?

    TGraphErrors *geK = new TGraphErrors(8);
    sprintf(name, "%s_%s_%s_k", Simulation.c_str(), FC.c_str(), key.c_str());
    geK->SetName(name);
    sprintf(name, "%s, %s, correlation loss; Deposit; #it{k}", FC.c_str(), Simulation.c_str());
    geK->SetTitle(name);
    TGraphErrors *geC = new TGraphErrors(8);
    sprintf(name, "%s_%s_%s_C", Simulation.c_str(), FC.c_str(), key.c_str());
    geC->SetName(name);
    sprintf(name, "%s, %s, correction factor; Deposit; #it{C}", FC.c_str(), Simulation.c_str());
    geC->SetTitle(name);
    TGraphErrors *geG = new TGraphErrors(8);
    sprintf(name, "%s_%s_Gate_%ins", FC.c_str(), Simulation.c_str(), RIGHT); // schema: PuFC_Geant4_Gate_15ns
    geG->SetName(name);
    sprintf(name, "%s, %s, ToF gating correction; Deposit; #it{G}", FC.c_str(), Simulation.c_str());
    geG->SetTitle(name);
    for (Int_t i = 0; i < 8; i++)
    {
        // neutron scattering correction factor
        sprintf(name, "Simulation/%s/%s_%s/EffToF/%s_ProjT_%s_%i", Simulation.c_str(), FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), i+1);
        TH1D *hProjReal = (TH1D*)fAna->Get(name); if (!hProjReal) cout << "Could not get " << name << endl;

        sprintf(name, "Simulation/%s/%s_%s/EffToF/Scattered/%s_ProjT_%s_Sc_%i", Simulation.c_str(), FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), i+1);
        TH1D *hProjScat = (TH1D*)fAna->Get(name); if (!hProjScat) cout << "Could not get " << name << endl;

        sprintf(name, "Simulation/%s/%s_ideal/EffToF/%s_ProjT_ideal_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH1D *hProjIdeal = (TH1D*)fAna->Get(name); if (!hProjIdeal) cout << "Could not get " << name << endl;
        Int_t tMax = (Int_t)hProjIdeal->GetBinCenter(hProjIdeal->GetMaximumBin());
        Int_t bl = hProjIdeal->FindBin(tMax - Left(FC));
        Int_t br = hProjIdeal->FindBin(tMax + Right(FC) + 1.0) - 1;
//        cout << bl << " " << br << endl;

        /// Calculate correction factor //////////////////////////////////////////////////////
        geC->SetPoint(i, i+1, hProjIdeal->Integral(bl, br) / hProjReal->Integral(bl, br)); ///
        //////////////////////////////////////////////////////////////////////////////////////
//        cout << "Channel " << i+1 << " \t" << hProjIdeal->Integral() << " / " << hProjReal->Integral() << " * " << fTarget << " = \t" << hProjIdeal->Integral() / hProjReal->Integral() * fTarget << endl;

        /// Calculate correlation loss
        Double_t nfTot, nfSc, DnfTot, DnfSc, k, Dk;
        nfTot = hProjReal->IntegralAndError(bl, br, DnfTot);
        nfSc = hProjScat->IntegralAndError(bl, br, DnfSc);
        k = 1. - nfSc / nfTot;
        Dk = nfSc / nfTot * sqrt(pow(DnfSc / nfSc, 2) + pow(DnfTot / nfTot, 2));
        geK->SetPoint(i, i+1, k);
        geK->SetPointError(i, 0, Dk);
        cout << "Channel " << i+1 << " \t" << nfSc << " (" << DnfSc << ") " << " / " << nfTot << " (" << DnfTot << ") " << " --> \t" << 1. - nfSc / nfTot << " (" << Dk << ") " << endl;

        Double_t C = hProjIdeal->Integral(bl, br) / hProjReal->Integral(bl, br);
        Double_t D = sqrt(1.0 / hProjIdeal->GetEntries() + 1.0 / hProjReal->GetEntries());
//        if (!strcmp(Simulation.c_str(), "MCNP"))
//        {
//            Double_t x, T;
//            geT->GetPoint(i, x, T);
//            C /= T;
//        }

        if (!strcmp(key.c_str(), "real")) {
            sprintf(name, "Simulation/%s/%s_real/FitToF/%s_FitT_real_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
            TH1D *hFitReal = (TH1D*)fAna->Get(name); if (!hFitReal) cout << "Could not get " << name << endl;
            Double_t G = GatingCorrection(hProjReal, hFitReal, i, FC);
            geG->SetPoint(i, i+1, G);
            cout << i+1 << "   " << C << " +- " << C*D << "   " << G << endl;
        } else
            cout << i+1 << "   " << C << " +- " << C*D << endl;
    }
    if (!strcmp(key.c_str(), "real")) {
        Save(fAna, FC+"/Correction", geG);
        Save(fAna, FC+"/Correction", geK);
    }
    Save(fAna, FC+"/Correction", geC);
    fAna->Save();
    fAna->Close();
}//*/

void SimSB(string Simulation, string FC = "PuFC", Int_t ch = 0)
{
    cout << "Analyzing " << Simulation << " " << FC << " Shadow Bar" << endl;
    char name[128] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    if (!fAna) cout << "Could not open " << "Analysis.root" << endl;

    sprintf(name, "Simulation/%s/%s_SB/EffToF/%s_ProjT_SB_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), ch + 1);
    TH1D *hProjSB = (TH1D*) fAna->Get(name); if (!hProjSB) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/%s/%s_real/EffToF/%s_ProjT_real_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), ch + 1);
    TH1D *hProjTot = (TH1D*) fAna->Get(name); if (!hProjTot) cout << "Could not get " << name << endl;

}

void AnaSim()
{
//    MCNPtoROOT();
//    AnaSim("Geant4", "PuFC", "real");
    AnaSim("Geant4", "UFC", "real");
//    AnaSim("Geant4", "UFC", "SB");
//    AnaSim("MCNP", "PuFC", "real");
//    AnaSim("MCNP", "UFC", "real");

//    PeakPosition("NIF", "MCNP");
//    PeakPosition("NIF", "Geant4");
//    PeakPosition("UFC_NIF", "Geant4");
//    PeakPosition("UFC_SB", "Geant4", "SB");
}
#endif
