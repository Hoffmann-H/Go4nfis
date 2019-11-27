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

/*TH1D* GetProjection(TH2 *pH2, TGraphErrors *pSigma, Int_t ch, string FC = "PuFC")
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
}//*/

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

//TH1D *h = 0;
//Double_t f(Double_t *x, Double_t *par)
//{ // return: histogram h, folded by gaussian (p0 ampl, p1 mean, p2 width), evaluated at x.
//    Double_t s = 0;
//    for (Int_t j = 1; j < h->GetNbinsX()+1; j++)
//    {
//        s += h->GetBinContent(j) * par[0] * exp(-pow((x[0] - h->GetBinCenter(j) - par[1]) / par[2], 2));
//    }
//    return s + par[3];
//}

//void FitPeakForm(TFile *fAna, string Run = "NIF", Int_t ch = 0, string Simulation = "Geant4", string key = "real", Bool_t Draw = 0)
//{ // Fit simulated ToF folding to data
//    string FC = Run[0] == 'U' ? "UFC" : "PuFC";
//    char name[128] = "";
////    sprintf(name, "Simulation/%s/%s_real/EffToF/%s_ProjT_real_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), ch+1);
////    h = (TH1D*)fAna->Get(name); if (!h) cout << "Could not get " << name << endl;
//    h = GetProjection(fAna, ch, "Simulation/"+Simulation+"/"+FC+"_"+key, FC, key);

//    // With Background
//    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
//    TFile *fExp = TFile::Open(name); if (!fExp) cout << "Could not open " << name << endl;
//    sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
//    TH1D *hExp = (TH1D*)fExp->Get(name); if (!hExp) cout << "Could not get " << name << endl;
//    hExp->Sumw2();
//    Double_t bg = (hExp->Integral(Gate_0(ch, FC), Gate_a(ch, FC)) + hExp->Integral(Gate_b(ch, FC), Gate_3(ch, FC))) / (Gate_a(ch, FC) - Gate_0(ch, FC) + Gate_3(ch, FC) - Gate_b(ch, FC));
//    // Background-subtracted
////    sprintf(name, "%s/ToF/Signal/%s/H1AnaHZDRDtG_%i", FC.c_str(), Run.c_str(), ch+1);
////    TH1D *hExp = (TH1D*)fAna->Get(name); if (!hExp) cout << "Could not get " << name << endl;
////    Double_t bg = 0;

//    Double_t Range[] = {hExp->GetBinLowEdge(Gate_0(ch, FC)), hExp->GetBinLowEdge(Gate_3(ch, FC))};
//    Double_t Par[] = {(hExp->GetMaximum() - bg) / h->Integral(), 0, 1, bg};
//    sprintf(name, "f%s_%i", Run.c_str(), ch+1);
//    TF1 *fit = new TF1(name, f, Range[0], Range[1], 4);
//    fit->SetParameters(Par[0], Par[1], Par[2], Par[3]);
//    fit->FixParameter(3, bg);
////    fit->FixParameter(1, 0);

//    hExp->Fit(name, "LR0Q");
////    for (Int_t p = 0; p < 4; p++)
////        cout << fit->GetParameter(p) << " +- " << fit->GetParError(p) << endl;
//    cout << FC << "\t" << ch+1 << " \t" << fit->GetParameter(2) << " +- " << fit->GetParError(2) << endl;

//    if (Draw) {
//        sprintf(name, "c%s_%i", Run.c_str(), ch+1);
//        new TCanvas(name, "Fit folding");
//        hExp->Draw("hist");
//        fit->SetNpx(1000);
//        fit->Draw("same");
//    } else { // Save
//        Int_t N = h->GetNbinsX();
//        Double_t xmin = h->GetBinLowEdge(1);
//        Double_t xmax = h->GetBinLowEdge(N+1);
//        sprintf(name, "%s_FoldT_Fit_%i", FC.c_str(), ch+1);
//        TH1D* hFit = new TH1D(name, name, N, xmin, xmax);
//        sprintf(name, "%s, Ch. %i, Fit simulated ToF spectrum", FC.c_str(), ch+1);
//        hFit->SetTitle(name);
//        for (Int_t bin = 1; bin < N+1; bin++)
//            hFit->SetBinContent(bin, fit->Eval(hFit->GetBinCenter(bin)) /*- bg*/);
//        Save(fAna, "Simulation/"+Simulation+"/"+FC+"_"+key+"/FitToF", hFit);
//    }
//}

//void PeakPosition(string Run = "NIF", string Simulation = "Geant4", string key = "real")
//{
//    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
////    Int_t i = 2;
//    for (Int_t i = 0; i < 4; i++)
//        FitPeakForm(fAna, Run, i, Simulation, key, 1);
//    fAna->Save();
//    fAna->Close();
//}

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

Double_t GatingCorrection(TH1D *pH1Peak, Int_t ch, string FC = "PuFC", Double_t width = 15.0)
{
    Double_t bg = pH1Peak->GetBinContent(pH1Peak->GetNbinsX());
    Double_t xc = pH1Peak->GetBinCenter(pH1Peak->GetMaximumBin());
    Double_t x0 = xc - LEFT; // - width;
    Double_t x1 = xc + RIGHT; // + width;
    Int_t bin0 = pH1Peak->FindBin(x0);
    Int_t bin1 = pH1Peak->FindBin(x1);
    Double_t w0 = (pH1Peak->GetBinLowEdge(bin0 + 1) - x0) / pH1Peak->GetBinWidth(bin0); // Constant bin interpolation
    Double_t w1 = (x1 - pH1Peak->GetBinLowEdge(bin1)) / pH1Peak->GetBinWidth(bin1);
    Double_t Integral = pH1Peak->Integral() - pH1Peak->GetNbinsX() * bg;
    Double_t GatedIntegral = w0 * pH1Peak->GetBinContent(bin0) + pH1Peak->Integral(bin0 + 1, bin1 - 1) + w1 * pH1Peak->GetBinContent(bin1) - (w0 + bin1 - bin0 + 1 + w1) * bg;
//    cout << x0 << " " << x1 << " " << bin0 << " " << bin1 << " " << w0 << " " << w1 << " " << bg << " " << Integral << " " << GatedIntegral << endl;
    return Integral / GatedIntegral;
}

Double_t TargetFactor(TFile *fAna)
{
    char name[32] = "";
    sprintf(name, "Simulation/Target/E_Dir");
    TH1D *hTargetDir = (TH1D*)fAna->Get(name);
    if (!hTargetDir) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/Target/E_Tot");
    TH1D *hTargetTot = (TH1D*)fAna->Get(name);
    if (!hTargetTot) cout << "Could not get " << name << endl;
    return hTargetDir->Integral() / hTargetTot->Integral();
}

void AnaSim(string Simulation, string FC = "PuFC", string key = "real")
{ // SimToF == 0: force SimToF recalculation. 1: Use SimToF if possible. 2: Use FitToF if possible. Otherwise create it.
    cout << "Analyzing " << FC << " " << Simulation << " simulation results" << endl;
    char name[128] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    if (!fAna) cout << "Could not open " << "Analysis.root" << endl;
    Double_t fTarget = TargetFactor(fAna);

    TGraphErrors *geC = new TGraphErrors(8);
    sprintf(name, "%s_%s_%s_C", Simulation.c_str(), FC.c_str(), key.c_str());
    geC->SetName(name);
    sprintf(name, "%s, %s, correction factor; Deposit; C", FC.c_str(), Simulation.c_str());
    geC->SetTitle(name);
    TGraphErrors *geG = new TGraphErrors(8);
    sprintf(name, "%s_%s_Gate_%ins", FC.c_str(), Simulation.c_str(), RIGHT); // schema: PuFC_Geant4_Gate_15ns
    geG->SetName(name);
    sprintf(name, "%s, %s, ToF gating correction; Deposit; G", FC.c_str(), Simulation.c_str());
    geG->SetTitle(name);
    for (Int_t i = 0; i < 8; i++)
    {
        // neutron scattering correction factor
        sprintf(name, "Simulation/%s/%s_%s/EffToF/%s_ProjT_%s_%i", Simulation.c_str(), FC.c_str(), key.c_str(), FC.c_str(), key.c_str(), i+1);
        TH1D *hProjReal = (TH1D*)fAna->Get(name); if (!hProjReal) cout << "Could not get " << name << endl;
        sprintf(name, "Simulation/%s/%s_ideal/EffToF/%s_ProjT_ideal_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
        TH1D *hProjIdeal = (TH1D*)fAna->Get(name); if (!hProjIdeal) cout << "Could not get " << name << endl;

        /// Calculate correction factor //////////////////////////////////////////////////////
        geC->SetPoint(i, i+1, hProjIdeal->Integral() / hProjReal->Integral() * fTarget); /////
        //////////////////////////////////////////////////////////////////////////////////////
//        cout << "Channel " << i+1 << " \t" << hProjIdeal->Integral() << " / " << hProjReal->Integral() << " * " << fTarget << " = \t" << hProjIdeal->Integral() / hProjReal->Integral() * fTarget << endl;

        if (!strcmp(key.c_str(), "real")) {
            sprintf(name, "Simulation/%s/%s_real/FitToF/%s_FitT_real_%i", Simulation.c_str(), FC.c_str(), FC.c_str(), i+1);
            TH1D *hFitReal = (TH1D*)fAna->Get(name); if (!hFitReal) cout << "Could not get " << name << endl;
            Double_t G = GatingCorrection(hFitReal, i, FC);
            geG->SetPoint(i, i+1, G);
        }
    }
    if (!strcmp(key.c_str(), "real"))
        Save(fAna, FC+"/Correction", geG);
    Save(fAna, FC+"/Correction", geC);
    fAna->Save();
    fAna->Close();
}//*/

void AnaSim()
{
//    MCNPtoROOT();
    AnaSim("Geant4", "PuFC", "real");
    AnaSim("Geant4", "UFC", "real");
    AnaSim("Geant4", "UFC", "SB");
//    AnaSim("MCNP", "PuFC", GATE, 1);
//    PeakPosition("NIF", "MCNP");
//    PeakPosition("NIF", "Geant4");
//    PeakPosition("UFC_NIF", "Geant4");
//    PeakPosition("UFC_SB", "Geant4", "SB");
}
#endif
