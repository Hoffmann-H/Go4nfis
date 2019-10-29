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
    else
        cout << "Opened " << h->GetName() << endl;
    return h;
}

TH1D* GetProjection(TFile *f, Int_t ch, string SimulationPath = "Simulation/Geant4/PuFC_real", string FC = "PuFC", string key = "real")
{
    char name[64] = "";
    TH2D *pH2 = GetSpectrum(f, ch, SimulationPath, FC, key);
    TH1D *h = (TH1D*)pH2->ProjectionX();
    sprintf(name, "%s_ProjT_%s_%i", FC.c_str(), key.c_str(), ch + 1);
    h->SetName(name);
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
    Double_t xc = pH1Peak->GetBinCenter(pH1Peak->GetMaximumBin());
    Double_t x0 = xc - width;
    Double_t x1 = xc + width;
    Int_t bin0 = pH1Peak->FindBin(x0);
    Int_t bin1 = pH1Peak->FindBin(x1);
    Double_t w0 = (pH1Peak->GetBinLowEdge(bin0 + 1) - x0) / pH1Peak->GetBinWidth(bin0); // Constant bin interpolation
    Double_t w1 = (x1 - pH1Peak->GetBinLowEdge(bin1)) / pH1Peak->GetBinWidth(bin1);
    Double_t Integral = pH1Peak->Integral();
    Double_t GatedIntegral = w0 * pH1Peak->GetBinContent(bin0) + pH1Peak->Integral(bin0 + 1, bin1 - 1) + w1 * pH1Peak->GetBinContent(bin1);
//    cout << x0 << " " << x1 << " " << bin0 << " " << bin1 << " " << w0 << " " << w1 << " " << Integral << " " << GatedIntegral << endl;
    return Integral / GatedIntegral;
}

void AnaSim(string Simulation, string FC = "PuFC", Double_t SimToF = 0)
{ // SimToF = Gate width. 0: exact.
    char name[128] = "";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    if (!fAna) cout << "Could not open " << "Analysis.root" << endl;
    sprintf(name, "Simulation/Target/Direct");
    TH1D *hTargetDir = (TH1D*)fAna->Get(name);
    if (!hTargetDir) cout << "Could not get " << name << endl;
    sprintf(name, "Simulation/Target/Total");
    TH1D *hTargetTot = (TH1D*)fAna->Get(name);
    if (!hTargetTot) cout << "Could not get " << name << endl;
    Double_t fTarget = hTargetDir->Integral() / hTargetTot->Integral();

    TGraphErrors *geC = new TGraphErrors(8);
    sprintf(name, "%s_%s_C", FC.c_str(), Simulation.c_str());
    geC->SetName(name);
    sprintf(name, "%s, %s, correction factor; Deposit; C", FC.c_str(), Simulation.c_str());
    geC->SetTitle(name);
    TGraphErrors *geG = new TGraphErrors(8);
    sprintf(name, "%s_%s_Gate_%.0fns", FC.c_str(), Simulation.c_str(), SimToF); // schema: PuFC_Geant4_Gate_15ns
    geG->SetName(name);
    sprintf(name, "%s, %s, ToF gating correction; Deposit; G", FC.c_str(), Simulation.c_str());
    geG->SetTitle(name);
    for (Int_t i = 0; i < 8; i++)
    {
        string SimulationPath = "Simulation/"+Simulation+"/"+FC+"_real";
        TH1D *pH1real = GetProjection(fAna, i, SimulationPath, FC, "real");
        Save(fAna, SimulationPath+"/EffToF/", pH1real);
        if (SimToF != 0)
        {
            // Fold projection
            TH1D *pH1Peak = GetPeakForm(pH1real, i, FC, 1.0);
            sprintf(name, "%s_FoldT_real_%i", FC.c_str(), i+1);
            pH1Peak->SetName(name);
            // Find Gating correction
            Double_t G = GatingCorrection(pH1Peak, i, FC, SimToF);
            geG->SetPoint(i, i+1, G);
            // Scale to experimental data
            if (!strcmp(FC.c_str(), "PuFC")) sprintf(name, "PuFC/ToF/Signal/NIF/H1AnaHZDRDtG_%i", i+1);
            else                          sprintf(name, "UFC/ToF/Signal/UFC_NIF/H1AnaHZDRDtG_%i", i+1);
            TH1F *pH1Exp = (TH1F*)fAna->Get(name);
            if (!pH1Exp) cout << "Could not get " << name << endl;
            Double_t Scale = pH1Exp->Integral(Gate_1(i, FC), Gate_2(i, FC)-1) * pH1Exp->GetBinWidth(1) / pH1Peak->GetBinWidth(1) * G;
            pH1Peak->Scale(Scale);
            // Save
            Save(fAna, SimulationPath+"/SimToF/", pH1Peak);

            sprintf(name, "fG_%i", i+1);
            TF1 *f = new TF1(name, "gaus");
            f->SetRange(0, 100);
            f->SetParameters(pH1Exp->GetMaximum(), 30, 5);
            pH1Exp->Fit(name, "LR0Q");
            cout << i+1 << " " << f->GetParameter(0);
            pH1Peak->Fit(name, "LR0Q");
            cout << " " << f->GetParameter(0) << endl;
//            sprintf(name, "fPol2_%i", i+1);
//            TF1 *f2 = new TF1(name, "[0]+[1]*(x-[2])*(x-[2])");
//            f2->SetRange(pH1Exp->GetBinCenter(pH1Exp->GetMaximumBin()) - 3, pH1Exp->GetBinCenter(pH1Exp->GetMaximumBin()) + 3);
//            f2->SetParameters(pH1Exp->GetMaximum(), -10, 30);
//            pH1Exp->Fit(name, "LR0Q");
//            cout << i+1 << " " << f2->GetParameter(2) << " " << pH1Peak->GetBinCenter(pH1Peak->GetMaximumBin()) << endl;
//            Save(fAna, SimulationPath+"/SimToF/", f2);
        }

        SimulationPath = "Simulation/"+Simulation+"/"+FC+"_ideal";
        TH1D *pH1ideal = GetProjection(fAna, i, SimulationPath, FC, "ideal");
        Save(fAna, SimulationPath+"/EffToF/", pH1ideal);

        geC->SetPoint(i, i+1, pH1ideal->Integral() / pH1real->Integral() * fTarget);
    }
    if (SimToF != 0)
        Save(fAna, FC+"/Correction", geG);
    Save(fAna, FC+"/Correction", geC);
    fAna->Save();
    fAna->Close();
}//*/

void AnaSim()
{
//    MCNPtoROOT();
//    AnaSim("Geant4", "PuFC", GATE);
//    AnaSim("Geant4", "UFC", GATE);
    AnaSim("MCNP", "PuFC", GATE);
}
#endif
