/// file effUmA.dat is used
#include "SaveToFile.C"

void NumberOfPuAtoms()
{
    Double_t eff = 0.986;
    char name[128] = "";
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TFile *fSF = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SF.root", "READ");
    if (!fSF) cout << "Could not open " << "/home/homma93/Programme/Go4nfis/offline/results/SF.root" << endl;
    Double_t PuSFT2 = 6.76E10 * 365.24*24*60*60;
    Double_t DPuSFT2 = 7E8 * 365.24*24*60*60;
    Double_t rSF_06[] = {4.3275, 3.8991, 3.4868, 3.3583, 3.5389, 3.5293, 3.2918, 4.2490};
    Double_t DrSF_06[] = {0.0018, 0.0017, 0.0016, 0.0016, 0.0017, 0.0017, 0.0016, 0.0018};
    Double_t rSF_11[] = {4.3106, 3.8920, 3.4781, 3.3610, 3.5363, 3.5257, 3.2965, 4.2455};
    Double_t DrSF_11[] = {0.0025, 0.0024, 0.0022, 0.0022, 0.0023, 0.0023, 0.0022, 0.0025};
    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
    Double_t DrSF_Toni[] = {0.0013, 0.0012, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0013};
    Double_t N_Toni[] = {1.35, 1.216, 1.088, 1.049, 1.105, 1.102, 1.029, 1.327};
    Double_t DN_Toni[] = {0.023, 0.021, 0.018, 0.018, 0.019, 0.019, 0.017, 0.022};
    TGraphErrors *geFisBg = new TGraphErrors(8);
    TGraphErrors *geBgRate = new TGraphErrors(8);
    TGraphErrors *geAtoms = new TGraphErrors(8);
    sprintf(name, "PuFC/ToF/Background/NIF/BackgroundRate");
    TGraphErrors *geNIF = (TGraphErrors*)fAna->Get(name);
    if (!geNIF) cout << "Could not get " << name << endl;
    sprintf(name, "PuFC/ToF/Background/SB/BackgroundRate");
    TGraphErrors *geSB = (TGraphErrors*)fAna->Get(name);
    if (!geSB) cout << "Could not get " << name << endl;
    sprintf(name, "Histograms/Raw/Scaler/Rates/H1RawRate_47");
    TH1D *hTlive = (TH1D*)fSF->Get(name);
    if (!hTlive) cout << "Could not get " << name << endl;
    Double_t t_live = hTlive->Integral();
    Double_t x, BgRateNIF, DBgRateNIF, BgRateSB, DBgRateSB, BgRateBeam, DBgRateBeam, BgRateSF, DBgRateSF, BgRate, DBgRate;
    for (Int_t i = 0; i < 8; i++)
    {
        // Get NIF, SB, SF background
        geNIF->GetPoint(i, x, BgRateNIF);
        DBgRateNIF = geNIF->GetErrorY(i);
        geSB->GetPoint(i, x, BgRateSB);
        DBgRateSB = geSB->GetErrorY(i);
        sprintf(name, "Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", i+1);
        TH1I *hDt = (TH1I*)fSF->Get(name);
        if (!hDt) cout << "Could not get " << name << endl;
        Double_t FissionBackground = hDt->Integral();
        BgRateSF = FissionBackground / t_live;
        DBgRateSF = sqrt(FissionBackground) / t_live;

        // Save SF
        geFisBg->SetPoint(i, i+1, FissionBackground);
        geFisBg->SetPointError(i, 0, sqrt(FissionBackground));
        geBgRate->SetPoint(i, i+1, BgRateSF);
        geBgRate->SetPointError(i, 0, DBgRateSF);

        // Print values and average
        BgRateBeam = (BgRateNIF/DBgRateNIF + BgRateSB/DBgRateSB) / (1/DBgRateNIF + 1/DBgRateSB);
        DBgRateBeam = sqrt(2) / (1/DBgRateNIF + 1/DBgRateSB);
        Double_t w = (1/DBgRateNIF + 1/DBgRateSB + 1/DBgRateSF + 1/DrSF_06[i] + 1/DrSF_11[i]);
        BgRate = (BgRateNIF/DBgRateNIF + BgRateSB/DBgRateSB + BgRateSF/DBgRateSF + rSF_06[i]/DrSF_06[i] + rSF_11[i]/DrSF_11[i]) / w;
        DBgRate = sqrt(5)/w;
        sprintf(name, "%.3f(%i) & %.4f(%i) & %.4f(%i) & %.4f(%i) & %.4f(%i)",
                BgRateBeam, (Int_t)(1000*DBgRateBeam),
                BgRateSF, (Int_t)(10000*DBgRateSF),
                rSF_06[i], (Int_t)(10000*DrSF_06[i]),
                rSF_11[i], (Int_t)(10000*DrSF_11[i]),
                BgRate, (Int_t)(10000*DBgRate));
//        cout << i+1 << " & " << name << " \\\\" << endl;

        Double_t nAtoms = BgRate * PuSFT2 / log(2.0);
        Double_t DnAtoms = DBgRate * PuSFT2 / log(2.0);
        sprintf(name, "%.3f(%i) & %.3f(%i)", N_Toni[i], (Int_t)(1000*DN_Toni[i]), nAtoms*1.E-19 / eff, (Int_t)(DnAtoms*1.E-15 / eff));
        cout << i+1 << " & " << name << " \\\\" << endl;
        geAtoms->SetPoint(i, i+1, nAtoms);
        geAtoms->SetPointError(i, 0, DnAtoms);
    }
    Save(fAna, "PuFC/nAtoms/SF", geFisBg, "FissionBackground");
    Save(fAna, "PuFC/nAtoms/SF", geBgRate, "BackgroundRate");
    Save(fAna, "PuFC/nAtoms", geAtoms, "PuFC_effN");
    fAna->Save();
    fAna->Close();
}

TGraphErrors* EffUmA(string file_name = "/home/hoffma93/Programme/ROOT/Data/effUmA.dat")
{ // return a TGraphErrors with eps_nELBE * m_A from calibration with H19
  // units: mg/cm^2
    TGraphErrors *ge = new TGraphErrors(file_name.c_str(), "%lg %lg %lg");
    if (ge == 0)
        cout << "Fehler beim Öffnen von " << file_name << endl;
    return ge;
}

void NumberOfUAtoms(string StrRes = "")
{
//    Double_t Deposit_radius = 3.70; // cm
//    Double_t Deposit_area = TMath::Pi() * pow(Deposit_radius, 2); // cm^2
//    Double_t Deposit_area[] = {42.7736, 41.5982, 43.3861, 43.4792, 43.2923, 43.1283, 42.7877, 43.6696}; // cm^2
    Double_t Deposit_area[] = {43, 43, 43, 43, 43, 43, 43, 43};
    Double_t DeltaDepositArea = 0.7 / sqrt(8); // cm^2
    Double_t MolarMass = 235.3175644086;
    Double_t u = 1.660539E-24; // [g]
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    string Name = "/home/hoffma93/Programme/ROOT/Data/effUmA" + StrRes + ".dat";
    TGraphErrors *gEffUmA = EffUmA(Name);
    TGraphErrors *gEffNU = new TGraphErrors(8);
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t x, y, yerr;
        gEffUmA->GetPoint(i, x, y); // y = effUmA [mg/cm^2]
        yerr = gEffUmA->GetErrorY(i);
        Double_t effNU = Deposit_area[i] * (y / 1000) / (MolarMass * u);
        Double_t DeffNU = effNU * sqrt( pow(DeltaDepositArea / Deposit_area[i], 2) + pow(yerr / y, 2) );
        cout << " " << i+1 << "   " << effNU << "+-" << DeffNU << endl;
        gEffNU->SetPoint(i, i+1, effNU);
        gEffNU->SetPointError(i, 0, DeffNU);
    }
    Name = "effUmA" + StrRes;
    Save(fAna, "UFC/nAtoms", gEffUmA, Name);
    Name = "UFC_effN" + StrRes;
    Save(fAna, "UFC/nAtoms", gEffNU, Name);
    fAna->Save();
    fAna->Close();
}

void NumberOfAtoms()
{
    NumberOfPuAtoms();
    NumberOfUAtoms("");
}
