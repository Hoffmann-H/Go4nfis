/// file effUmA.dat is used
#include "SaveToFile.C"

void NumberOfPuAtoms()
{
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    Double_t PuSFT2 = 6.76E10 * 365.24*24*60*60;
    Double_t DPuSFT2 = 7E8 * 365.24*24*60*60;
    Double_t rSF_Toni[] = {4.3234, 3.8945, 3.4859, 3.3611, 3.5388, 3.5293, 3.2948, 4.2508};
    Double_t DrSF_Toni[] = {0.0013, 0.0012, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011, 0.0013};
    TGraphErrors *ge = new TGraphErrors(8);
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t nAtoms = rSF_Toni[i] * PuSFT2 / log(2.0);
        Double_t DnAtoms = DrSF_Toni[i] * PuSFT2 / log(2.0);
        cout << " " << i+1 << " " << nAtoms << "+-" << DnAtoms << endl;
        ge->SetPoint(i, i+1, nAtoms);
        ge->SetPointError(i, 0, DnAtoms);
    }
    Save(fAna, "PuFC/nAtoms", ge, "PuFC_effN");
    fAna->Save();
    fAna->Close();
}

TGraphErrors* EffUmA()
{ // return a TGraphErrors with eps_nELBE * m_A from calibration with H19
  // units: mg/cm^2
    string file_name = "/home/hoffma93/Programme/ROOT/Data/effUmA.dat";
    TGraphErrors *ge = new TGraphErrors(file_name.c_str(), "%lg %lg %lg");
    if (ge == 0)
        cout << "Fehler beim Öffnen von " << file_name << endl;
    return ge;
}

void NumberOfUAtoms()
{
    Double_t Deposit_radius = 3.70; // cm
    Double_t Deposit_area = TMath::Pi() * pow(Deposit_radius, 2); // cm^2
    Double_t MolarMass = 235.3175644086;
    Double_t u = 1.660539E-24; // [g]
    TFile* fAna = TFile::Open("~/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    TGraphErrors *gEffUmA = EffUmA();
    TGraphErrors *gEffNU = new TGraphErrors(8);
    for (Int_t i = 0; i < 8; i++)
    {
        Double_t x, y, yerr;
        gEffUmA->GetPoint(i, x, y); // y = effUmA [mg/cm^2]
        yerr = gEffUmA->GetErrorY(i);
        Double_t effNU = Deposit_area * (y / 1000) / (MolarMass * u);
        Double_t DeffNU = Deposit_area * (yerr / 1000) / (MolarMass * u);
//        cout << " " << i+1 << "   " << effNU << "+-" << DeffNU << endl;
        gEffNU->SetPoint(i, i+1, effNU);
        gEffNU->SetPointError(i, 0, DeffNU);
    }
    Save(fAna, "UFC/nAtoms", gEffUmA, "effUmA");
    Save(fAna, "UFC/nAtoms", gEffNU, "UFC_effN");
    fAna->Save();
    fAna->Close();
}

void NumberOfAtoms()
{
    NumberOfPuAtoms();
    NumberOfUAtoms();
}
