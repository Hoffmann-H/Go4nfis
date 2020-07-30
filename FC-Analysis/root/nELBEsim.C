#include "SaveToFile.C"
#include "/gpfs/home/hoffma93/StyleSheets/StyleSheet.C"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TRandom2.h"

TGraphErrors* GetCrossSectionData(string StrElement, Int_t Z, string StrEval)
{
    string DataPath = "/gpfs/home/hoffma93/Programme/ROOT/Data";
    char name[64] = "";
    sprintf(name, "%s/%s%i_nf_%s.dat", DataPath.c_str(), StrElement.c_str(), Z, StrEval.c_str());
    TGraphErrors *gr = new TGraphErrors(name, "%lg %lg");
    return gr;
}

TGraphErrors* GetCrossSectionToni(string format)
{
    string FilePath = "/gpfs/home/hoffma93/Programme/ROOT/Data/CS_U_Toni.dat";
    TGraphErrors *gr = new TGraphErrors(FilePath.c_str(), format.c_str());
    return gr;
}

void CheckCrossSections(string format, Int_t Z)
{
    TGraphErrors *gJEFF32 = GetCrossSectionToni(format.c_str());
    TGraphErrors *gJEFF33 = GetCrossSectionData("U", Z, "JEFF_3.3");
    new TCanvas("U233");
    gJEFF32->Draw();
    gJEFF33->SetLineColor(kRed);
    gJEFF33->Draw("same");
}

void WeightCrossSectionH19(string StrFileOutPath, string StrEval)
{//                      234U       235U      236U       238U
    Double_t IsoVec[] = {0.0003620, 0.999183, 0.0000940, 0.0003610};
    Double_t norm = IsoVec[0] + IsoVec[1] + IsoVec[2] + IsoVec[3];
    cout << "Normierung: " << norm << endl;
//    TGraphErrors *gU233 = GetCrossSectionData("U", 233, StrEval);
    TGraphErrors *gU234 = GetCrossSectionData("U", 234, StrEval);
    TGraphErrors *gU235 = GetCrossSectionData("U", 235, StrEval);
    TGraphErrors *gU236 = GetCrossSectionData("U", 236, StrEval);
    TGraphErrors *gU238 = GetCrossSectionData("U", 238, StrEval);
//    TGraphErrors *gU233 = GetCrossSectionToni("%lg %lg %lg");
//    TGraphErrors *gU234 = GetCrossSectionToni("%lg %*lg %*lg %lg %lg");
//    TGraphErrors *gU235 = GetCrossSectionToni("%lg %*lg %*lg %*lg %*lg %lg %lg");
//    TGraphErrors *gU236 = GetCrossSectionToni("%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
//    TGraphErrors *gU238 = GetCrossSectionToni("%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
    TGraphErrors *gH19 = (TGraphErrors*)gU235->Clone("H19Target");
//    gH19->SetName("H19Target");
    Int_t N = gU235->GetN();
    Double_t x, y;
    for (Int_t i = 0; i < N; i++)
    {
        gU235->GetPoint(i, x, y);
        Double_t cs = IsoVec[0] * gU234->Eval(x) +
                      IsoVec[1] * y +
                      IsoVec[2] * gU236->Eval(x) +
                      IsoVec[3] * gU238->Eval(x);
        gH19->SetPoint(i, x, cs / norm);
    }
    TFile *f = TFile::Open(StrFileOutPath.c_str(), "UPDATE");
    SaveToFile(f, StrEval, gH19);
    f->Save(); f->Close();
}

void WeightCrossSectionUFC(string StrFileOutPath, string StrEval)
{//                      233U       234U     235U   236U     238U
    Double_t IsoVec[] = {0.0000011, 0.00459, 0.904, 0.00401, 0.0912};
    Double_t norm = IsoVec[0] + IsoVec[1] + IsoVec[2] + IsoVec[3] + IsoVec[4];
    cout << "Normierung: " << norm << endl;
    TGraphErrors *gU233 = GetCrossSectionData("U", 233, StrEval);
    TGraphErrors *gU234 = GetCrossSectionData("U", 234, StrEval);
    TGraphErrors *gU235 = GetCrossSectionData("U", 235, StrEval);
    TGraphErrors *gU236 = GetCrossSectionData("U", 236, StrEval);
    TGraphErrors *gU238 = GetCrossSectionData("U", 238, StrEval);
//    TGraphErrors *gU233 = GetCrossSectionToni("%lg %lg %lg");
//    TGraphErrors *gU234 = GetCrossSectionToni("%lg %*lg %*lg %lg %lg");
//    TGraphErrors *gU235 = GetCrossSectionToni("%lg %*lg %*lg %*lg %*lg %lg %lg");
//    TGraphErrors *gU236 = GetCrossSectionToni("%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
//    TGraphErrors *gU238 = GetCrossSectionToni("%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
    TGraphErrors *gUFC = (TGraphErrors*)gU235->Clone();
    gUFC->SetName("UFCTarget");
    Int_t N = gU235->GetN();
    Double_t x, y;
    for (Int_t i = 0; i < N; i++)
    {
        gU235->GetPoint(i, x, y);
        Double_t cs = IsoVec[0] * gU233->Eval(x) +
                      IsoVec[1] * gU234->Eval(x) +
                      IsoVec[2] * y +
                      IsoVec[3] * gU236->Eval(x) +
                      IsoVec[4] * gU238->Eval(x);
        gUFC->SetPoint(i, x, cs / norm);
    }
    TFile *f = TFile::Open(StrFileOutPath.c_str(), "UPDATE");
    SaveToFile(f, StrEval, gUFC);
    f->Save(); f->Close();
}

void GetJEFF()
{/// Create and save deposits' cross section graphs
    string StrFileOutPath = "/home/hoffma93/TrackLength/FissionXS.root";
    WeightCrossSectionH19(StrFileOutPath, "JEFF_3.3");
    WeightCrossSectionUFC(StrFileOutPath, "JEFF_3.3");
}
/*
void WeightG4withXS(TFile *fSpec, string StrEval = "JEFF_3.3")
{
    char name[64] = "";
    TH1F *pH1src = (TH1F*) fSpec->Get("Source/nEnergy/Source_Ekin");
    TFile *fXS = TFile::Open("~/TrackLength/FissionXS.root", "READ");
    if (!fXS) cout << "Could not get " << "Fission cross section data" << endl;
    string FCs[] = {"UFC", "H19"};
    Int_t NumCh[] = {8, 10};
    for (Int_t iFC = 0; iFC < 2; iFC++)
    {
	string FC = FCs[iFC];
        // UFC: r = 3,7cm, H19: r = 3,8cm
        Double_t r = strcmp(FC.c_str(), "H19") ? 3.7 : 3.8;
        Double_t A = TMath::Pi() * r * r; // [A] = cm^2
	//cout << FC << ". A = " << A << endl;
	TGraphErrors *pGrXS;
	sprintf(name, "%s/%sTarget", StrEval.c_str(), FC.c_str());
        fXS->GetObject(name, pGrXS);
        if (!pGrXS) cout << "Could not get " << name << endl;
	Int_t NumBinX, NumBinY;
        Double_t Ekin, XS, Content, EContent;
	for (Int_t i = 0; i < NumCh[iFC]; i++)
	{
            sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
	    cout << "Weight " << name << " by " << pGrXS->GetName() << endl;
	    TH2D *pH2In = (TH2D*)fSpec->Get(name); if (!pH2In) cout << "Could not get " << name << endl;
	    sprintf(name, "%s/EToFvsEkin/Scattered/%s_EToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
	    TH2D *pH2ScIn = (TH2D*)fSpec->Get(name); if (!pH2ScIn) cout << "Could not get " << name << endl;

	    //////////////////////////////
	    TH2F *pH2 = (TH2F*)pH2In->Clone();
            sprintf(name, "%s_EToFvsEkin_wXS_Ch.%i", FC.c_str(), i+1);
            pH2->SetName(name);
            pH2->SetDirectory(0);
	    TH2F *pH2Sc = (TH2F*)pH2ScIn->Clone();
            sprintf(name, "%s_EToFvsEkin_wXS_Sc_Ch.%i", FC.c_str(), i+1);
            pH2Sc->SetName(name);
            pH2Sc->SetDirectory(0);
	    NumBinX = pH2In->GetNbinsX();       //E(t) bins
            NumBinY = pH2In->GetNbinsY();       //Ekin bins
            for (int y_i=1; y_i<=NumBinY; y_i++)
            {   //get neutron kinetic energy of corresponding bin y_i
                Ekin    = pH2In->GetYaxis()->GetBinCenter(y_i);
                //get cross section @ E=Ekin
                XS      = pGrXS->Eval(Ekin);
                //loop over all Time-of-Flight bins
                for (int x_i=1; x_i<=NumBinX; x_i++)
                {   //multiply bin content with xs section to get fission rate
                    Content = pH2In->GetBinContent(x_i, y_i) * XS;
                    EContent= pH2In->GetBinError(x_i, y_i) * XS;
                    pH2->SetBinContent(x_i, y_i, Content);
                    pH2->SetBinError(x_i, y_i, EContent);
		    Content = pH2ScIn->GetBinContent(x_i, y_i) * XS;
                    EContent= pH2ScIn->GetBinError(x_i, y_i) * XS;
                    pH2Sc->SetBinContent(x_i, y_i, Content);
                    pH2Sc->SetBinError(x_i, y_i, EContent);
                }
            }
	    //cout << pH2->Integral() << " / " << A << " / " << pH1src->Integral() << endl;
	    pH2->Scale( 1. / A / pH1src->Integral() );
	    //cout << pH2->Integral() << endl;
	    pH2Sc->Scale( 1. / A / pH1src->Integral() );
            SaveToFile(fSpec, FC + "/EToFvsEkin_wXS", pH2);
	    SaveToFile(fSpec, FC + "/EToFvsEkin_wXS/Scattered", pH2Sc);
	    //////////////////////////////
	}
	cout << "Weight " << pH1src->GetName() << " by " << pGrXS->GetName() << endl;
	TH1F *pH1vac = (TH1F*) pH1src->Clone();
	sprintf(name, "%s_Ekin", FC.c_str());
	pH1vac->SetName(name);
	pH1vac->SetDirectory(0);
	NumBinX = pH1src->GetNbinsX();
	for (int x_i = 1; x_i < NumBinX; x_i++)
        {
            Ekin = pH1src->GetXaxis()->GetBinCenter(x_i);
	    XS = pGrXS->Eval(Ekin);
	    Content = pH1src->GetBinContent(x_i) * XS;
            EContent = pH1src->GetBinError(x_i) * XS;
            pH1vac->SetBinContent(x_i, Content);
	    pH1vac->SetBinError(x_i, EContent);
	}
	pH1vac->Scale( 1. / A / pH1src->Integral() );
	SaveToFile(fSpec, "Source/nfEnergy", pH1vac);
//	fSpec->Save();
    }
    fSpec->Save();
}

void DoWeightG4withXS(string FileName, string StrEval = "JEFF_3.3")
{
    TFile *f = TFile::Open(FileName.c_str(), "UPDATE");
    if (!f) cout << "Could not open " << FileName << endl;
    WeightG4withXS(f, StrEval);
    f->Save();
    f->Close();
}
*/
TH2F* GetGeant4Data(TFile *f, string FC, bool scat, Int_t Channel)
{   //function to get E_ToF vs. Ekin correlation histograms from a G4 simulation
    string InputPath;

    if (scat)
    {   if (FC=="H19")
            InputPath = FC + "/EToFvsEkin/Scattered/" + FC + "_EToFvsEkin_Sc_Ch."
                        + std::to_string(Channel+1);
        else
            InputPath = FC + "/EToFvsEkin/Scattered/" + FC + "_EToFvsEkin_Sc_Ch."
                        + std::to_string(Channel+1);
    } else {
        if (FC=="H19")
            InputPath = FC + "/EToFvsEkin/" + FC + "_EToFvsEkin_Ch." +
                        std::to_string(Channel+1);
        else
            InputPath = FC + "/EToFvsEkin/" + FC + "_EToFvsEkin_Ch." +
                        std::to_string(Channel+1);
    }
    cout << InputPath << endl;

    TH2F *pH2 = (TH2F*) f->Get(InputPath.c_str())->Clone();
        pH2->SetDirectory(0);
    return pH2;
}

TH1F* GetSourceSpectrum(TFile *f)
{
    string InputPath = "Source/nEnergy/Source_Ekin";
    TH1F *pH1source = (TH1F*)f->Get(InputPath.c_str())->Clone();
    pH1source->SetDirectory(0);
    return pH1source;
}

TH2F* WeightWithXS(TH2F *pH2In, TGraphErrors *pGrXS)
{   //function to weight the E_ToF vs Ekin histograms with the fission cross
    //section to get a unnormalized fission rate
    //-> weight content in ToF bin by the cross section at energy bin
    Int_t NumBinX = pH2In->GetNbinsX();       //E(t) bins
    Int_t NumBinY = pH2In->GetNbinsY();       //Ekin bins
    Double_t Ekin, XS, Content, EContent;

    TH2F *pH2 = (TH2F*)pH2In->Clone(); pH2->Clear();

    for (int y_i=1; y_i<=NumBinY; y_i++)
    {   //get neutron kinetic energy of corresponding bin y_i
        Ekin    = pH2In->GetYaxis()->GetBinCenter(y_i);
        //get cross section @ E=Ekin
        XS      = pGrXS->Eval(Ekin);
        //loop over all Time-of-Flight bins
        for (int x_i=1; x_i<=NumBinX; x_i++)
        {   //multiply bin content with xs section to get fission rate
            Content = pH2In->GetBinContent(x_i, y_i) * XS;
            EContent= pH2In->GetBinError(x_i, y_i) * XS;
            pH2->SetBinContent(x_i, y_i, Content);
            pH2->SetBinError(x_i, y_i, EContent);
        }
    }
    return pH2;
}

TH1F* WeightTH1withXS(TH1F *pH1In, TGraphErrors *pGrXS)
{   //function to weight a Ekin histogram with the fission cross
    //section to get a unnormalized fission rate

    Int_t NumBinX = pH1In->GetNbinsX();       //Ekin bins
    Double_t Ekin, XS, Content, EContent;

    TH1F *pH1 = (TH1F*)pH1In->Clone(); pH1->Clear();

    for (int x_i=1; x_i<=NumBinX; x_i++)
    {   //get neutron kinetic energy of corresponding bin x_i
        Ekin    = pH1In->GetXaxis()->GetBinCenter(x_i);
        //get cross section @ E=Ekin
        XS      = pGrXS->Eval(Ekin);
        //multiply bin content with cross section to get fission rate
        Content = pH1In->GetBinContent(x_i) * XS;
        EContent= pH1In->GetBinError(x_i) * XS;
        pH1->SetBinContent(x_i, Content);
        pH1->SetBinError(x_i, EContent);
    }
    return pH1;
}

TGraphErrors* GetXS(string StrXSFileInPath, string StrFC, string StrEval)
{
    TFile *f = new TFile(StrXSFileInPath.c_str(), "READ");
    TGraphErrors *pGrXS;
    f->GetObject((StrEval+"/"+StrFC+"Target").c_str(), pGrXS); if (!pGrXS) cout << "Could not get " << StrEval << "/" << StrFC << "Target" << endl;
    f->Close();
    return pGrXS;
}

TH1F* GetCorrection(TFile *f, string FC, uint ch)
{   //method to determine a Correction factor C from geometric and vacuum simulation
    //C(E(t)) = vacuum fission rate at E(t) / total detected fission rate at E(t)
    cout << "Calculate correction factor for " << FC << " ch." << ch+1 << endl;
    char name[64] = "";

    //get Energy-ToF-Correlation histogram of all G4 simulated neutrons
    sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2F* pH2Tot = (TH2F*)f->Get(name); if (!pH2Tot) cout << "Could not get " << name << endl;
    //get source spectrum
    sprintf(name, "Source/nEnergy/%s_Ekin", FC.c_str());
    TH1F *pH1src = (TH1F*)f->Get(name); if (!pH1src) cout << "Could not get " << name << endl;
    //get projection to the E(t) axis
    sprintf(name, "G4Total_%s_Ch.%i_Hist", FC.c_str(), ch+1);
    TH1F* pFCTot = (TH1F*)pH2Tot->ProjectionX(name, 0, -1, "e");
    //calculate Correction factor C(E(t))
    sprintf(name, "C_%s_%i", FC.c_str(), ch+1); // This is the final name of the correction factor histogram
    TH1F* pFCCorr = (TH1F*)pFCTot->Clone(name); //pFCCorr->Clear();
    pFCCorr->Divide(pH1src, pFCTot, 1.0, 1.0);
    pFCCorr->SetDirectory(0);
    return pFCCorr;
}

TH1F* GetCorrLoss(TFile *f, string FC, uint ch)
{
    // k = unscattered neutron fission rate / total detected neutron fission rate
    cout << "Calculate correlation loss factor for " << FC << " ch." << ch+1 << endl;
    char name[64] = "";

    // get E(t)-Ekin-Correlation histogram of G4 simulated total fission
    sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2F* pH2Tot = (TH2F*)f->Get(name); if (!pH2Tot) cout << "Could not get " << name << endl;
    // get E(t)-Ekin-Correlation histogram of G4 simulated scattered fission
    sprintf(name, "%s/EToFvsEkin/Scattered/%s_EToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2F* pH2Scat = (TH2F*)f->Get(name); if (!pH2Scat) cout << "Could not get " << name << endl;
    // get E(t)-Ekin-Correlation histogram of G4 un-scattered fission
    TH2F* pH2UnScat = (TH2F*)pH2Tot->Clone("UnScat");
        pH2UnScat->Add(pH2Tot, pH2Scat, 1., -1.);
    // get projections to the E(t) axis
    TH1F *pFCUnScat, *pFCTot, *pFCCorr;
    sprintf(name, "G4UnScattered_%s_Ch.%i_Hist", FC.c_str(), ch+1);
    pFCUnScat     = (TH1F*)pH2UnScat->ProjectionX(name, 0, -1, "e");
    sprintf(name, "G4Total_%s_Ch.%i_Hist", FC.c_str(), ch+1);
    pFCTot      = (TH1F*)pH2Tot->ProjectionX(name, 0, -1, "e");
    // divide
    pFCCorr = (TH1F*)pFCUnScat->Clone("CorrLoss");
        pFCCorr->Divide(pFCUnScat, pFCTot, 1.0, 1.0);
    return pFCCorr;
}

TH1F* GetTransm(TFile *f, string FC, uint ch)
{
    //T = direct neutron rate / source neutron rate
    cout << "Calculate neutron transmission factor for " << FC << " ch." << ch+1 << endl;
    char name[64] = "";

    //get Energy-ToF-Correlation histogram of all G4 simulated neutrons
    sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2F* pH2Tot = (TH2F*)f->Get(name); if (!pH2Tot) cout << "Could not get " << name << endl;
    //get Energy-ToF-Correlation histogram of all scattered neutrons
    sprintf(name, "%s/EToFvsEkin/Scattered/%s_EToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2F* pH2Scat = (TH2F*)f->Get(name); if (!pH2Scat) cout << "Could not get " << name << endl;
    //get Energy-ToF-Correlation histogram of all un-scattered neutrons
    TH2F* pH2UnScat = (TH2F*)pH2Tot->Clone("UnScat");
        pH2UnScat->Add(pH2Tot, pH2Scat, 1., -1.);
    //get neutron source spectrum
    sprintf(name, "Source/nEnergy/%s_Ekin", FC.c_str());
    TH1F *pFCSource = (TH1F*)f->Get(name); if (!pFCSource) cout << "Could not get " << name << endl;
    //get projection to the E(t) axis
    TH1F *pFCUnScat, *pFCTransm;
    sprintf(name, "G4UnScattered_%s_Ch.%i_Hist", FC.c_str(), ch+1);
    pFCUnScat = (TH1F*)pH2UnScat->ProjectionY(name, 0, -1, "e");
    pFCTransm = (TH1F*)pFCUnScat->Clone("CorrLoss");
    pFCTransm->Divide(pFCUnScat, pFCSource, 1., 1.);
    return pFCTransm;
}

void CreateColorGradient()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
/*
TH1F* GetCorrLoss(string StrSimuFile, string StrFC, uint ch, Bool_t Draw = 0)
{   //method to determine a Correction factor k concerning the loss of kinetic
    //energy and ToF of scattered neutrons
    //k = unscattered neutron fission rate / total detected neutron fission rate

    //get Energy-ToF-Correlation histogram of all G4 simulated neutrons
    TH2F* h2Tot    = GetGeant4Data(StrSimuFile, StrFC, false, ch);
    //get Energy-ToF-Correlation histogram of all scattered neutrons
    TH2F* h2Scat   = GetGeant4Data(StrSimuFile, StrFC, true,  ch);
    //get Energy-ToF-Correlation histogram of all un-scattered neutrons
    TH2F* h2UnScat = (TH2F*)h2Tot->Clone("UnScat");
        h2UnScat->Add(h2Tot, h2Scat, 1., -1.);
//    cout << pH2Tot->Integral() << " - " << pH2Scat->Integral() << " -> " << pH2UnScat->Integral() << endl;

    //convert total and scattered histogram to fission rate histograms
    TH2F* pH2Tot      = WeightWithXS(h2Tot,   "~/TrackLength/FissionXS.root", StrFC, "JEFF_3.3"); pH2Tot->SetDirectory(0); h2Tot->Delete();
    TH2F* pH2Scat     = WeightWithXS(h2Scat,  "~/TrackLength/FissionXS.root", StrFC, "JEFF_3.3"); pH2Scat->SetDirectory(0); h2Scat->Delete();
    TH2F* pH2UnScat   = WeightWithXS(h2UnScat,"~/TrackLength/FissionXS.root", StrFC, "JEFF_3.3"); pH2UnScat->SetDirectory(0); h2UnScat->Delete();

    //get projection to the E(t) axis
    TH1F *pFCUnScat, *pFCTot, *pFCCorr;
    string ProjName;
    ProjName = "G4UnScattered_"+StrFC+"_Ch."+std::to_string(ch+1)+"_Hist";
    pFCUnScat     = (TH1F*)pH2UnScat->ProjectionX(ProjName.c_str(), 0, -1, "e");
    ProjName = "G4Total_"+StrFC+"_Ch."+std::to_string(ch+1)+"_Hist";
    pFCTot      = (TH1F*)pH2Tot->ProjectionX(ProjName.c_str(), 0, -1, "e");
    pFCCorr = (TH1F*)pFCUnScat->Clone("CorrLoss");
        pFCCorr->Divide(pFCUnScat, pFCTot, 1.0, 1.0);

    if (!Draw)
        return pFCCorr;

    TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(60);
    TLatex T2; T2.SetTextFont(42); T2.SetTextSize(0.07);

    CreateColorGradient();

    TGaxis *axis[4];

    TCanvas *c1 = new TCanvas();
    c1->Divide(2,1);
    c1->cd(1); c1->cd(1)->SetLogz();
    c1->cd(1)->SetMargin(0.18, 0.05, 0.14, 0.035);
        pH2Scat->GetXaxis()->SetRangeUser(0., 10.);
        pH2Scat->GetYaxis()->SetRangeUser(0., 10.);
        pH2Scat->GetXaxis()->SetTitle("");
        pH2Scat->GetYaxis()->SetTitle("");
        pH2Scat->SetMinimum(0.5);
        pH2Scat->SetMaximum(3e5);
        pH2Scat->SetStats(0);
        pH2Scat->DrawCopy("A COL");
        Tl.DrawText(1, 9, "(a) scattered");

        gPad->Modified();
        gPad->Update();

        axis[0] = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
              gPad->GetUxmax(), gPad->GetUymin(),0,gPad->GetUxmax(),504,"L");

        axis[1] = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
              gPad->GetUxmin(), gPad->GetUymax(),0,gPad->GetUymax(),504);

    c1->cd(2); c1->cd(2)->SetLogz();
    c1->cd(2)->SetMargin(0.05, 0.18, 0.14, 0.035);
        pH2Tot->GetXaxis()->SetRangeUser(0., 10.);
        pH2Tot->GetYaxis()->SetRangeUser(0., 10.);
        pH2Tot->GetXaxis()->SetTitle("");
        pH2Tot->GetYaxis()->SetTitle("");
        pH2Tot->GetZaxis()->SetLabelSize(0.07);
        pH2Tot->SetMinimum(0.5);
        pH2Tot->SetMaximum(3e5);
        pH2Tot->SetStats(0);
        pH2Tot->DrawCopy("A COLZ");
        Tl.DrawText(1, 9, "(b) total");

        gPad->Modified();
        gPad->Update();

        axis[2] = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
              gPad->GetUxmax(), gPad->GetUymin(),0,gPad->GetUxmax(),504,"L");
        axis[3] = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
              gPad->GetUxmin(), gPad->GetUymax(),0,gPad->GetUymax(),504, "U");

        for (int i=0; i<4; i++)
        {
            axis[i]->SetLabelSize(0.08);
            axis[i]->SetTitleSize(0.09);
            axis[i]->SetTextFont(132);
            axis[i]->SetLabelFont(42);


            if (i%2==0)
            {   axis[i]->SetTitle("#it{E}_{n}(#it{t}) / MeV");
                axis[i]->SetTitleOffset(0.75);
            }
            if (i==1)
            {    axis[i]->SetTitle("#it{E}_{n} / MeV");
                 axis[i]->SetTitleOffset(1.0);
                 axis[i]->SetLabelOffset(0.015);
            }
        }

        c1->cd(1);
            axis[0]->Draw();
            axis[1]->Draw();
        c1->cd(2);
            axis[2]->Draw();
            axis[3]->Draw();


    pFCCorr->Divide(pFCCorr, pFCTot, 1.0, 1.0);
//    pFCCorr->DrawCopy("HIST");
    pH2Tot->SetMaximum(3e5);
    pH2Tot->GetXaxis()->SetRangeUser(0., 10.);
    pH2Tot->GetYaxis()->SetRangeUser(0., 10.);
//    pH2Tot->GetYaxis()->SetLabelOffset(1.1);
    pH2Tot->GetXaxis()->SetTitle("#it{E}_{n}(#it{t}) / MeV");
    pH2Tot->GetYaxis()->SetTitle("#it{E}_{n} / MeV");
//    pH2Tot->DrawCopy("COLZ");

    string glStrSimPreFix = "G4";
    string ObjName = glStrSimPreFix +
                     "CorrLoss_"+StrFC+"_Ch." + std::to_string(ch+1);
    string ObjTitle= glStrSimPreFix +
                    " correlation loss factor #it{k}, " + StrFC +
                     "Ch. " + std::to_string(ch+1);

    pFCCorr->SetName(ObjName.c_str());
    pFCCorr->SetTitle(ObjTitle.c_str());

    delete pH2Scat, pH2UnScat, pH2Tot;

    return pFCCorr;
}
*/
TGraphErrors *TH1toTGraphError(TH1F *pH)
{
    TGraphErrors *pGrE = new TGraphErrors(pH);
    char name[64] = "";
    sprintf(name, "Gr_%s", pH->GetName());
    pGrE->SetName(name);

    for (int i=1; i<pH->GetNbinsX(); i++)
    {   pGrE->SetPoint(i-1, pH->GetBinCenter(i), pH->GetBinContent(i));
        pGrE->SetPointError(i-1, 0, pH->GetBinError(i));
    }
    return pGrE;
}

void CorrectionRefH19(TFile *f)
{
    string Name = "";
    TH1F *pH1_H19[10], *pH1_sumH19, *pH1_UFC[8], *pH1_Ref[8];
    TGraphErrors *ge_sumH19, *ge_Ref[8];
    for (Int_t channel = 0; channel < 10; channel++)
    {
        Name = "H19/Correction/C_H19_"+std::to_string(channel+1);
        pH1_H19[channel] = (TH1F*)f->Get(Name.c_str());
        if (!pH1_H19[channel]) cout << "Could not get " << Name << endl;
        if (channel==0) {
            pH1_sumH19 = (TH1F*)pH1_H19[0]->Clone("C_H19");
            if (!pH1_sumH19) cout << "Could not get " << Name << endl;
            pH1_sumH19->SetDirectory(0);
        } else
            pH1_sumH19->Add(pH1_H19[channel], 1.0);
    }
    pH1_sumH19->Scale(0.1);
    Name = "H19/Correction";
    SaveToFile(f, Name, pH1_sumH19);
    ge_sumH19 = TH1toTGraphError(pH1_sumH19);
    Name = "H19/Correction/Graph";
    SaveToFile(f, Name, ge_sumH19);
    for (Int_t channel = 0; channel < 8; channel++)
    {
        Name = "UFC/Correction/C_UFC_"+std::to_string(channel+1);
        pH1_UFC[channel] = (TH1F*)f->Get(Name.c_str());
        if (!pH1_UFC[channel]) cout << "Could not get " << Name << endl;
        Name = "C_UFC_RefH19_"+std::to_string(channel+1);
        pH1_Ref[channel] = (TH1F*)pH1_UFC[channel]->Clone(Name.c_str());
        pH1_Ref[channel]->SetDirectory(0);
        pH1_Ref[channel]->Divide(pH1_UFC[channel], pH1_sumH19, 1.0, 1.0);
        if (!pH1_Ref[channel]) cout << "Could not get " << Name << endl;
        SaveToFile(f, "UFC_RefH19", pH1_Ref[channel]);
	ge_Ref[channel] = TH1toTGraphError(pH1_Ref[channel]);
	SaveToFile(f, "UFC_RefH19/Graph", ge_Ref[channel]);
    }
}

void DoCorrectionRefH19(string Filename)
{
    TFile *f = TFile::Open(Filename.c_str(), "UPDATE");
    if (!f) cout << "Could not open " << Filename << endl;
    CorrectionRefH19(f);
    f->Save();
    f->Close();
}

void TransmCorrLoss(TFile *f)
{
    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        for (Int_t channel = 0; channel < NumCh[i_FC]; channel++)
        {
            TH1F *pk = GetCorrLoss(f, FC, channel);
            pk->SetDirectory(0);
            pk->SetName(("k_"+FC+"_"+std::to_string(channel+1)).c_str());
            SaveToFile(f, FC+"/CorrLoss", pk);
            TGraphErrors *gk = TH1toTGraphError(pk);
            SaveToFile(f, FC+"/CorrLoss/Graph", gk);

            TH1F *pT = GetTransm(f, FC, channel);
            pT->SetDirectory(0);
            pT->SetName(("T_"+FC+"_"+std::to_string(channel+1)).c_str());
            SaveToFile(f, FC+"/Transmission", pT);
            TGraphErrors *gT = TH1toTGraphError(pT);
            SaveToFile(f, FC+"/Transmission/Graph", gT);
        }
    }
    f->Save();
}

void DoTransmCorrLoss(string StrHistoFileName)
{ // file name with .root ending
    TFile *f = TFile::Open((StrHistoFileName).c_str(), "UPDATE"); // histogrammed simulation data
    if (!f) cout << "Could not open " << StrHistoFileName << endl;
    TransmCorrLoss(f);
    f->Save();
    f->Close();
}

void Correction(TFile *f)
{
    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        for (Int_t channel = 0; channel < NumCh[i_FC]; channel++)
        {
            TH1F *pC = GetCorrection(f, FC, channel);
            pC->SetDirectory(0);
            pC->SetName(("C_"+FC+"_"+std::to_string(channel+1)).c_str());
            SaveToFile(f, FC+"/Correction", pC);
	    TGraphErrors *gC = TH1toTGraphError(pC);
        SaveToFile(f, FC+"/Correction/Graph", gC);
        }
    }
    f->Save();
}

void DoCorrection(string StrHistoFileName)
{ // file name with .root ending
    TFile *f = TFile::Open((StrHistoFileName).c_str(), "UPDATE"); // histogrammed simulation data
    if (!f) cout << "Could not open " << StrHistoFileName << ".root" << endl;
    Correction(f);
    f->Save();
    f->Close();
}

void avWeight(string FileHit, string FileTL)
{
    TFile *fHit = TFile::Open(FileHit.c_str(), "READ"); if (!fHit) cout << "Could not get " << FileHit << endl;
    TFile *fTL = TFile::Open(FileTL.c_str(), "UPDATE"); if (!fTL) cout << "Could not get " << FileTL << endl;

    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        for (Int_t channel = 0; channel < NumCh[i_FC]; channel++)
        {
            for (Int_t sc = 0; sc < 2; sc++)
            {
                string Name = "";
                // get E-E(t) histograms: pH2_w weighted by 1/cos(theta), pH2 unweighted
                if (sc) Name = FC + "/ToFvsEkin/Scattered/" + FC + "_ToFvsEkin_Sc_Ch." + std::to_string(channel+1);
                else    Name = FC + "/ToFvsEkin/" + FC + "_ToFvsEkin_Ch." + std::to_string(channel+1);
                TH2F *pH2_w = (TH2F*) fTL->Get(Name.c_str()); if (!pH2_w) cout << "Could not get " << Name << endl;
                TH2F *pH2 = (TH2F*) fHit->Get(Name.c_str()); if (!pH2) cout << "Could not get " << Name << endl;

                // make projections
                TH1F *px, *px_w, *py, *py_w;
                Name = "pE(t)_w_"+FC+std::to_string(channel+1);
                px_w = (TH1F*)pH2_w->ProjectionX(Name.c_str(), 0, -1, "e");
                Name = "pE(t)_"+FC+std::to_string(channel+1);
                px = (TH1F*)pH2->ProjectionX(Name.c_str(), 0, -1, "e");
                Name = "pEkin_w_"+FC+std::to_string(channel+1);
                py_w = (TH1F*)pH2_w->ProjectionY(Name.c_str(), 0, -1, "e");
                Name = "pEkin_"+FC+std::to_string(channel+1);
                py = (TH1F*)pH2->ProjectionY(Name.c_str(), 0, -1, "e");

                // divide x projections
                if (sc) Name = FC+"_WvsToF_Sc_Ch."+std::to_string(channel+1);
                else    Name = FC+"_WvsToF_Ch."+std::to_string(channel+1);
                TH1F *pWx  = (TH1F*)px_w->Clone(Name.c_str());
                pWx->SetDirectory(0);
                pWx->Divide(px_w, px, 1.0, 1.0);
                pWx->GetYaxis()->SetTitle("<#font[12]{W}>");
                if (sc) Name = "Weight/Scattered";
                else    Name = "Weight";
                SaveToFile(fTL, Name, pWx);
                pWx->Delete();

                // divide y projections
                if (sc) Name = FC+"_WvsEkin_Sc_Ch."+std::to_string(channel+1);
                else    Name = FC+"_WvsEkin_Ch."+std::to_string(channel+1);
                TH1F *pWy = (TH1F*)py_w->Clone(Name.c_str());
                pWy->SetDirectory(0);
                pWy->Divide(py_w, py, 1.0, 1.0);
                pWy->GetYaxis()->SetTitle("<#font[12]{W}>");
                if (sc) Name = "Weight/Scattered";
                else    Name = "Weight";
                SaveToFile(fTL, Name, pWy);
                pWy->Delete();
            }
        }
    }
    fTL->Save();
    fTL->Close();
}//*/

void script(string StrFileInPath = "~/TrackLength/G4UFCvsH19_1E9_hist_w.root")
{ /// Copy the 9 FC channels' Correction factor histograms
  /// from StrFileInPath into the 2016/06 analysis root file

    TFile *f = TFile::Open(StrFileInPath.c_str()); if (!f) cout << "Could not open " << StrFileInPath << endl;
    TH1F *h[9];
    for (Int_t i = 0; i < 8; i++)
    {
        h[i] = (TH1F*)f->Get(("UFC/Correction/C_UFC_"+std::to_string(i+1)).c_str());
        h[i]->SetDirectory(0);
        h[i]->SetName(("G4Correction_UFC_Ch."+std::to_string(i+1)).c_str());
    }
    h[8] = (TH1F*)f->Get("H19/Correction/C_H19");
    h[8]->SetDirectory(0);
    h[8]->SetName("G4Correction_H19");

    string StrFileOutPath = "~/FC-Analysis/2016.06/results/nfis_2016.06.root";
    TFile *g = TFile::Open(StrFileOutPath.c_str(), "UPDATE"); if (!f) cout << "Could not open " << StrFileOutPath << endl;
    for (Int_t i = 0; i < 9; i++)
    {
        SaveToFile(g, "Analysis/Simulation/Geant4/Correction", h[i]);
        cout << "Saved " << h[i]->GetName() << endl;
    }
    g->Save();
    g->Close();
}

TH2D* MakeEvsT(string file_to_read, Bool_t draw, string FC = "UFC", string key = "real", Int_t ch = 0)
{ // Get MCNP 2D histogram
    cout << "Reading MCNP result " << file_to_read << endl;
    char name[64] = "";
    string PrototypeFile = "~/TrackLength/G4UFCvsH19_1E8_hist.root";

    // Open tabular as graph
    TGraphErrors *gTE = new TGraphErrors(file_to_read.c_str(), "%lg %lg %lg %lg");
    Int_t N = gTE->GetN();
    cout << N << " entries" << endl;

    // Open a 2D histogram as axis prototype
    TFile *f = TFile::Open(PrototypeFile.c_str()); if (!f) cout << "Could not open " << PrototypeFile << endl;
    sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH2F *pProto = (TH2F*)f->Get(name); if (!pProto) cout << "Could not get prototype histogram file " << name << endl;
    sprintf(name, "%s_ToFvsEkin_%s_Ch.%i", FC.c_str(), key.c_str(), ch + 1);
    TH2D *pHist = (TH2D*)pProto->Clone(name);
    sprintf(name, "%s, ToF vs Ekin, %s, ch.%i", FC.c_str(), key.c_str(), ch + 1);
    pHist->SetTitle(name);
    pHist->GetXaxis()->SetTitle("#font[12]{t} / ns");
    pHist->GetYaxis()->SetTitle("#font[12]{E} / MeV");

    Double_t T, E, binwx, binwy, bT, bE, C, DC;
    Int_t binT, binE;
    TRandom2 *rand = new TRandom2();

    pHist->Scale(0.0);

    for (Int_t i = 0; i < N; i++)
    {
        // read from graph
        gTE->GetPoint(i, T, E);
        T *= 10; // unit: 1 ns
        C = gTE->GetErrorX(i);
        DC = gTE->GetErrorY(i);

        // skip emty entries
        if (C == 0.0) continue;

        binwx = pHist->GetXaxis()->GetBinWidth(
                        pHist->GetXaxis()->FindBin(T));
        binwy = pHist->GetYaxis()->GetBinWidth(
                        pHist->GetYaxis()->FindBin(E));
        bT = T + 0.1*binwx;
        bE = E + 0.1*binwy;
        binT = pHist->GetXaxis()->FindBin(bT);
        binE = pHist->GetYaxis()->FindBin(bE);

        pHist->SetBinContent(binT, binE, C);
        pHist->SetBinError(binT, binE, DC * C);
    }
    if (draw)
    {
        sprintf(name, "%s_EvsT_%i", FC.c_str(), ch + 1);
        TCanvas *c1 = new TCanvas(name, name, 200, 10, 700, 500);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        pHist->SetStats(0);
        pHist->Draw("colz");
    }
//    cout << pHist->Integral() << endl;
    return pHist;
}

TH2F* GetMCNPbins(string in_file_path, string FC = "H19", Int_t channel = 0)
{ /// Find MCNP's time and energy binning.
    string TabularFile = in_file_path + std::to_string(channel+1) + "4";
    // Open tabular as filestream
    std::ifstream input(TabularFile.c_str());
    Bool_t stop = 0;
    string line, word, HistName, HistTitle;
    Double_t E, T, C[10];
    vector<float> *ToFBins = new vector<float>();
    vector<float> *EnergyBins = new vector<float>();
    // Loop over all time headlines
    do {
        // Find start of time headline
        do {
            input >> word;
        } while (strcmp(word.c_str(), "time:"));
        // Read head line
        Int_t i = 0; // i: time columns
        for (i = 0; i < 5; i++) {
            input >> word;
            if (!strcmp(word.c_str(), "energy")) break;
            if (!strcmp(word.c_str(), "total")) break;
            T = std::stod(word); // unit in file: [t] = 10ns
            //Save T
            ToFBins->push_back(10*T); // unit: [t] = 1ns
        }
        // Check if last table block is reached
        if (!strcmp(word.c_str(), "total"))
            stop = 1;
        // Move on until "energy"
        while (strcmp(word.c_str(), "energy"))
            input >> word;
    } while (!stop);
    input >> word;
    while (strcmp(word.c_str(), "total")) {
        E = std::stod(word);
        EnergyBins->push_back(E);
        getline(input, line);
        input >> word;
    }
    int nBinsX = ToFBins->size()-1;
    int nBinsY = EnergyBins->size()-1;
    cout << FC << " ch." << channel+1 << " MCNP bins:" << endl;
    cout << "\tToF:\t" << nBinsX << " bins \t" << ToFBins->at(0) << " - " << ToFBins->at(nBinsX) << endl;
    cout << "\tE:\t"   << nBinsY << " bins \t" << EnergyBins->at(0) << " - " << EnergyBins->at(nBinsY) << endl;

    HistName = FC + "_TvsE_MCNPbins_Ch." + std::to_string(channel+1);
    HistTitle = "ToF vs Ekin, " + FC + " ch." + std::to_string(channel+1) + ", MCNP bins";
    TH2F *h2 = new TH2F( HistName.c_str(), HistTitle.c_str(),
                         nBinsX, &(ToFBins->at(0)),
                         nBinsY, &(EnergyBins->at(0)));
    return h2;
}

Double_t CalculateNeutronEnergy(Float_t FlightPath, Float_t ToF)
{   //calculates the kinetic energy of the neutron in MeV, [L] = m, [ToF] = ns!
    //E_kin = m_n c^2 (gamma-1)
    Double_t E_kin   = 0;
    Double_t L       = FlightPath;       //n flight path to the i-th deposit [L] = m
    Double_t t       = ToF*1e-9;         //neutron time of flight [t] = s
    Double_t c       = 299792458.0;      //speed of light [c] = m/s
    Double_t m_nc2   = 939.56537921;     //neutron rest mass [MeV/c²]*c²

    Double_t beta = L/(c*t);
    Double_t gamma= 1./sqrt(1.-pow(beta,2));

    E_kin = m_nc2*(gamma - 1.0);

    return E_kin;
}

void HistoMCNP(string in_file_path, string FC, Int_t channel, TFile* f)
{ // Create TvsE and EvsE histograms from MCNP simulation results, save to file f
    char name[64] = "";
    TH2F *hM = GetMCNPbins(in_file_path, FC, channel);
    // FCscat_c_tally_2*4_xyz_0.dmp
    TGraphErrors *gTE = new TGraphErrors((in_file_path+std::to_string(channel+1)+"4_xyz_0.dmp").c_str(), "%lg %lg %lg %lg");
    Int_t N = gTE->GetN();
    cout << "N=" << N << endl;
    Double_t T, E, C, DC;
    Int_t binT, binE;
    for (Int_t i = 0; i < N; i++)
    {
        gTE->GetPoint(i, T, E); // unit [T]: 10 ns
        T *= 10; // unit: 1 ns
        C = gTE->GetErrorX(i);
        DC = gTE->GetErrorY(i);
        if (C == 0.0) continue;
        binT = hM->GetXaxis()->FindBin(T);
        binE = hM->GetYaxis()->FindBin(E);
        if(T > hM->GetXaxis()->GetBinCenter(binT))
            binT++;
        if(E > hM->GetYaxis()->GetBinCenter(binE))
            binE++;
//        cout << binT << " " << binE << " " << C << " " << DC * C << endl;
        hM->SetBinContent(binT, binE, C);
        hM->SetBinError(binT, binE, DC * C);
    }
//    new TCanvas();
//    gPad->SetLogz();
//    hM->Draw("colz");
//    cout << "hM: " << hM->Integral() << endl;

    // Open a 2D histogram as axis prototype
    string PrototypeFile = "~/TrackLength/G4UFCvsH19_1E8_hist.root";
    TFile *p = TFile::Open(PrototypeFile.c_str()); if (!p) cout << "Could not open " << PrototypeFile << endl;
    sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), channel+1);
    TH2F *protoTvsE = (TH2F*)p->Get(name); if (!protoTvsE) cout << "Could not get prototype histogram file " << name << endl;
    sprintf(name, "%s_ToFvsEkin_prototype_Ch.%i", FC.c_str(), channel+1);
    protoTvsE->SetName(name);
    sprintf(name, "%s_ToFvsEkin_Ch.%i", FC.c_str(), channel+1);
    TH2D *hTvsE = (TH2D*)protoTvsE->Clone(name);
    sprintf(name, "%s, ToF vs Ekin, ch.%i", FC.c_str(), channel+1);
    hTvsE->SetTitle(name);
    hTvsE->GetXaxis()->SetTitle("#font[12]{t} / ns");
    hTvsE->GetYaxis()->SetTitle("#font[12]{E} / MeV");
    hTvsE->Scale(0.0);

    sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), channel+1);
    TH2F *protoEvsE = (TH2F*)p->Get(name); if (!protoEvsE) cout << "Could not get prototype histogram file " << name << endl;
    sprintf(name, "%s_EToFvsEkin_prototype_Ch.%i", FC.c_str(), channel+1);
    protoEvsE->SetName(name);
    sprintf(name, "%s_EToFvsEkin_Ch.%i", FC.c_str(), channel+1);
    TH2D *hEvsE = (TH2D*)protoEvsE->Clone(name);
    sprintf(name, "%s, E(ToF) vs Ekin, ch.%i", FC.c_str(), channel+1);
    hEvsE->SetTitle(name);
    hEvsE->GetXaxis()->SetTitle("#font[12]{E(t)} / MeV");
    hEvsE->GetYaxis()->SetTitle("#font[12]{E} / MeV");
    hEvsE->Scale(0.0);

    TRandom2 *rnd = new TRandom2();
    Int_t n = 2000000;
    Double_t step = hM->Integral() / n;
    Double_t H19_FlightPath[10]  = { 5.976, 5.977, 5.986, 5.987, 5.996,
                                       5.997, 6.006, 6.007, 6.016, 6.017};
    Double_t UFC_FlightPath[8]   = { 6.262, 6.251, 6.241, 6.230,
                                       6.220, 6.209, 6.199, 6.188};
    Double_t FlightPath = strcmp(FC.c_str(), "H19") ? UFC_FlightPath[channel] : H19_FlightPath[channel];
    for (binT = 1; binT <= hM->GetNbinsX(); binT++)
        for (binE = 1; binE <= hM->GetNbinsY(); binE++)
        {
            C = hM->GetBinContent(binT, binE);
            while (C > 0)
            {
                T = hM->GetXaxis()->GetBinLowEdge(binT) + rnd->Uniform(0, hM->GetXaxis()->GetBinWidth(binT));
                E = hM->GetYaxis()->GetBinLowEdge(binE) + rnd->Uniform(0, hM->GetYaxis()->GetBinWidth(binE));
                DC = (C > step) ? step : C;
                hTvsE->Fill(T, E, DC);
                hEvsE->Fill(CalculateNeutronEnergy(FlightPath, T), E, DC);
                C -= step;
            }
        }
//    cout << "hTvsE: " << hTvsE->Integral() << endl;
    string Path = FC+"/ToFvsEkin";
    SaveToFile(f, Path, hTvsE);
    Path = FC+"/EToFvsEkin";
    SaveToFile(f, Path, hEvsE);
//    hEvsE->Draw("colz");
}

void ReadMCNP(string InFilePath, string OutFileName)
{
    cout << "Reading MCNP results from .dmp files" << endl
         << "\tInput paths: " << InFilePath << "*" << "4_xyz_0.dmp" << endl
         << "\tHistogram output file: " << OutFileName << endl;
    TFile *f = TFile::Open(OutFileName.c_str(), "UPDATE");
    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 1; i_FC < 2; i_FC++) // No H19 results yet!
    {
        string FC = FCs[i_FC];
        for (Int_t ch = 0; ch < NumCh[i_FC]; ch++)
        {
            HistoMCNP(InFilePath, FC, ch, f);
        }
    }
    f->Save();
    f->Close();
}

void SourceMCNP(string StrFileFilled, string StrFileVoid)
{
// Add fission cross section weighted source spectrum to the Filled MCNP simulation root file.
// The weighted source spectrum is equal to a weighted deposit's Void neutron energy distribution
    string FCs[] = {"UFC", "H19"};
    char name[64] = "";
    // Open Filled and Void MCNP root files
    TFile *f = TFile::Open(StrFileFilled.c_str(), "UPDATE");
    if (!f) cout << "Could not open " << StrFileFilled << endl;
    TFile *fVoid = TFile::Open(StrFileVoid.c_str(), "READ");
    if (!fVoid) cout << "Could not open " << StrFileVoid << endl;

    for (Int_t i_FC = 0; i_FC < 1; i_FC++) // No MCNP H19 data yet!
    { // For both fission chambers:
        string FC = FCs[i_FC];
        // Get a fission chamber deposit's E(t) vs Ekin histogram from Void simulation
        sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.1", FC.c_str(), FC.c_str());
        TH2D *h2EE = (TH2D*)fVoid->Get(name);
        if (!h2EE) cout << "Could not get " << name << endl;
        // Projection on Ekin axis
        TH1D *h1E = (TH1D*)h2EE->ProjectionY();
        sprintf(name, "%s_Ekin", FC.c_str());
        h1E->SetName(name);
        SaveToFile(f, "Source/nEnergy", h1E);
    }
    f->Save();
    f->Close();
}

void MCNP1()
{
    string InFileFilled = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_filled_2/tally/FCscat_b_tally-2";
    string OutFileFilled = "/gpfs/home/hoffma93/TrackLength/MCNP_0213_filled.root";
    string InFileVoid = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_void/tally/FCscat_c_tally-2";
    string OutFileVoid = "/gpfs/home/hoffma93/TrackLength/MCNP_0213_void.root";
    ReadMCNP(InFileFilled, OutFileFilled);
    ReadMCNP(InFileVoid, OutFileVoid);
    SourceMCNP(OutFileFilled, OutFileVoid);
//    DoCorrection(OutFileFilled);
}

void MCNP2()
{
    string InFileFilled = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_filled_2/tally/FCscat_c_tally-2";
    string OutFileFilled = "/gpfs/home/hoffma93/TrackLength/MCNP_0313_filled.root";
    string InFileVoid = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_void/tally/FCscat_d_tally-2";
    string OutFileVoid = "/gpfs/home/hoffma93/TrackLength/MCNP_0313_void.root";
    ReadMCNP(InFileFilled, OutFileFilled);
    ReadMCNP(InFileVoid, OutFileVoid);
    SourceMCNP(OutFileFilled, OutFileVoid);
//    DoCorrection(OutFileFilled);
}

void StatUncertainty(TH2F *pH2, TGraphErrors *pGrXS, Double_t A, Int_t N)
{ // pH2: Ekin on y-axis
  // XS: Cross section used by G4TTree
  // A: Deposit area
  // N: Simulated projectiles
/// Uncertainty for n' = n * CrossSection / A' / N, (n counting statistics)
/// Delta n' = sqrt(n' * CrossSection / A' / N)
/// A' >= A, but most direct with A' = A
    for (uint biny = 0; biny < pH2->GetNbinsY() + 2; biny++)
    {
        Double_t CrossSection = pGrXS->Eval(pH2->GetYaxis()->GetBinCenter(biny));
        for (uint binx = 0; binx < pH2->GetNbinsX() + 2; binx++)
        {
            Double_t BinContent = pH2->GetBinContent(binx, biny);
            Double_t BinError = sqrt(BinContent * CrossSection / A / N);
            pH2->SetBinError(binx, biny, BinError);
        }
    }
}

void DoStatUncertainty(string FileName, Int_t N)
{
    char name[64] = "";
    // Open cross section file
    string StrEval = "JEFF_3.3";
    TFile *fXS = TFile::Open("~/TrackLength/FissionXS.root", "READ");
    if (!fXS) cout << "Could not get " << "Fission cross section data" << endl;
    // Open histograms file
    TFile *f = TFile::Open(FileName.c_str(), "UPDATE");
    if (!f) cout << "Could not open " << FileName << endl;
    // Loop over FCs
    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};//{1,1};
    for (Int_t iFC = 0; iFC < 2; iFC++)
    {
        string FC = FCs[iFC];
        // UFC: r = 3,7cm, H19: r = 3,8cm
        Double_t r = strcmp(FC.c_str(), "H19") ? 3.7 : 3.8;
        Double_t A = TMath::Pi() * r * r; // [A] = cm^2
        // Get Cross section graph
        TGraphErrors *pGrXS;
        sprintf(name, "%s/%sTarget", StrEval.c_str(), FC.c_str());
        fXS->GetObject(name, pGrXS);
        if (!pGrXS) cout << "Could not get " << name << endl;
        for (Int_t ch = 0; ch < NumCh[iFC]; ch++)
        {
            sprintf(name, "%s/EToFvsEkin/%s_EToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
            cout << name << endl;
            TH2F* pH2EE = (TH2F*)f->Get(name); if (!pH2EE) cout << "Could not get " << name << endl;
            StatUncertainty(pH2EE, pGrXS, A, N);

//            cout << "XS / b: " << pGrXS->Eval(pH2EE->GetYaxis()->GetBinCenter(4900)) << endl;
//            cout << "A / cm^2: " << A << endl;
//            cout << "N: " << N << endl;
//            cout << pH2EE->GetBinContent(4900, 4900) << " +- " << pH2EE->GetBinError(4900, 4900) << endl;

            sprintf(name, "%s/EToFvsEkin", FC.c_str());
            f->GetDirectory(name)->cd();
            pH2EE->Write(pH2EE->GetName(), TObject::kOverwrite);
        }
    }
    f->Save();
    f->Close();
}

void PrintEE()
{
    TFile *f = TFile::Open("/home/hoffma93/TrackLength/Backup/G4UFCvsH19_1E9_hist_w.root", "READ");
    string FC = "H19";
    string Name = FC + "/EToFvsEkin/" + FC + "_EToFvsEkin_Ch." + "1";
    TH2F *h = (TH2F*)f->Get(Name.c_str());
    cout << h->GetBinContent(4900, 4900) << " +- " << h->GetBinError(4900, 4900) << endl;
    f->Close();
}

// CompSim("~/Programme/Geant4-Work/results/PuFC_real_FG_1E7.root", "~/Programme/Geant4-Work/results/PuFC_real_FG_TL_1E7.root", "PuFC", 0)
void CompSim(string StrSimFile1, string StrSimFile2, string FC = "PuFC", Int_t ch = 0)
{// Compare the ToF spectrum of a deposit for two simulation runs; norm to sum.
    /// open root files
    TFile *f1 = TFile::Open(StrSimFile1.c_str(), "READ"); if (!f1) cout << "Could not open " << StrSimFile1 << endl;
    TFile *f2 = TFile::Open(StrSimFile2.c_str(), "READ"); if (!f2) cout << "Could not open " << StrSimFile2 << endl;

    /// Get histograms
    char name[128] = "";
    sprintf(name, "%s/ToF/%s_ToF_Ch.%i", FC.c_str(), FC.c_str(), ch+1);
    TH1F *h1 = (TH1F*) f1->Get(name); if (!h1) cout << "Could not get " << name << endl;
    TH1F *h2 = (TH1F*) f2->Get(name); if (!h2) cout << "Could not get " << name << endl;

    /// Normalization
    Double_t scale = h2->Integral() / h1->Integral();
    cout << "scale " << h2->Integral() << "\t/ " << h1->Integral() << "\t= " << scale << endl;
    h1->Scale(scale);

    TCanvas *c = new TCanvas();
    h1->SetLineColor(kBlue);
    h1->SetStats(0);
    h1->Draw("hist");
    h2->SetLineColor(kRed);
    h2->Draw("same hist");

    TLegend *l = new TLegend(0.3, 0.3, 0.7, 0.4, "n weight");
    l->AddEntry(h1, "1/cos(#Theta)");
    l->AddEntry(h2, "#it{l / V}");
    l->Draw();
}

TH2D* GetEvsT(string file_to_read, string FC = "PuFC", Int_t ch = 0)
{ // Get MCNP 2D histogram
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
    Tmin = 0.0, Tmax = 1500.0, Emin = 0.0, Emax = 16.0;
    Tbins = 15000, Ebins = 1600;
    sprintf(name, "%s_ToFvsEkin_Ch.%i", FC.c_str(), ch + 1);
    TH2D *pHist = new TH2D(name, name,
                           Tbins, Tmin, Tmax, Ebins, Emin, Emax);
    sprintf(name, "%s, ToF vs Ekin, ch.%i", FC.c_str(), ch + 1);
    pHist->SetTitle(name);
    pHist->GetXaxis()->SetTitle("#font[12]{t} / ns");
    pHist->GetYaxis()->SetTitle("#font[12]{E} / MeV");

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
//    cout << pHist->Integral() << endl;
    return pHist;
}

// /net/cns/projects/NTOF/Hemera/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_void/tally/FCscat_b_tally-2*4_xyz_0.dmp
// /net/cns/projects/NTOF/Hemera/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_filled_2/tally/FCscat_b_tally-2*4_xyz_0.dmp
void ConvertMCNP()
{
    string FC = "PuFC";
    string DirName = "/net/cns/projects/NTOF/Hemera/MCNP/FissionChamberScattering";
//    string Path = "FCscat_H19_realspec_tracklength_void/tally/FCscat_b";
//    string SimRootFile = "~/TrackLength/MCNP_0724_void.root";
    string Path = "FCscat_H19_realspec_tracklength_filled_2/tally/FCscat_b";
    string SimRootFile = "~/TrackLength/MCNP_0724_filled.root";

    char name[128] = "";
    TH2D *h[8];

    TFile *fSim = TFile::Open(SimRootFile.c_str(), "UPDATE");

    // Loop over files
    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "%s/%s_tally-2%i4_xyz_0.dmp", DirName.c_str(), Path.c_str(), 8 - i);
        cout << name << endl;

        // Create 2D histogram
        h[i] = GetEvsT(name, FC, i);

        Save(fSim, FC+"/ToFvsEkin", h[i]);
//        TH1D *hProj = TimeProjection(h[i], i, FC, key);
//        Save(fSim, "Simulation/MCNP/"+FC+"_"+key+"/EffToF", hProj);

    }
    fSim->Save();
    fSim->Close();
}

void CorrectionVacuum()
{
    string FC = "PuFC";
    string SimTool = "MCNP";
    string SimToolShortcut = SimTool;
    if (!strcmp(SimTool.c_str(), "Geant4")) SimToolShortcut = "G4";

    string VoidSimHistFile = "~/TrackLength/MCNP_0724_void.root";
    string FilledSimHistFile = "~/TrackLength/MCNP_0724_filled.root";
//    string FCAnalysisFile = "~/Programme/FC-Analysis/2016.06/nfis_2020.root";
    string FCAnalysisFile = "~/TrackLength/nfis_2020.root";

    cout << "**** Correction factors using cross-section weighted E-t-histograms" << endl;
    cout << "\t Fission chamber: " << FC << endl;
    cout << "\t Simulation: " << SimTool << endl;
//    cout << "\t Shortcut: " << SimToolShortcut << endl;
    cout << "\t Void simulation histograms: " << VoidSimHistFile << endl;
    cout << "\t Filled simulation histograms: " << FilledSimHistFile << endl;


    char name[64];

    TH2D *pH2Void;
    TH2D *pH2Filled;
    TH1D *pH1Void;
    TH1D *pH1Filled;
    TH1D *pH1Correction[8];
    TGraphErrors *pGrCorrection[8];

    TFile *fVoid = TFile::Open(VoidSimHistFile.c_str(), "READ"); if (!fVoid) cout << "Could not open " << VoidSimHistFile << endl;
    TFile *fFilled = TFile::Open(FilledSimHistFile.c_str(), "READ"); if (!fFilled) cout << "Could not open " << FilledSimHistFile << endl;
    TFile *fAna = TFile::Open(FCAnalysisFile.c_str(), "UPDATE"); if (!fAna) cout << "Could not open " << FCAnalysisFile << endl;

    for (Int_t i = 0; i < 1; i++)
    {
        // Get E vs t histograms
        sprintf(name, "%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        pH2Void = (TH2D*) fVoid->Get(name); if (!pH2Void) cout << "Could not get " << name << " from " << fVoid->GetName() << endl;
        pH2Filled = (TH2D*) fFilled->Get(name); if (!pH2Filled) cout << "Could not get " << name << " from " << fFilled->GetName() << endl;

        // Make projections on the time axis
        sprintf(name, "%s_Void_ToF_Ch.%i", FC.c_str(), i+1);
        pH1Void = (TH1D*) pH2Void->ProjectionX(name, 0, -1, "e");
        sprintf(name, "%s_Filled_ToF_Ch.%i", FC.c_str(), i+1);
        pH1Filled = (TH1D*) pH2Filled->ProjectionX(name, 0, -1, "e");

//        pH1Void->Rebin(10);
//        pH1Filled->Rebin(10);

        // Divide projections to get correction factor
//        pH1Correction[i] = (TH1D*) pH1Void->Clone("Correction");
//        pH1Correction[i]->Divide(pH1Void, pH1Filled, 1., 1.);

        // Name and format
//        sprintf(name, "%sCorrection_%s_Ch.%i", SimToolShortcut.c_str(), FC.c_str(), i+1);
//        pH1Correction[i]->SetName(name);
        pGrCorrection[i] = new TGraphErrors();
        sprintf(name, "GrCorrection_%s_Ch.%i", FC.c_str(), i+1);
        pGrCorrection[i]->SetName(name);

        Int_t nPoint = 0;
        for (Int_t bin = 1; bin <= pH1Void->GetNbinsX(); bin++)
        {
            if (pH1Void->GetBinContent(bin) != 0 && pH1Filled->GetBinContent(bin) != 0)
                pGrCorrection[i]->SetPoint(nPoint++, pH1Void->GetXaxis()->GetBinCenter(bin), pH1Void->GetBinContent(bin) / pH1Filled->GetBinContent(bin));
        }


        // save
        sprintf(name, "Analysis/Simulation/%s/Correction/", SimToolShortcut.c_str());
//        SaveToFile(fAna, name, pH1Correction[i]);
        SaveToFile(fAna, name, pGrCorrection[i]);
        fAna->Save();

        delete pH2Void, pH2Filled, pH1Void, pH1Filled;
    }

//    fAna->Save();
    fAna->Close();
}

void nELBEsim(string FileName = "G4UFCvsH19_1E7_hist.root")
{
//    string FileName = "G4UFCvsH19_23082016_hist.root"; // 2016 original
//    string FileName = "G4UFCvsH19_v10.2_hist.root";    // G4 10.2 no TL
//    string FileName = "G4UFCvsH19_v10.2_hist_w.root";  // G4 10.2 with TL
//    string FileName = "G4UFCvsH19_1E9_hist.root";      // G4 10.5 no TL
//    string FileName = "G4UFCvsH19_1E9_hist_w.root";    // G4 10.5 with TL
    string File1 = "/home/hoffma93/TrackLength/" + FileName;
    TFile *f = TFile::Open(File1.c_str(), "UPDATE");
    Correction(f);
    CorrectionRefH19(f);
    TransmCorrLoss(f);
    f->Save();
    f->Close();
}
