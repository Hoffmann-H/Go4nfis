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

//TH1F* WeightTH1withXS(TH1F *pH1In, string StrXSFileInPath, string StrFC, string StrEval)
//{   //function to weight a Ekin histogram with the fission cross
//    //section to get a unnormalized fission rate

//    Int_t NumBinX = pH1In->GetNbinsX();       //Ekin bins
//    Double_t Ekin, XS, Content, EContent;
//    TGraphErrors *pGrXS;

//    TH1F *pH1 = (TH1F*)pH1In->Clone(); pH1->Clear();

//    //open cross section file and get it
//    TFile *f = new TFile(StrXSFileInPath.c_str(), "READ");
//        f->GetObject((StrEval+"/"+StrFC+"Target").c_str(), pGrXS); if (!pGrXS) cout << "Could not get " << StrEval << "/" << StrFC << "Target" << endl;

//    for (int x_i=1; x_i<=NumBinX; x_i++)
//    {   //get neutron kinetic energy of corresponding bin x_i
//        Ekin    = pH1In->GetXaxis()->GetBinCenter(x_i);
//        //get cross section @ E=Ekin
//        XS      = pGrXS->Eval(Ekin);
//        //multiply bin content with cross section to get fission rate
//        Content = pH1In->GetBinContent(x_i) * XS;
//        EContent= pH1In->GetBinError(x_i) * XS;
//        pH1->SetBinContent(x_i, Content);
//        pH1->SetBinError(x_i, EContent);
//    }
//    f->Close();
//    return pH1;
//}

TGraphErrors* GetXS(string StrXSFileInPath, string StrFC, string StrEval)
{
    TFile *f = new TFile(StrXSFileInPath.c_str(), "READ");
    TGraphErrors *pGrXS;
    f->GetObject((StrEval+"/"+StrFC+"Target").c_str(), pGrXS); if (!pGrXS) cout << "Could not get " << StrEval << "/" << StrFC << "Target" << endl;
    f->Close();
    return pGrXS;
}

//TH1F* GetCorrection(TFile *f, string StrFC, uint ch, TGraphErrors *pGrXS)
//{   //method to determine a correction factor C from geometric and vacuum simulation
//    //C(E(t)) = vacuum fission rate at E(t) / total detected fission rate at E(t)
//    string ProjName, CorrName;

//    //get Energy-ToF-Correlation histogram of all G4 simulated neutrons
//    TH2F* pH2Tot = GetGeant4Data(f, StrFC, false, ch);
//    //get source spectrum
//    TH1F *pFCSource = (TH1F*)GetSourceSpectrum(f);

//    //convert total and source histogram to fission rate histograms
//    TH2F *pH2TotW = WeightWithXS(pH2Tot, pGrXS);
//    TH1F *pFCSourceW = WeightTH1withXS(pFCSource, pGrXS);
//    // ...or even dont
////    TH2F *pH2TotW = pH2Tot->Clone();
////    TH1F *pFCSourceW = pFCSource->Clone();
//    //get projection to the E(t) axis
//    ProjName = "G4Total_"+StrFC+"_Ch."+std::to_string(ch+1)+"_Hist";
//    TH1F* pFCTot = (TH1F*)pH2TotW->ProjectionX(ProjName.c_str(), 0, -1, "e");
//    //calculate correction factor C(E(t))
//    CorrName = "C_"+StrFC+"_"+std::to_string(ch+1);
//    TH1F* pFCCorr = (TH1F*)pFCTot->Clone(CorrName.c_str()); //pFCCorr->Clear();
//          pFCCorr->Divide(pFCSourceW, pFCTot, 1.0, 1.0);
//    pFCCorr->SetDirectory(0);
//    return pFCCorr;
//}

TH1F* GetCorrection(TFile *f, string StrFC, uint ch, TGraphErrors *pGrXS)
{   //method to determine a correction factor C from geometric and vacuum simulation
    //C(E(t)) = vacuum fission rate at E(t) / total detected fission rate at E(t)
    string SumName, ProjName, CorrName;

    //get Energy-ToF-Correlation histogram of all G4 simulated neutrons
    TH2F* pH2Tot = GetGeant4Data(f, StrFC, false, ch);
    //get Energy-ToF-Correlation histogram of scattered neutrons
    TH2F* pH2Sc = GetGeant4Data(f, StrFC, true, ch);
    //get Energy-ToF-Correlation histogram of direct neutrons
    SumName = StrFC + "_EToFvsEkin_Dir_Ch." + std::to_string(ch+1);
    TH2F* pH2Dir = (TH2F*)pH2Tot->Clone(SumName.c_str());
    pH2Dir->Add(pH2Tot, pH2Sc, 1.0, -1.0);
    //get Energy distribution of direct neutrons
    ProjName = StrFC + "_Ekin_Dir_Ch." + std::to_string(ch+1);
    TH1F* pH1Dir = (TH1F*)pH2Dir->ProjectionY(ProjName.c_str(), 0, -1, "e");
        // x and y projections should be equal,
        // but they aren't, probably due to binning.
    //get source spectrum
    TH1F *pFCSource = (TH1F*)GetSourceSpectrum(f);

    //convert direct, scattered and source histogram to fission rate histograms
    TH1F *pFCDir = WeightTH1withXS(pH1Dir, pGrXS);
    TH2F *pH2ScW = WeightWithXS(pH2Sc, pGrXS);
    TH1F *pFCSourceW = WeightTH1withXS(pFCSource, pGrXS);
    // ...or even dont
//    TH1F *pFCDir = pH1Dir->Clone();
//    TH2F *pH2ScW = (TH2F*)pH2Sc->Clone();
//    TH1F *pFCSourceW = (TH1F*)pFCSource->Clone();

    //get projection to the E(t) axis
    ProjName = "G4Sc_"+StrFC+"_Ch."+std::to_string(ch+1)+"_Hist";
    TH1F* pFCSc = (TH1F*)pH2ScW->ProjectionX(ProjName.c_str(), 0, -1, "e");
    //add direct and scattered fission rate
    SumName = "G4Total_"+StrFC+"_Ch."+std::to_string(ch+1)+"_Hist";
    TH1F* pFCTot = (TH1F*)pFCDir->Clone(SumName.c_str());
    pFCTot->Add(pFCDir, pFCSc, 1.0, 1.0);
    CorrName = "C_"+StrFC+"_"+std::to_string(ch+1);
    TH1F* pFCCorr = (TH1F*)pFCTot->Clone(CorrName.c_str()); //pFCCorr->Clear();
    pFCCorr->Divide(pFCSourceW, pFCTot, 1.0, 1.0);
    pFCCorr->SetDirectory(0);
    //if correction factor >> 1, probably normation 1/A_dep missing
    if (pFCCorr->Integral() / pFCCorr->GetNbinsX() > 5.0)
    {
        // UFC: r = 3,7cm, H19: r = 3,8cm
        Double_t r = strcmp(StrFC.c_str(), "H19") ? 3.7 : 3.8;
        Double_t A = TMath::Pi() * r * r; // [A] = cm^2
        cout << "Scaling source spectrum with 1/" << A << endl;
        pFCCorr->Scale(1.0/A);
    }
    return pFCCorr;
}

TH1F* GetTransm(TFile *f, string StrFC, uint ch)
{
    //T = direct neutron rate / source neutron rate

    //get Energy-ToF-Correlation histogram of all G4 simulated neutrons
    TH2F* pH2Tot    = GetGeant4Data(f, StrFC, false, ch);
    //get Energy-ToF-Correlation histogram of all scattered neutrons
    TH2F* pH2Scat   = GetGeant4Data(f, StrFC, true,  ch);
    //get Energy-ToF-Correlation histogram of all un-scattered neutrons
    TH2F* pH2UnScat = (TH2F*)pH2Tot->Clone("UnScat");
        pH2UnScat->Add(pH2Tot, pH2Scat, 1., -1.);

    //get projection to the E(t) axis
    TH1F *pFCUnScat, *pFCSource, *pFCTransm;
    string ProjName;
    ProjName = "G4UnScattered_"+StrFC+"_Ch."+std::to_string(ch+1)+"_Hist";
    pFCUnScat = (TH1F*)pH2UnScat->ProjectionY(ProjName.c_str(), 0, -1, "e");
    pFCSource = (TH1F*)GetSourceSpectrum(f);
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
{   //method to determine a correction factor k concerning the loss of kinetic
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


void CorrectionRefH19(TFile *f)
{
    string Name = "";
    TH1F *pH1_H19[10], *pH1_sumH19, *pH1_UFC[8], *pH1_Ref[8];
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
    SaveToFile(f, Name.c_str(), pH1_sumH19);
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
    }
}

void DoCorrection(string StrHistoFileName)
{ // file name without .root ending
    TFile *f = TFile::Open((StrHistoFileName+".root").c_str(), "UPDATE"); // histogrammed simulation data
    if (!f) cout << "Could not open " << StrHistoFileName << ".root" << endl;

    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        TGraphErrors *pGrXS = GetXS("~/TrackLength/FissionXS.root", FC, "JEFF_3.3");
        for (Int_t channel = 0; channel < NumCh[i_FC]; channel++)
        {
//            TH1F *pk_H = GetCorrLoss((StrHistoFileName+".root").c_str(), FC, channel);
//            pk_H->SetName(("H_k_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "Hit/CorrLoss", pk_H);
//            TH1F *pT_H = GetTransm((StrHistoFileName+".root").c_str(), FC, channel);
//            pT_H->SetName(("H_T_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "Hit/Transmission", pT_H);
            TH1F *pC = GetCorrection(f, FC, channel, pGrXS);
            pC->SetDirectory(0);
            pC->SetName(("C_"+FC+"_"+std::to_string(channel+1)).c_str());
            SaveToFile(f, FC+"/Correction", pC);
        }
    }
    f->Save();
    CorrectionRefH19(f);
    f->Save();
    f->Close();
}
/*
void DoCorrectionW(string StrHistoFileName)
{ // file name without _w.root ending
    TFile *fW = TFile::Open((StrHistoFileName+"_w.root").c_str(), "READ"); // histogrammed simulation data, with track length weighting
    TFile *fout = TFile::Open((StrHistoFileName+"_result.root").c_str(), "UPDATE");

    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        TGraphErrors *pGrXS = GetXS("~/TrackLength/FissionXS.root", FC, "JEFF_3.3");
        for (Int_t channel = 0; channel < NumCh[i_FC]; channel++)
        {
//            TH1F *pk_W = GetCorrLoss((StrHistoFileName+"_w.root").c_str(), FC, channel);
//            pk_W->SetName(("W_k_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(fout, "TrackLength/CorrLoss", pk_W);
//            TH1F *pT_W = GetTransm((StrHistoFileName+"_w.root").c_str(), FC, channel);
//            pT_W->SetName(("W_T_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(fout, "TrackLength/Transmisssion", pT_W);
            TH1F *pC_W = GetCorrection(fW, FC, channel, pGrXS);
            pC_W->SetName(("W_C_"+FC+"_"+std::to_string(channel+1)).c_str());
            SaveToFile(fout, "TrackLength/Correction", pC_W);
        }
    }
    fW->Close();
    CorrectionRefH19(fout, 1);
    fout->Save();
    fout->Close();
}//*/

/*void avWeight(string StrHistoFileName)
{
    TFile *f = TFile::Open((StrHistoFileName+"_result.root").c_str(), "UPDATE");

    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        for (Int_t channel = 0; channel < NumCh[i_FC]; channel++)
        {
            for (Int_t sc = 0; sc < 2; sc++)
            {
                // get E-E(t) histograms: pH2_w weighted by 1/cos(theta), pH2 unweighted
                TH2F *pH2_w = GetGeant4Data(StrHistoFileName+"_w.root", FC, sc, channel);
                TH2F *pH2 = GetGeant4Data(StrHistoFileName+".root", FC, sc, channel);

                // make projections
                TH1F *pH1_w, *pH1, *pHavW;
                string Name;
                Name = "pE(t)_w_"+FC+std::to_string(channel+1);
                pH1_w = (TH1F*)pH2_w->ProjectionX(Name.c_str(), 0, -1, "e");
                Name = "pE(t)_"+FC+std::to_string(channel+1);
                pH1 = (TH1F*)pH2->ProjectionX(Name.c_str(), 0, -1, "e");

                // divide projections
                if (sc) Name = FC+"_WvsEToF_sc_Ch."+std::to_string(channel+1);
                else    Name = FC+"_WvsEToF_Ch."+std::to_string(channel+1);
                pHavW = (TH1F*)pH1_w->Clone(Name.c_str());
                pHavW->SetDirectory(0);
                pHavW->Divide(pH1_w, pH1, 1.0, 1.0);
                pHavW->GetYaxis()->SetTitle("<#font[12]{W}>");

                if (sc) Name = "TrackLength/Weight/Scattered";
                else    Name = "TrackLength/Weight";
                SaveToFile(f, Name, pHavW);
                pHavW->Delete();
            }
        }
    }
    f->Save();
    f->Close();
}*/

void scriptW(string StrFileInPath = "~/TrackLength/G4UFCvsH19_1E8_hist_result.root")
{ /// Copy the 9 FC channels' correction factor histograms
  /// from StrFileInPath into the 2016/06 analysis root file

    TFile *f = TFile::Open(StrFileInPath.c_str()); if (!f) cout << "Could not open " << StrFileInPath << endl;
    TH1F *h[9];
    for (Int_t i = 0; i < 8; i++)
    {
        h[i] = (TH1F*)f->Get(("TrackLength/Correction/W_C_UFC_"+std::to_string(i+1)).c_str());
        h[i]->SetDirectory(0);
        h[i]->SetName(("G4Correction_UFC_Ch."+std::to_string(i+1)).c_str());
    }
    h[8] = (TH1F*)f->Get("TrackLength/Correction/W_C_H19");
    h[8]->SetDirectory(0);
    h[8]->SetName("G4Correction_H19");

    string StrFileOutPath = "~/FC-Analysis/2016.06/results/nfis_2016.06.root";
    TFile *g = TFile::Open(StrFileOutPath.c_str(), "UPDATE"); if (!f) cout << "Could not open " << StrFileOutPath << endl;
    for (Int_t i = 0; i < 9; i++)
    {
        SaveToFile(g, "Analysis/Simulation/Geant4/Correction", h[i]);
    }
    g->Save();
    g->Close();
}

void scriptH(string StrFileInPath = "~/TrackLength/G4UFCvsH19_1E9_hist_result.root")
{ /// Copy the 9 FC channels' correction factor histograms
  /// from StrFileInPath into the 2016/06 analysis root file

    TFile *f = TFile::Open(StrFileInPath.c_str()); if (!f) cout << "Could not open " << StrFileInPath << endl;
    TH1F *h[9];
    for (Int_t i = 0; i < 8; i++)
    {
        h[i] = (TH1F*)f->Get(("Hit/Correction/H_C_UFC_"+std::to_string(i+1)).c_str());
        h[i]->SetDirectory(0);
        h[i]->SetName(("G4Correction_UFC_Ch."+std::to_string(i+1)).c_str());
    }
    h[8] = (TH1F*)f->Get("Hit/Correction/H_C_H19");
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

void nELBEsim()
{
    string File1 = "/home/hoffma93/TrackLength/G4UFCvsH19_1E7_hist";
    string File2 = "/home/hoffma93/TrackLength/G4UFCvsH19_23082016_hist_repr";
    string File3 = "/home/hoffma93/TrackLength/G4UFCvsH19_v10.2_hist";
    string File4 = "/home/hoffma93/TrackLength/G4UFCvsH19_";
    DoCorrection(File1);
//    avWeight(StrFileInPath);
//    scriptW(File1+"_result.root");
//    scriptH(File1+"_result.root");
//    DoCorrectionH(File2);
//    scriptH(File2+"_result.root");

//    string InFilePath = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_filled_2/tally/FCscat_b_tally-2";
//    string OutFileName = "/gpfs/home/hoffma93/TrackLength/MCNP_0213_filled.root";
//    string InFilePath = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_void/tally/FCscat_c_tally-2";
//    string OutFileName = "/gpfs/home/hoffma93/TrackLength/MCNP_0213_void.root";
//    string InFilePath = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_filled_2/tally/FCscat_c_tally-2";
//    string OutFileName = "/gpfs/home/hoffma93/TrackLength/MCNP_0313_filled.root";
//    string InFilePath = "/net/cns/projects/NTOF/Hypnos/MCNP/FissionChamberScattering/FCscat_H19_realspec_tracklength_void/tally/FCscat_d_tally-2";
//    string OutFileName = "/gpfs/home/hoffma93/TrackLength/MCNP_0313_void.root";
//    ReadMCNP(InFilePath, OutFileName);
}
