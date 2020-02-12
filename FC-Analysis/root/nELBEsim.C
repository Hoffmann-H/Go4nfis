#include "SaveToFile.C"
#include "/gpfs/home/hoffma93/StyleSheets/StyleSheet.C"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TLatex.h"

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
{
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
void CompareWeighting(string StrHistoFileName)
{ // file name without .root ending
    TFile *fH = TFile::Open((StrHistoFileName+".root").c_str(), "READ"); // histogrammed simulation data, without track length weighting
//    TFile *fW = TFile::Open((StrHistoFileName+"_w.root").c_str(), "READ");
    TFile *fout = TFile::Open((StrHistoFileName+"_result.root").c_str(), "UPDATE");

    string FCs[] = {"H19", "UFC"};
    Int_t NumCh[] = {10, 8};
    for (Int_t i_FC = 0; i_FC < 2; i_FC++)
    {
        string FC = FCs[i_FC];
        TGraphErrors *pGrXS = GetXS("~/TrackLength/FissionXS.root", FC, "JEFF_3.3");
        for (Int_t channel = 0; channel < 2/*NumCh[i_FC]*/; channel++)
        {
//            TH1F *pk_H = GetCorrLoss((StrHistoFileName+".root").c_str(), FC, channel);
//            pk_H->SetName(("H_k_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "Hit/CorrLoss", pk_H);
//            TH1F *pk_W = GetCorrLoss((StrHistoFileName+"_w.root").c_str(), FC, channel);
//            pk_W->SetName(("W_k_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "TrackLength/CorrLoss", pk_W);

//            TH1F *pT_H = GetTransm((StrHistoFileName+".root").c_str(), FC, channel);
//            pT_H->SetName(("H_T_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "Hit/Transmission", pT_H);
//            TH1F *pT_W = GetTransm((StrHistoFileName+"_w.root").c_str(), FC, channel);
//            pT_W->SetName(("W_T_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "TrackLength/Transmisssion", pT_W);

            TH1F *pC_H = GetCorrection(fH, FC, channel, pGrXS);
            pC_H->SetDirectory(0);
            pC_H->SetName(("H_C_"+FC+"_"+std::to_string(channel+1)).c_str());
            SaveToFile(fout, "Hit/Correction", pC_H);
//            TH1F *pC_W = GetCorrection((StrHistoFileName+"_w.root").c_str(), FC, channel);
//            pC_W->SetName(("W_C_"+FC+"_"+std::to_string(channel+1)).c_str());
//            SaveToFile(f, "TrackLength/Correction", pC_W);
        }
    }
    fH->Close();
    fout->Save();
    fout->Close();

//    new TCanvas((FC+"_"+std::to_string(channel+1)).c_str());
//    pCorrH->Draw();
//    pCorrW->SetLineColor(kRed);
//    pCorrW->Draw("same");
//    TLegend *l = new TLegend(0.2, 0.5, 0.6, 0.8);
//    l->AddEntry(pCorrH, "Hit");
//    l->AddEntry(pCorrW, "Track length");
//    l->Draw();
}

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

void CorrectionRefH19(string StrHistoFileName, Bool_t w = 0)
{
    TFile *f = TFile::Open((StrHistoFileName+"_result.root").c_str(), "UPDATE");
    string Name = "";
    TH1F *pH1_H19[10], *pH1_sumH19, *pH1_UFC[8], *pH1_Ref[8];
    for (Int_t channel = 0; channel < 10; channel++)
    {
        if (w) Name = "TrackLength/Correction/W_C_H19_"+std::to_string(channel+1);
        else   Name = "Hit/Correction/H_C_H19_"+std::to_string(channel+1);
        pH1_H19[channel] = (TH1F*)f->Get(Name.c_str());
        if (!pH1_H19[channel]) cout << "Could not get " << Name << endl;
        if (channel==0) {
            pH1_sumH19 = (TH1F*)pH1_H19[0]->Clone(w?"W_C_H19":"H_C_H19");
            if (!pH1_sumH19) cout << "Could not get " << Name << endl;
            pH1_sumH19->SetDirectory(0);
        } else
            pH1_sumH19->Add(pH1_H19[channel], 1.0);
    }
    pH1_sumH19->Scale(0.1);
    if (w) Name = "TrackLength/Correction";
    else   Name = "Hit/Correction";
    SaveToFile(f, Name.c_str(), pH1_sumH19);
    for (Int_t channel = 0; channel < 8; channel++)
    {
        if (w) Name = "TrackLength/Correction/W_C_UFC_"+std::to_string(channel+1);
        else   Name = "Hit/Correction/H_C_UFC_"+std::to_string(channel+1);
        pH1_UFC[channel] = (TH1F*)f->Get(Name.c_str());
        if (!pH1_UFC[channel]) cout << "Could not get " << Name << endl;
        if (w) Name = "W_C_UFC_RefH19_"+std::to_string(channel+1);
        else   Name = "H_C_UFC_RefH19_"+std::to_string(channel+1);
        pH1_Ref[channel] = (TH1F*)pH1_UFC[channel]->Clone(Name.c_str());
        pH1_Ref[channel]->SetDirectory(0);
        pH1_Ref[channel]->Divide(pH1_UFC[channel], pH1_sumH19, 1.0, 1.0);
        if (!pH1_Ref[channel]) cout << "Could not get " << Name << endl;
        if (w) Name = "TrackLength/Correction/RefH19";
        else   Name = "Hit/Correction/RefH19";
        SaveToFile(f, Name.c_str(), pH1_Ref[channel]);
    }
    f->Save();
    f->Close();
}

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

void nELBEsim()
{
//    GetCorrLoss("~/TrackLength/G4UFCvsH19_1E8_hist.root", "UFC", 0, 1);

    string File1 = "/home/hoffma93/TrackLength/G4UFCvsH19_1E9_hist";
    string File2 = "/home/hoffma93/TrackLength/G4UFCvsH19_23082016_hist_repr";
    string File3 = "/home/hoffma93/TrackLength/G4UFCvsH19_1E6_fast_Phi";
    string File4 = "/home/hoffma93/TrackLength/G4UFCvsH19_1E6_fast_Theta";
    CompareWeighting(File1);
//    avWeight(StrFileInPath);
//    scriptW(File1+"_result.root");
//    scriptH(File1+"_result.root");
//    CompareWeighting(File2);
//    CorrectionRefH19(File2);
//    scriptH(File2+"_result.root");
}
