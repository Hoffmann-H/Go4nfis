#ifndef TOF_H
#define TOF_H
#include "FC.C"
#include "Runs.C"
#include "SaveToFile.C"

TH1F* GetToF(TFile *f, Int_t ch, string FC = "PuFC")
{
    char name[64] = "";
    sprintf(name, "/Histograms/Analysis/FC/TimeDiff/PH-Gated/H1AnaHZDRDtG_%i", ch+1);
    TH1I* pH1I= (TH1I*) f->Get(name);
    Int_t nbins = pH1I->GetNbinsX();
    Int_t xmin = pH1I->GetBinLowEdge(1);
    Int_t xmax = pH1I->GetBinLowEdge(nbins + 1);
    TH1F* pH1F = new TH1F(pH1I->GetName(), "ToF; t / ns; counts", nbins, xmin, xmax);
    pH1F->Add(pH1I, 1);
    pH1F->Sumw2();
    return pH1F;
}

Double_t ToF_BG_min[]   = {25., 850.};               //ns
Double_t ToF_BG_max[]   = {95.,2450.};               //ns
Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{   //compute the function only in selected sub-ranges
    if (reject)
    {   //only include points in sub-range into the fit
        if ((x[0] >= ToF_BG_min[0] && x[0] <= ToF_BG_max[0]) ||
            (x[0] >= ToF_BG_min[1] && x[0] <= ToF_BG_max[1]))
            return par[0];
        else
        {   TF1::RejectPoint();
            return 0;
        }
    }
    //compute the function values over the whole range
    else
        return par[0];
}

TF1* FitLeft(TH1F *pH, Int_t i, string FC, string Run)
{
    char name[64] = "";
    ToF_BG_min[0] = pH->GetBinLowEdge(Gate_0(i, FC));
    ToF_BG_max[0] = pH->GetBinLowEdge(Gate_a(i, FC));
    ToF_BG_min[1] = pH->GetBinLowEdge(Gate_b(i, FC));
    ToF_BG_max[1] = pH->GetBinLowEdge(Gate_3(i, FC));
    Double_t xmin = ToF_BG_min[0];
    Double_t xmax = ToF_BG_max[0];
//    cout << ToF_BG_min[0] << " " << ToF_BG_max[0] << " " << ToF_BG_min[1] << " " << ToF_BG_max[1] << " " << xmin << " " << xmax << endl;

    sprintf(name, "%s_fL_%i", Run.c_str(), i+1);
    TF1 *fLeft = new TF1(name, fline, xmin, xmax, 1);
    reject = kTRUE;
    pH->Fit(name, "LR0Q");
    return fLeft;
}

TF1* FitRight(TH1F *pH, Int_t i, string FC, string Run)
{
    char name[64] = "";
    ToF_BG_min[0] = pH->GetBinLowEdge(Gate_0(i, FC));
    ToF_BG_max[0] = pH->GetBinLowEdge(Gate_a(i, FC));
    ToF_BG_min[1] = pH->GetBinLowEdge(Gate_b(i, FC));
    ToF_BG_max[1] = pH->GetBinLowEdge(Gate_3(i, FC));
    Double_t xmin = ToF_BG_min[1];
    Double_t xmax = ToF_BG_max[1];
//    cout << ToF_BG_min[0] << " " << ToF_BG_max[0] << " " << ToF_BG_min[1] << " " << ToF_BG_max[1] << " " << xmin << " " << xmax << endl;

    sprintf(name, "%s_fR_%i", Run.c_str(), i+1);
    TF1 *fRight = new TF1(name, fline, xmin, xmax, 1);
    reject = kTRUE;
    pH->Fit(name, "LR0Q");
    return fRight;
}

TF1* FitTotal(TH1F *pH, Int_t i, string FC, string Run)
{
    char name[64] = "";
    ToF_BG_min[0] = pH->GetBinLowEdge(Gate_0(i, FC));
    ToF_BG_max[0] = pH->GetBinLowEdge(Gate_a(i, FC));
    ToF_BG_min[1] = pH->GetBinLowEdge(Gate_b(i, FC));
    ToF_BG_max[1] = pH->GetBinLowEdge(Gate_3(i, FC));
    Double_t xmin = pH->GetBinLowEdge(1);
    Double_t xmax = pH->GetBinLowEdge(pH->GetNbinsX() + 1);

    sprintf(name, "%s_fT_%i", Run.c_str(), i+1);
    TF1 *fTotal = new TF1(name, fline, xmin, xmax, 1);
    reject = kTRUE;
    pH->Fit(name, "LR0Q");
    reject = kFALSE;
    return fTotal;
}

string SubtractBackground(string Run, Bool_t print = 0)
{
    string FC = (Run[0] == 'U') ? "UFC" : "PuFC";
    TFile *fAna = TFile::Open("/home/hoffma93/Programme/Go4nfis/FC-Analysis/results/Analysis.root", "UPDATE");
    char name[128] = "";
    sprintf(name, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", Run.c_str());
    TFile *f = TFile::Open(name, "READ"); if (!f) cout << "Error opening " << name << endl;

    // live time
    TH1D *pHt = (TH1D*)f->Get("Histograms/Raw/Scaler/Rates/H1RawRate_47");
    Double_t t_live = pHt->Integral();

    // Create graphs
    TGraphErrors *geNIF = new TGraphErrors(8);
    TGraphErrors *geBG = new TGraphErrors(8);
    TGraphErrors *geNIFrate = new TGraphErrors(8);
    TGraphErrors *geBGrate = new TGraphErrors(8);
    TGraph *geChi2dofLeft = new TGraphErrors(8);
    geChi2dofLeft->SetNameTitle("chi2dofLeft", "Reduced #chi^{2}; Deposit; #chi^{2}/#font[12]{dof}");
    TGraph *geChi2dofRight = new TGraphErrors(8);
    geChi2dofRight->SetNameTitle("chi2dofRight", "Reduced #chi^{2}; Deposit; #chi^{2}/#font[12]{dof}");
    TGraph *geChi2dofTotal = new TGraphErrors(8);
    geChi2dofTotal->SetNameTitle("chi2dofTotal", "Reduced #chi^{2}; Deposit; #chi^{2}/#font[12]{dof}");

    // Initialize returned string
    std::stringstream line;
    line << Run;

    for (Int_t i = 0; i < 8; i++)
    {
        // Open ToF spectrum
        TH1F *pH = GetToF(f, i, FC);
        Save(fAna, FC+"/ToF/Total/"+Run, pH);
        Double_t C_total, DC_total;
        C_total = pH->IntegralAndError(0, -1, DC_total);
        // Fit
        TF1 *fLeft = FitLeft(pH, i, FC, Run);
        Save(fAna, FC+"/ToF/Background/"+Run+"/Left", fLeft);
        TF1 *fRight = FitRight(pH, i, FC, Run);
        Save(fAna, FC+"/ToF/Background/"+Run+"/Right", fRight);
        TF1 *fTotal = FitTotal(pH, i, FC, Run);
        Save(fAna, FC+"/ToF/Background/"+Run+"/Total", fTotal);
        // Subtract
        Double_t x1, x2;
        fTotal->GetRange(x1, x2);
//        cout << x1 << " " << x2 << " " << pH->FindBin(x1) << " " << pH->FindBin(x2) << endl;
//        if (Run[0] == 'U') // UFC only uses left background interval
//            pH->Add(fLeft, -1);
//        else
            pH->Add(fTotal, -1);
        Save(fAna, FC+"/ToF/Signal/"+Run, pH);

        // chi2 / dof
        geChi2dofLeft->SetPoint(i, i+1, fLeft->GetChisquare() / fLeft->GetNDF());
        geChi2dofRight->SetPoint(i, i+1, fRight->GetChisquare() / fRight->GetNDF());
        geChi2dofTotal->SetPoint(i, i+1, fTotal->GetChisquare() / fTotal->GetNDF());
        if (print)
        {
            if (Run[0] == 'U') // UFC
                sprintf(name, "%.2f(%i) & %.2f & %.2f(%i) & %.2f & %.2f(%i)",
                        fLeft->GetParameter(0), (Int_t)(100*fLeft->GetParError(0)),
                        fLeft->GetChisquare() / fLeft->GetNDF(),
                        fRight->GetParameter(0), (Int_t)(100*fRight->GetParError(0)),
                        fRight->GetChisquare() / fRight->GetNDF(),
                        fTotal->GetParameter(0), (Int_t)(100*fTotal->GetParError(0)));
            else
                sprintf(name, "%.1f(%i) & %.2f & %.1f(%i) & %.2f & %.1f(%i)",
                        fLeft->GetParameter(0), (Int_t)(10*fLeft->GetParError(0)),
                        fLeft->GetChisquare() / fLeft->GetNDF(),
                        fRight->GetParameter(0), (Int_t)(10*fRight->GetParError(0)),
                        fRight->GetChisquare() / fRight->GetNDF(),
                        fTotal->GetParameter(0), (Int_t)(10*fTotal->GetParError(0)));
            cout << FC << " " << i+1 << " & " << name << " \\\\" << endl;
        }

        // Integrate
        Double_t C_nif, DC_nif;
        C_nif = pH->IntegralAndError(Gate_1(i, FC), Gate_2(i, FC), DC_nif);
        Double_t C_bg, DC_bg;
        C_bg = C_total - C_nif;
        DC_bg = sqrt(pow(DC_total , 2) + pow(DC_nif, 2));

//        cout << DC_nif << "  " << sqrt(C_total + C_sf) << endl;
        geBG->SetPoint(i, i+1, C_bg);
        geBG->SetPointError(i, 0, DC_bg);
        geNIF->SetPoint(i, i+1, C_nif);
        geNIF->SetPointError(i, 0, DC_nif);
        geBGrate->SetPoint(i, i+1, C_bg / t_live);
        geBGrate->SetPointError(i, 0, DC_bg / t_live);
        geNIFrate->SetPoint(i, i+1, C_nif / t_live);
        geNIFrate->SetPointError(i, 0, DC_nif / t_live);

//        cout << DC_nif / C_nif << endl;
        line << " " << C_nif << " " << DC_nif << " " << C_nif / t_live << " " << DC_nif / t_live;
    }
    Save(fAna, FC+"/ToF/Signal/"+Run, geNIF, "InducedFission");
    Save(fAna, FC+"/ToF/Background/"+Run, geBG, "FissionBackground");
    Save(fAna, FC+"/ToF/Signal/"+Run, geNIFrate, "FissionRate");
    Save(fAna, FC+"/ToF/Background/"+Run, geBGrate, "BackgroundRate");
    Save(fAna, FC+"/ToF/Background/"+Run+"/Left", geChi2dofLeft);
    Save(fAna, FC+"/ToF/Background/"+Run+"/Right", geChi2dofRight);
    Save(fAna, FC+"/ToF/Background/"+Run+"/Total", geChi2dofTotal);

//    TCanvas *c1 = new TCanvas();
//    pH->Draw("hist");
//    fTotal->Draw("same");
//    fLeft->Draw("same");
//    fRight->Draw("same");
//    c1->Draw();

    fAna->Save();
    fAna->Close();
    return line.str();
}

void DoToF(string FC)
{
    string FileName = "InducedFission_"+FC+".txt";
    ofstream output("../results/"+FileName);
    output << "# nr run {indfis unc indfis[1/s] unc[1/s]}" << endl;
    Int_t nRuns = GetnRuns(FC);
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(FC, j);
        output << j+1 << " " << SubtractBackground(Run) << endl;
    }
    output.close();
}

void DoToF()
{
    string FileName = "InducedFission.txt";
    ofstream output("../results/"+FileName);
    output << "# nr run {indfis unc indfis[1/s] unc[1/s]}" << endl;
    Int_t nRuns = GetnRuns();
    for (Int_t j = 0; j < nRuns; j++)
    {
        string Run = GetRunName(j);
        output << j+1 << " " << SubtractBackground(Run) << endl;
    }
    output.close();
}

void ToF()
{
    DoToF("UFC");
    DoToF("UFC_FG");
    DoToF("UFC_BG");
    cout << SubtractBackground("UFC_NIF") << endl;
    SubtractBackground("UFC_SB");

    DoToF("PuFC");
    DoToF("PuFC_FG");
    DoToF("PuFC_BG");
    cout << SubtractBackground("NIF") << endl;
    SubtractBackground("SB");
    //*/
    /*SubtractBackground("PuFC_FG_MS4", "PuFC");
    SubtractBackground("PuFC_FG_MS5", "PuFC");
    SubtractBackground("PuFC_FG_MS6", "PuFC");
    SubtractBackground("PuFC_FG_MS7", "PuFC");
    SubtractBackground("PuFC_BG_MS9", "PuFC");
    SubtractBackground("PuFC_BG_MS10", "PuFC");
    SubtractBackground("PuFC_BG_MS11", "PuFC");
    SubtractBackground("NIF", "PuFC");
    SubtractBackground("SB", "PuFC");
    SubtractBackground("UFC_FG_MS20_2", "UFC");
    SubtractBackground("UFC_FG_MS20_3", "UFC");
    SubtractBackground("UFC_FG_MS20_4", "UFC");
    SubtractBackground("UFC_FG_MS21_2", "UFC");
    SubtractBackground("UFC_FG_MS21_3", "UFC");
    SubtractBackground("UFC_BG_MS20_5", "UFC");
    SubtractBackground("UFC_BG_MS21_4", "UFC");
    SubtractBackground("UFC_NIF", "UFC");
    SubtractBackground("UFC_SB", "UFC");//*/

}
#endif
