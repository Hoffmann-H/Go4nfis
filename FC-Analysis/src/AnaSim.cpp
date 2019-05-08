#include "AnaSim.h"
using namespace std;

AnaSim::AnaSim(string file_name, string fc, string setup, Plot *p)
{
    cout << endl << "Creating simulation analysis instance for " << fc << ", " << setup << endl;
//    FilePath = "";
    FileName = file_name;
    FC = fc;
    PuFC = strcmp(fc.c_str(), "UFC");
    Setup = setup;
    // File name syntax: "PxxxxO..." PuFC, Open
    //                   "Pxxxxx..." PuFC, SB
    //                   "UxxxO..." UFC, Open
    //                   "Uxxxx..." UFC, SB
//    if (FileName[0] == 'P') {
//        FC = "PuFC";
//        PuFC = kTRUE;
//        if (FileName[5] == 'O')
//            Setup = "Open";
//        else
//            Setup = "SB";
//    } else if(FileName[0] == 'U') {
//        FC = "UFC";
//        PuFC = kFALSE;
//        if (FileName[4] == 'O')
//            Setup = "Open";
//        else
//            Setup = "SB";
//    } else {

    DrawSingle = kTRUE;
    DrawMulti = kTRUE;
    if (p != 0)
        plot = p;
    else {
        DrawSingle = kFALSE;
        DrawMulti = kFALSE; }
    CommentFlag = kTRUE;
    DoneHist = kFALSE;
    DoneProjection = kFALSE;
    DoneDistance = kFALSE;
    DoneToFwidth = kFALSE;
    DoneSigma = kFALSE;
    DoneFis = kFALSE;
    DoneDirect = kFALSE;
    DoneCorrections = kFALSE;

    cout << " Opening " << FileName << endl;
    f = TFile::Open(file_name.c_str(), "UPDATE");
    cout << " Opened successfully!" << endl;
    OpenHists();
    SetDistances();
    SetToFwidths();
    OpenSigma("/home/hoffma93/Programme/ROOT/Data");
    UisoVec[0] = 0.904, DUisoVec[0] = 0.005; // U-235 portion
    UisoVec[1] = 0.0912, DUisoVec[1] = 0.0006; // U-238 portion

    cout << "Done: Create simulation analysis instance for " << FC << ", " << Setup << endl;
}


AnaSim::~AnaSim()
{

}


void AnaSim::OpenHists()
{
    cout << "Open histograms..." << endl;
    char name[64] = "";
    for (Int_t i = 0; i < NumCh; i++)
    {
        // Open Ekin vs ToF
        sprintf(name, "/%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        cout << " " << name << endl;
        pH2TvsE[i] = (TH2F*)f->Get(name);
    }
    DoneHist = kTRUE;
    cout << "Done: Open histograms" << endl;
}


void AnaSim::Projections()
{
    if (!DoneHist)
        OpenHists();
    if (!DoneToFwidth)
        SetToFwidths();
    cout << endl << "Projections..." << endl;

    char name[64] = "";
    for (Int_t i = 0; i < NumCh; i++)
    {
        // y projection
        sprintf(name, "proj_E_%i", i+1);
        pH1Eproj[i] = (TH1F*)pH2TvsE[i]->ProjectionY(name, binToFmin[i], binToFmax[i]);
        sprintf(name, "%s, time-gated neutron spectrum, ch. %i", FC.c_str(), i+1);
        pH1Eproj[i]->SetTitle(name);
        sprintf(name, "proj_y_%i", i+1);
        if (DrawMulti)
            plot->simply(pH1Eproj[i], name, "xTitle", "yTitle");
    }
    DoneProjection = kTRUE;
    cout << "Done: Projections" << endl;
}


void AnaSim::SetDistances()
{// Define distances from front window to deposit in m
    cout << "Define distances source-deposit..." << endl;
    if (strcmp(FC.c_str(), "UFC"))
        // PuFC
        for (int i = 0; i < NumCh; i++)
            Distance[i] = 1.6990 - i * 0.0205;
    else
        // UFC
        for (int i = 0; i < NumCh; i++)
            Distance[i] = 1.6299 - i * 0.0108;
    DoneDistance = kTRUE;
    cout << "Done: Distances for " << FC << endl;
}


void AnaSim::SetToFwidths()
{// Define ToF and energy range
    if (!DoneHist)
        OpenHists();
    if (!DoneDistance)
        SetDistances();

    char name[64] = "/home/hoffma93/Programme/ROOT/Data/Source_E.dat";
    cout << endl << "ToF range and energy spectrum using " << name << endl;
    TGraph *gE = new TGraph(name);
    GetEwidth(gE);
    if (DrawSingle)
        plot->Source_E(gE, En, DEn);
    cout << "Check 1" << endl;
    Double_t dummy = 0;
    gE->GetPoint(0, Emin, dummy);
    gE->GetPoint(gE->GetN()-1, Emax, dummy);
    cout << "Check 2" << endl;
    Double_t binOffset = pH2TvsE[0]->GetYaxis()->GetBinLowEdge(0);
    Double_t binWidth = pH2TvsE[0]->GetYaxis()->GetBinWidth(0);
    cout << "Check 3" << endl;
    binEmin = (Emin - binOffset) / binWidth;
    binEmax = (Emax - binOffset) / binWidth;
    cout << " Bin nr: " << binEmin << "-" << binEmax << ", " << Emin << " < E < " << Emax << endl;
    cout << "Ch   bins   ToF" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
//        cout << "Check " << i << endl;
        ToFmin[i] = ToF(Emax, Distance[i]);
        ToFmax[i] = 200;//ToF(Emin, sqrt( pow(Distance[i], 2) + pow(0.037, 2) ));
        Double_t binOffset = pH2TvsE[i]->GetXaxis()->GetBinLowEdge(0);
        Double_t binWidth = pH2TvsE[i]->GetXaxis()->GetBinWidth(0);
        binToFmin[i] = (ToFmin[i] - binOffset) / binWidth;
        binToFmax[i] = (ToFmax[i] - binOffset) / binWidth;
        cout << " " << i+1 << "   " << binToFmin[i] << "-" << binToFmax[i] << "   " << ToFmin[i] << " < ToF < " << ToFmax[i] << endl;
    }
    DoneToFwidth = kTRUE;
    cout << "Done: ToF and energy range" << endl;
}


void AnaSim::GetEwidth(TGraph *gE)
{
    cout << " Calculate source energy standard deviation..." << endl;
    Int_t N = gE->GetN();
    Double_t E, w; // Energy, weight read-out vars
    Double_t sum = 0, sumw = 0, var = 0;
    for (Int_t i = 0; i < N; i++)
    {
        gE->GetPoint(i, E, w);
        sum += w * E;
        sumw += w;
    }
    En = sum / sumw;
    for (Int_t i = 0; i < N; i++)
    {
        gE->GetPoint(i, E, w);
        var += w / sumw * pow(E - En, 2); // variance = sum_i(h_i * (x_i - <x>)^2)
    }
    DEn = sqrt(var);
    cout << "  E_n = " << En << " +- " << DEn << " MeV" << endl;
    cout << " Done: Energy standard deviation" << endl;
}


void AnaSim::OpenSigma(string path)
{
    char name[64] = "";

    sprintf(name, "%s/Pu242.dat", path.c_str());
    gPu242 = new TGraph(name);
    DgPu242 = new TGraph(name, "%lg %*lg %lg");
    sprintf(name, "%s/U235.dat", path.c_str());
    gU235 = new TGraph(name);
    DgU235 = new TGraph(name, "%lg %*lg %lg");
    sprintf(name, "%s/U238.dat", path.c_str());
    gU238 = new TGraph(name);
    DgU238 = new TGraph(name, "%lg %*lg %lg");

    DoneSigma = kTRUE;

    if (!DrawSingle)
        return;
    DrawSigma();
}


void AnaSim::DrawSigma()
{
    if (!DoneSigma)
        cerr << "No cross section data! " << endl;

    TCanvas *pC1 = new TCanvas("Pu242", "Pu-242 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gPu242->SetTitle("Pu-242 (n,f); E / MeV; cross section / barn");
    gPu242->Draw();

    TCanvas *pC2 = new TCanvas("U235", "U-235 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gU235->SetTitle("U-235 (n,f); E / MeV; cross section / barn");
    gU235->Draw();

    TCanvas *pC3 = new TCanvas("U238", "U-238 cross sections", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    gU238->SetTitle("U-238 (n,f); E / MeV; cross section / barn");
    gU238->Draw();
}


void AnaSim::SaveToFile(string RootPath, TObject *pObj)
{   //saves a TObject into member TFile *f
    f->ReOpen("UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder path already exists, otherwise create it
    if (f->Get(RootPath.c_str())!=0)
        EvalDir = (TDirectory*) f->Get(RootPath.c_str());
    else
    {   cout << " Creating root directory " << RootPath << endl;
        f->mkdir(RootPath.c_str(), "Folder");
        EvalDir = f->GetDirectory(RootPath.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    f->Save(); //file->Close();
    cout << " Saved " << GraphName << endl;
}


TH1D* AnaSim::CopyRange(TH1D* pH, char* name, Double_t x0, Double_t x1, Double_t yoffset)
{
    Double_t offset = pH->GetBinLowEdge(0);
    Double_t ChPerBin = pH->GetBinWidth(0);
    int b0 = (x0 - offset) / ChPerBin;
    int b1 = (x1 - offset) / ChPerBin + 1;
    TH1D* pH2 = new TH1D(name, name, b1 - b0, pH->GetBinLowEdge(b0), pH->GetBinLowEdge(b1));
    for (int bin = 0; bin < pH2->GetNbinsX(); bin++)
        pH2->SetBinContent(bin + 1, pH->GetBinContent(bin + b0) + yoffset);
    return pH2;
}


Double_t AnaSim::Et(Double_t ToF, Double_t FlightPath)
{// FlightPath in m, ToF in ns
    Double_t beta = FlightPath / ToF * 1.E9 / c;
    Double_t g = 1.0 / sqrt(1.0 - beta * beta);
    Double_t Ekin = (g - 1.0) * m0; // MeV
//    cout << ToF << " ns --> " << Ekin << " MeV" << endl;
    return Ekin;
}


Double_t AnaSim::ToF(Double_t Ekin, Double_t FlightPath)
{
    Double_t g = Ekin / m0 + 1.0;
    Double_t beta = sqrt(1.0 - 1.0 / (g*g));
    Double_t ToF = FlightPath / beta * 1.E9 / c;
//    cout << "E " << Ekin << ", Path " << FlightPath << ", ToF " << ToF << endl;
    return ToF;
}


void AnaSim::nProjectiles()
{// return the number of incident neutrons assuming 100% transmission
    if (!DoneDistance)
        SetDistances();
    cout << endl << "Get number of projectiles towards deposits..." << endl;
    char name[64] = "";
    cout << "Ch   projetiles" << endl;
    if (f->Get("Source/Source_Theta_Ch.1") != 0) // If single angular hists are available...
    {
        for (int i = 0; i < NumCh; i++)
        {
            sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
            TH1F *pH1Ang = (TH1F*) f->Get(name);
            nProj[i] = pH1Ang->Integral();
            cout << " " << i+1 << "   " << nProj[i] << endl;
        }
    } else { // No angular hists. Calculate...
        TH1F *pH1Ekin = (TH1F*) f->Get("Source/nEnergy/Source_Ekin");
        Double_t N = pH1Ekin->Integral();
        Double_t theta_max = 0.074641914;
        for (int i = 0; i < NumHist; i++)
        {
            Double_t theta_dep = atan(0.037 / Distance[i]);
            nProj[i] = (1 - pow(cos(theta_dep), 2)) / (1 - pow(cos(theta_max), 2)) * N;
            cout << " " << i+1 << "   " << nProj[i] << endl;
        }
    }
    DoneProjectiles = kTRUE;
    cout << "Done: projectiles" << endl;
    return;
}


void AnaSim::relSigma(Double_t E, Double_t &w, Double_t &Dw)
{
    if(PuFC)
    {
        w = gPu242->Eval(E) / gPu242->Eval(15.0);
        Dw = sqrt( pow(DgPu242->Eval(E) / gPu242->Eval(15.0), 2) +
                        pow(DgPu242->Eval(15.0) * gPu242->Eval(E) / pow(gPu242->Eval(15.0), 2), 2) );
    } else { // UFC
        Double_t numerator = gU235->Eval(E) * UisoVec[0] + gU238->Eval(E) * UisoVec[1];
        Double_t denominator = gU235->Eval(15.0) * UisoVec[0] + gU238->Eval(15.0) * UisoVec[1];
        w = numerator / denominator;
        Dw = sqrt( pow(DgU235->Eval(E)    * UisoVec[0] / denominator                          , 2) + // U-235 sigma(E)
                        pow(DgU235->Eval(15.0) * UisoVec[0] * numerator / (denominator*denominator), 2) + // U-235 sigma(15MeV)
                        pow(DUisoVec[0]        * (gU235->Eval(E) * gU238->Eval(15.0) - gU235->Eval(15.0) * gU238->Eval(E))
                                               * UisoVec[1] / (denominator*denominator)            , 2) + // U-235 portion
                        pow(DgU238->Eval(E)    * UisoVec[1] / denominator, 2) + // U-238 sigma(E)
                        pow(DgU238->Eval(15.0) * UisoVec[1] * numerator / (denominator*denominator), 2) + // U-238 sigma(15MeV)
                        pow(DUisoVec[1]        * (gU238->Eval(E) * gU235->Eval(15.0) - gU238->Eval(15.0) * gU235->Eval(E))
                                               * UisoVec[0] / (denominator*denominator)            , 2) ); // U-238 portion
    }
}


void AnaSim::nFissions()
{// Calculate the number of n-induced fission events within the analysis' time gate. Neglect pulse duration
    if (!DoneHist) // Open ToFvsEkin histogram
        OpenHists();
    if (!DoneSigma)
        cerr << "No cross section data! " << endl;
    if (!DoneDistance)
        SetDistances();
    if (!DoneToFwidth)
        SetToFwidths();
    if (!DoneProjection)
        Projections();
    // ...
    cout << endl << "Calculating effective number of fissioning neutrons..." << endl;

    char name[64] = "";
    char title[64] = "";
    Int_t nbinsx;
    Double_t xmin, xmax;

    for (int i = 0; i < NumCh; i++)
    {
        // Draw projection
        if (DrawMulti)
        {
            sprintf(name, "%s_%s_EkinG_%i", FC.c_str(), Setup.c_str(), i+1);
            plot->simply(pH1Eproj[i], (string)name, "Ekin / MeV", "N");
        }

        // Calculate number of sc + unsc effective neutrons for fission
        nbinsx = pH1Eproj[i]->GetNbinsX();
        xmin = pH1Eproj[i]->GetXaxis()->GetBinLowEdge(1);
        xmax = pH1Eproj[i]->GetXaxis()->GetBinUpEdge(nbinsx);
        sprintf(name, "%s_%s_Y_%i", FC.c_str(), Setup.c_str(), i+1);
        sprintf(title, "%s, %s, fission Yield in ToF range, ch %i", FC.c_str(), Setup.c_str(), i+1);
        TH1D* pH1YG = new TH1D(name, title, nbinsx, xmin, xmax);
        Double_t D2nFis = 0;
        for (Int_t bin = 1; bin <= pH1Eproj[i]->GetNbinsX(); bin++)
        {// Iterate over projection's energy bins
            Double_t E = pH1Eproj[i]->GetBinCenter(bin);
            Double_t nNeutrons = pH1Eproj[i]->GetBinContent(bin);
            Double_t DnNeutrons = sqrt(nNeutrons); // Poisson statistics
            Double_t weight, Dweight;

            relSigma(E, weight, Dweight); //// Get relative cross section

            pH1YG->SetBinContent(bin, weight * nNeutrons);
            D2nFis += pow(Dweight * nNeutrons, 2) + pow(weight * DnNeutrons, 2);
        }
        nFis[i] = pH1YG->Integral(); // Fission count as estimated by analysis
        DnFis[i] = sqrt(D2nFis);
        cout << "Ch " << i+1 << ", eff.n " << nFis[i] << " +- " << DnFis[i] << endl;
        if (DrawMulti)
            plot->TvsE(i, pH2TvsE[i], binToFmin[i], binToFmax[i], binEmin, binEmax);
        if (DrawMulti)
            PeakForm(i);
    }
    DoneFis = kTRUE;
    cout << "Done: effective fissioning neutrons" << endl;
}


void AnaSim::PeakForm(Int_t ch)
{
    if (!DoneHist)
        OpenHists();
    if (!DoneToFwidth)
        SetToFwidths();
    if (!DoneProjectiles)
        nProjectiles();
//    if (!DoneDirect)
//        DirectN();
    char name[64] = "";
    sprintf(name, "proj_ToF_%i", ch+1);
    TH1F *pH1X = (TH1F*)pH2TvsE[ch]->ProjectionX(name);
    Int_t Nx = pH1X->GetNbinsX();
    Int_t Ny = pH2TvsE[ch]->GetNbinsY();
    Double_t x0 = pH1X->GetBinLowEdge(1);
    Double_t x1 = pH1X->GetBinLowEdge(Nx+2);
    sprintf(name, "t_Sc_%i", ch+1);
    TH1F *pH1Sc = new TH1F(name, "Effective incident neutrons", Nx, x0, x1);
    sprintf(name, "t_Unsc_%i", ch+1);
    TH1F *pH1Unsc = new TH1F(name, "Effective incident neutrons", Nx, x0, x1);

    for (Int_t j = 0; j < Nx+2; j++)
    {
        if (pH1X->GetBinContent(j) != 0)
        {
            Double_t sumSc = 0;
            Double_t sumUnsc = 0;
            for (Int_t k = 1; k < Ny+1; k++)
            {
                Double_t E = pH1X->GetBinCenter(k);
                Double_t weight = 0, Dweight = 0; // weighting by (n,f) cross section with respect to 15 MeV
                if(PuFC)
                    weight = gPu242->Eval(E) / gPu242->Eval(15.0);
                else { // UFC
                    Double_t numerator = gU235->Eval(E) * UisoVec[0] + gU238->Eval(E) * UisoVec[1];
                    Double_t denominator = gU235->Eval(15.0) * UisoVec[0] + gU238->Eval(15.0) * UisoVec[1];
                    weight = numerator / denominator; }
                sumSc += (E < Emin) ? weight * pH2TvsE[ch]->GetBinContent(j, k) : 0;
                sumUnsc += (E >= Emin) ? weight * pH2TvsE[ch]->GetBinContent(j, k) : 0;
            } // for (k)
            pH1Sc->SetBinContent(j, sumSc);
            pH1Unsc->SetBinContent(j, sumUnsc);
        } // if (!=0)
    } // for (j)
    plot->DtPeakForm(ch, pH1Sc, pH1Unsc, Emin, nProj[ch]);
}


void AnaSim::DirectN()
{
    if (!DoneHist)
        OpenHists();
    if (!DoneToFwidth)
        SetToFwidths();
    if (!DoneProjection)
        Projections();
    cout << endl << "Getting number of direct neutrons" << endl;
    char name[64] = "";
    for (int i = 0; i < NumCh; i++)
    {
        nFullE[i] = pH2TvsE[i]->Integral(binToFmin[i], binToFmax[i], binEmin, binEmax);
        DnFullE[i] = sqrt(nFullE[i]); // use Poisson statistics because nFullE is counted
        Double_t E, N;
        Double_t w = 0, Dw = 0;
        Double_t effN = 0, D2effN = 0;
        for (Int_t bin = binEmin; bin <= binEmax; bin++)
        {
            E = pH1Eproj[i]->GetBinCenter(bin);
            N = pH1Eproj[i]->GetBinContent(bin);
            relSigma(E, w, Dw);
            effN += w * N;
            D2effN += /*pow(Dw * N, 2) +*/ w * w * N;
        }
        sprintf(name, "/%s/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        TH2F *pH2TvsE_Sc = (TH2F*) f->Get(name);
        effFullE[i] = effN;
        DeffFullE[i] = sqrt(D2effN);
        nDirect[i] = pH2TvsE[i]->Integral() - pH2TvsE_Sc->Integral();
        DnDirect[i] = sqrt(nDirect[i]);
        if (CommentFlag)
            cout << " Ch " << i+1 << ", full-energy neutrons " << nFullE[i] << "+-" << DnFullE[i] <<
                    ", effective full-E neutrons " << effFullE[i] << "+-" << DeffFullE[i] <<
                    ", direct neutrons " << nDirect[i] << "+-" << DnDirect[i] << endl;
    }
    DoneDirect = kTRUE;
    cout << "Done: Direct neutrons" << endl;
}


void AnaSim::Corrections()
{// Calculate correction factors for data analysis from simulation.
    // T := Direct incident n / emitted neutrons towards deposit
    // S := fissions induced by direct / induced fissions (both in time gate)
    //    = direct fissions / ( direct fissions + scattered fissions )
    if (!DoneProjectiles)
        nProjectiles();
    if (!DoneDirect)
        DirectN();

    /////////////////////////////////////////////////////////////////////////////
    /// Get effective number of scattered neutrons
    if (!DoneFis)
        nFissions();
    /////////////////////////////////////////////////////////////////////////////

    cout << endl << "Correction factors..." <<
            endl << "Ch   T   S   F" << endl;

    for (int i = 0; i < NumCh; i++)
    {
//        cout << 0 / nProj[i] << " " << DnDirect[i] / nDirect[i] << " " << DnFis[i] / nFis[i] << " " << DnFullE[i] / nFullE[i] << " " << DeffFullE[i] / effFullE[i] << endl;
        T[i] = nDirect[i] / nProj[i];
        DT[i] = DnDirect[i] / nProj[i];
        S[i] = nDirect[i] / nFis[i] * effFullE[i] / nFullE[i];
        DS[i] = S[i] * sqrt( pow(DnDirect[i] / nDirect[i], 2) +
                             pow(DnFis[i] / nFis[i], 2) +
                             pow(DeffFullE[i] / effFullE[i], 2) +
                             pow(DnFullE[i] / nFullE[i], 2) );
        F[i] = nProj[i] / nFis[i] * effFullE[i] / nFullE[i];
        DF[i] = F[i] * sqrt( pow(DnFis[i] / nFis[i], 2) +
                             pow(DeffFullE[i] / effFullE[i], 2) +
                             pow(DnFullE[i] / nFullE[i], 2) );
        cout << " " << i+1 << "  " << T[i] << "+-" << DT[i] << "  " << S[i] << "+-" << DS[i] << "  " << F[i] << "+-" << DF[i] << endl;
//        cout << " " << nDirect[i] << " " << nProj[i] << endl;
    }
    cout << "Uncertainty influences on Scattering correction factor" << endl
         << " Ch   nDirect   nFis   effFullE   nFullE" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        Double_t unc[] = {S[i] * DnDirect[i] / nDirect[i],
                          S[i] * DnFis[i] / nFis[i],
                          S[i] * DeffFullE[i] / effFullE[i],
                          S[i] * DnFullE[i] / nFullE[i]};
        cout << " " << i+1 << "   " << unc[0] << "   " << unc[1] << "   " << unc[2] << "   " << unc[3] << endl;
    }
    DoneCorrections = kTRUE;
    cout << "Done: Corrections" << endl;
    if (!DrawSingle)
        return;
    plot->SimT(T, DT, S, DS, F, DF);
}

//*/
