#include "Sim.h"
using namespace std;

Sim::Sim(string file_name, string fc, string setup, Bool_t use_track_id, Bool_t draw_flag)
{
    cout << endl << "Starting " << fc << " " << setup << " simulation analysis" << endl;
    // Input settings
    FileName = file_name;
    FC = fc;
    PuFC = strcmp(fc.c_str(), "UFC");
    Setup = setup;

    // processing settings
    tID = use_track_id;
    DoneHist = kFALSE;
    DoneDistance = kFALSE;
    DoneSigma = kFALSE;
    DoneProjectiles = kFALSE;
    DoneEwidth = kFALSE;
    DoneTwidth = kFALSE;
    DoneDirect = kFALSE;
    DoneCalc = kFALSE;

    // output settings
    CommentFlag = kFALSE; // set manually
    DrawSingle = draw_flag;
    DrawMulti = kTRUE;

    cout << " Opening " << file_name << endl;
    f = TFile::Open(FileName.c_str(), "UPDATE");
    if (f == 0)
        cout << " Error while opening " << FileName << endl;
    if (DrawSingle || DrawMulti)
        plot = new Plot(FC, Setup, f);
    OpenHists();
    SetDistances();
    nProjectiles();
    GetSigma("/home/hoffma93/Programme/ROOT/Data");
    GetEwidth();
    GetTwidth();

    UisoVec[0] = 0.904, DUisoVec[0] = 0.005; // U-235 portion
    UisoVec[1] = 0.0912, DUisoVec[1] = 0.0006; // U-238 portion
    cout << "Started simulation analysis" << endl;
}


Sim::~Sim()
{

}


void Sim::OpenHists()
{
    cout << endl << "Open histograms..." << endl;
    char name[64] = "";
    for (Int_t i = 0; i < NumCh; i++)
    {
        // Open Ekin vs ToF
        sprintf(name, "/%s/ToFvsEkin/%s_ToFvsEkin_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        cout << " " << name << endl;
        pH2TvsE[i][0] = (TH2F*) f->Get(name);
        sprintf(name, "/%s/ToFvsEkin/Scattered/%s_ToFvsEkin_Sc_Ch.%i", FC.c_str(), FC.c_str(), i+1);
        pH2TvsE[i][1] = (TH2F*) f->Get(name);
    }
    DoneHist = kTRUE;
    cout << "Done: Open histograms" << endl;
}


void Sim::SetDistances()
{// Define distances from front window to deposit in m
    cout << endl << "Define distances source-deposit..." << endl;
    if (strcmp(FC.c_str(), "UFC"))
        // PuFC
        for (int i = 0; i < NumCh; i++)
            Distance[i] = 1.7092 - i * 0.0205;
    else
        // UFC
        for (int i = 0; i < NumCh; i++)
            Distance[i] = 1.6342 - i * 0.0105;
    cout << " Ch   Distance/m" << endl;
    for (int i = 0; i < NumCh; i++)
        cout << " " << i+1 << "   " << Distance[i] << endl;
    DoneDistance = kTRUE;
    cout << "Done: Distances for " << FC << endl;
}


void Sim::GetSigma(string path)
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
    plot->Sigma(gPu242, gU235, gU238);
}


void Sim::relSigma(Double_t E, Double_t &w, Double_t &Dw)
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
        Dw = 0;
    }
}


void Sim::GetEwidth()
{
    cout << endl << "Getting source's energy width..." << endl;
    char name[64] = "/home/hoffma93/Programme/ROOT/Data/Source_E.dat";
//    char name[64] = "/home/hoffma93/Programme/ROOT/Data/nELBE_E.dat";
    TGraph *gE = new TGraph(name);
    cout << "Source energy spectrum at " << name << endl;

    // Energy spectrum's standard deviation
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
    cout << " E_n = " << En << " +- " << DEn << " MeV" << endl;
    if (DrawSingle)
        plot->Source_E(gE, En, DEn, 1);

    // Energy spectrum's edges
    Double_t dummy = 0;
    gE->GetPoint(0, Emin, dummy);
    gE->GetPoint(gE->GetN()-1, Emax, dummy);
    Double_t binOffset = pH2TvsE[0][0]->GetYaxis()->GetBinLowEdge(0);
    Double_t binWidth = pH2TvsE[0][0]->GetYaxis()->GetBinWidth(0);
    binEmin = (Emin - binOffset) / binWidth;
    binEmax = (Emax - binOffset) / binWidth;

    DoneEwidth = kTRUE;
    cout << "Done: Source energy width" << endl;
}


void Sim::GetTwidth()
{
    if (!DoneEwidth)
        GetEwidth();
    cout << endl << "Getting direct neutrons' time width..." << endl;
    cout << " Ch      range      bin nr.s" << endl;
    cout << " E   " << Emin << "-" << Emax << "   " << binEmin << "-" << binEmax << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        ToFmin[i] = ToF(Emax, Distance[i]);
        ToFmax[i] = ToF(Emin, sqrt( pow(Distance[i], 2) + pow(0.037, 2) ));
        Double_t binOffset = pH2TvsE[i][0]->GetXaxis()->GetBinLowEdge(0);
        Double_t binWidth = pH2TvsE[i][0]->GetXaxis()->GetBinWidth(0);
        binToFmin[i] = (ToFmin[i] - binOffset) / binWidth;
//        binToFmax[i] = (ToFmax[i] - binOffset) / binWidth;
        binToFmax[i] = pH2TvsE[i][0]->GetNbinsX();
        cout << " " << i+1 << "   " << ToFmin[i] << "-" << ToFmax[i] << "   " << binToFmin[i] << "-" << binToFmax[i] << endl;
    }
    cout << "Done: Time width" << endl;
}


void Sim::Projections()
{
    if (!DoneProjectiles)
        nProjectiles();
    if (!DoneSigma)
        cout << "No Cross section data!" << endl;
    cout << endl << "Projections..." << endl;
    char name[128] = "";
    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "%s_%s_ProjE_%i", FC.c_str(), Setup.c_str(), i+1);
        pH1Eproj[i][0] = (TH1F*) pH2TvsE[i][0]->ProjectionY(name, binToFmin[i], binToFmax[i]);
        sprintf(name, "%s, %s, neutron energy spectrum, ch. %i", FC.c_str(), Setup.c_str(), i+1);
        pH1Eproj[i][0]->SetTitle(name);
        sprintf(name, "%s_%s_ProjE_Sc_%i", FC.c_str(), Setup.c_str(), i+1);
        pH1Eproj[i][1] = (TH1F*) pH2TvsE[i][1]->ProjectionY(name, binToFmin[i], binToFmax[i]);
        sprintf(name, "%s, %s, scattered neutron energy spectrum, ch. %i", FC.c_str(), Setup.c_str(), i+1);
        pH1Eproj[i][1]->SetTitle(name);
        if (DrawMulti)
        { // Calculate weighted time spectrum. Could take some time.
            if (PuFC)
                plot->Eproj(i, pH1Eproj[i][1], pH1Eproj[i][0], nProj[i]);//, gPu242);
            else
                plot->Eproj(i, pH1Eproj[i][1], pH1Eproj[i][0], nProj[i]);
            Int_t NbinsX = pH2TvsE[i][0]->GetNbinsX();
            Int_t NbinsY = pH2TvsE[i][0]->GetNbinsY();
            Double_t xmin = pH2TvsE[i][0]->GetXaxis()->GetBinLowEdge(1);
            Double_t xmax = pH2TvsE[i][0]->GetXaxis()->GetBinLowEdge(NbinsX + 1);
            Double_t ymin = pH2TvsE[i][0]->GetYaxis()->GetBinLowEdge(1);
            Double_t ymax = pH2TvsE[i][0]->GetYaxis()->GetBinLowEdge(NbinsX + 1);
            sprintf(name, "%s_%s_ProjT_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Tproj[i][0] = new TH1F(name, name, NbinsX, xmin, xmax);
            sprintf(name, "%s, %s, Time profile, Ch.%i; #font[12]{t} / ns; Effective neutrons", FC.c_str(), Setup.c_str(), i+1);
            pH1Tproj[i][0]->SetTitle(name);
            sprintf(name, "%s_%s_ProjT_Sc_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Tproj[i][1] = new TH1F(name, name, NbinsX, xmin, xmax);
            sprintf(name, "%s, %s, Time profile, scattered, Ch.%i; #font[12]{t} / ns; Effective neutrons", FC.c_str(), Setup.c_str(), i+1);
            pH1Tproj[i][1]->SetTitle(name);
            sprintf(name, "%s_%s_ProjEeff_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Eeff[i][0] = new TH1F(name, name, NbinsY, ymin, ymax);
            sprintf(name, "%s, %s, Effective neutron spectrum, Ch.%i; #font[12]{t} / ns; Effective neutrons", FC.c_str(), Setup.c_str(), i+1);
            pH1Eeff[i][0]->SetTitle(name);
            sprintf(name, "%s_%s_ProjEeff_Sc_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Eeff[i][1] = new TH1F(name, name, NbinsY, ymin, ymax);
            sprintf(name, "%s, %s, Effective scattered neutron spectrum, Ch.%i; #font[12]{t} / ns; Effective neutrons", FC.c_str(), Setup.c_str(), i+1);
            pH1Eeff[i][1]->SetTitle(name);

//            cout << "E   w   Total   Scattered" << endl;

            for (Int_t binE = 0; binE <= binEmax; binE++)
            {
                Double_t E = pH1Eproj[i][0]->GetBinCenter(binE);
                Double_t w = 0, Dw = 0;
                relSigma(E, w, Dw);
                pH1Eeff[i][0]->AddBinContent(binE, w * pH1Eproj[i][0]->GetBinContent(binE));
                pH1Eeff[i][1]->AddBinContent(binE, w * pH1Eproj[i][1]->GetBinContent(binE));

//                Double_t Total = 0, Scattered = 0;

                for (Int_t binT = 0; binT <= NbinsY + 1; binT ++)
                {
//                    Total     += pH2TvsE[i][0]->GetBinContent(binT, binE);
//                    Scattered += pH2TvsE[i][1]->GetBinContent(binT, binE);
                    pH1Tproj[i][0]->AddBinContent(binT, w * pH2TvsE[i][0]->GetBinContent(binT, binE));
                    if (tID)
                        pH1Tproj[i][1]->AddBinContent(binT, w * pH2TvsE[i][1]->GetBinContent(binT, binE));
                    else if (binE < binEmin)
                        pH1Tproj[i][1]->AddBinContent(binT, w * pH2TvsE[i][0]->GetBinContent(binT, binE));
                }

//                cout << E << "  " << w << "  " << Total << "  " << Scattered << endl;

            } // for(binE)
//            cout << endl << "=== Total:     " << pH1Tproj[i][0]->Integral() << "\t===" << endl << endl;
//            cout << endl << "=== Scattered: " << pH1Tproj[i][1]->Integral() << "\t===" << endl << endl;

            SaveToFile("Analysis/EffEnergy", pH1Eeff[i][0]);
            SaveToFile("Analysis/EffEnergy/Scattered", pH1Eeff[i][1]);
            SaveToFile("Analysis/EffToF", pH1Tproj[i][0]);
            SaveToFile("Analysis/EffToF/Scattered", pH1Tproj[i][1]);
            plot->DtPeakForm(i, pH1Tproj[i][1], pH1Tproj[i][0], Emin, nProj[i]);
//            plot->Eproj(i, pH1Eeff[i][1], pH1Eeff[i][0], nProj[i]);
        } else {
            cout << "Error loading hists? Try DrawMulti = kTRUE" << endl;
            sprintf(name, "Analysis/EffToF/%s_%s_ProjT_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Tproj[i][0] = (TH1F*) f->Get(name);
            sprintf(name, "Analysis/EffToF/Scattered/%s_%s_ProjT_Sc_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Tproj[i][1] = (TH1F*) f->Get(name);
            sprintf(name, "Analysis/EffEnergy/%s_%s_ProjEeff_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Eeff[i][0] = (TH1F*) f->Get(name);
            sprintf(name, "Analysis/EffEnergy/Scattered/%s_%s_ProjEeff_Sc_%i", FC.c_str(), Setup.c_str(), i+1);
            pH1Eeff[i][1] = (TH1F*) f->Get(name);
        } // if(Draw)
    } // for(i)
    cout << "Done: Projections" << endl;
    DoneProjections = kTRUE;
}


Double_t Sim::ToF(Double_t Ekin, Double_t FlightPath)
{
    Double_t g = Ekin / m0 + 1.0;
    Double_t beta = sqrt(1.0 - 1.0 / (g*g));
    Double_t ToF = FlightPath / beta * 1.E9 / c;
//    cout << "E " << Ekin << ", Path " << FlightPath << ", ToF " << ToF << endl;
    return ToF;
}


void Sim::nProjectiles()
{
    cout << endl << "Get number of projectiles towards deposits..." << endl;
    char name[64] = "";
    cout << " Ch  projectiles" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        sprintf(name, "Source/Source_Theta_Ch.%i", i+1);
        TH1F *pH1Ang = (TH1F*) f->Get(name);
        nProj[i] = pH1Ang->Integral();
        cout << " " << i+1 << "   " << nProj[i] << endl;
    }
    cout << "Done: projectiles" << endl;
    DoneProjectiles = kTRUE;
}


void Sim::DirectN()
{
    cout << endl << "Getting number of direct neutrons..." << endl;
    cout << " Ch  nDirect" << endl;

    for (Int_t i = 0; i < NumCh; i++)
    {
        nDirect[i] = tID ? pH2TvsE[i][0]->Integral() - pH2TvsE[i][1]->Integral()
                      : nN(i, binEmin, binEmax, 0);
        DnDirect[i] = sqrt(nDirect[i]);
        cout << " " << i+1 << "   " << nDirect[i] << "+-" << DnDirect[i] << endl;
        if (DrawMulti)
            plot->TvsE(i, pH2TvsE[i][0], binToFmin[i], binToFmax[i], nDirect[i]);
    }
    cout << "Done: Direct neutrons" << endl;
    DoneDirect = kTRUE;
}


void Sim::DirectEff()
{
    if (!DoneDirect)
        DirectN();
    cout << endl << "Effective direct neutrons..." << endl;
    cout << " Ch  effDirect" << endl;

    Double_t E, w = 0, Dw = 0, N, D2N;
    for (Int_t i = 0; i < NumCh; i++)
    {
        N = 0;
        D2N = 0;
        for (Int_t binE = binEmin; binE <= binEmax; binE++)
        {
            E = pH1Eproj[i][0]->GetBinCenter(binE);
            relSigma(E, w, Dw);
            Double_t n = tID ? pH1Eproj[i][0]->GetBinContent(binE) - pH1Eproj[i][1]->GetBinContent(binE)
                             : pH1Eproj[i][0]->GetBinContent(binE);
            N += w * n;
            D2N += w * w * n;
        } // for(binE)
        effDirect[i] = N;
        DeffDirect[i] = sqrt(D2N);
        cout << " " << i+1 << "   " << effDirect[i] << "+-" << DeffDirect[i] << endl;
    } // for(i)
    cout << "Done: Effective direct neutrons" << endl;
}


Int_t Sim::nN(Int_t ch, Int_t binElow, Int_t binEup, Int_t Sc)
{
    return pH1Eproj[ch][Sc]->Integral(binElow, binEup);
}


Double_t Sim::effN(Int_t ch, Int_t binElow, Int_t binEup, Int_t Sc)
{
    Double_t E, w = 0, Dw = 0, N = 0;
    for (Int_t binE = binElow; binE <= binEup; binE++)
    {
        E = pH1Eproj[ch][Sc]->GetBinCenter(binE);
        relSigma(E, w, Dw);
        N += w * pH1Eproj[ch][Sc]->GetBinContent(binE);
    }
    return N;
}


Double_t Sim::DeffN(Int_t ch, Int_t binElow, Int_t binEup, Int_t Sc)
{
    Double_t E, w = 0, Dw = 0, D2N = 0;
    for (Int_t binE = binElow; binE <= binEup; binE++)
    {
        E = pH1Eproj[ch][Sc]->GetBinCenter(binE);
        relSigma(E, w, Dw);
        D2N += w * w * pH1Eproj[ch][Sc]->GetBinContent(binE);
    }
    return sqrt(D2N);
}


void Sim::ScatN()
{
    cout << endl << "Getting number of scattered neutrons..." << endl;
    cout << " Ch  nScat" << endl;

    for (Int_t i = 0; i < NumCh; i++)
    {
        nScat[i] = tID ? pH1Eproj[i][1]->Integral()
                       : nN(i, 0, binEmin - 1, 0);
        DnScat[i] = sqrt(nScat[i]);
        cout << " " << i+1 << "   " << nScat[i] << "+-" << DnScat[i] << endl;
    }
    cout << "Done: Scattered neutrons" << endl;
}


void Sim::ScatEff()
{
    cout << endl << "Getting effective number of scattered neutrons..." << endl;
    cout << " Ch  effScat" << endl;

    for (Int_t i = 0; i < NumCh; i++)
    {
        effScat[i] = tID ? effN(i, 0, binEmax, 1)
                         : effN(i, 0, binEmin - 1, 0);
        DeffScat[i] = tID ? DeffN(i, 0, binEmax, 1)
                          : DeffN(i, 0, binEmin - 1, 0);
        cout << " " << i+1 << "   " << effScat[i] << "+-" << DeffScat[i] << endl;
    } // for(i)
    cout << "Done: Effective scattered neutrons" << endl;
}


void Sim::SigmaEff()
{
    cout << endl << "Efficient cross sections..." << endl;
    cout << " Ch  Direct/b   Scattered/b" << endl;
    Double_t sigma15 = gPu242->Eval(15.0);
    for (Int_t i = 0; i < NumCh; i++)
    {
        CsDir[i] = sigma15 * effDirect[i] / nDirect[i];
        CsSc[i]  = sigma15 * effScat[i] / nScat[i];
        cout << " " << i+1 << "   " << CsDir[i] << "   " << CsSc[i] << endl;
    }
    cout << "Done: efficient cross sections" << endl;
}


void Sim::Calculate()
{
    Projections();
    DirectN();
    DirectEff();
    ScatN();
    ScatEff();
    SigmaEff();
    DoneCalc = kTRUE;
}


void Sim::PrintEffN()
{
    if (!DoneProjections)
        Projections();
    if (!DoneSigma)
        cout << "No cross section data!" << endl;
    cout << "E   w   {dir   sc}" << endl;
    for (Int_t binE = 1; binE <= pH1Eproj[0][0]->GetNbinsX(); binE++)
    {
        Double_t E = pH1Eproj[0][0]->GetBinCenter(binE), w, Dw;
        relSigma(E, w, Dw);
        cout << E << " " << w << " ";
        for (Int_t i = 0; i < NumCh; i++)
        {
            if (tID)
            {
                cout << pH1Eproj[i][0]->GetBinContent(binE) - pH1Eproj[i][1]->GetBinContent(binE) << " "
                     << pH1Eproj[i][1]->GetBinContent(binE) << " ";
            } else {
                cout << (binE < binEmin ? 0 : pH1Eproj[i][0]->GetBinContent(binE)) << " "
                     << (binE < binEmin ? pH1Eproj[i][0]->GetBinContent(binE) : 0) << " ";
            }
        }
        cout << endl;
    }
}


void Sim::SimToF()
{
    Projections();
    char name[64] = "";
    char title[128] = "";

    /// simulation properties
    Double_t AccPulseLength = 7.0; // ns
    Double_t TimeResolution = 2.5; // ns
    Int_t NbinsProj = pH1Tproj[0][0]->GetNbinsX();
    Double_t ProjBinWidth = pH1Tproj[0][1]->GetXaxis()->GetBinWidth(1);
    Int_t AccPulseBins = AccPulseLength / ProjBinWidth;
    cout << endl << "Simulating ToF spectra..." << endl
         << " Acc. pulse length: " << AccPulseLength << " ns" << endl
         << " Time resolution: " << TimeResolution << " ns" << endl
         << " Original bins: " << NbinsProj << endl
         << " Original bin width: " << ProjBinWidth << endl
         << " Acc. pulse bins: " << AccPulseBins << endl;

    /// experimental properties - depending on FC
    TFile* fExp; // experimental data analysis file for right time offset
    Double_t t0 = 42;
    Double_t tm[NumCh];
    Double_t t3 = 402;
    if (PuFC)
    {
        if (!strcmp(Setup.c_str(), "FG") || !strcmp(Setup.c_str(), "Open"))
        {
            fExp = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/NIF.root");
            cout << "PuFC, Open" << endl;
        }
        else
            fExp = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/SB.root");
        Double_t m[] = {139.7705, 128.2225, 125.879, 129.663, 129.5655, 128.589, 127.466, 127.49};
        for (Int_t i = 0; i < NumCh; i++)
            tm[i] = m[i];
    } else {
        if (!strcmp(Setup.c_str(), "FG") || !strcmp(Setup.c_str(), "Open"))
            fExp = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_NIF.root");
        else
            fExp = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/UFC_SB.root");
        Double_t m[] = {273.4865, 268.945, 260.5955, 265.2345, 265.21, 264.3555, 263.501, 263.086};
        for (Int_t i = 0; i < NumCh; i++)
            tm[i] = m[i];
    }
    TH1I* pH1Exp[NumCh];
    Double_t scale[NumCh];
    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "Analysis/ToF/H1AnaHZDRDtG_%i", i+1);
        pH1Exp[i] = (TH1I*) fExp->Get(name);
        cout << pH1Exp[i]->GetName() << endl;
        scale[i] = pH1Exp[i]->Integral(tm[i] - 25, tm[i] + 25) / pH1Tproj[i][0]->Integral();
//        cout << pH1Exp[i]->Integral(tm[i] - 25, tm[i] + 25) << "  " << pH1Tproj[i][0]->Integral() << "  " << scale[i] << endl;
    }
    Double_t minToF = pH1Exp[0]->GetXaxis()->GetBinLowEdge(1);
    Double_t maxToF = pH1Exp[0]->GetXaxis()->GetBinLowEdge(pH1Exp[0]->GetNbinsX() + 1);
    Int_t NbinsToF = (maxToF - minToF) / ProjBinWidth; // use Tproj's bin width
//    Int_t NbinsToF = pH1Exp[0]->GetNbinsX(); // use exp. data's bin width
//    cout << " Sim. bins: " << NbinsToF << endl
//         << " Sim. ToF range: " << minToF << "-" << maxToF << " ns" << endl;

    /// Prepare simulated ToF histograms
    for (Int_t i = 0; i < NumCh; i++)
    {
        sprintf(name, "%s_Dt_Ch.%i", FC.c_str(), i+1);
        sprintf(title, "%s, Ch. %i, simulated ToF spectrum; #font[12]{t} [ns]; Counts", FC.c_str(), i+1);
        pH1ToF[i][0] = new TH1F(name, title,  NbinsToF, minToF, maxToF);
    }
    cout << " Sim. bins: " << pH1ToF[0][0]->GetNbinsX() << endl
         << " Sim. ToF range: " << pH1ToF[0][0]->GetBinLowEdge(1) << "-" << pH1ToF[0][0]->GetBinLowEdge(pH1ToF[0][0]->GetNbinsX()+1) << " ns" << endl;

    /// create folding function
    Double_t ampl = 1.0 / (Double_t)AccPulseBins / sqrt(2 * TMath::Pi()) * pH1ToF[0][0]->GetBinWidth(1) / TimeResolution;
    cout << " Ampl: " << ampl << endl;
    TF1* g = new TF1("fG", "gaus", -440, 440);
    g->SetParameters(ampl, 0, TimeResolution);

    /// fold
    for (Int_t i = 0; i < 1/*NumCh*/; i++)
    {
        Int_t ProjMaxBin = pH1Tproj[i][0]->GetMaximumBin();
        Double_t ProjMax = pH1Tproj[i][0]->GetBinCenter(ProjMaxBin);
//        cout << ProjMax << " -> " << tm[i] << endl;
        for (Int_t j = 1; j <= NbinsToF; j++)
        {
            Double_t t0 = pH1ToF[i][0]->GetBinCenter(j); // ToF sampling point
//            cout << t0 << endl;
            Double_t sum = 0;
            for (Int_t k = 1; k < NbinsProj - AccPulseBins; k++)
            {
                Double_t t1 = 0.5 * (pH1Tproj[i][0]->GetBinCenter(k) + pH1Tproj[i][0]->GetBinCenter(k + AccPulseBins));
//                cout << " " << t1;
                sum += g->Eval(t1 - t0 + tm[i] - ProjMax) * pH1Tproj[i][0]->Integral(k, k + AccPulseBins - 1);
            }
//            cout << j << "  " << t0 << "  " << sum << endl;
            pH1ToF[i][0]->SetBinContent(j, scale[i] * sum);
        }
//        cout << NbinsToF << endl;
//        SaveToFile("Analysis/TimeDiff", pH1ToF[i][0]);
    }

    /// Draw
    cout << pH1Exp[0]->GetName() << ", " << pH1ToF[0][0]->GetName() << endl;
    TLegend *l = new TLegend(0.6, 0.6, 0.8, 0.8);
    l->AddEntry(pH1Exp[0], "Experiment");
    l->AddEntry(pH1ToF[0][0], "Simulation");
    TCanvas *c1 = new TCanvas("cSimToF", "Simulated ToF spectrum", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    pH1Exp[0]->SetStats(0);
    sprintf(title, "%s, %s, #Delta#font[12]{t} ch.%i; #Delta#font[12]{t} [ns]; counts/ns", FC.c_str(), Setup.c_str(), 1);
    pH1Exp[0]->SetTitle(title);
    pH1Exp[0]->GetXaxis()->SetLabelSize(0.06);
    pH1Exp[0]->GetXaxis()->SetTitleSize(0.07);
    pH1Exp[0]->GetXaxis()->SetTitleOffset(0.7);
    pH1Exp[0]->GetYaxis()->SetLabelSize(0.06);
    pH1Exp[0]->GetYaxis()->SetTitleSize(0.07);
    pH1Exp[0]->GetYaxis()->SetTitleOffset(0.6);
    pH1Exp[0]->Draw("hist");
    pH1ToF[0][0]->Scale(10);
    pH1ToF[0][0]->SetLineColor(kRed);
    pH1ToF[0][0]->Draw("same hist");
    l->Draw();
    c1->Modified();
    c1->Update();
    c1->Draw();
    cout << pH1Tproj[0][0]->Integral() << "  " << pH1ToF[0][0]->Integral() << endl;
}


void Sim::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    cout << "Saving " << pObj->GetName() << " to " << path << endl;
    TFile* f = TFile::Open(FileName.c_str(), "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder "Analysis" already exists, otherwise create it
    if (f->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) f->Get(path.c_str());
    else
    {   cout << "Creating " << path << endl;
        f->mkdir(path.c_str(), "Folder containing offline Analysis objects");
        EvalDir = f->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    f->Save(); f->Close();
}
