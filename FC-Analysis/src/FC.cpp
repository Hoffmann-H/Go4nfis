//
#include "FC.h"

using namespace std;

FC::FC()//Double_t nFieldFG, Double_t DnFieldFG, Double_t nFieldBG, Double_t DnFieldBG)
{

}


FC::~FC()
{

}


void FC::InitVar(Bool_t draw)
{
    cout << endl << "Initializing common variables..." << endl;
    DrawSingle = draw;
    DrawMulti = kFALSE;
    if (DrawSingle || DrawMulti)
        plot = new Plot(Name, "Open");
    FgRuns = 0;
    BgRuns = 0;
    DoneQDC = kFALSE;
    DoneThresholds = kFALSE;
    DoneLimits = kFALSE;
    DoneNatoms = kFALSE;
    DoneDtBG = kFALSE;
    DoneDt = kFALSE;
    DoneScatCorr = kFALSE;
    DoneSimFg = kFALSE;
    DoneIso = kFALSE;
    DoneRawCS = kFALSE;
    DoneCorrections = kFALSE;
    DoneTransmission = kFALSE;

    u = 1.660539E-24; // [g]
    Area = 4300;
    DArea = 120;
    cout << "Done: common variables" << endl;
}


void FC::UseHists(Int_t start, Int_t stop, string setup, Int_t run)
{
    char name[16] = "";
    char file[128] = "";
    for (Int_t num = start; num <= stop; num++)
    {
        // Add zeros
        if (num < 10)
            sprintf(name, "lmd/00%i", num);
        else if (num < 100)
            sprintf(name, "lmd/0%i", num);
        else
            sprintf(name, "lmd/%i", num);
        // Define file name
        sprintf(file, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", name);
        // Create Hist instance
        pH[nHist] = new Hist(file, setup, Name, name, run);
        nHist++;
    }
}


void FC::UseHist(string file_name, string setup, Int_t run)
{
    char file[128] = "";
    sprintf(file, "/home/hoffma93/Programme/Go4nfis/offline/results/%s.root", file_name.c_str());
    pH[nHist] = new Hist(file, setup, Name, file_name, run);
    nHist++;
}


void FC::RegisterHists()
{
    cout << endl << "Conclusion: Runs and Hists" << endl;
    Int_t run = -1;
    tFG = 0;
    tBG = 0;
    for (Int_t j = 0; j < nHist; j++)
    {
        if (pH[j]->RunNr != run)
            cout << " Run " << pR[++run]->Name << endl;
        pR[run]->GetHist(pH[j]);
        cout << "  Hist " << pH[j]->Name << endl;
        if (pH[j]->IsForeground())
            tFG += pH[j]->t_live;
        else
            tBG += pH[j]->t_live;
    }
    if (DrawSingle)
        Stability();
}


void FC::Stability()
{
    char name[64] = "";
    sprintf(name, "%s (n,f)-Rate; #font[12]{t}; Rate [1/s]", Name.c_str());
    TMultiGraph *mg1 = new TMultiGraph("mgStabNIF", name);
    sprintf(name, "%s Kanal", Name.c_str());
    TLegend *l1 = new TLegend(0.9, 0.5, 1.0, 1.0, name);
    sprintf(name, "%s Spontanspaltrate; #font[12]{t}; Rate [1/s]", Name.c_str());
    TMultiGraph *mg2 = new TMultiGraph("mgStabSF", name);
    sprintf(name, "%s Kanal", Name.c_str());
    TLegend *l2 = new TLegend(0.9, 0.5, 1.0, 1.0, name);
    TGraphErrors *g1[NumCh];
    TGraphErrors *g2[NumCh];
    Double_t X[nHist];
    Double_t Xerr[nHist];
    for (Int_t j = 0; j < nHist; j++)
    {
        X[j] = 0.5*((Double_t)pH[j]->t_start + (Double_t)pH[j]->t_stop)/* - 1.4011E+09*/; // TODO: Insert time here
        Xerr[j] = 0.5*((Double_t)pH[j]->t_stop - (Double_t)pH[j]->t_start);
        if (CommentFlag)
            cout << pH[j]->t_start << ", " << pH[j]->t_stop << " -> " << X[j] << "+-" << Xerr[j] << endl;
    }
//    if (CommentFlag)
//        cout << " Ch   Hist   rNIF   rSF" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        Double_t rNIF[nHist];
        Double_t DrNIF[nHist];
        Double_t rSF[nHist];
        Double_t DrSF[nHist];
        for (Int_t j = 0; j < nHist; j++)
        {
            rNIF[j] = pH[j]->nNIF[i] / pH[j]->t_live;
            DrNIF[j] = pH[j]->DnNIF[i] / pH[j]->t_live;
            rSF[j] = pH[j]->nSF[i] / pH[j]->t_live;
            DrSF[j] = pH[j]->DnSF[i] / pH[j]->t_live;
//            if (CommentFlag)
//                cout << " " << i+1 << "   " << pH[j]->Name << "   " << rNIF[j] << "+-" << DrNIF[j] << "   " << rSF[j] << "+-" << DrSF[j] << endl;
        }
        g1[i] = new TGraphErrors(nHist, X, rNIF, Xerr, DrNIF);
        sprintf(name, "gStabNIF_%i", i+1);
        g1[i]->SetName(name);
        g1[i]->SetLineColor(i);
        g1[i]->SetLineWidth(2);
        mg1->Add(g1[i]);
        sprintf(name, "%i", i+1);
        l1->AddEntry(g1[i], name);
        g2[i] = new TGraphErrors(nHist, X, rSF, Xerr, DrSF);
        sprintf(name, "gStabSF_%i", i+1);
        g2[i]->SetName(name);
        g2[i]->SetLineColor(i);
        g2[i]->SetLineWidth(2);
        mg2->Add(g2[i]);
        sprintf(name, "%i", i+1);
        l2->AddEntry(g2[i], name);
    }
    TCanvas *c1 = new TCanvas("cStabNIF", "Stability", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg1->Draw("AP");
    mg1->GetXaxis()->SetTimeDisplay(1);
    mg1->GetXaxis()->SetTimeFormat("%d.%m. %H:%M%F1970-01-01 00:00:00s0");
    l1->Draw();
    c1->Modified();
    c1->Update();
    TCanvas *c2 = new TCanvas("cStabSF", "Stability", 200, 10, 700, 500);
    gPad->SetTicks(1, 1);
    mg2->Draw("AP");
    mg2->GetXaxis()->SetTimeDisplay(1);
    mg2->GetXaxis()->SetTimeFormat("%d.%m. %H:%M%F1970-01-01 00:00:00s0");
    l2->Draw();
    c2->Modified();
    c2->Update();
}


//void FC::GetHistPointers()
//{
//    Int_t start = 0;
//    for (Int_t k = 0; k < FgRuns + BgRuns; k++)
//    {
//        for (Int_t j = 0; j < pR[k]->nHist; j++)
//            pH[start + j] = pR[k]->pH[j];
//        start += pR[k]->nHist;
//    }
//}


//void FC::AnalyzeDt()
//{   // Extract (n,f) counts from Dt spectra
//    // and sum up neutron flux
//    if (!DoneLimits)
//        GetLimits();
//    if (!DoneDtBG)
//        this->AnalyzeDtBG(); // Estimate background
//    cout << endl << "Analyzing ToF peaks..." << endl;

//    // Monitor rates
//    Int_t n = 0;
//    Int_t div[FgRuns + BgRuns];
//    Double_t mon[FgRuns + BgRuns];
//    for (Int_t j = 0; j < FgRuns + BgRuns; j++)
//    {
//        n += pR[j]->nHist;
//        div[j] = n;
//        mon[j] = pR[j]->GetMonitor();
//    }

//    for (int i = 0; i < NumCh; i++)
//    {
//        for (Int_t j = 0; j < FgRuns + BgRuns; j++)
//            pR[j]->AnalyzeDt(i, lim[0][i], lim[1][i], lim[2][i], lim[3][i]);
//        for (Int_t j = 0; j < nHist; j++)
//        {
//            Double_t rNIF[nHist];
//            Double_t DrNIF[nHist];
//            Double_t rSF[nHist];
//            Double_t DrSF[nHist];
//            rNIF[j] = pH[j]->nNIF[i] / pH[j]->t_live;
//            DrNIF[j] = pH[j]->DnNIF[i] / pH[j]->t_live;
//            rSF[j] = pH[j]->nSF[i] / pH[j]->t_live;
//            DrSF[j] = pH[j]->DnSF[i] / pH[j]->t_live;
//            if (DrawSingle)
//                plot->Stability(i, nHist, rNIF, DrNIF);
//        }
//    }
//    cout << "Ch   nFG   tFG   nBG   tBG" << endl;
//    for (int i = 0; i < NumCh; i++)
//        cout << " " << i+1 << "   " << nFg1[i] << "+-" << DnFg1[i] << "   " << tFG <<
//                              "   " << nBg1[i] << "+-" << DnBg1[i] << "   " << tBG << endl;
//    DoneDt = kTRUE;
//    cout << "Done: ToF" << endl;
//}


//void FC::ScatCorr(Int_t i)
//{
//    pDirect[i][0] = 1 - pHFG->NeutronFlux[i] / pHBG->NeutronFlux[i] * nBg1[i] / nFg1[i];
//    DpDirect[i][0] = sqrt( pow(pHFG->DNeutronFlux[i] / pHBG->NeutronFlux[i] * nBg1[i] / nFg1[i], 2) +
//                           pow(pHFG->NeutronFlux[i] * pHBG->DNeutronFlux[i] / pow(pHBG->NeutronFlux[i], 2) * nBg1[i] / nFg1[i], 2) +
//                           pow(pHFG->NeutronFlux[i] / pHBG->NeutronFlux[i] * DnBg1[i] / nFg1[i], 2) +
//                           pow(pHFG->NeutronFlux[i] / pHBG->NeutronFlux[i] * nBg1[i] * DnFg1[i] / pow(nFg1[i], 2), 2) );
//}


//void FC::ScatCorrDiff()
//{
//    if (!DoneDt)
//        AnalyzeDt();
//    if (!DoneRawCS)
//        CrossSection();
//    cout << endl << "Experimental Scattering correction: Shadow Cone" << endl;
//    cout << "Ch   pDirect" << endl;
//    for (int i = 0; i < NumCh; i++)
//    {
////        ScatCorr(i);
//        cout << " " << i+1 << "   " << pDirect[i][0] << " +- " << DpDirect[i][0] << endl;
//    }
//    DoneScatCorr = kTRUE;
//    cout << "Done: Scattering correction" << endl;
//}


//void FC::ScatCorrFit()
//{
//    if (!DoneLimits)
//        GetLimits();
//    if (!DoneDtBG)
//        AnalyzeDtBG();
//    cout << endl << "Experimental Scattering correction: Fit" << endl;

//    char name[32] = "";
//    Double_t pi = 3.1415926535898;
//    Double_t ChPerBin = pHFg1[0]->pHDtg1[0]->GetBinWidth(1);
//    Double_t ampl, Dampl, center, width, Dwidth; // center and width in units of QDC channels
//    Double_t low = pHFg1[0]->pHDtg1[0]->GetBinCenter(lim[0][0]);
//    Double_t up = pHFg1[0]->pHDtg1[0]->GetBinCenter(lim[3][0]);
//    for (int i = 0; i < NumCh; i++)
//    {
//        width = ChPerBin * (lim[2][i]-lim[1][i]) / 6;
//        center = 0.5 * ( pHFg1[0]->pHDtg1[0]->GetBinCenter(lim[1][i]) +
//                         pHFg1[0]->pHDtg1[0]->GetBinCenter(lim[2][i]) );
//        sprintf(name, "f%sNIFDtPeak_%i", Name.c_str(), i+1);
//        TF1* fFG = new TF1(name, func_peak, low, up, 4);
//        fFG->SetNpx(1000);
//        fFG->SetRange(low, up);
//        fFG->SetParameters(100, center, width, pHFG->t_live*avBg1[i]);
//        pHFG->pHDtg1[i]->Fit(name, "LQR");
//        ampl = (Double_t)fFG->GetParameter(0);
//        Dampl = (Double_t)fFG->GetParError(0);
//        width = (Double_t)fFG->GetParameter(2);
//        Dwidth = (Double_t)fFG->GetParError(2);
//        nFg1[i] = sqrt(pi) * ampl * width / ChPerBin / pHFG->t_real;
//        DnFg1[i] = sqrt(pi) / ChPerBin / pHFG->t_real * sqrt( pow(Dampl*width, 2) + pow(ampl*Dwidth, 2) );
////        cout << "  FG " << fFG->GetParameter(0) << "  " << fFG->GetParameter(1) << " " << fFG->GetParameter(2) << " " << fFG->GetParameter(3) << " " << nFg1[i] << "+-" << DnFg1[i] << endl;

//        sprintf(name, "f%sSBDtPeak_%i", Name.c_str(), i+1);
//        TF1* fBG = new TF1(name, func_peak, low, up, 4);
//        fBG->SetNpx(1000);
//        fBG->SetRange(low, up);
//        fBG->SetParameters(100, center, width, pHBG->t_live*avBg1[i]);
//        fBG->FixParameter(2, width);
//        pHBG->pHDtg1[i]->Fit(name, "LQR");
//        ampl = (Double_t)fBG->GetParameter(0);
//        Dampl = (Double_t)fBG->GetParError(0);
//        nBg1[i] = sqrt(pi) * ampl * width / ChPerBin / pHBG->t_real;
//        DnBg1[i] = sqrt(pi) / ChPerBin  / pHBG->t_real* sqrt( pow(Dampl*width, 2) + pow(ampl*Dwidth, 2) );
////        cout << "  BG " << fBG->GetParameter(0) << "  " << fBG->GetParameter(1) << " " << fBG->GetParameter(2) << " " << fBG->GetParameter(3) << " " << nBg1[i] << "+-" << DnBg1[i] << endl;

//        ScatCorr(i);
//    }

//    DoneScatCorr = kTRUE;
//    cout << "Done: Scattering correction" << endl;
//}


void FC::ScatCorrSim()
{
//    if (!DoneDt)
//        AnalyzeDt();
    if (!DoneSimFg)
        GetSimFg();
    cout << endl << "Scattering correction: Simulaion" << endl;

    DoneScatCorr = kTRUE;
    cout << "Done: Scattering correction" << endl;
}


void FC::CrossSection()
{// Calculate cross section from uncorrected numbers.
 // Need: fluence of incident neutrons
 //       efficient number of target atoms
 //       number of (n,f) events
//    if (!DoneDt)
//        AnalyzeDt();
    if (!DoneNatoms)
        GetNatoms();
    cout << endl << "Uncorrected cross section..." << endl;
    Double_t avCS = 0; // average cross section
    Double_t D2avCS = 0;
    Double_t avFG;
    Double_t D2avFG;
    Double_t avBG;
    Double_t D2avBG;
    for (int i = 0; i < NumCh; i++)
    {
        cout << " Chanel " << i+1 << endl;
        if (CommentFlag)
            cout << "  " << nAtoms[i] << "+-" << DnAtoms[i] << " atoms" << endl
                 << "  Run   t_live[s]   (n,f)-events   n-Flux[mm^-2 s^-1]   sigma[b]" << endl;
        avFG = 0;
        D2avFG = 0;
        avBG = 0;
        D2avBG = 0;

        // Cross section ////////////////////////////////////////////////////////////
        //  weighting with live times. Exact would be weighting with monitor counts
        for (Int_t k = 0; k < FgRuns + BgRuns; k++)
        {
            pR[k]->CrossSection(i);
            if (CommentFlag)
                cout << "  " << pR[k]->Name << "   " << pR[k]->t_live << "   " << pR[k]->nNIF[i] << "+-" << pR[k]->DnNIF[i] << "   " << pR[k]->NeutronFlux[i] << "+-" << pR[k]->DstatNeutronFlux[i] << "   " << pR[k]->uncCS[i] << "+-" << pR[k]->DuncCS[i] << endl;
            if (pR[k]->IsForeground())
            {
                avFG += pR[k]->t_live / tFG * pR[k]->uncCS[i];
                D2avFG += pow(pR[k]->t_live / tFG * pR[k]->DuncCS[i], 2);
//                cout << "   " << pR[k]->DuncCS[i] << "  " << pR[k]->t_live / tFG * pR[k]->DuncCS[i] << "  " << D2avFG << endl;
            }
            else
            {
                avBG += pR[k]->t_live / tBG * pR[k]->uncCS[i];
                D2avBG += pow(pR[k]->t_live / tBG * pR[k]->DuncCS[i], 2);
            }
        }
        avCS += avFG / NumCh;
        D2avCS += D2avFG / NumCh / NumCh;
        pDirect[i][0] = 1.0 - avBG / avFG;
        DpDirect[i][0] = sqrt(D2avBG) / avFG;
//        if (CommentFlag)
//            cout << "   nFlux " << pHFG->NeutronFlux[i] << "+-" << pHFG->DNeutronFlux[i] << "mm^-2 s^-1" << endl <<
//                    "   nAtoms " << nAtoms[i] << "+-" << DnAtoms[i] << endl <<
//                    "   Foreground " << nFg1[i] << "+-" << DnFg1[i]  << " s^-1" << endl;
        uCS[i] = avFG;
        DuCS[i] = sqrt(D2avFG);
        cout << "  raw CS = " << uCS[i] << "+-" << DuCS[i] << " barn" << endl;
    }
    cout << "Average: sigma = " << avCS << "+-" << sqrt(D2avCS) << endl;
    DoneRawCS = kTRUE;
    cout << "Done: raw cross section" << endl;
}


void FC::GetSimFg()
{
    cout << endl << "Simulation results requested. Starting " << Name << " foreground simulation analysis..." << endl;
    sim = new AnaSim(Name, 1, plot);
    sim->Corrections();
    for (Int_t i = 0; i < NumCh; i++)
    {
        fTS[i] = sim->F[i];
        DfTS[i] = sim->DF[i];
        pDirect[i][1] = sim->S[i];
        DpDirect[i][1] = sim->DS[i];
        SimT[i] = sim->T[i];
        DSimT[i] = sim->DT[i];
    }
    DoneSimFg = kTRUE;
    cout << "Done: Got foreground simulation" << endl;
}


//void FC::GetSimBg()
//{
//    if (!DoneSimFg)
//        GetSimFg();
//    cout << endl << "Simulation results requested. Starting " << Name << " background simulation analysis..." << endl;
//    sim->ShadowCone();
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        pDirect[i][2] = sim->SC[i];
//        DpDirect[i][2] = sim->DSC[i];
//    }
//    DoneSimBg = kTRUE;
//}


void FC::Corrections()
{
    if (!DoneRawCS)
        CrossSection();
    if (!DoneScatCorr)
        ScatCorrSim();
    if (!DoneIso)
        IsoVec();
    cout << endl << "Corrections..." << endl;
//    cout << "Ch   F   Isotopes" << endl;
//    for (int i = 0; i < NumCh; i++)
//    {
//        cout << i+1 << "   ";
//        cout << 1 / Transm[i] << "   ";
//        cout << InScat[i] << "   ";
//        cout << fIsoVec - sIsoVec / InScat[i] * Transm[i] / uCS[i] << "" << endl;
//    }
    if (strcmp(Name.c_str(), "UFC"))
        cout << "eval CS(14.97 MeV) = " << sim->Fg->gPu242->Eval(14.97) << "+-" << sim->Fg->DgPu242->Eval(14.97) << endl;
    else
        cout << "eval CS(14.97 MeV) = " << sim->Fg->gU235->Eval(14.97) << "+-" << sim->Fg->DgU235->Eval(14.97) << endl;
    cout << "Ch   uncorrected CS   corrected CS" << endl;
    Double_t sUCS = 0, sCS = 0, D2UCS = 0, D2CS = 0;
    for (int i = 0; i < NumCh; i++)
    {
        CS[i] = fTS[i] * uCS[i] * fIsoVec - sIsoVec;
        Double_t uncert[] = {DfIsoVec * fTS[i] * uCS[i],
                             fIsoVec * DfTS[i] * uCS[i],
                             fIsoVec * fTS[i] * DuCS[i],
                             DsIsoVec
                            };
        DCS[i] = sqrt( pow(uncert[0], 2) +
                       pow(uncert[1], 2) +
                       pow(uncert[2], 2) +
                       pow(uncert[3], 2) );
        cout << " " << i+1 << ",  " << uCS[i] << "+-" << DuCS[i] << " barn,  " << CS[i] << "+-" << DCS[i] << " barn" << endl;
//        cout << "  unc: f " << uncert[0] << ", TS " << uncert[1] << ", u " << uncert[2] << ", s " << uncert[3] << endl;
        sUCS += uCS[i];
        D2UCS += pow(DuCS[i], 2);
        sCS += CS[i];
        D2CS += pow(DCS[i], 2);
    }
    cout << "Average: " << sUCS / 8 << "+-" << sqrt(D2UCS) / 8 << "    " << sCS / 8 << "+-" << sqrt(D2CS) / 8 << endl;
    DoneCorrections = kTRUE;
    cout << "Done: Corrections" << endl;
//    if (!DrawSingle)
//        return;
//    plot->Result(uCS, DuCS, CS, DCS);
}


//void FC::CompareShadowCone()
//{
//    if (!DoneScatCorr)
//        ScatCorrDiff();
//    if (!DoneSimBg)
//        GetSimBg();

//    cout << endl << "Comparing methods for In-scattering contribution..." << endl;
//    cout << "Foreground" << endl <<
//            " " << sim->Fg->nProj[0] << " Projectiles" << endl <<
//            " " << sim->Fg->nDirect[0] << " direct neutrons" << endl <<
//            " " << sim->Fg->effDirect[0] << " effective direct neutrons" << endl <<
//            " " << sim->Fg->effScat[0] << " effective scattered neutrons" << endl;
//    cout << "Background" << endl <<
//            " " << sim->Bg->nProj[0] << " Projectiles" << endl <<
//            " " << sim->Bg->nDirect[0] << " direct neutrons" << endl <<
//            " " << sim->Bg->effDirect[0] << " effective direct neutrons" << endl <<
//            " " << sim->Bg->effDirect[0] + sim->Bg->effScat[0] << " effective neutrons" << endl;
//    cout << " Ch   Exp   Sim(1)   Sim(2)" << endl;
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        cout << " " << i+1 << "  " << pDirect[i][0] << " +- " << DpDirect[i][0]
//                           << "  " << pDirect[i][1] << " +- " << DpDirect[i][1]
//                           << "  " << pDirect[i][2] << " +- " << DpDirect[i][2] << endl;
//    }

//    if (!DrawSingle)
//        return;
//    Double_t pExp[NumCh];
//    Double_t DpExp[NumCh];
//    Double_t pSim1[NumCh];
//    Double_t DpSim1[NumCh];
//    Double_t pSim2[NumCh];
//    Double_t DpSim2[NumCh];
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        pExp[i] = pDirect[i][0];
//        DpExp[i] = DpDirect[i][0];
//        pSim1[i] = pDirect[i][1];
//        DpSim1[i] = DpDirect[i][1];
//        pSim2[i] = pDirect[i][2];
//        DpSim2[i] = DpDirect[i][2];
//    }

//    plot->CompSc(&(pExp[0]), &(DpExp[0]), &(pSim1[0]), &(DpSim1[0]), &(pSim2[0]), &(DpSim2[0]));
//}


//void FC::CompareTransmission()
//{
//    if (!DrawSingle)
//        return;
//    if (!DoneSimFg)
//        GetSimFg();
//    if (!DoneTransmission)
//        ExpTrans();
//    plot->CompT(ExpT, DExpT, SimT, DSimT);
//}


//void FC::ExpTrans()
//{
//    if (!DoneDt)
//        AnalyzeDt();
//    if (!DoneDtBG)
//        AnalyzeDtBG();
//    GetExpT();
//    if (!DoneSimFg)
//        GetSimFg();
//    cout << endl << "Experimental transmission" << endl;

//    cout << " Ch   expT   simT   anaT" << endl;
//    for (int i = 0; i < NumCh; i++)
//    {
//        ExpT[i] = uT[i] * pow(sd[i] / sd[7], 2);
//        DExpT[i] = DuT[i] * pow(sd[i] / sd[7], 2);
//        if (CommentFlag)
//            cout << " " << i+1 << "   " << ExpT[i] << "+-" << DExpT[i] << "   " << SimT[i] << "+-" << DSimT[i] << endl;
//    }
//    //
//    cout << "Done: Experimental transmission" << endl;
//    DoneTransmission = kTRUE;
//    if (!DrawSingle)
//        return;
//    plot->ExpT(ExpT, DExpT, SimT, DSimT);
//}


//void FC::SaveToFile(string path, TObject *pObj)
//{   //saves a TObject into the selected file fname
//    TFile* file = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/Evaluation.root", "UPDATE");
//    TDirectory *EvalDir;
//    TObject *pGraph;
//    pGraph = (TObject*) pObj;
//    string GraphName = pGraph->GetName();
//    //check if folder path already exists, otherwise create it
//    if (file->Get(path.c_str())!=0)
//        EvalDir = (TDirectory*) file->Get(path.c_str());
//    else
//    {   file->mkdir(path.c_str(), "Folder containing offline Analysis objects");
//        EvalDir = file->GetDirectory(path.c_str());
//    }
//    EvalDir->cd();
//    if (EvalDir->Get(GraphName.c_str())!=0)
//        EvalDir->Delete((GraphName+";*").c_str());
//    pGraph->Clone()->Write();
//    file->Save(); file->Close();
//}


//Double_t FC::func_peak(Double_t *x, Double_t *par)
//{
//    return par[0] * exp( - pow((x[0] - par[1]) / par[2], 2)) + par[3];
//}
