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
    DtBgLow = 63600; // QDC channels
    DtBgUp = 78500;
    Area = 4300;
    DArea = 120;
    cout << "Done: common variables" << endl;
}


void FC::GetLimits(Double_t n)
{   // Set integration limits to n sigma
    cout << endl << "Opening peak integration limits..." << endl;

    Double_t ChPerBin = pHFG[0]->pHDtG[0]->GetBinWidth(1);
    Double_t BinOffset = pHFG[0]->pHDtG[0]->GetBinCenter(1);
    if (CommentFlag)
        cout << " ChPerBin: " << ChPerBin << ", BinOffset: " << BinOffset << endl
             << "Ch   Peak   bins" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        Double_t ToF_low = pHFG[0]->GetPeakLow(i);
        Double_t ToF_up = pHFG[0]->GetPeakUp(i); // Limits for background are the same
        Double_t center = 0.5 * (ToF_low + ToF_up);
        Double_t width = (ToF_up - ToF_low) / 6.0;
        DtPeakLow[i] = center - 500;//n * width;
        DtPeakUp[i] = center + 1500;//n * width;

        lim[0][i] = (DtBgLow - BinOffset) / ChPerBin;
        lim[1][i] = (DtPeakLow[i] - BinOffset) / ChPerBin;
        lim[2][i] = (DtPeakUp[i] - BinOffset) / ChPerBin;
        lim[3][i] = (DtBgUp - BinOffset) / ChPerBin;
        if (CommentFlag)
            cout << " " << i+1 << "   " << DtPeakLow[i] << " " << DtPeakUp[i] << "   " << lim[0][i] << " " << lim[1][i] << " " << lim[2][i] << " " << lim[3][i] << " " << endl;
    }
    DoneLimits = kTRUE;
    cout << "Done: peak integration limits" << endl;
}


void FC::AnalyzeDt()
{   // Extract (n,f) counts from Dt spectra
    // and sum up neutron flux
    if (!DoneLimits)
        GetLimits();
//    if (!DoneDtBG)
//        this->AnalyzeDtBG(); // Estimate background
    cout << endl << "Analyzing ToF peaks..." << endl;

    for (int i = 0; i < NumHist; i++)
    {
        nFG[i] = 0;
        Double_t D2nFG = 0;
        tFG = 0;
        nBG[i] = 0;
        Double_t D2nBG = 0;
        tBG = 0;
        Double_t nFluence = 0;
        Double_t D2nFluence = 0;
        Double_t tMon = 0;
        for (Int_t j = 0; j < FgRuns; j++)
        {
            pHFG[j]->AnalyzeDtPeak(i, lim[0][i], lim[1][i], lim[2][i], lim[3][i]);
//            cout << "DEBUG: " << j+1 << " " << pHFG[j]->nNIF[i] << "+-" << pHFG[j]->DnNIF[i] << " " << pHFG[j]->t_live << " " << pHFG[j]->t_mon << endl;
            nFG[i] += pHFG[j]->nNIF[i];
            D2nFG += pow(pHFG[j]->DnNIF[i], 2);
            tFG += pHFG[j]->t_live;
            nFluence += pHFG[j]->NeutronFlux[i] * pHFG[j]->t_mon;
            D2nFluence += pow(pHFG[j]->DstatNeutronFlux[i] * pHFG[j]->t_mon, 2);
            tMon += pHFG[j]->t_mon;
        }
        DnFG[i] = sqrt(D2nFG);
        nFlux[i] = nFluence / tMon;
        DnFlux[i] = sqrt(D2nFluence) / tMon;
        for (Int_t j = 0; j < BgRuns; j++)
        {
            pHBG[j]->AnalyzeDtPeak(i, lim[0][i], lim[1][i], lim[2][i], lim[3][i]);
            nBG[i] += pHBG[j]->nNIF[i];
            D2nBG += pow(pHBG[j]->DnNIF[i], 2);
            tBG += pHBG[j]->t_live;
//            cout << "DEBUG: " << nBG[i] << "+-" << sqrt(D2nBG) << "   " << tBG << endl;
        }
        DnBG[i] = sqrt(D2nBG);
        if (DrawMulti)
        {
            plot->Dt(i, pHFG[0]->pHDtG[i], nFG[i], DnFG[i],
                     lim[0][i], lim[1][i], lim[2][i], lim[3][i], "FG");
        }
    }
    cout << "Ch   nFG   tFG   nBG   tBG" << endl;
    for (int i = 0; i < NumHist; i++)
        cout << " " << i+1 << "   " << nFG[i] << "+-" << DnFG[i] << "   " << tFG <<
                              "   " << nBG[i] << "+-" << DnBG[i] << "   " << tBG << endl;

    DoneDt = kTRUE;
    cout << "Done: ToF" << endl;
}


//void FC::ScatCorr(Int_t i)
//{
//    pDirect[i][0] = 1 - pHFG->NeutronFlux[i] / pHBG->NeutronFlux[i] * nBG[i] / nFG[i];
//    DpDirect[i][0] = sqrt( pow(pHFG->DNeutronFlux[i] / pHBG->NeutronFlux[i] * nBG[i] / nFG[i], 2) +
//                           pow(pHFG->NeutronFlux[i] * pHBG->DNeutronFlux[i] / pow(pHBG->NeutronFlux[i], 2) * nBG[i] / nFG[i], 2) +
//                           pow(pHFG->NeutronFlux[i] / pHBG->NeutronFlux[i] * DnBG[i] / nFG[i], 2) +
//                           pow(pHFG->NeutronFlux[i] / pHBG->NeutronFlux[i] * nBG[i] * DnFG[i] / pow(nFG[i], 2), 2) );
//}


void FC::ScatCorrDiff()
{
    if (!DoneDt)
        AnalyzeDt();
    if (!DoneRawCS)
        CrossSection();
    cout << endl << "Experimental Scattering correction: Shadow Cone" << endl;
    cout << "Ch   pDirect" << endl;
    for (int i = 0; i < NumCh; i++)
    {
//        ScatCorr(i);
        cout << " " << i+1 << "   " << pDirect[i][0] << " +- " << DpDirect[i][0] << endl;
    }
    DoneScatCorr = kTRUE;
    cout << "Done: Scattering correction" << endl;
}


//void FC::ScatCorrFit()
//{
//    if (!DoneLimits)
//        GetLimits();
//    if (!DoneDtBG)
//        AnalyzeDtBG();
//    cout << endl << "Experimental Scattering correction: Fit" << endl;

//    char name[32] = "";
//    Double_t pi = 3.1415926535898;
//    Double_t ChPerBin = pHFG[0]->pHDtG[0]->GetBinWidth(1);
//    Double_t ampl, Dampl, center, width, Dwidth; // center and width in units of QDC channels
//    Double_t low = pHFG[0]->pHDtG[0]->GetBinCenter(lim[0][0]);
//    Double_t up = pHFG[0]->pHDtG[0]->GetBinCenter(lim[3][0]);
//    for (int i = 0; i < NumHist; i++)
//    {
//        width = ChPerBin * (lim[2][i]-lim[1][i]) / 6;
//        center = 0.5 * ( pHFG[0]->pHDtG[0]->GetBinCenter(lim[1][i]) +
//                         pHFG[0]->pHDtG[0]->GetBinCenter(lim[2][i]) );
//        sprintf(name, "f%sNIFDtPeak_%i", Name.c_str(), i+1);
//        TF1* fFG = new TF1(name, func_peak, low, up, 4);
//        fFG->SetNpx(1000);
//        fFG->SetRange(low, up);
//        fFG->SetParameters(100, center, width, pHFG->t_live*avBg[i]);
//        pHFG->pHDtG[i]->Fit(name, "LQR");
//        ampl = (Double_t)fFG->GetParameter(0);
//        Dampl = (Double_t)fFG->GetParError(0);
//        width = (Double_t)fFG->GetParameter(2);
//        Dwidth = (Double_t)fFG->GetParError(2);
//        nFG[i] = sqrt(pi) * ampl * width / ChPerBin / pHFG->t_real;
//        DnFG[i] = sqrt(pi) / ChPerBin / pHFG->t_real * sqrt( pow(Dampl*width, 2) + pow(ampl*Dwidth, 2) );
////        cout << "  FG " << fFG->GetParameter(0) << "  " << fFG->GetParameter(1) << " " << fFG->GetParameter(2) << " " << fFG->GetParameter(3) << " " << nFG[i] << "+-" << DnFG[i] << endl;

//        sprintf(name, "f%sSBDtPeak_%i", Name.c_str(), i+1);
//        TF1* fBG = new TF1(name, func_peak, low, up, 4);
//        fBG->SetNpx(1000);
//        fBG->SetRange(low, up);
//        fBG->SetParameters(100, center, width, pHBG->t_live*avBg[i]);
//        fBG->FixParameter(2, width);
//        pHBG->pHDtG[i]->Fit(name, "LQR");
//        ampl = (Double_t)fBG->GetParameter(0);
//        Dampl = (Double_t)fBG->GetParError(0);
//        nBG[i] = sqrt(pi) * ampl * width / ChPerBin / pHBG->t_real;
//        DnBG[i] = sqrt(pi) / ChPerBin  / pHBG->t_real* sqrt( pow(Dampl*width, 2) + pow(ampl*Dwidth, 2) );
////        cout << "  BG " << fBG->GetParameter(0) << "  " << fBG->GetParameter(1) << " " << fBG->GetParameter(2) << " " << fBG->GetParameter(3) << " " << nBG[i] << "+-" << DnBG[i] << endl;

//        ScatCorr(i);
//    }

//    DoneScatCorr = kTRUE;
//    cout << "Done: Scattering correction" << endl;
//}


void FC::ScatCorrSim()
{
    if (!DoneDt)
        AnalyzeDt();
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
    if (!DoneDt)
        AnalyzeDt();
    if (!DoneNatoms)
        GetNatoms();
    cout << endl << "Uncorrected cross section..." << endl;
    Double_t avCS; // average cross section
    Double_t D2avCS;
    for (int i = 0; i < NumHist; i++)
    {
        cout << " Chanel " << i+1 << endl;
        if (CommentFlag)
            cout << " " << nAtoms[i] << "+-" << DnAtoms[i] << " atoms" << endl
                 << "  Run   t_live[s]   (n,f)-events   n-Flux[mm^-2 s^-1]   sigma[b]" << endl;
        Double_t avFG = 0;
        Double_t D2avFG = 0;
        Double_t avBG = 0;
        Double_t D2avBG = 0;

        /// Cross section ////////////////////////////////////////////////////////////
        // weighting with live times. Exact would be weighting with monitor counts
        for (Int_t j = 0; j < FgRuns; j++)
        {
            pHFG[j]->UncorrectedCrossSection(i);
            if (CommentFlag)
                cout << "  Fg " << j+1 << "   " << pHFG[j]->t_live << "   " << pHFG[j]->nNIF[i] << "+-" << pHFG[j]->DnNIF[i] << "   " << pHFG[j]->NeutronFlux[i] << "+-" << pHFG[j]->DstatNeutronFlux[i] << "   " << pHFG[j]->uncCS[i] << "+-" << pHFG[j]->DuncCS[i] << endl;
            avFG += pHFG[j]->t_live / tFG * pHFG[j]->uncCS[i];
            D2avFG += pow(pHFG[j]->t_live / tFG * pHFG[j]->DuncCS[j], 2);
        }
        avCS += avFG / NumHist;
        D2avCS += D2avFG / pow(NumHist, 2);
        for (Int_t j = 0; j < BgRuns; j++)
        {
            pHBG[j]->UncorrectedCrossSection(i);
            if (CommentFlag)
                cout << "  Bg " << j+1 << "   " << pHBG[j]->t_live << "   " << pHBG[j]->nNIF[i] << "+-" << pHBG[j]->DnNIF[i] << "   " << pHBG[j]->NeutronFlux[i] << "+-" << pHBG[j]->DstatNeutronFlux[i] << "   " << pHBG[j]->uncCS[i] << "+-" << pHBG[j]->DuncCS[i] << endl;
            avBG += pHBG[j]->t_live / tBG * pHBG[j]->uncCS[i];
            D2avBG += pow(pHBG[j]->t_live / tBG * pHBG[j]->DuncCS[j], 2);
        }
        pDirect[i][0] = 1.0 - avBG / avFG;
        DpDirect[i][0] = sqrt(D2avBG) / avFG;
//        if (CommentFlag)
//            cout << "   nFlux " << pHFG->NeutronFlux[i] << "+-" << pHFG->DNeutronFlux[i] << "mm^-2 s^-1" << endl <<
//                    "   nAtoms " << nAtoms[i] << "+-" << DnAtoms[i] << endl <<
//                    "   Foreground " << nFG[i] << "+-" << DnFG[i]  << " s^-1" << endl;
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


void FC::GetSimBg()
{
    if (!DoneSimFg)
        GetSimFg();
    cout << endl << "Simulation results requested. Starting " << Name << " background simulation analysis..." << endl;
    sim->ShadowCone();
    for (Int_t i = 0; i < NumCh; i++)
    {
        pDirect[i][2] = sim->SC[i];
        DpDirect[i][2] = sim->DSC[i];
    }
    DoneSimBg = kTRUE;
}


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
    if (!DrawSingle)
        return;
    plot->Result(uCS, DuCS, CS, DCS);
}


void FC::CompareShadowCone()
{
    if (!DoneScatCorr)
        ScatCorrDiff();
    if (!DoneSimBg)
        GetSimBg();

    cout << endl << "Comparing methods for In-scattering contribution..." << endl;
    cout << "Foreground" << endl <<
            " " << sim->Fg->nProj[0] << " Projectiles" << endl <<
            " " << sim->Fg->nDirect[0] << " direct neutrons" << endl <<
            " " << sim->Fg->effDirect[0] << " effective direct neutrons" << endl <<
            " " << sim->Fg->effScat[0] << " effective scattered neutrons" << endl;
    cout << "Background" << endl <<
            " " << sim->Bg->nProj[0] << " Projectiles" << endl <<
            " " << sim->Bg->nDirect[0] << " direct neutrons" << endl <<
            " " << sim->Bg->effDirect[0] << " effective direct neutrons" << endl <<
            " " << sim->Bg->effDirect[0] + sim->Bg->effScat[0] << " effective neutrons" << endl;
    cout << " Ch   Exp   Sim(1)   Sim(2)" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
        cout << " " << i+1 << "  " << pDirect[i][0] << " +- " << DpDirect[i][0]
                           << "  " << pDirect[i][1] << " +- " << DpDirect[i][1]
                           << "  " << pDirect[i][2] << " +- " << DpDirect[i][2] << endl;
    }

    if (!DrawSingle)
        return;
    Double_t pExp[NumCh];
    Double_t DpExp[NumCh];
    Double_t pSim1[NumCh];
    Double_t DpSim1[NumCh];
    Double_t pSim2[NumCh];
    Double_t DpSim2[NumCh];
    for (Int_t i = 0; i < NumCh; i++)
    {
        pExp[i] = pDirect[i][0];
        DpExp[i] = DpDirect[i][0];
        pSim1[i] = pDirect[i][1];
        DpSim1[i] = DpDirect[i][1];
        pSim2[i] = pDirect[i][2];
        DpSim2[i] = DpDirect[i][2];
    }

    plot->CompSc(&(pExp[0]), &(DpExp[0]), &(pSim1[0]), &(DpSim1[0]), &(pSim2[0]), &(DpSim2[0]));
}


void FC::CompareTransmission()
{
    if (!DrawSingle)
        return;
    if (!DoneSimFg)
        GetSimFg();
    if (!DoneTransmission)
        ExpTrans();
    plot->CompT(ExpT, DExpT, SimT, DSimT);
}


void FC::ExpTrans()
{
    if (!DoneDt)
        AnalyzeDt();
    if (!DoneDtBG)
        AnalyzeDtBG();
    GetExpT();
    if (!DoneSimFg)
        GetSimFg();
    cout << endl << "Experimental transmission" << endl;

    cout << " Ch   expT   simT   anaT" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        ExpT[i] = uT[i] * pow(sd[i] / sd[7], 2);
        DExpT[i] = DuT[i] * pow(sd[i] / sd[7], 2);
        if (CommentFlag)
            cout << " " << i+1 << "   " << ExpT[i] << "+-" << DExpT[i] << "   " << SimT[i] << "+-" << DSimT[i] << endl;
    }
    //
    cout << "Done: Experimental transmission" << endl;
    DoneTransmission = kTRUE;
    if (!DrawSingle)
        return;
    plot->ExpT(ExpT, DExpT, SimT, DSimT);
}


void FC::SaveToFile(string path, TObject *pObj)
{   //saves a TObject into the selected file fname
    TFile* file = TFile::Open("/home/hoffma93/Programme/Go4nfis/offline/results/Evaluation.root", "UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder path already exists, otherwise create it
    if (file->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) file->Get(path.c_str());
    else
    {   file->mkdir(path.c_str(), "Folder containing offline Analysis objects");
        EvalDir = file->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    file->Save(); file->Close();
}


Double_t FC::func_peak(Double_t *x, Double_t *par)
{
    return par[0] * exp( - pow((x[0] - par[1]) / par[2], 2)) + par[3];
}
