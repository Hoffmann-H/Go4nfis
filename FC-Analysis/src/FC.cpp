//
#include "FC.h"

using namespace std;

FC::FC()//Double_t nFieldFG, Double_t DnFieldFG, Double_t nFieldBG, Double_t DnFieldBG)
{

}


FC::~FC()
{

}


void FC::InitVar()
{
    cout << endl << "Initializing common variables..." << endl;
    Draw = kFALSE;
    DoneQDC = kFALSE;
    DoneThresholds = kFALSE;
    DoneLimits = kFALSE;
    DoneNatoms = kFALSE;
    DoneDtBG = kFALSE;
    DoneDt = kFALSE;
    DoneScatCorr = kFALSE;
    DoneSim = kFALSE;
    DoneIso = kFALSE;
    DoneRawCS = kFALSE;

    u = 1.660539E-24;
    DtBgLow = 63600; // QDC channels
    DtBgUp = 78500;
    Area = 4300;
    DArea = 120;
    Yield = 2.15882E4;
    DYield = 543.3;
    cout << "Done: common variables" << endl;
}


void FC::GetLimits(Double_t n)
{   // Set integration limits to n sigma
    cout << endl << "Opening peak integration limits..." << endl;

    Double_t ChPerBin = pHFG->pHDtG[0]->GetBinWidth(1);
    Double_t BinOffset = pHFG->pHDtG[0]->GetBinCenter(1);
    if (CommentFlag)
        cout << " ChPerBin: " << ChPerBin << ", BinOffset: " << BinOffset << endl;
    cout << "Ch   Peak   bins" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        Double_t ToF_low = pHFG->GetPeakLow(i);
        Double_t ToF_up = pHFG->GetPeakUp(i); // Limits for background are the same
        Double_t center = 0.5 * (ToF_low + ToF_up);
        Double_t width = (ToF_up - ToF_low) / 6.0;
        DtPeakLow[i] = center - n * width;
        DtPeakUp[i] = center + n * width;

        lim[0][i] = (DtBgLow - BinOffset) / ChPerBin;
        lim[1][i] = (DtPeakLow[i] - BinOffset) / ChPerBin;
        lim[2][i] = (DtPeakUp[i] - BinOffset) / ChPerBin;
        lim[3][i] = (DtBgUp - BinOffset) / ChPerBin;
        cout << " " << i+1 << "   " << DtPeakLow[i] << " " << DtPeakUp[i] << "   " << lim[0][i] << " " << lim[1][i] << " " << lim[2][i] << " " << lim[3][i] << " " << endl;
    }
    DoneLimits = kTRUE;
    cout << "Done: peak integration limits" << endl;
}


void FC::AnalyzeDt()
{   // Extract (n,f) counts from Dt spectra
    if (!DoneLimits)
        GetLimits();
    if (!DoneDtBG)
        this->AnalyzeDtBG(); // Estimate background
    cout << endl << "Analyzing TimeDiff peaks..." << endl;

    Double_t PeakInt;
    cout << "Ch   nFG   nBG" << endl;
    for (int i = 0; i < NumHist; i++)
    {
        // Integrate peaks: FG
        PeakInt = pHFG->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
//        cout << "DEBUG: " << PeakInt << " " << pHFG->t_live << " " << avBg[i] << endl;
        nFG[i] = PeakInt - (lim[2][i] - lim[1][i] + 1) * pHFG->t_live * avBg[i];
        DnFG[i] = sqrt(PeakInt + pow((lim[2][i] - lim[1][i] + 1) * pHFG->t_live * DavBg[i], 2));

        PeakInt = pHBG->pHDtG[i]->Integral(lim[1][i], lim[2][i]);
        nBG[i] = PeakInt - (lim[2][i] - lim[1][i] + 1) * pHBG->t_live * avBg[i];
        DnBG[i] = sqrt(PeakInt + pow((lim[2][i] - lim[1][i] + 1) * pHBG->t_live * DavBg[i], 2));
        cout << " " << i+1 << "   " << nFG[i] << "+-" << DnFG[i] << "   " << nBG[i] << "+-" << DnBG[i] << endl;
        if (Draw)
        {
            plot->Dt(i, pHFG->pHDtG[i], nFG[i], DnFG[i],
                     lim[0][i], lim[1][i], lim[2][i], lim[3][i], pHFG->t_live * avBg[i], "FG");

        }
    }

    DoneDt = kTRUE;
    cout << "Done: TimeDiff" << endl;
}


void FC::ScatCorr(Int_t i)
{
    pDirect[i][0] = 1 - pHFG->Monitor / pHBG->Monitor * nBG[i] / nFG[i];
    DpDirect[i][0] = sqrt( pow(pHFG->DMonitor / pHBG->Monitor * nBG[i] / nFG[i], 2) +
                           pow(pHFG->Monitor * pHBG->DMonitor / pow(pHBG->Monitor, 2) * nBG[i] / nFG[i], 2) +
                           pow(pHFG->Monitor / pHBG->Monitor * DnBG[i] / nFG[i], 2) +
                           pow(pHFG->Monitor / pHBG->Monitor * nBG[i] * DnFG[i] / pow(nFG[i], 2), 2) );
}


void FC::ScatCorrDiff()
{
    if (!DoneDt)
        AnalyzeDt();
    cout << endl << "Experimental Scattering correction: Shadow Cone" << endl;
    cout << "Ch   pDirect" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        ScatCorr(i);
        cout << " " << i+1 << "   " << pDirect[i][0] << " +- " << DpDirect[i][0] << endl;
    }
    DoneScatCorr = kTRUE;
    cout << "Done: Scattering correction" << endl;
}


void FC::ScatCorrFit()
{
    if (!DoneLimits)
        GetLimits();
    if (!DoneDtBG)
        AnalyzeDtBG();
    cout << endl << "Experimental Scattering correction: Fit" << endl;

    char name[32] = "";
    Double_t pi = 3.1415926535898;
    Double_t ChPerBin = pHFG->pHDtG[0]->GetBinWidth(1);
    Double_t ampl, Dampl, center, width, Dwidth; // center and width in units of QDC channels
    Double_t low = pHFG->pHDtG[0]->GetBinCenter(lim[0][0]);
    Double_t up = pHFG->pHDtG[0]->GetBinCenter(lim[3][0]);
    for (int i = 0; i < NumHist; i++)
    {
        width = ChPerBin * (lim[2][i]-lim[1][i]) / 6;
        center = 0.5 * ( pHFG->pHDtG[0]->GetBinCenter(lim[1][i]) +
                         pHFG->pHDtG[0]->GetBinCenter(lim[2][i]) );
        sprintf(name, "f%sNIFDtPeak_%i", Name.c_str(), i+1);
        TF1* fFG = new TF1(name, func_peak, low, up, 4);
        fFG->SetNpx(1000);
        fFG->SetRange(low, up);
        fFG->SetParameters(100, center, width, pHFG->t_live*avBg[i]);
        pHFG->pHDtG[i]->Fit(name, "LQR");
        ampl = (Double_t)fFG->GetParameter(0);
        Dampl = (Double_t)fFG->GetParError(0);
        width = (Double_t)fFG->GetParameter(2);
        Dwidth = (Double_t)fFG->GetParError(2);
        nFG[i] = sqrt(pi) * ampl * width / ChPerBin;
        DnFG[i] = sqrt(pi) / ChPerBin * sqrt( pow(Dampl*width, 2) + pow(ampl*Dwidth, 2) );
//        cout << "  FG " << fFG->GetParameter(0) << "  " << fFG->GetParameter(1) << " " << fFG->GetParameter(2) << " " << fFG->GetParameter(3) << " " << nFG[i] << "+-" << DnFG[i] << endl;

        sprintf(name, "f%sSBDtPeak_%i", Name.c_str(), i+1);
        TF1* fBG = new TF1(name, func_peak, low, up, 4);
        fBG->SetNpx(1000);
        fBG->SetRange(low, up);
        fBG->SetParameters(100, center, width, pHBG->t_live*avBg[i]);
        fBG->FixParameter(2, width);
        pHBG->pHDtG[i]->Fit(name, "LQR");
        ampl = (Double_t)fBG->GetParameter(0);
        Dampl = (Double_t)fBG->GetParError(0);
        nBG[i] = sqrt(pi) * ampl * width / ChPerBin;
        DnBG[i] = sqrt(pi) / ChPerBin * sqrt( pow(Dampl*width, 2) + pow(ampl*Dwidth, 2) );
//        cout << "  BG " << fBG->GetParameter(0) << "  " << fBG->GetParameter(1) << " " << fBG->GetParameter(2) << " " << fBG->GetParameter(3) << " " << nBG[i] << "+-" << DnBG[i] << endl;

        ScatCorr(i);
    }

    DoneScatCorr = kTRUE;
    cout << "Done: Scattering correction" << endl;
}


void FC::ScatCorrSim()
{
    if (!DoneDt)
        AnalyzeDt();
    if (!DoneSim)
        GetSimRes();
    cout << endl << "Scattering correction: Simulaion" << endl;
//    for (int i = 0; i < NumHist; i++)
//    {
////        nfDirect[i] = nFG[i] * (1.0 - InScat[i]);
////        DnfDirect[i] = sqrt( pow(DnFG[i] * (1.0 - InScat[i]), 2) +
////                            pow(nFG[i] * DInScat[i], 2) );
//        cout << " ch " << i+1 << " " << S[i] << " +- " << DInScat[i] << endl;
//    }
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
//    Double_t sumCS = 0; // cross section all channels' sum
//    Double_t D2sumCS = 0;
//    Double_t avCS; // average cross section
//    Double_t DavCS;
    for (int i = 0; i < NumHist; i++)
    {
        /// Cross section ////////////////////////////////////////////////////////////
        uCS[i] = nFG[i] / (nAtoms[i] * pHFG->NeutronFluence[i]); ///////////////////
        //////////////////////////////////////////////////////////////////////////////

        DuCS[i] = sqrt( pow(DnFG[i] / (nAtoms[i] * pHFG->NeutronFluence[i]), 2) +
                       pow(nFG[i] * DnAtoms[i] / (nAtoms[i] * nAtoms[i] * pHFG->NeutronFluence[i]), 2) +
                       pow(nFG[i] * pHFG->DNeutronFluence[i] / (nAtoms[i] * pHFG->NeutronFluence[i]*pHFG->NeutronFluence[i]), 2) );
//        sumCS += uCS[i];
//        D2sumCS += DuCS[i] * DuCS[i];
        cout << " Ch " << i+1 << endl;
        if (CommentFlag)
            cout << "   nFluence " << pHFG->NeutronFluence[i] << "+-" << pHFG->DNeutronFluence[i] << endl <<
                    "   nAtoms " << nAtoms[i] << "+-" << DnAtoms[i] << endl <<
                    "   Foreground " << nFG[i] << "+-" << DnFG[i] << endl;
        cout << "  raw CS = " << uCS[i] * 1.E22 << "+-" << DuCS[i] * 1.E22 << " barn" << endl;
    }
//    avCS = sumCS / NumHist;
//    DavCS = sqrt(D2sumCS) / NumHist;
//    cout << " Average: sigma = " << avCS << "+-" << DavCS << endl;
    DoneRawCS = kTRUE;
    cout << "Done: raw cross section" << endl;
}


void FC::GetSimRes()
{
    cout << endl << "Simulation results requested. Starting " << Name << " analysis..." << endl;
    sim = new AnaSim(SimPath, Name, "Open", plot);
    sim->Corrections();
    for (Int_t i = 0; i < NumCh; i++)
    {
        fTS[i] = sim->F[i];
        DfTS[i] = sim->DF[i];
        pDirect[i][1] = sim->S[i]; // simplified identity. Assuming equal spectra for Direct and Full-E neutrons.
        DpDirect[i][1] = sim->DS[i];
    }
    DoneSim = kTRUE;
    cout << "Done: Got simulation results" << endl;
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
    cout << "Ch   uncorrected CS   corrected CS" << endl;
    for (int i = 0; i < NumCh; i++)
    {
        CS[i] = fIsoVec * fTS[i] * uCS[i] - sIsoVec;
        Double_t uncert[] = {DfIsoVec * fTS[i] * uCS[i],
                             fIsoVec * DfTS[i] * uCS[i],
                             fIsoVec * fTS[i] * DuCS[i],
                             DsIsoVec };
        DCS[i] = sqrt( pow(uncert[0], 2) +
                       pow(uncert[1], 2) +
                       pow(uncert[2], 2) +
                       pow(uncert[3], 2) );
        cout << " " << i+1 << ",  " << uCS[i] * 1.E22 << "+-" << DuCS[i] * 1.E22 << " barn,  " << CS[i] * 1.E22 << "+-" << DCS[i] * 1.E22 << " barn" << endl;
        cout << "  unc: f " << uncert[0] << ", TS " << uncert[1] << ", u " << uncert[2] << ", s " << uncert[3] << endl;
    }
    DoneCorrections = kTRUE;
    cout << "Done: Corrections" << endl;
}


void FC::CompareShadowCone(string BgPath)
{
    if (!DoneScatCorr)
        ScatCorrDiff();
    if (!DoneSim)
        GetSimRes();
    cout << "FC::CompareShadowCone(): sim->nProj = " << sim->nProj[0] << endl;
    AnaSim *simBg = new AnaSim(BgPath, Name, "SB", plot);
    simBg->Corrections();
    cout << endl << "Comparing methods for In-scattering contribution..." << endl;
    cout << "Foreground" << endl <<
            " " << sim->nProj[0] << " Projectiles" << endl <<
            " " << sim->nDirect[0] << " direct neutrons" << endl <<
            " " << sim->effFullE[0] << " effective Full-E neutrons" << endl <<
            " " << sim->effSc[0] << " effective scattered neutrons" << endl;
    cout << "Background" << endl <<
            " " << simBg->nProj[0] << " Projectiles" << endl <<
            " " << simBg->nDirect[0] << " direct neutrons" << endl <<
            " " << simBg->effFullE[0] << " effective Full-E neutrons" << endl <<
            " " << simBg->effFullE[0] + simBg->effSc[0] << " effective neutrons" << endl;
    cout << " Ch   Exp   Sim(1)   Sim(2)" << endl;
    for (Int_t i = 0; i < NumCh; i++)
    {
//        cout << "   " << simBg->effSc[i] << "  " << sim->effSc[i] << "  " << sim->nProj[i] << "  " << simBg->nProj[i] << endl;
        pDirect[i][2] = 1.0 - (simBg->effFullE[i] + simBg->effSc[i]) / (sim->effFullE[i] + sim->effSc[i]) * sim->nProj[i] / simBg->nProj[i];
        DpDirect[i][2] = sqrt( pow(sim->DeffFullE[i] * (simBg->effFullE[i] + simBg->effSc[i]), 2) / pow(sim->effFullE[i] + sim->effSc[i], 4) +
                               pow(sim->DeffSc[i] * (simBg->effFullE[i] + simBg->effSc[i]), 2) / pow(sim->effFullE[i] + sim->effSc[i], 4) +
                               pow(simBg->DeffFullE[i] / (sim->effFullE[i] + sim->effSc[i]), 2) +
                               pow(simBg->DeffSc[i] / (sim->effFullE[i] + sim->effSc[i]), 2) ); // nProj: exact
        cout << " " << i+1 << "  " << pDirect[i][0] << "+-" << DpDirect[i][0]
                           << "  " << pDirect[i][1] << "+-" << DpDirect[i][1]
                           << "  " << pDirect[i][2] << "+-" << DpDirect[i][2] << endl;
    }
//    cout << "Uncertainty influences on shadow cone simulation" << endl
//         << " Ch   nFis(Fg)   nFis(Bg)" << endl;
//    for (Int_t i = 0; i < NumCh; i++)
//    {
//        Double_t unc[] = {sim->DnFis[i] * simBg->nFis[i] / sim->nFis[i] / sim->nFis[i] * sim->nProj[i] / simBg->nProj[i],
//                          simBg->DnFis[i] / sim->nFis[i] * sim->nProj[i] / simBg->nProj[i]};
//        cout << " " << i+1 << "   " << unc[0] << "   " << unc[1] << endl;
//    }
    if (!Draw)
//        return;
        plot = new Plot(Name, "Open");
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
