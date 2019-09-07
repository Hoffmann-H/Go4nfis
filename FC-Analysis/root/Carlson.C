//#include "SaveToFile.C"
#include "NumberOfAtoms.C"

TGraph* Spectrum(string fac = "nELBE")
{
    string file_name;
    if (!strcmp(fac.c_str(), "nELBE"))
        file_name = "/home/hoffma93/Programme/ROOT/Data/nELBE_E.dat";
    else
        file_name = "/home/hoffma93/Programme/ROOT/Data/Source_E.dat";
    TGraph *g = new TGraphErrors(file_name.c_str());
    if (g == 0)
        cout << "Fehler beim Öffnen von " << file_name << endl;
    return g;
}

TGraphErrors* Result(Int_t isotope, string entry = "eta")
{
    // Read format: E eta(235) W(235) a2(235) eta(238) W(238) a2(238) DE Deta(235) DW(235) Da2(235) Deta(238) DW(238) Da2(238)
    string file_name = "/home/hoffma93/Experiment/Carlson-Korrektur/results/Carlson_nELBE.dat";
    string format;
    if (isotope == 235)
    {
        if (!strcmp(entry.c_str(), "eta"))
            format = "%lg %lg %*lg %*lg %*lg %*lg %*lg %lg %lg";
        if (!strcmp(entry.c_str(), "W"))
            format = "%lg %*lg %lg %*lg %*lg %*lg %*lg %lg %*lg %lg";
        if (!strcmp(entry.c_str(), "a2"))
            format = "%lg %*lg %*lg %lg %*lg %*lg %*lg %lg %*lg %*lg %lg";
    } else if (isotope == 238) {
        if (!strcmp(entry.c_str(), "eta"))
            format = "%lg %*lg %*lg %*lg %lg %*lg %*lg %lg %*lg %*lg %*lg %lg";
        if (!strcmp(entry.c_str(), "W"))
            format = "%lg %*lg %*lg %*lg %*lg %lg %*lg %lg %*lg %*lg %*lg %*lg %lg";
        if (!strcmp(entry.c_str(), "a2"))
            format = "%lg %*lg %*lg %*lg %*lg %*lg %lg %lg %*lg %*lg %*lg %*lg %*lg %lg";
    } else
        cout << "Unknown isotope " << isotope << endl;

    TGraphErrors *ge = new TGraphErrors(file_name.c_str(), format.c_str());
    if (ge == 0)
        cout << "Fehler beim Öffnen von " << file_name << " im Format " << format << endl;
    return ge;
}

void Carlson()
{
    TGraphErrors *geU235eta = Result(235, "eta");
    TGraphErrors *geU235W   = Result(235, "W");
    TGraphErrors *geU235a2  = Result(235, "a2");
    TGraphErrors *geU238eta = Result(238, "eta");
    TGraphErrors *geU238W   = Result(238, "W");
    TGraphErrors *geU238a2  = Result(238, "a2");
//    TGraph *g = Spectrum("nELBE");
    TGraphErrors *geEffUmA = EffUmA();
    new TCanvas();
    geEffUmA->Draw();
    TFile *f = TFile::Open("/home/hoffma93/Experiment/Carlson-Korrektur/results/Carlson.root", "UPDATE");
    TDirectory *pDir;
    pDir = Prepare(f, "nELBE/U-235");
    Save(pDir, geU235eta, "geU235eta");
    Save(pDir, geU235W, "geU235W");
    Save(pDir, geU235a2, "geU235a2");
    pDir = Prepare(f, "nELBE/U-238");
    Save(pDir, geU238eta, "geU238eta");
    Save(pDir, geU238W, "geU238W");
    Save(pDir, geU238a2, "geU238a2");
    f->Save();
    f->Close();
}
