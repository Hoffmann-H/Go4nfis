// Class for loading and accessing histograms from a root file. 
// This is meant to be created once for each exp. setup: NIF, SB. 
// Open histograms: Hist *pHist = new Hist("path/to/file.root")

#include "Hist.h"
#include <TH1I.h>
#include "TFile.h"
#include <string.h>
#include <iostream>
#include <sstream>

using namespace std;

Hist::Hist(const char *file_path)
{
    // Constructor
    CommentFlag = kFALSE;
    if(CommentFlag)
        cout << "Create instance Hist using " << file_path << endl;
    FilePath = file_path;
}

Hist::~Hist()
{
    // Destructor
    cout << "Destroy" << endl;
}

TH1I* Hist::GetTH1I(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile *f = TFile::Open(FilePath);
    if (f==0) //check if root file exists
    {   cerr << "File not found!" << endl;
        return 0;
    }
    else
    {   TH1I* histo = (TH1I*) f->Get(hname)->Clone();
        if (histo==0)
        {   cerr << "Histogram not found!" << endl;
            return 0;
        }
        else
            return histo;
    }
    f->Close();
}

TH1F* Hist::GetTH1F(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile *f = TFile::Open(FilePath);
    if (f==0) //check if root file exists
    {   cerr << "File not found!" << endl;
        return 0;
    }
    else
    {   TH1F* histo = (TH1F*) f->Get(hname)->Clone();
        if (histo==0)
        {   cerr << "Histogram not found!" << endl;
            return 0;
        }
        else
            return histo;
    }
    f->Close();
}

TH1D* Hist::GetTH1D(const char* hname)
{   //opens a histogram in rootfile fname with (root path and name) hname
    TFile *f = TFile::Open(FilePath);
    if (f==0) //check if root file exists
    {   cerr << "File not found!" << endl;
        return 0;
    }
    else
    {   TH1D* histo = (TH1D*) f->Get(hname)->Clone();
        if (histo==0)
        {   cerr << "Histogram not found!" << endl;
            return 0;
        }
        else
            return histo;
    }
    f->Close();
}

void Hist::Save(TH1I* pHtoSave, string hpath, string hname)
{
    // Open root file
    TFile *f = TFile::Open(FilePath, "UPDATE");
    // Check if folder exists. If not, create.
    TDirectory *SaveDir;
    if (f->Get(hpath.c_str())!=0)
        SaveDir = (TDirectory*) f->Get(hpath.c_str());
    else
    {
        cout << "Create dir " << hpath << endl;
        SaveDir = f->mkdir(hpath.c_str(), "subdir");
    }
    string dir_name = "/"+hpath;
    SaveDir->cd(dir_name.c_str());
    // Delete old data
    if (SaveDir->Get(hname.c_str())!=0) {
        string del_name = hname+";*";
        SaveDir->Delete(del_name.c_str());
        if (CommentFlag)
            cout << "Delete old " << hname << endl;
    }
    else if (CommentFlag)
        cout << "Hist not existing yet." << endl;
    // Save
    pHtoSave->Clone()->Write(hname.c_str());
    f->Save();
    f->Close();
}

void Hist::SaveCanvas(TCanvas* pCtoSave, string hpath, string hname)
{
    // Open root file
    TFile *f = TFile::Open(FilePath, "UPDATE");
    // Check if folder exists. If not, create.
    TDirectory *SaveDir;
    if (f->Get(hpath.c_str())!=0)
        SaveDir = (TDirectory*) f->Get(hpath.c_str());
    else
    {
        cout << "Create dir " << hpath << endl;
        SaveDir = f->mkdir(hpath.c_str(), "subdir");
    }
    string dir_name = "/"+hpath;
    SaveDir->cd(dir_name.c_str());
    // Delete old data
    if (SaveDir->Get(hname.c_str())!=0) {
        string del_name = hname+";*";
        SaveDir->Delete(del_name.c_str());
        if (CommentFlag)
            cout << "Delete old " << hname << endl;
    }
    else if (CommentFlag)
        cout << "Hist not existing yet." << endl;
    // Save
    pCtoSave->Clone()->Write(hname.c_str());
    f->Save();
    f->Close();
}
