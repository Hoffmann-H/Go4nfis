#ifndef SAVE_TO_FILE
#define SAVE_TO_FILE

#include "TF1.h"

TDirectory* Prepare(TFile *f, string path)
{
    TDirectory *pDir = f->GetDirectory(path.c_str());
    if (pDir == 0)
    {
        cout << " Creating root directory " << path << endl;
        f->mkdir(path.c_str(), "Folder");
        pDir = f->GetDirectory(path.c_str());
    }
    pDir->cd();
    return pDir;
}

void Save(TFile *f, string path, TH1 *pObj, string setName = "")
{
    TDirectory *pDir = Prepare(f, path);
    pObj->Write(setName.c_str(), TObject::kOverwrite);
    return;
}

void Save(TFile *f, string path, TGraph *pObj, string setName = "")
{
    TDirectory *pDir = Prepare(f, path);
    pObj->Write(setName.c_str(), TObject::kOverwrite);
    return;
}

void Save(TFile *f, string path, TF1 *pObj, string setName = "")
{
    TDirectory *pDir = Prepare(f, path);
    pObj->Write(setName.c_str(), TObject::kOverwrite);
    return;
}

//// Usage:
////    TFile *f = TFile::Open("Path/To/Root/File", "UPDATE");
////    string path = "Root/Tree";
////    TDirectory *pDir = Prepare(f, path);
////
////    TH1F *h = new TH1F("NameDummy", ...);
////    Fill...
////    Save(pDir, h, "Name");
////
////    f->Save();
////    f->Close();
void Save(TDirectory *pDir, TH1 *pObj, string name)
{
    if (pDir->Get(name.c_str()) != 0)
    {
        cout << name << " already exists!" << endl;
        pDir->Delete((name+";*").c_str());
    } else {
        cout << name << " did not exist yet!" << endl;
    }
    pObj->SetName(name.c_str());
    pDir->cd();//evtl
    pObj->Write();
}


void Save(TDirectory *pDir, TGraph *pObj, string name)
{
    if (pDir->Get(name.c_str()) != 0)
    {
        cout << name << " already exists!" << endl;
        pDir->Delete((name+";*").c_str());
    } else {
        cout << name << " did not exist yet!" << endl;
    }
    pObj->SetName(name.c_str());
    pDir->cd();//evtl
    pObj->Write();
}
void Save(TDirectory *pDir, TF1 *pObj, string name)
{
    if (pDir->Get(name.c_str()) != 0)
    {
        cout << name << " already exists!" << endl;
        pDir->Delete((name+";*").c_str());
    } else {
        cout << name << " did not exist yet!" << endl;
    }
    pObj->SetName(name.c_str());
    pDir->cd();//evtl
    pObj->Write();
}

//// Usage:
////    TFile *f = TFile::Open("Path/To/Root/File", "UPDATE");
////    string path = "Root/Tree";
////    TDirectory *pDir = Prepare(f, path);
////
////    string name = "...";
////    Clear(pDir, name);
////    TH1F *h = new TH1F(name, ...);
////    Fill...
////    h->Write();
////
////    f->Save();
////    f->Close();
void Clear(TDirectory *pDir, string name)
{ // Only works properly if new histogram with same name is not created yet
  // Call this before histogram creation
    if (pDir->Get(name.c_str()) != 0)
    {
        cout << name << " already exists!" << endl;
        pDir->Delete((name+";*").c_str());
    } else {
        cout << name << " did not exist yet!" << endl;
    }
}

//void TestStF()
//{
//    TFile *f = TFile::Open("~/Programme/ROOT/TestStF.root", "UPDATE");
//    string path = "Folder1/SubFolder1";
//    TDirectory *pDir = Prepare(f, path);

//    string name = "h1";
//    TH1F *h1 = new TH1F(name.c_str(), name.c_str(), 100, -5, 5);
//    h1->FillRandom("gaus", 500);
//    Save(pDir, h1, "H1");

//    name = "h2";
//    TH1F *h2 = new TH1F(name.c_str(), name.c_str(), 100, -5, 5);
//    h2->FillRandom("gaus", 500);
//    Save(pDir, h2, "H2");

//    f->Save();
//    f->Close();
//}

void SaveToFile(TFile *f, string path, TObject *pObj)
{
    f->ReOpen("UPDATE");
    TDirectory *EvalDir;
    TObject *pGraph;
    pGraph = (TObject*) pObj;
    string GraphName = pGraph->GetName();
    //check if folder path already exists, otherwise create it
    if (f->Get(path.c_str())!=0)
        EvalDir = (TDirectory*) f->Get(path.c_str());
    else
    {   cout << " Creating root directory " << path << endl;
        f->mkdir(path.c_str(), "Folder");
        EvalDir = f->GetDirectory(path.c_str());
    }
    EvalDir->cd();
    if (EvalDir->Get(GraphName.c_str())!=0)
        EvalDir->Delete((GraphName+";*").c_str());
    pGraph->Clone()->Write();
    f->Save(); //file->Close();
//    cout << " Saved " << GraphName << endl;
}
#endif
