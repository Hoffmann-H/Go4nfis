#ifndef HIST_H
#define HIST_H
#include "TH1F.h"
#include "TH1I.h"
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#define NumHist 8

class Hist : public TObject
{
public:
    Hist(const char *file_path = "");  // Constructor
//     Hist(); // Standard constructor, needed for object i/o
    ~Hist(); // Destructor
    
    const char* FilePath = "";
    
    Bool_t CommentFlag;

    TH1I* GetTH1I(const char *hname);
    TH1F* GetTH1F(const char *hname);
    TH1D* GetTH1D(const char *hname);

    void Save(TH1I* pHtoSave, string hpath, string hname);
    void SaveCanvas(TCanvas* pCtoSave, string hpath, string hname);
    
//    ClassDef(Hist, 1);
};
#endif
