#include "Hist.cxx"
#include "Xsection.cxx"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TFile.h>

using namespace std;

void run()
{
    Xsection *pXs = new Xsection();
    pXs -> DoAnalyze();
    return;
}
