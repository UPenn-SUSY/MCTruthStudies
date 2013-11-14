#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"

#include "include/TruthNtupleLooper.h"

// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  std::cout << "testing basic truth ntuple looper\n";

  std::cout << "input file name: " << argv[1] << "\n";

  TFile* f = new TFile(argv[1]);
  TTree* t = static_cast<TTree*>(f->Get("truth"));

  TruthNtuple::TruthNtupleLooper tnl(t);
  tnl.Loop();

  delete f;
}
