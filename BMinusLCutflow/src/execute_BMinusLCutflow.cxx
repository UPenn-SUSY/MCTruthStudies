#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"

#include "include/BMinusLCutflow.h"


// -----------------------------------------------------------------------------
void help()
{
  std::cout << "usage:\n\t./Cutflow <INPUT FILE NAME>\n\n";
}

// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  std::cout << "testing basic truth ntuple looper\n";

  if (argc < 2) {
    help();
    return 0;
  }

  std::cout << "Number of input files: " << argc-1 << "\n";
  std::cout << "\n";

  TChain t("truth");
  bool isSignal = true;
  std::cout << "input files:\n";
  for (int it = 1; it != argc; ++it) {
    std::cout << "\t" << it << " -- " << argv[it] << "\n";
    if (strcmp(argv[it],"-s") == 0) {
    std::cout << strcmp(argv[it],"-s") << "\n";
      isSignal = true;
    }
    else if (strcmp(argv[it],"-b") == 0) {
    std::cout << strcmp(argv[it],"-b") << "\n";
      isSignal = false;
    }
    else {
      std::cout << "\t" << it << " -- " << argv[it] << "\n";
    t.Add(argv[it]);
    }
  }
  std::cout << "\n";

  std::cout << "Retrieved tree from input files. Constructing Cutflow object\n";
  BMinusL::Cutflow bmlcf(&t, isSignal);

  std::cout << "Preparing to loop over events!\n";
  bmlcf.Loop();

  std::cout << "Done Looping! Writing to file!\n";

  bmlcf.writeToFile();

  std::cout << "\n\n";
}
