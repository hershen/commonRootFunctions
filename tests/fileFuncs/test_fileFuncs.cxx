#include "catch.hpp"

#include "fileFuncs.h"

#include "TChain.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

using namespace myFuncs;

SCENARIO("Copy ttree from file to file", "[copyTTreeFromFileToFile]") {
  GIVEN("A TFile with a TTree with some info") {
    const std::string fromFilename = "/home/hershen/temp/fromFile.root";
    TFile fromFile(fromFilename.c_str(), "RECREATE");
    TTree tree("someTreeName", "");
    TRandom3 r(0);
    double a = r.Uniform();
    double b = r.Uniform();
    double c = r.Uniform();
    tree.Branch("a", &a);
    tree.Branch("b", &b);
    tree.Branch("c", &c);
    tree.Fill();
    fromFile.Write();
    fromFile.Close();

    WHEN("Copying the ttree to another file") {
      const std::string toFilename = "/home/hershen/toFile.root";
      copyTTreeFromFileToFile("someTreeName", fromFilename, toFilename);
      THEN("The new file should contain the tree") {
        std::vector<std::string> branches{"a", "b", "c"};
        double newA = 0;
        double newB = 0;
        double newC = 0;
        std::vector<void *> pointers{&newA, &newB, &newC};
        auto chain = myFuncs::openChain_setBranch(toFilename, "someTreeName", branches, pointers);
        CHECK(chain->GetEntries() == 1);

        chain->GetEntry(0);

        CHECK(a == Approx(newA));
        CHECK(b == Approx(newB));
        CHECK(c == Approx(newC));

        delete chain;

      }
    }
  }
}
