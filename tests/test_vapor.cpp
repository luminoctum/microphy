#include <iostream>
#include <debugger.hpp>
#include <mesh/mesh.hpp>
#include <cantera/thermo/IdealGasPhase.h>
#include <vapor.hpp>
#include <condensates.hpp>
//#include <cantera/base/Solution.h>

int main() {
  // no MPI
  Globals::my_rank = 0;
  Globals::nranks  = 1;

  pdebug = std::make_unique<Debugger>(1);

  IOWrapper infile;
  infile.Open("test_vapor.inp", IOWrapper::FileMode::read);

  ParameterInput *pinput = new ParameterInput();
  pinput->LoadFromFile(infile);
  infile.Close();

  // set up mesh
  int restart = false;
  int mesh_only = false;

  Mesh *pmesh = new Mesh(pinput, mesh_only);

  // initialize mesh
  pmesh->Initialize(restart, pinput);

  // load chemistry
  //std::string kinetics_file = pinput->GetString("chemistry", "jupiter_h2o.yaml");
  //auto sol = Cantera::newSolution(kinetics_file, "vapor");

  //Vapor *pvapor = new Vapor(pmesh->my_blocks(0), pinput);
  auto pvapor = new Vapor<Cantera::IdealGasPhase>(pinput);
  //auto ptracer = new Tracer<Cantera::IdealGasPhase>(pinput);

  std::cout << "Number of species = ";
  std::cout << pvapor->nSpecies() << std::endl;

  std::vector<Real> wghts(pvapor->nSpecies());
  pvapor->getMolecularWeights(wghts.data());

  for (int i = 0; i < pvapor->nSpecies(); ++i) {
    std::cout << pvapor->speciesName(i) << " ";
    std::cout << wghts[i] << " ";
    std::cout << pvapor->moleFraction(i) << std::endl;
  }
  std::cout << "Number of dry species = ";
  std::cout << pvapor->nDryMixtures() << std::endl;

  delete pinput;
  delete pmesh;
  delete pvapor;
}
