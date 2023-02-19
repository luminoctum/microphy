#include <iostream>
#include <debugger.hpp>
#include <mesh/mesh.hpp>
#include <vapor.hpp>
#include <suspension.hpp>
#include <condensates.hpp>
#include <cantera/thermo.h>
#include <cantera/kinetics.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/thermo/IdealSolidSolnPhase.h>
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

  //auto pcloud = new Condensates<1>(pmesh->my_blocks(0), pinput);

  // kinetics manager
  std::string kinfile = pinput->GetString("chemistry", "kinetics_file");
  auto surface = Cantera::newPhase(kinfile, "vapor:cloud");
  auto vapor = new Vapor<Cantera::IdealGasPhase>(pinput);
  auto cloud = new Suspension<Cantera::IdealGasPhase,
    Cantera::IdealSolidSolnPhase>(pinput, vapor, "cloud");
  //auto cloud = Cantera::newPhase(kinfile, "cloud");
  auto kincond = Cantera::newKinetics({surface, vapor, cloud}, kinfile, "vapor:cloud");

  vapor->setState_TP(200., 1.E5);

  // set small amount of trial condensates
  Real ym[2] = {1., 2.};
  cloud->setPressure(1.E5);
  cloud->setRelativeMole(1.E-6);
  cloud->setMassFractions(ym);

  std::cout << "Total species = "
            << kincond->nTotalSpecies() << std::endl;
  std::cout << "Number of Phases = " << kincond->nPhases() << std::endl;
  std::cout << "Species names = ";
  for (int i = 0; i < kincond->nTotalSpecies(); ++i)
    std::cout << kincond->kineticsSpeciesName(i) << " ";
  std::cout << std::endl;

  std::cout << "Number of kinetics reactions = "
            << kincond->nReactions() << std::endl;
  for (int i = 0; i < kincond->nReactions(); ++i) {
    std::cout << kincond->reactionString(i) << std::endl;
  }

  vapor->setSaturationIndicator(*kincond);

  for (int i = 1; i <= NumVapors; ++i)
    std::cout << vapor->testSaturation(i) << std::endl;

  delete pinput;
  delete pmesh;
  delete surface;
  delete vapor;
  delete cloud;
}
