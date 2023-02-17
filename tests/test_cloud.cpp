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

  auto pvapor = new Vapor<Cantera::IdealGasPhase>(pinput);
  auto pcloud = new Condensates<1>(pmesh->my_blocks(0), pinput);
  //auto ptracer = new Tracer<Cantera::IdealGasPhase>(pinput);

  std::cout << "Scalar indices = ";
  std::cout << pcloud->scalarIndexFirst() << " - ";
  std::cout << pcloud->scalarIndexLast() << std::endl;

  std::cout << "Number of solution groups = ";
  std::cout << pcloud->nSolutions() << std::endl;

  for (int i = 0; i < pcloud->nSolutions(); ++i) {
    auto thermo = pcloud->solution(i)->thermo();
    std::cout << thermo->name() << std::endl;;
    for (int j = 0; j < thermo->nSpecies(); ++j) {
      std::cout << thermo->speciesName(j) << std::endl;
    }
    std::cout << std::endl;
    //std::cout << wghts[i] << " ";
    //std::cout << pcloud->moleFraction(i) << std::endl;
  }

  std::cout << "cloud species start at = ";
  std::cout << pcloud->speciesScalarIndexStart("cloud") << std::endl;

  std::cout << "precipitation species start at = ";
  std::cout << pcloud->speciesScalarIndexStart("precipitation") << std::endl;

  std::cout << "Index H2O(l) = ";
  std::cout << pcloud->speciesScalarIndex("H2O(l)") << std::endl;

  std::cout << "Index H2O(s) = ";
  std::cout << pcloud->speciesScalarIndex("H2O(s)") << std::endl;

  std::cout << "Index H2O(l)-rain = ";
  std::cout << pcloud->speciesScalarIndex("H2O(l)-rain") << std::endl;

  std::cout << "Index H2O(s)-snow = ";
  std::cout << pcloud->speciesScalarIndex("H2O(s)-snow") << std::endl;

  std::cout << "Index H2O(l)-snow = ";
  std::cout << pcloud->speciesScalarIndex("H2O(l)-snow") << std::endl;

  /*std::vector<Real> wghts(pcloud->nSpecies());
  pcloud->getMolecularWeights(wghts.data());

  for (int i = 0; i < pcloud->nSpecies(); ++i) {
    std::cout << pcloud->speciesName(i) << " ";
    std::cout << wghts[i] << " ";
    std::cout << pcloud->moleFraction(i) << std::endl;
  }*/

  delete pinput;
  delete pmesh;
  delete pvapor;
  delete pcloud;
}
