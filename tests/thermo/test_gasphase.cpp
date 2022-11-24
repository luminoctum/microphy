#include <iostream>
#include <gasphase/gasphase.hpp>
#include <debugger.hpp>

int main() {
  // no MPI
  Globals::my_rank = 0;
  Globals::nranks  = 1;

  pdebug = std::make_unique<Debugger>(1);

  IOWrapper infile;
  infile.Open("test_gasphase.inp", IOWrapper::FileMode::read);

  ParameterInput *pinput = new ParameterInput();
  pinput->LoadFromFile(infile);
  infile.Close();

  // set up mesh
  int restart = false;
  int mesh_only = false;

  Mesh *pmesh = new Mesh(pinput, mesh_only);

  // initialize mesh
  pmesh->Initialize(restart, pinput);

  GasPhase *pthermo = new GasPhase(pmesh->my_blocks(0), pinput);

  std::cout << *pthermo << std::endl;

  for (int i = 1; i <= NumVapors; ++i)
    std::cout << pthermo->getMoist(i) << std::endl;

  delete pinput;
  delete pmesh;
  delete pthermo;
}
