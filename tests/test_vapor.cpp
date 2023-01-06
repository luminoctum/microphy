#include <iostream>
#include <vapor.hpp>
#include <debugger.hpp>

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

  Vapor *pvapor = new Vapor(pmesh->my_blocks(0), pinput);

  std::cout << *pvapor << std::endl;

  for (int i = 1; i <= NumVapors; ++i)
    std::cout << pvapor->getMolecule(i) << std::endl;

  delete pinput;
  delete pmesh;
  delete pvapor;
}
