// C/C++ headers
#include <iostream>
#include <fstream>

// debugger
#include <debugger.hpp>

// microphy headers
#include <microphy/vapor.hpp>
#include <microphy/suspension.hpp>
#include <microphy/construct_atmosphere.hpp>
//#include <condensates.hpp>

// cantera headers
#include <cantera/thermo.h>
#include <cantera/kinetics.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/thermo/IdealSolidSolnPhase.h>
//#include <cantera/base/Solution.h>

// Athena++ headers
#include <mesh/mesh.hpp>
#include <hydro/hydro.hpp>
#include <coordinates/coordinates.hpp>

// climath headers
extern "C" {
  #include <interpolation.h>
}

// cliutils
#include <cliutils/ndarrays.hpp>
#include <cliutils/stride_iterator.hpp>

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
  auto kincond = Cantera::newKinetics({surface, vapor, cloud}, kinfile, "vapor:cloud");

  vapor->setState_TP(200., 1.E5);
  vapor->updateThermo();


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

  vapor->setLatentHeats(kincond.get());

  for (int i = 1; i <= NumVapors; ++i)
    std::cout << vapor->testSaturation(i) << std::endl;

  // construct adiabatic atmosphere
  // mesh limits
  Real x1min = pmesh->mesh_size.x1min;
  Real x1max = pmesh->mesh_size.x1max;
  Real H0 = pinput->GetReal("mesh", "pscale_height");

  Real dlnp, **w1, *z1;
  dlnp = (x1max - x1min)/pmesh->mesh_size.nx1/(2.*H0);
  int nx1 = (int)((x1max - x1min)/(H0*dlnp));
  NewCArray(w1, nx1, NumHydros);
  std::fill(w1[0], w1[0] + nx1*NumHydros, 0.);
  z1 = new Real [nx1];
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + H0*dlnp;

  std::ofstream outfile("adiabatic_atm.temp-vapors.out");

  // read in parameters
  Real Ts = pinput->GetReal("problem", "Ts");
  Real Ps = pinput->GetReal("problem", "Ps");

  vapor->setState_TP(Ts, Ps);

  AtmConstructionMethod method;
  method.type = "pseudo";
  method.grav = -pmesh->my_blocks(0)->phydro->hsrc.GetG1();
  method.dlnp = -dlnp;

  for (int b = 0; b < pmesh->nblocal; ++b) {
    MeshBlock *pmb = pmesh->my_blocks(b);
    auto phydro = pmb->phydro;
    auto pcoord = pmb->pcoord;

    std::cout << "I'm now constructing atmosphere" << std::endl << std::endl;
    ConstructAtmosphere(w1, nx1, method, vapor, cloud, kincond.get());

    // Change to log quantity
    for (int i = 0; i < nx1; ++i)
      for (int n = 0; n < NumHydros; ++n) {
        if (n == iv1 || n == iv2 || n == iv3)
          w1[i][n] = 0.;
        else
          w1[i][n] = log(w1[i][n]);
      }

    // interpolate to model grid
    int is = pmb->is, js = pmb->js, ks = pmb->ks;
    int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real buf[NumHydros];

          // set gas concentration and velocity
          interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NumHydros);
          for (int n = 0; n < NumHydros; ++n)
            for (int k = ks; k <= ke; ++k)
              for (int j = js; j <= je; ++j) {
                if (n == iv1 || n == iv2 || n == iv3)
                  phydro->w(n,k,j,i) = 0.;
                else {
                  phydro->w(n,k,j,i) = exp(buf[n]);
                }
              }
        }

    // print temperature and vapor
    int iH2O = vapor->vaporIndex("H2O");
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          outfile << vapor->getTemp(phydro->w.at(k,j,i)) << " "
                  << phydro->w(iH2O,k,j,i) << std::endl;
        }
  }

  FreeCArray(w1);
  delete[] z1;

  delete pinput;
  delete pmesh;
  delete surface;
  delete vapor;
  delete cloud;
}
