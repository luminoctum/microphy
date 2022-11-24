#ifndef CREATE_PHASES_HPP
#define CREATE_PHASES_HPP

#include "molecules.hpp"

template<template<int> class T, int N = 0>
void create_phases(Molecule *species, int i)
{
  species[i] = T<N>();
  create_phases<T, N+1>(species, i + NumVapors);
}

template<>
void create_phases<Hydrogen,1>(Molecule*,int) {}

template<>
void create_phases<Helium,1>(Molecule*,int) {}

template<>
void create_phases<Water,NumPhases>(Molecule*,int) {}

template<>
void create_phases<Ammonia,NumPhases>(Molecule*,int) {}

template<>
void create_phases<Methane,NumPhases>(Molecule*,int) {}

template<>
void create_phases<HydrogenSulfide,NumPhases>(Molecule*,int) {}

#endif
