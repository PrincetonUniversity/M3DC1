#include "fusion_io.h"
#include <sstream>

fio_species::fio_species()
{
  nucleons = protons = electrons = 0;
}

fio_species::fio_species(const fio_species& s)
{
  nucleons = s.nucleons;
  protons = s.protons;
  electrons = s.electrons;
}

fio_species::fio_species(const int m, const int p, const int e)
{
  nucleons = m;
  protons = p;
  electrons = e;
}

fio_species::fio_species(const int i)
{
  int n = i / 65536; 
  protons = (i-n*65536) / 256;
  electrons = (i-n*65536-protons*256);
  nucleons = protons + n;
}

fio_species::operator int() const
{
  return electrons + protons*256 + (nucleons-protons)*65536;
}

int fio_species::charge() const
{
  return protons - electrons;
}

int fio_species::mass() const
{
  return nucleons;
}

bool fio_species::operator==(const fio_species& s) const
{
  return (s.nucleons == nucleons) &&
    (s.protons == protons) && 
    (s.electrons == electrons);
}

std::string fio_species::name() const
{
  std::ostringstream ss;

  if(protons==0) {
    if(electrons==1) return "e";
  } else if (protons==1) {
    switch(nucleons) {
    case(1): ss << "H"; break;
    case(2): ss << "D"; break;
    case(3): ss << "T"; break;
    }
  } else {
    if(protons >= 0 && protons < 19) {
      ss << species_name[protons];
    } else {
      ss << protons << "?";
    }
  }

  // append charge state
  ss << nucleons << "_" << charge();

  return ss.str();
}

const std::string fio_species::species_name[19] = 
  {"e", 
   "H", "He", 
   "Li", "Be",  "B",  "C", "N", "O",  "F", "Ne",
   "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"};
