#ifndef FUSION_IO_SPECIES_H
#define FUSION_IO_SPECIES_H

#include <string>

class fio_species {
  static const std::string species_name[19];
  int nucleons, protons, electrons;

 public:
  fio_species();
  fio_species(const fio_species&);
  fio_species(const int);
  fio_species(const int, const int, const int);
  int charge() const;
  int mass() const;
  std::string name() const;
  bool operator==(const fio_species&) const;

  operator int() const;
};

static fio_species fio_electron(0, 0, 1);
static fio_species fio_proton(1, 1, 0);
static fio_species fio_alpha(4, 2, 0);

#endif
