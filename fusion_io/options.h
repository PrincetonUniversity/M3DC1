#ifndef FIO_OPTIONS_H
#define FIO_OPTIONS_H

#include <map>
#include <string>

struct fio_option {
  std::string name;
};

struct fio_int_option : public fio_option {
  int data;
  fio_int_option(const char* n, int d)
    : data(d)
  { name = n; }
};

struct fio_real_option : public fio_option {
  double data;
  fio_real_option(const char* n, double d)
    : data(d)
  { name = n; }
};

struct fio_str_option : public fio_option {
  std::string data;
  fio_str_option(const char* n, const char* d)
    : data(d)
  { name = n; }
};

struct fio_option_list {
  std::map<int, fio_option*> opts;

  fio_option_list();
  ~fio_option_list();
  int clear();
  int set_option(const int, const int);
  int set_option(const int, const double);
  int set_option(const int, const char*);
};


#endif
