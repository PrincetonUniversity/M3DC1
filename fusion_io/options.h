#ifndef FIO_OPTIONS_H
#define FIO_OPTIONS_H

#include <map>
#include <string>

struct fio_option_list : public std::map<int, void*> {
  fio_option_list();
  virtual ~fio_option_list();
  int clear();

  template<typename T>
  int add_option(const int, const T&);

  template<typename T>
  int set_option(const int, const T&);

  template<typename T>
  int get_option(const int, T*) const;

  template<typename T>
  static bool is_type(const int);
};


#endif
