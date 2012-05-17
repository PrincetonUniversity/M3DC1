#include "options.h"
#include <iostream>

fio_option_list::fio_option_list()
{
}

fio_option_list::~fio_option_list()
{
  clear();
}

int fio_option_list::clear()
{
  std::map<int, fio_option*>::iterator i = opts.begin();

  while(i != opts.end()) {
    delete(i->second);
    i++;
  }

  opts.clear();
  return 0;
}


int fio_option_list::set_option(const int i, const int v)
{
  std::map<int, fio_option*>::iterator o = opts.find(i);

  if(o==opts.end()) {
    std::cerr << "Error: option " << i << " not supported in this context" << std::endl;
    return 1;
  }
  
  if(typeid(*(o->second)) != typeid(fio_int_option)) {
    std::cerr << "Error: option " << i << " is not integer-valued" << std::endl;
    return 2;
  }

  ((fio_int_option*)o->second)->data = v;

  return 0;
}

int fio_option_list::set_option(const int i, const double v)
{
  std::map<int, fio_option*>::iterator o = opts.find(i);

  if(o==opts.end()) {
    std::cerr << "Error: option " << i << " not supported in this context" << std::endl;
    return 1;
  }
  
  if(typeid(*(o->second)) != typeid(fio_real_option)) {
    std::cerr << "Error: option " << i << " is not real-valued" << std::endl;
    return 2;
  }

  ((fio_real_option*)o->second)->data = v;

  return 0;
}

int fio_option_list::set_option(const int i, const char* v)
{
  std::map<int, fio_option*>::iterator o = opts.find(i);

  if(o==opts.end()) {
    std::cerr << "Error: option " << i << " not supported in this context" << std::endl;
    return 1;
  }
  
  if(typeid(*(o->second)) != typeid(fio_str_option)) {
    std::cerr << "Error: option " << i << " is not string-valued" << std::endl;
    return 2;
  }

  ((fio_str_option*)o->second)->data = v;

  return 0;
}
