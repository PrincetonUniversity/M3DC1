#include "fusion_io.h"
#include <iostream>
#include <typeinfo>

fio_option_list::fio_option_list()
{
}

fio_option_list::~fio_option_list()
{
  clear();
}

template<> bool fio_option_list::is_type<int>(const int i)
{
  return FIO_IS_INT_OPT(i);
}
template<> bool fio_option_list::is_type<double>(const int i)
{
  return FIO_IS_REAL_OPT(i);
}
template<> bool fio_option_list::is_type<std::string>(const int i)
{
  return FIO_IS_STR_OPT(i);
}


int fio_option_list::clear()
{
  iterator i = begin();

  while(i != end()) {
    if(FIO_IS_INT_OPT(i->first)) {
      delete((int*)i->second);
    } else if(FIO_IS_REAL_OPT(i->first)) {
      delete((double*)i->second);
    } else if(FIO_IS_STR_OPT(i->first)) {
      delete((std::string*)i->second);
    } else {
      std::cerr << "Error: encountered option of unknown type" << std::endl;
    }
    i++;
  }

  std::map<int, void*>::clear();
  return FIO_SUCCESS;
}

template<typename T>
int fio_option_list::add_option(const int i, const T& v)
{
  iterator o = find(i);

  if(o!=end()) {
    std::cerr << "Error: option " << i << " already present" << std::endl;
    return 1;
  }
  
  if(!is_type<T>(i)) {
    std::cerr << "Error: option " << i << " is not of type" 
	      << typeid(T).name() << std::endl;
    return 2;
  }

  (*this)[i] = new T(v);

  return FIO_SUCCESS;
}

template<typename T>
int fio_option_list::set_option(const int i, const T& v)
{
  iterator o = find(i);

  if(o==end()) {
    std::cerr << "Error: option " << i 
	      << " not supported in this context" << std::endl;
    return 1;
  }
  
  if(!is_type<T>(i)) {
    std::cerr << "Error: option " << i << " is not of type"
	      << typeid(T).name() << std::endl;
    return 2;
  }

  *((T*)(o->second)) = v;

  return FIO_SUCCESS;
}

template<typename T>
int fio_option_list::get_option(const int i, T* v) const
{
  const_iterator o = find(i);

  if(o==end()) {
    std::cerr << "Error: option " << i 
	      << " not supported in this context" << std::endl;
    return 1;
  }
  
  if(!is_type<T>(i)) {
    std::cerr << "Error: option " << i << " is not of type "
	      << typeid(T).name() << std::endl;
    return 2;
  }

  *v = *((T*)(o->second));

  return FIO_SUCCESS;
}

template int fio_option_list::get_option<int>(const int, int*) const;
template int fio_option_list::get_option<double>(const int, double*) const;
template int fio_option_list::get_option<std::string>(const int, 
						      std::string*) const;
template int fio_option_list::set_option<int>(const int, const int&);
template int fio_option_list::set_option<double>(const int, const double&);
template int fio_option_list::set_option<std::string>(const int, 
						      const std::string&);
template int fio_option_list::add_option<int>(const int, const int&);
template int fio_option_list::add_option<double>(const int, const double&);
template int fio_option_list::add_option<std::string>(const int, 
						      const std::string&);
