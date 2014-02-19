#include "auxiliary_functions.h"
#include <sstream>
#include <stdexcept>


//-------------------------------------------------------
//
// d2s - convert data to string
//
//-------------------------------------------------------

std::string d2s(double x, bool scientific, int precision)
{
  std::ostringstream o;
  if (scientific)
  {
    o.setf(std::ios::scientific);
    o.precision(precision);
  }
  if (!(o << x))
    throw std::runtime_error("Bad conversion from double to string!");
  return o.str();
}


std::string d2s(int x)
{
  std::ostringstream o;
  if (!(o << x))
    throw std::runtime_error("Bad conversion from int to string!");
  return o.str();
}


std::string d2s(unsigned int x)
{
  std::ostringstream o;
  if (!(o << x))
    throw std::runtime_error("Bad conversion from unsigned int to string!");
  return o.str();
}


#if defined(HAVE_64BIT_SIZE_T)
std::string d2s(size_t x)
{
  std::ostringstream o;
  if (!(o << x))
    throw std::runtime_error("Bad conversion from size_t to string!");
  return o.str();
}
#endif



//-------------------------------------------------------
//
// expect and require
//
//-------------------------------------------------------

void requirement_fails(const char *file,
                       unsigned int line,
                       std::string message)
{
  std::string exc = "Exception:\nfile = " + std::string(file) +
                    "\nline = " + d2s(line) +
                    "\nmessage = " + message + "\n";
  throw std::runtime_error(exc);
}
