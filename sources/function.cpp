#include "function.h"

NAMESPACE_FEM_OPEN


ConstantFunction::ConstantFunction(double value)
  : _value(value)
{ }



ConstantFunction::ConstantFunction(const ConstantFunction &cf)
  : _value(cf._value)
{ }



ConstantFunction& ConstantFunction::operator =(const ConstantFunction &cf)
{
  _value = cf._value;
  return *this;
}



double ConstantFunction::value(const Point &point, const double time) const
{
  return _value;
}


NAMESPACE_FEM_CLOSE
