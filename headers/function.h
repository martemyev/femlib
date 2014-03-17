#ifndef FEM_FUNCTION_H
#define FEM_FUNCTION_H

#include "config.h"

NAMESPACE_FEM_OPEN


class Point;

class Function
{
public:
  virtual double value(const Point &point,
                       const double time = 0) const = 0;
};



// =======================================================

class ConstantFunction : public Function
{
public:
  ConstantFunction(double value);
  ConstantFunction(const ConstantFunction &cf);
  ConstantFunction& operator =(const ConstantFunction &cf);

  virtual double value(const Point &point,
                       const double time = 0) const;

private:
  double _value;
};



NAMESPACE_FEM_CLOSE

#endif // FEM_FUNCTION_H
