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


NAMESPACE_FEM_CLOSE

#endif // FEM_FUNCTION_H
