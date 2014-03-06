#ifndef FEM_FUNCTION_H
#define FEM_FUNCTION_H

class Point;

class Function
{
public:
  virtual double value(const Point &point,
                       const double time = 0) const = 0;
};


#endif // FEM_FUNCTION_H
