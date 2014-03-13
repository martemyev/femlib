#include "math_functions.h"
#include "auxiliary_functions.h"

NAMESPACE_FEM_OPEN



double rel_error(const Vec &vec1, const Vec &vec2)
{
  Vec diff;
  VecDuplicate(vec1, &diff); // allocate the memory for diff
  VecCopy(vec1, diff); // copy vec1 to diff
  VecAXPY(diff, -1, vec2); // diff = vec1 - vec2
  double norm_diff, norm_vec2;
  VecNorm(diff, NORM_2, &norm_diff);
  VecNorm(vec2, NORM_2, &norm_vec2);
  VecDestroy(&diff); // free the memory
  return norm_diff / norm_vec2;
}



double dot_product(const Point &p1, const Point &p2)
{
  double dp = 0.;
  for (unsigned int i = 0; i < Point::n_coord; ++i)
    dp += p1.coord(i) * p2.coord(i);
  return dp;
}



void normalize(Point &vec)
{
  const double norm = sqrt(dot_product(vec, vec)); // norm of the vector
  expect(norm > 1e-14, "Vector has very small length (" + d2s(norm) + "), so it's dangerous to normalize it");
  vec /= norm;
}



double norm(const Point &vec)
{
  return sqrt(dot_product(vec, vec));
}


NAMESPACE_FEM_CLOSE
