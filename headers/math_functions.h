#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "petscvec.h"
#include "point.h"


/**
 * Relative error in 2-norm between two vectors
 * rel_error = norm2(vec1-vec2) / norm2(vec2)
 */
double rel_error(const Vec &vec1, const Vec &vec2);


/**
 * Dot (inner) product of two 3D vectors
 * @param p1 - one vector (represented as object of class Point)
 * @param p2 - another vector (represented as object of class Point)
 */
double dot_product(const Point &p1, const Point &p2);


/**
 * Vector normalization. The vector is changed
 * @param vec - vector (represented as object of class Point)
 * @return
 */
void normalize(Point &vec);



/**
 * Norm of a vector represented by a Point object
 */
double norm(const Point &vec);



#endif // MATH_FUNCTIONS_H
