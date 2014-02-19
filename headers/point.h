#ifndef POINT_H
#define POINT_H


/**
 * Point in 3-dimensional space.
 * It's used for 2-dimensional triangles as well, and
 * in this case one of the coordinates is 0 (usually it's z-coordinate).
 */
class Point
{
public:
            /**
             * The number of Cartesian coordinates, that describe the point.
             * Here we always use 3 coordinates to describe a point.
             */
  static const unsigned int n_coord = 3;

            /**
             * Default constructor.
             * Coordinates are initialized by 0.
             */
  Point();

            /**
             * Constructor with parameter.
             * Coordinates are initialized by array of numbers.
             * @param coordinates - array of point coordinates
             */
  Point(const double coordinates[]);

            /**
             * Constructor with parameters.
             * Coordinates are initialized by numbers.
             * @param x_coord - x-coordinate of the point
             * @param y_coord - y-coordinate of the point
             * @param z_coord - z-coordinate of the point
             */
  Point(const double x_coord,
        const double y_coord = 0,
        const double z_coord = 0);

            /**
             * Copy constructor
             */
  Point(const Point &p);

            /**
             * Copy assignment operator
             */
  Point& operator =(const Point &p);

  Point& operator /=(double d);

  friend Point operator -(const Point &p1, const Point &p2);

            /**
             * Get the coordinate of the point
             * @param number - the serial number of coordinate [0, n_coord)
             */
  double coord(unsigned int number) const;

            /**
             * Set the value of specific coordinate
             * @param number - the number of coordinate that we want to set
             * @param value - new value of coordinate
             */
  void coord(unsigned int number, double value);

private:
            /**
             * Cartesian coordinates of the point
             */
  double _coord[n_coord];
};



#endif // POINT_H
