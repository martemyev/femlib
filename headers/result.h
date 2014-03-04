#ifndef RESULT_H
#define RESULT_H

#include "petscvec.h"
#include <string>

class Parameters;
class DoFHandler;


/**
 * This class is designed to provide a user with tools for export any results in appropriate formats
 */
class Result
{
public:
            /**
             * Constructor
             * @param dof_handler - pointer to the handler of degrees of freedom
             */
  Result(const Parameters *param, const DoFHandler *dof_handler);

            /**
             * Destructor
             */
  ~Result();

            /**
             * Write the results to file in vtk(vtu) format to work then in Paraview
             * @param filename - the name of output file
             */
  void write_vtu(const std::string &filename, const Vec &solution, const Vec &exact_solution = 0) const;
  void write_vts(const std::string &filename, const Vec &solution, const Vec &exact_solution = 0) const;


private: //========================= PRIVATE ===================
            /**
             * Parameters of the problem
             */
  const Parameters *_param;

            /**
             * Constant pointer to dof_handler
             */
  const DoFHandler *_dof_handler;

            /**
             * Copy constructor is private to prevent copying of an object of the Result class
             */
  Result(const Result&);

            /**
             * Copy assignment operator is private for the same reason as copy constructor is.
             */
  Result& operator=(const Result&);
};


#endif // RESULT_H
