#ifndef RESULT_H
#define RESULT_H

#include "config.h"
#include "petscvec.h"
#include <string>
#include <vector>

NAMESPACE_FEM_OPEN



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
  Result(const DoFHandler *dof_handler);

            /**
             * Destructor
             */
  ~Result();

            /**
             * Write the results to file in vtk(vtu) unstructured format to work then in Paraview
             * @param filename - the name of output file
             */
  void write_vtu(const std::string &filename,
                 const Vec &solution,
                 const Vec &exact_solution = 0) const;

            /**
             * Write the results to file in vtk(vts) structured format to work then in Paraview
             * @param filename - the name of output file
             */
  void write_vts(const std::string &filename,
                 unsigned int N_FINE_X,
                 unsigned int N_FINE_Y,
                 const Vec &solution,
                 const Vec &exact_solution = 0,
                 const std::vector<double> &coef_alpha = std::vector<double>(),
                 const std::vector<double> &coef_beta = std::vector<double>()) const;


private: //========================= PRIVATE ===================
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


NAMESPACE_FEM_CLOSE

#endif // RESULT_H
