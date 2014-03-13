#ifndef SYM_CSR_PATTERN_H
#define SYM_CSR_PATTERN_H

#include "config.h"
#include "csr_pattern.h"

NAMESPACE_FEM_OPEN



class DoFHandler;
class FineMesh;


/**
 * Symmetric compressed sparse row (CSR) pattern
 */
class SymCSRPattern : public CSRPattern
{
public:
            /**
             * Constructor
             */
  SymCSRPattern();

            /**
             * Destructor
             */
  ~SymCSRPattern();

  void make_sparse_format(const DoFHandler &dof_handler);

  void make_sparse_format(const FineMesh &fmesh);

  void make_sparse_format(const std::vector<Triangle> &triangles, unsigned int order, CONNECTIVITY connectivity);

  int find(unsigned int num_row, unsigned num_col) const;


private:
            /**
             * Copy constructor. It's private, since we don't want to copy SymCSRPattern objects
             */
  SymCSRPattern(const SymCSRPattern&);

            /**
             * The same is valid as for copy constructor
             */
  SymCSRPattern& operator=(const SymCSRPattern&);

};


NAMESPACE_FEM_CLOSE

#endif // SYM_CSR_PATTERN_H
