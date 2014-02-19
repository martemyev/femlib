#ifndef SYM_CSR_PATTERN_H
#define SYM_CSR_PATTERN_H

#include "csr_pattern.h"
//#include "dof_handler.h"
//#include "fine_mesh.h"

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



#endif // SYM_CSR_PATTERN_H
