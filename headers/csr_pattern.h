#ifndef SPARSE_FORMAT_H
#define SPARSE_FORMAT_H

#include <vector>
#include <set>
#include "dof_handler.h"



enum CONNECTIVITY
{
  VERTICES, // the sparse pattern is based on connectivity between vertices
  DOFS      // the sparse pattern is based on connectivity between degrees of freedom
};



/**
 * Compressed sparse row (CSR) pattern
 */
class CSRPattern
{
public:
            /**
             * Constructor
             */
  CSRPattern();

            /**
             * Destructor
             */
  virtual ~CSRPattern();

            /**
             * Make sparse format
             * @param dof_handler - handler of degrees of freedom
             */
  void make_sparse_format(const DoFHandler &dof_handler, COUPLING coupling);

            /**
             * Get the serial number of nonzero element in the pattern
             * @param num_row - the number of the row
             * @param num_col - the number of the col
             */
  int find(unsigned int num_row, unsigned num_col) const;

            /**
             * Get matrix order
             */
  unsigned int order() const;

            /**
             * Get an element from row array
             * @param number - the serial number of the element
             */
  unsigned int row(unsigned int number) const;

            /**
             * Get an element from col array
             * @param number - the serial number of the element
             */
  unsigned int col(unsigned int number) const;

  const int* nnz() const;


protected: // ==================== PROTECTED ===================
            /**
             * Matrix order (matrix is square).
             */
  unsigned int _order;
            /**
             * The number of nonzero elements in each row. Starts from 0.
             * Its dimension = the number of matrix rows + 1.
             * Its last element means the number of nonzero elements in the matrix.
             */
  std::vector<unsigned int> _row;

            /**
             * The numbers of columns where nonzero elements take place.
             * Its dimension is equal to the amount of nonzero elements in the matrix.
             */
  std::vector<unsigned int> _col;

            /**
             * To free some memory
             */
  void clear();

  void make_cg_sparse_format(const DoFHandler &dof_handler);
  void make_dg_sparse_format(const DoFHandler &dof_handler);

            /**
             * Initialization of the pattern arrays
             * @param connect - connection between dofs
             */
  void pattern_initialization(const std::set<unsigned int> *connect);


private: // ============================= PRIVATE ==============================

            /**
             * Copy constructor. It's private, since we don't want to copy CSRPattern objects
             */
  CSRPattern(const CSRPattern&);

            /**
             * The same is valid as for copy constructor
             */
  CSRPattern& operator=(const CSRPattern&);
};

#endif // SPARSE_FORMAT_H
