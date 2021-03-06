#ifndef FEM_FINITE_ELEMENT_H
#define FEM_FINITE_ELEMENT_H

#include "config.h"

NAMESPACE_FEM_OPEN



class FiniteElement
{
public:
  FiniteElement(unsigned int order);

  unsigned int order() const;

private: // ======================== PRIVATE =======================
            /**
             * Order of basis functions
             */
  unsigned int _order;

  FiniteElement(const FiniteElement &fe);
  FiniteElement& operator =(const FiniteElement &fe);
};


NAMESPACE_FEM_CLOSE

#endif // FEM_FINITE_ELEMENT_H
