#include "finite_element.h"

NAMESPACE_FEM_OPEN



FiniteElement::FiniteElement(unsigned int order)
  : _order(order)
{ }



unsigned int FiniteElement::order() const
{
  return _order;
}


NAMESPACE_FEM_CLOSE
