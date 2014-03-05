#include "finite_element.h"



FiniteElement::FiniteElement(unsigned int order)
  : _order(order)
{ }



unsigned int FiniteElement::order() const
{
  return _order;
}
