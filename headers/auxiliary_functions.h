#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

#include "config.h"
#include <string>


//-------------------------------------------------------
//
// d2s - convert data to string
//
//-------------------------------------------------------
            /**
             * Convert float number to string
             * @param x - the number in double format
             * @param scientific - use scientific format (e.g. 1e-10), or not
             * @param precision - if scientific format is used, we can change the precision
             * @return data in string format
             */
std::string d2s(double x, bool scientific = false, int precision = 6);


            /**
             * Convert integer data to string
             * @param x - the integer number
             * @return data in string format
             */
std::string d2s(int x);


            /**
             * Convert unsigned integer data to string
             * @param x - unsigned integer number
             * @return data in string format
             */
std::string d2s(unsigned int x);


#if defined(HAVE_64BIT_SIZE_T)
            /**
             * Convert size_t data to string
             * @param x - size_t number
             * @return data in string format
             */
std::string d2s(size_t x);
#endif



//-------------------------------------------------------
//
// expect and require
//
//-------------------------------------------------------

#if defined(DEBUG)
  #define expect(condition, message)  \
    if (!(condition))                 \
      requirement_fails(__FILE__,     \
                        __LINE__,     \
                        message)
#else
  // in release (or release-like) versions
  // nothing happens
  #define expect(condition, message) { }
#endif // DEBUG

#define require(condition, message) \
  if (!(condition))                 \
    requirement_fails(__FILE__,     \
                      __LINE__,     \
                      message)


            /**
             * Throw an informative exception,
             * if requirement or expectation fails
             */
void requirement_fails(const char *file,
                       unsigned int line,
                       std::string message);


#endif // AUXILIARY_FUNCTIONS_H
