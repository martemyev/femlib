#ifndef FEM_CONFIG_H
#define FEM_CONFIG_H

#include <string>

// options defined in cmake
/* #undef DEBUG */
/* #undef HAVE_64BIT_SIZE_T */

// the following two definitions are used
// to wrap the code in namespace 'fem'.
// thus all classes will be used like 'fem::NameOfClass'
#define NAMESPACE_FEM_OPEN namespace fem {
#define NAMESPACE_FEM_CLOSE }

// constants from cmake
const std::string HOME_DIRECTORY = "/home/artemiev";

#endif // FEM_CONFIG_H
