#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/python.hpp>
#pragma GCC diagnostic pop

#include "include/histFuncs.h"
#include "src/histFuncs.cxx"
using namespace boost::python;


BOOST_PYTHON_MODULE(histFuncs)
{
    def("makeCanvas", &myFuncs::makeCanvas, return_value_policy<manage_new_object>());
}
