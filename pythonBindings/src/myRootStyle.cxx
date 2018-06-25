#include <boost/python.hpp>

// #include "myRootStyle.h"
#include "src/myRootStyle.cxx"
using namespace boost::python;

BOOST_PYTHON_MODULE(myRootStyle)
{
    def("setMyRootStyle", &myFuncs::setMyRootStyle);
}
