#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/python.hpp>
#pragma GCC diagnostic pop

#include "src/myRootStyle.cxx"
using namespace boost::python;


BOOST_PYTHON_MODULE(myRootStyle)
{
    def("setMyRootStyle", &myFuncs::setMyRootStyle);
}
