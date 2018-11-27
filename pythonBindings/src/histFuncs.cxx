#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/python.hpp>
#pragma GCC diagnostic pop

#include "TLegend.h"

#include "include/histFuncs.h"
#include "src/histFuncs.cxx"
using namespace boost::python;

// int f1(int x) { return f(x); }
// TLegendEntry* (myFuncs::Legend::*AddEntry1)(const TObject*) = &myFuncs::Legend::AddEntry;

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(addEntryOverloads, myFuncs::Legend::AddEntry, 1, 3)

BOOST_PYTHON_MODULE(histFuncs) {
  def("makeCanvas", &myFuncs::makeCanvas, return_value_policy<manage_new_object>());

  class_<myFuncs::Legend>("Legend", init<double, double, double, double>())
      // .def("AddEntry", &myFuncs::Legend::AddEntry, addEntryOverloads(args("obj", "label=", "Option_t")))
      .def("Draw", &myFuncs::Legend::Draw);
}
