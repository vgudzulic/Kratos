//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "includes/table.h"
#include "add_table_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

typedef Table<double,double> DoubleTableType;

double TableGetNearestValue(DoubleTableType& ThisTable, double X)
{
    return ThisTable.GetNearestValue(X);
}


void  AddTableToPython(pybind11::module& m)
{
    using namespace pybind11;
    
    class_<DoubleTableType, DoubleTableType::Pointer>(m,"PiecewiseLinearTable")
    .def(init<>())
    .def(init<Matrix const&>())
    .def("GetValue", &DoubleTableType::GetValue)
    .def("GetDerivative",&DoubleTableType::GetDerivative)
    .def("GetNearestValue", TableGetNearestValue)
    .def("AddRow", &DoubleTableType::PushBack)
    .def("__repr__", &DoubleTableType::Info)
    ;
}

}  // namespace Python.

} // Namespace Kratos

