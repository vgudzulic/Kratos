//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#if !defined(KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"
#include "includes/checks.h"

namespace Kratos
{
// Degrees of freedom
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, VELOCITY_POTENTIAL)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, AUXILIARY_VELOCITY_POTENTIAL)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, PSI)

//Embedded variables
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, GEOMETRY_DISTANCE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, ROTATION_ANGLE)

//Wake variables
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, WAKE_DISTANCE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, Vector, WAKE_ELEMENTAL_DISTANCES)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, Vector, WAKE_ORIGIN)

// Adjoint variables
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, ADJOINT_VELOCITY_POTENTIAL)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, ADJOINT_AUXILIARY_VELOCITY_POTENTIAL)

// Flow field magnitudes
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, VELOCITY_LOWER)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, PRESSURE_LOWER)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, POTENTIAL_JUMP)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, ENERGY_NORM_REFERENCE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, POTENTIAL_ENERGY_REFERENCE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, HEAT_CAPACITY_RATIO)

// Free stream magnitudes
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, FREE_STREAM_VELOCITY)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, FREE_STREAM_DENSITY)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, FREE_STREAM_MACH)

// Integral magnitudes
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, LIFT_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, DRAG_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, MOMENT_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, LIFT_COEFFICIENT_JUMP)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, LIFT_COEFFICIENT_FAR_FIELD)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, DRAG_COEFFICIENT_FAR_FIELD)

// Geometrical variables
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, double, REFERENCE_CHORD)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, WAKE_NORMAL)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, WING_SPAN_DIRECTION)

// Markers
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, WAKE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, KUTTA)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, WING_TIP)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, bool, TRAILING_EDGE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, bool, UPPER_SURFACE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, bool, LOWER_SURFACE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, bool, UPPER_WAKE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, bool, LOWER_WAKE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, AIRFOIL)

// To be removed
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, TRAILING_EDGE_ELEMENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, DECOUPLED_TRAILING_EDGE_ELEMENT)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, DEACTIVATED_WAKE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, ALL_TRAILING_EDGE)
KRATOS_DEFINE_APPLICATION_VARIABLE(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION, int, ZERO_VELOCITY_CONDITION)
} // namespace Kratos

#endif /* KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_VARIABLES_H_INCLUDED */
