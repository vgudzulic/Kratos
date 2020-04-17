// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//         -        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//         -        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//                  in the documentation and/or other materials provided with the distribution.
//         -        All advertising materials mentioning features or use of this software must display the following acknowledgement:
//                         This product includes Kratos Multi-Physics technology.
//         -        Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "includes/define.h"
#include "fem_to_dem_application.h"
#include "fem_to_dem_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

  using namespace pybind11;



  PYBIND11_MODULE(KratosFemToDemApplication, m)
  {


    class_<KratosFemToDemApplication, KratosFemToDemApplication::Pointer, KratosApplication>(m, "KratosFemToDemApplication")
        .def(init<>());
        ;

    AddCustomStrategiesToPython (m);
    AddCustomUtilitiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,BACKUP_LAST_STRUCTURAL_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,BACKUP_LAST_STRUCTURAL_DISPLACEMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,SMOOTHED_STRUCTURAL_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,OLD_RELAXED_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,RELAXED_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,FSI_INTERFACE_RESIDUAL)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,ACCELERATION_BACKUP);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,DISPLACEMENT_BACKUP);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HARDENING_MODULUS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_INITIAL_VOLUME);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_VOLUME);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FRAGILE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VOLUME_COUNTED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ERASED_VOLUME);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,COHESION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,RECOMPUTE_NEIGHBOURS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DEMFEM_CONTACT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,GENERATE_DEM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DAMAGE_ELEMENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TIME_UNIT_CONVERTER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISPLACEMENT_INCREMENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YIELD_STRESS_C);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YIELD_STRESS_T);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FRAC_ENERGY_T)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FRAC_ENERGY_C)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,INTERNAL_PRESSURE_ITERATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PFEM_PRESSURE_ITERATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YIELD_SURFACE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRESS_VECTOR_INTEGRATED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TANGENT_CONSTITUTIVE_TENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SMOOTHING);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,IS_DAMAGED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,RECONSTRUCT_PRESSURE_LOAD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,IS_DYNAMIC);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRESS_THRESHOLD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,INITIAL_THRESHOLD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,INTEGRATION_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MAPPING_PROCEDURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,IS_DEM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,IS_SKIN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DEM_RADIUS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DEM_GENERATED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,INACTIVE_NODE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NUMBER_OF_ACTIVE_ELEMENTS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_FORCE_APPLIED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_FORCE_X);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_FORCE_Y);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_FORCE_Z);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_DAMAGE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_EXPANDED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EQUIVALENT_NODAL_STRESS);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EQUIVALENT_NODAL_STRESS_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,AUXILIAR_GRADIENT);

    // 3D case
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRAIN_TENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRESS_TENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STRESS_TENSOR_INTEGRATED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HENCKY_STRAIN_VECTOR);
    
    // Composite calculations
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MATRIX_STRESS_TENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FIBER_STRESS_TENSOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MATRIX_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FIBER_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YOUNG_MODULUS_FIBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DENSITY_FIBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,POISSON_RATIO_FIBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FIBER_VOLUMETRIC_PART);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MATRIX_STRESS_TENSOR_INTEGRATED);
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YIELD_STRESS_C_FIBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YIELD_STRESS_T_FIBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FRACTURE_ENERGY_FIBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ACUMULATED_PLASTIC_STRAIN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EQUIVALENT_STRESS_VM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HARDENING_LAW);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MAXIMUM_STRESS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MAXIMUM_STRESS_POSITION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,IS_TAKEN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_ID);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PLASTIC_UNIAXIAL_STRESS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PLASTIC_STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,LINE_LOAD);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,SURFACE_LOAD);
  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
