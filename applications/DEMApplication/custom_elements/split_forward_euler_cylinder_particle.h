//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//
//

#if !defined(KRATOS_SPLIT_FORWARD_EULER_CYLINDER_PARTICLE_H_INCLUDED )
#define  KRATOS_SPLIT_FORWARD_EULER_CYLINDER_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "split_forward_euler_spheric_particle.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_conditions/RigidFace.h"
#include "../custom_conditions/dem_wall.h"
#include "../custom_strategies/schemes/dem_integration_scheme.h"
#include "includes/kratos_export_api.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

class DEMWall;

class KRATOS_API(DEM_APPLICATION) SplitForwardEulerCylinderParticle : public SplitForwardEulerSphericParticle
{

public:

/// Pointer definition of SplitForwardEulerCylinderParticle
KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SplitForwardEulerCylinderParticle);

typedef GlobalPointersVector<Condition> ConditionWeakVectorType;
typedef GlobalPointersVector<Condition >::iterator ConditionWeakIteratorType;

typedef GlobalPointersVector<Element> ParticleWeakVectorType;
typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
/// Default constructor.
SplitForwardEulerCylinderParticle();
SplitForwardEulerCylinderParticle( IndexType NewId, GeometryType::Pointer pGeometry );
SplitForwardEulerCylinderParticle( IndexType NewId, NodesArrayType const& ThisNodes);
SplitForwardEulerCylinderParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~SplitForwardEulerCylinderParticle();

SplitForwardEulerCylinderParticle& operator=(const SplitForwardEulerCylinderParticle& rOther);

double CalculateVolume() override;
double CalculateMomentOfInertia() override;

void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) override;
void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info) override;
void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info) override;
void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info) override;


/// Turn back information as a string.
virtual std::string Info() const override
{
std::stringstream buffer;
buffer << "SplitForwardEulerCylinderParticle" ;
return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "SplitForwardEulerCylinderParticle";}

protected:

double CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness) override;

double CalculateDampingNorm(const double& r_normal_damping_coeff, const double& r_tangential_damping_coeff) override;

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SplitForwardEulerSphericParticle );

    // public members

    // protected members
}

virtual void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SplitForwardEulerSphericParticle );

    // public members

    // protected members
}

}; // Class SplitForwardEulerCylinderParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            SplitForwardEulerCylinderParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const SplitForwardEulerCylinderParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_SPLIT_FORWARD_EULER_CYLINDER_PARTICLE_H_INCLUDED  defined
