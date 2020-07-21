// System includes
#include <string>
#include <iostream>
#include <cmath>

#include <fstream>

// External includes

// Project includes
#include "split_forward_euler_cylinder_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_strategies/schemes/glued_to_wall_scheme.h"


namespace Kratos
{
// using namespace GeometryFunctions;

SplitForwardEulerCylinderParticle::SplitForwardEulerCylinderParticle()
    : SplitForwardEulerSphericParticle() {
}

SplitForwardEulerCylinderParticle::SplitForwardEulerCylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SplitForwardEulerSphericParticle(NewId, pGeometry) {
}

SplitForwardEulerCylinderParticle::SplitForwardEulerCylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SplitForwardEulerSphericParticle(NewId, pGeometry, pProperties) {
}

SplitForwardEulerCylinderParticle::SplitForwardEulerCylinderParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SplitForwardEulerSphericParticle(NewId, ThisNodes) {
}

Element::Pointer SplitForwardEulerCylinderParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SplitForwardEulerCylinderParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/// Destructor.
SplitForwardEulerCylinderParticle::~SplitForwardEulerCylinderParticle(){
}

SplitForwardEulerCylinderParticle& SplitForwardEulerCylinderParticle::operator=(const SplitForwardEulerCylinderParticle& rOther) {
    SplitForwardEulerSphericParticle::operator=(rOther);

    return *this;
}

double SplitForwardEulerCylinderParticle::CalculateVolume(){
    return Globals::Pi * GetRadius() * GetRadius();
}

double SplitForwardEulerCylinderParticle::CalculateMomentOfInertia() {
    return 0.5 * GetMass() * GetRadius() * GetRadius();
}

void SplitForwardEulerCylinderParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info){}
void SplitForwardEulerCylinderParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info){}
void SplitForwardEulerCylinderParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
void SplitForwardEulerCylinderParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}

double SplitForwardEulerCylinderParticle::CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness) {
    return std::sqrt(r_normal_stiffness*r_normal_stiffness+r_tangential_stiffness*r_tangential_stiffness);
}

double SplitForwardEulerCylinderParticle::CalculateDampingNorm(const double& r_normal_damping_coeff, const double& r_tangential_damping_coeff) {
    return std::sqrt(r_normal_damping_coeff*r_normal_damping_coeff+r_tangential_damping_coeff*r_tangential_damping_coeff);
}

}  // namespace Kratos.
