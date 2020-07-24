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

#if !defined(KRATOS_SPLIT_FORWARD_EULER_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_SPLIT_FORWARD_EULER_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
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

class KRATOS_API(DEM_APPLICATION) SplitForwardEulerSphericParticle : public SphericParticle
{

public:

/// Pointer definition of SplitForwardEulerSphericParticle
KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SplitForwardEulerSphericParticle);

/// Default constructor.
SplitForwardEulerSphericParticle();
SplitForwardEulerSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
SplitForwardEulerSphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
SplitForwardEulerSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~SplitForwardEulerSphericParticle();

SplitForwardEulerSphericParticle& operator=(const SplitForwardEulerSphericParticle& rOther);

typedef typename SphericParticle::ParticleDataBuffer ParticleDataBuffer;

// typedef std::unique_ptr<ParticleDataBuffer> BufferPointerType;

void Initialize(const ProcessInfo& r_process_info) override;
void CalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity, int search_control) override;

/// Turn back information as a string.
virtual std::string Info() const override
{
std::stringstream buffer;
buffer << "SplitForwardEulerSphericParticle" ;
return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "SplitForwardEulerSphericParticle";}

protected:

virtual void ComputeBallToBallStiffnessAndDamping(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                            unsigned int& r_num_sum,
                                                                            double& r_nodal_stiffness,
                                                                            double& r_nodal_damping,
                                                                            double& r_nodal_rotational_stiffness,
                                                                            double& r_nodal_rotational_damping);

virtual void ComputeBallToRigidFaceStiffnessAndDamping(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                                unsigned int& r_num_sum,
                                                                                double& r_nodal_stiffness,
                                                                                double& r_nodal_damping,
                                                                                double& r_nodal_rotational_stiffness,
                                                                                double& r_nodal_rotational_damping);

virtual double CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness);

virtual double CalculateDampingNorm(const double& r_normal_damping_coeff, const double& r_tangential_damping_coeff);

void AddUpForcesAndProject(double OldCoordSystem[3][3],
                                   double LocalCoordSystem[3][3],
                                   double LocalContactForce[3],
                                   double LocalElasticContactForce[3],
                                   double LocalElasticExtraContactForce[3],
                                   double GlobalContactForce[3],
                                   double GlobalElasticContactForce[3],
                                   double GlobalElasticExtraContactForce[3],
                                   double TotalGlobalElasticContactForce[3],
                                   double ViscoDampingLocalContactForce[3],
                                   const double cohesive_force,
                                   array_1d<double, 3>& other_ball_to_ball_forces,
                                   array_1d<double, 3>& rElasticForce,
                                   array_1d<double, 3>& rContactForce,
                                   const unsigned int i_neighbour_count,
                                   ProcessInfo& r_process_info) override;

virtual void AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                                      double LocalContactForce[3],
                                      double LocalElasticContactForce[3],
                                      double GlobalContactForce[3],
                                      double GlobalElasticContactForce[3],
                                      double ViscoDampingLocalContactForce[3],
                                      const double cohesive_force,
                                      array_1d<double, 3>& rElasticForce,
                                      array_1d<double, 3>& rContactForce,
                                      array_1d<double, 3>& elastic_force_backup,
                                      array_1d<double, 3>& total_force_backup) override;

virtual void AddUpMomentsAndProject(double LocalCoordSystem[3][3],
                                    double ElasticLocalRotationalMoment[3],
                                    double ViscoLocalRotationalMoment[3]) override;

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );

    // public members

    // protected members
}

virtual void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );

    // public members

    // protected members
}

}; // Class SplitForwardEulerSphericParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            SplitForwardEulerSphericParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const SplitForwardEulerSphericParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_SPLIT_FORWARD_EULER_SPHERIC_PARTICLE_H_INCLUDED  defined
