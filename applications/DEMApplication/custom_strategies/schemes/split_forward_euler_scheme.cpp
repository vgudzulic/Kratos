// Project includes
#include "split_forward_euler_scheme.h"

namespace Kratos {

    void SplitForwardEulerScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning SplitForwardEulerScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void SplitForwardEulerScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning SplitForwardEulerScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void SplitForwardEulerScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 >& i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]) {

        const double nodal_damping = i.GetValue(NODAL_DAMPING);

        if (nodal_damping > std::numeric_limits<double>::epsilon()) {
            const double nodal_aux_mass = i.GetValue(NODAL_AUX_MASS);
            array_1d<double, 3 >& aux_displacement = i.FastGetSolutionStepValue(AUX_DISPLACEMENT);
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    aux_displacement[k] += delta_t * force[k] / nodal_damping;
                    delta_displ[k] = (delta_t * aux_displacement[k] + nodal_aux_mass * displ[k]) / (delta_t + nodal_aux_mass) - displ[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                    vel[k] = delta_displ[k] / delta_t;
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        } else {
            // TODO: is this right?
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = delta_t * vel[k] + 0.5 * force[k] / mass * delta_t*delta_t;
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                    vel[k] = delta_displ[k] / delta_t;
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
    }

    void SplitForwardEulerScheme::CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node < 3 >& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        const double nodal_rotational_damping = i.GetValue(NODAL_ROTATIONAL_DAMPING);

        if (nodal_rotational_damping > std::numeric_limits<double>::epsilon()) {
            const double nodal_aux_inertia = i.GetValue(NODAL_AUX_INERTIA);
            array_1d<double, 3 >& aux_rotated_angle = i.FastGetSolutionStepValue(AUX_ROTATION);
            for (int k = 0; k < 3; k++) {
                if (Fix_Ang_vel[k] == false) {
                    aux_rotated_angle[k] += delta_t * torque[k] / nodal_rotational_damping;
                    delta_rotation[k] = (delta_t * aux_rotated_angle[k] + nodal_aux_inertia * rotated_angle[k]) / (delta_t + nodal_aux_inertia) - rotated_angle[k];
                    rotated_angle[k] += delta_rotation[k];
                    angular_velocity[k] = delta_rotation[k] / delta_t;
                } else {
                    delta_rotation[k] = delta_t * angular_velocity[k];
                    rotated_angle[k] += delta_rotation[k];
                }
            }
        } else {
            // TODO: is this right?
            for (int k = 0; k < 3; k++) {
                if (Fix_Ang_vel[k] == false) {
                    delta_rotation[k] = delta_t * angular_velocity[k] + 0.5 * torque[k] / moment_of_inertia * delta_t*delta_t;
                    rotated_angle[k] += delta_rotation[k];
                    angular_velocity[k] = delta_rotation[k] / delta_t;
                } else {
                    delta_rotation[k] = delta_t * angular_velocity[k];
                    rotated_angle[k] += delta_rotation[k];
                }
            }
        }
    }

    void SplitForwardEulerScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 > moments_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        array_1d<double, 3 >& local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);

        array_1d<double, 3 > local_angular_acceleration, local_torque, angular_acceleration;

        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
        CalculateLocalAngularAccelerationByEulerEquations(local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);
        GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_acceleration, angular_acceleration);

        // TODO: this is from symplectic euler...
        UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);

        double ang = DEM_INNER_PRODUCT_3(delta_rotation, delta_rotation);

        if (ang) {
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        } //if ang
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
    }

    void SplitForwardEulerScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        for (int k = 0; k < 3; k++) {
            if (Fix_Ang_vel[k] == false) {
                angular_velocity[k] += delta_t * angular_acceleration[k];
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }

    void SplitForwardEulerScheme::CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration) {

        for (int j = 0; j < 3; j++) {
            local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
            local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
        }
    }
} //namespace Kratos
