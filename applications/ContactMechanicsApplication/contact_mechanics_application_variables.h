//-------------------------------------------------------------
//          ___         _           _   
//  KRATOS / __|___ _ _| |_ __ _ __| |_ 
//        | (__/ _ \ ' \  _/ _` / _|  _|
//         \___\___/_||_\__\__,_\__|\__|MECHANICS
//                                            
//  License:(BSD)    ContactMechanicsApplication/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   ...
//
//-------------------------------------------------------------
//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//


#if !defined(KRATOS_CONTACT_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CONTACT_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
KRATOS_DEFINE_VARIABLE( int, NUMBER_OF_ACTIVE_CONTACTS )
KRATOS_DEFINE_VARIABLE( int, NUMBER_OF_STICK_CONTACTS )
KRATOS_DEFINE_VARIABLE( int, NUMBER_OF_SLIP_CONTACTS )
KRATOS_DEFINE_VARIABLE( bool, FRICTION_ACTIVE )
KRATOS_DEFINE_VARIABLE( double, PENALTY_PARAMETER )
KRATOS_DEFINE_VARIABLE( double, LAGRANGE_MULTIPLIER_NORMAL )
KRATOS_DEFINE_VARIABLE( double, LAGRANGE_MULTIPLIER_NORMAL_REACTION )
KRATOS_DEFINE_VARIABLE( double, LAGRANGE_MULTIPLIER_TANGENTIAL )
KRATOS_DEFINE_VARIABLE( double, LAGRANGE_MULTIPLIER_TANGENTIAL_REACTION )
KRATOS_DEFINE_VARIABLE( double, TAU_STAB )
KRATOS_DEFINE_VARIABLE( double, MU_STATIC )
KRATOS_DEFINE_VARIABLE( double, MU_DYNAMIC )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )
KRATOS_DEFINE_VARIABLE( double, CONTACT_ADHESION )
KRATOS_DEFINE_VARIABLE( double, CONTACT_FRICTION_ANGLE )
KRATOS_DEFINE_VARIABLE( double, TANGENTIAL_PENALTY_RATIO )
KRATOS_DEFINE_VARIABLE( double, CONTACT_PLASTIC_SLIP )

}

#endif	/* KRATOS_CONTACT_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
