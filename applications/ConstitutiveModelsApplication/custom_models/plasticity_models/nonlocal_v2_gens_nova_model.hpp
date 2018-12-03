//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                MCiantia $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NONLOCAL_V2_GENS_NOVA_MODEL_H_INCLUDED )
#define  KRATOS_NONLOCAL_V2_GENS_NOVA_MODEL_H_INCLUDED

// System includes

// External includes
#include <iostream>
#include <fstream>

// Project includes
#include "custom_models/plasticity_models/non_associative_plasticity_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/gens_nova_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/gens_nova_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"
#include "custom_models/elasticity_models/tamagnini_model.hpp"



//***** the hardening law associated to this Model has ... variables
// 0. Plastic multiplier
// 1. Plastic Volumetric deformation
// 2. Plastic Deviatoric deformation
// 3. ps (mechanical)
// 4. pt (ageing)
// 5. pcSTAR = ps + (1+k) p_t
// 6. Plastic Volumetric deformation Absolut Value
// 7. NonLocal Plastic Vol Def
// 8. NonLocal Plastic Dev Def
// 9. NonLocal Plastic Vol Def ABS
// ... (the number now is then..., xD)

namespace Kratos
{
   ///@addtogroup ConstitutiveModelsApplication
   ///@{

   ///@name Kratos Globals
   ///@{

   ///@}
   ///@name Type Definitions
   ///@{

   ///@}
   ///@name  Enum's
   ///@{

   ///@}
   ///@name  Functions
   ///@{

   ///@}
   ///@name Kratos Classes
   ///@{

   /// Short class definition.
   /** Detail class definition.
    */
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonlocalV2GensNovaModel : public NonAssociativePlasticityModel<TamagniniModel, GensNovaYieldSurface<GensNovaHardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         //typedef BorjaModel                                     ElasticityModelType;
         typedef TamagniniModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef GensNovaHardeningRule                             HardeningRuleType;
         typedef GensNovaYieldSurface<HardeningRuleType>    YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

         //base type
         typedef NonAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of NonlocalV2GensNovaModel
         KRATOS_CLASS_POINTER_DEFINITION( NonlocalV2GensNovaModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         NonlocalV2GensNovaModel() : BaseType() { mInitialized = false; }

         /// Copy constructor.
         NonlocalV2GensNovaModel(NonlocalV2GensNovaModel const& rOther) : BaseType(rOther), mInitialized(rOther.mInitialized) {}

         /// Assignment operator.
         NonlocalV2GensNovaModel& operator=(NonlocalV2GensNovaModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return ( NonlocalV2GensNovaModel::Pointer(new NonlocalV2GensNovaModel(*this)) );
         }

         /// Destructor.
         virtual ~NonlocalV2GensNovaModel() {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{

         
         /**
          * Initialize member data
          */    
         void InitializeModel(ModelDataType& rValues) override
         {
            KRATOS_TRY

            if (mInitialized == false) {
               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               const ModelDataType & rModelData = Variables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();

               double k = rMaterialProperties[KSIM];

               double & rPS     = Variables.Internal.Variables[3];
               double & rPT     = Variables.Internal.Variables[4];
               double & rPCstar = Variables.Internal.Variables[5];

               rPS = -rMaterialProperties[PS];
               rPT = -rMaterialProperties[PT];
               rPCstar = rPS + (1.0+k)*rPT;

               MatrixType Stress;
               this->UpdateInternalVariables(rValues, Variables, Stress);


               mInitialized = true;
            }
            mElasticityModel.InitializeModel( rValues );

            KRATOS_CATCH("")
         }

         /**
          * Check
          */    
         virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
         {
            KRATOS_TRY

            //LMV: to be implemented. but should not enter in the base one

            return 0;

            KRATOS_CATCH("")
         }

         ///@}
         ///@name Access
         ///@{

         /**
          * Has Values
          */   
         virtual bool Has(const Variable<double>& rThisVariable) override
         {
            if(rThisVariable == PLASTIC_STRAIN || rThisVariable == DELTA_PLASTIC_STRAIN )
               return true;

            return false;
         }


         /**
          * Get Values
          */
         void SetValue(const Variable<double>& rVariable,
               const double& rValue,
               const ProcessInfo& rCurrentProcessInfo) override 
         {
            KRATOS_TRY

            if ( rVariable == NONLOCAL_PLASTIC_VOL_DEF) {
               mInternal.Variables[7] = rValue;
            }
            else if ( rVariable == NONLOCAL_PLASTIC_DEV_DEF) {
               mInternal.Variables[8] = rValue;
            }
            else if ( rVariable == NONLOCAL_PLASTIC_VOL_DEF_ABS) {
               mInternal.Variables[9] = rValue;
            }

            KRATOS_CATCH("")
         }


         /**
          * Get Values
          */
         virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
         {
            KRATOS_TRY

            rValue=0;

            if (rThisVariable==PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0];
            }
            else if (rThisVariable==DELTA_PLASTIC_STRAIN)
            {
               rValue = this->mInternal.Variables[0]-mPreviousInternal.Variables[0];
            }
            else if ( rThisVariable == PS)
            {
               rValue = this->mInternal.Variables[3];
            }
            else if ( rThisVariable == PT)
            {
               rValue = this->mInternal.Variables[4];
            }
            else if ( rThisVariable == PM)
            {
               rValue = this->mInternal.Variables[5];
            } 
            else if ( rThisVariable == PLASTIC_VOL_DEF)
            {
               rValue = this->mInternal.Variables[1];
            }
            else if ( rThisVariable == PLASTIC_DEV_DEF)
            {
               rValue = this->mInternal.Variables[2];
            }
            else if ( rThisVariable == PLASTIC_VOL_DEF_ABS)
            {
               rValue = this->mInternal.Variables[6];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_VOL_DEF)
            {
               rValue = this->mInternal.Variables[7];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_DEV_DEF)
            {
               rValue = this->mInternal.Variables[8];
            }
            else if ( rThisVariable == NONLOCAL_PLASTIC_VOL_DEF_ABS)
            {
               rValue = this->mInternal.Variables[9];
            }
            else {
               rValue = NonAssociativePlasticityModel::GetValue( rThisVariable, rValue);
            }
            return rValue;

            KRATOS_CATCH("")
         }

         ///@}
         ///@name Inquiry
         ///@{


         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         virtual std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "NonlocalV2GensNovaModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "NonlocalV2GensNovaModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "NonlocalV2GensNovaModel Data";
         }


         ///@}
         ///@name Friends
         ///@{


         ///@}

      protected:
         ///@name Protected static Member Variables
         ///@{


         ///@}
         ///@name Protected member Variables
         ///@{

         bool mInitialized; 

         ///@}
         ///@name Protected Operators
         ///@{


         ///@}
         ///@name Protected Operations
         ///@{
            // Calculate Stress and constitutive tensor
            void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               { // modify the internal variables to make it nonlocal
               }

               double LocalPlasticVolStrain = mInternal.Variables[1];
               double NonLocalPlasticVolStrain = mInternal.Variables[7];
               mInternal.Variables[1] = mInternal.Variables[7];

               double LocalPlasticDevStrain = mInternal.Variables[2];
               double NonLocalPlasticDevStrain = mInternal.Variables[8];
               mInternal.Variables[2] = mInternal.Variables[8];

               double LocalPlasticVolStrainAbs = mInternal.Variables[6];
               double NonLocalPlasticVolStrainAbs = mInternal.Variables[9];
               mInternal.Variables[6] = mInternal.Variables[9];

               // integrate "analytically" ps and pt from plastic variables. Then update the internal variables.

               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               const ModelDataType & rModelData = Variables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();

               const double & rPs0 = rMaterialProperties[PS];
               const double & rPt0 = rMaterialProperties[PT];
               const double & rChis = rMaterialProperties[CHIS];
               const double & rChit = rMaterialProperties[CHIT];
               const double & rhos = rMaterialProperties[RHOS];
               const double & rhot = rMaterialProperties[RHOT];
               const double & k = rMaterialProperties[KSIM];

               const double & rPlasticVolDef = Variables.Internal.Variables[1]; 
               const double & rPlasticDevDef = Variables.Internal.Variables[2];
               const double & rPlasticVolDefAbs = Variables.Internal.Variables[6];

               double sq2_3 = sqrt(2.0/3.0);

               double ps;
               ps = rPlasticVolDef + sq2_3 * rChis * rPlasticDevDef; 
               ps = (-rPs0) * std::exp( -rhos*ps);

               double pt;
               pt = rPlasticVolDefAbs + sq2_3 * rChit * rPlasticDevDef; 
               pt = (-rPt0) * std::exp( rhot*pt);

               double pm;
               pm = ps + (1.0+k)*pt;


               mInternal.Variables[3] = ps;
               mInternal.Variables[4] = pt;
               mInternal.Variables[5] = pm;

               NonAssociativePlasticityModel::CalculateStressAndConstitutiveTensors( rValues, rStressMatrix, rConstitutiveMatrix);

               if ( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES) ) {
                  mInternal.Variables[1] = LocalPlasticVolStrain + ( mInternal.Variables[1] -  NonLocalPlasticVolStrain);
                  mInternal.Variables[2] = LocalPlasticDevStrain + ( mInternal.Variables[2] -  NonLocalPlasticDevStrain);
                  mInternal.Variables[6] = LocalPlasticVolStrainAbs + ( mInternal.Variables[6] -  NonLocalPlasticVolStrainAbs);
               } else {
                  mInternal.Variables[1] = LocalPlasticVolStrain;
                  mInternal.Variables[2] = LocalPlasticDevStrain;
                  mInternal.Variables[6] = LocalPlasticVolStrainAbs;
               }



               KRATOS_CATCH("")
            }
            //***************************************************************************************
            //***************************************************************************************
            // Compute Elasto Plastic Matrix
            void ComputeElastoPlasticTangentMatrix( ModelDataType & rValues, PlasticDataType & rVariables, Matrix & rEPMatrix) override
            {

               KRATOS_TRY

               // evaluate constitutive matrix and plastic flow
               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV


               MatrixType PlasticPotDerTensor;
               PlasticPotDerTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialDerivative, PlasticPotDerTensor);
               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerTensor);

               VectorType AuxF = prod( trans(DeltaStressYieldCondition), rEPMatrix);
               VectorType AuxG = prod( rEPMatrix, PlasticPotentialDerivative);

               Matrix PlasticUpdateMatrix(6,6);
               noalias(PlasticUpdateMatrix) = ZeroMatrix(6,6);
               double denom = 0;
               for (unsigned int i = 0; i < 6; i++) {
                  denom += AuxF(i)*PlasticPotentialDerivative(i);
                  for (unsigned int j = 0; j < 6; j++) {
                     PlasticUpdateMatrix(i,j) = AuxF(i) * AuxG(j);
                  }
               }

               rEPMatrix -= PlasticUpdateMatrix / ( H + denom);

               KRATOS_CATCH("")
            }
            //***********************************************************************************
            //***********************************************************************************
            // Compute one step of the elasto-plastic problem
            void ComputeOneStepElastoPlasticProblem( ModelDataType & rValues, PlasticDataType & rVariables, const MatrixType & rDeltaDeformationMatrix) override
            {
               KRATOS_TRY
            

               /*{ // debugging, or whatever
                  const ModelDataType & rModelData = rVariables.GetModelData();
                  const Properties & rMaterialProperties = rModelData.GetProperties();
                  const double & rPs0 = rMaterialProperties[PS];
                  const double & rPt0 = rMaterialProperties[PT];
                  const double & rChis = rMaterialProperties[CHIS];
               const double & rChit = rMaterialProperties[CHIT];
                  const double & rhos = rMaterialProperties[RHOS];
               const double & rhot = rMaterialProperties[RHOT];

                  double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
                  double & rPlasticDevDef = rVariables.Internal.Variables[2];
                  double & rPS     = rVariables.Internal.Variables[3];
                  double & rPT     = rVariables.Internal.Variables[4];
                  double & rPlasticVolDefAbs = rVariables.Internal.Variables[6];

                  double ps;
                  ps = rPlasticVolDef + rChis * rPlasticDevDef; 
                  ps = (-rPs0) * std::exp( -rhos*ps);

                  double pt;
                  pt = rPlasticVolDefAbs + rChit * rPlasticDevDef; 
                  pt = (-rPt0) * std::exp( rhot*pt);
                  

               }*/

               const ModelDataType & rModelData = rVariables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();
               const double & rhos = rMaterialProperties[RHOS];
               const double & rhot = rMaterialProperties[RHOT];
               double k =  rMaterialProperties[KSIM];
      
               const double & rChis = rMaterialProperties[CHIS];
               const double & rChit = rMaterialProperties[CHIT];

               MatrixType StressMatrix;
               // evaluate constitutive matrix and plastic flow
               double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
               double & rPlasticMultiplier = rVariables.Internal.Variables[0];
               double & rPlasticDevDef = rVariables.Internal.Variables[2];
               double & rPS     = rVariables.Internal.Variables[3];
               double & rPT     = rVariables.Internal.Variables[4];
               double & rPCstar = rVariables.Internal.Variables[5];

               Matrix ElasticMatrix(6,6);
               noalias(ElasticMatrix) = ZeroMatrix(6,6);
               this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

               VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
               VectorType PlasticPotentialDerivative;
               PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

               MatrixType PlasticPotDerTensor;
               PlasticPotDerTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialDerivative, PlasticPotDerTensor);
               double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerTensor);

               MatrixType StrainMatrix = prod( rDeltaDeformationMatrix, trans( rDeltaDeformationMatrix) );
               VectorType StrainVector; 
               ConvertCauchyGreenTensorToHenckyVector( StrainMatrix, StrainVector);

               VectorType AuxVector;
               AuxVector = prod( ElasticMatrix, StrainVector);
               double DeltaGamma;
               DeltaGamma = MathUtils<double>::Dot( AuxVector, DeltaStressYieldCondition);

               double Denominador = H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) );

               DeltaGamma /= Denominador;

               if ( DeltaGamma < 0)
                  DeltaGamma = 0;

               MatrixType UpdateMatrix;
               ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);
               UpdateMatrix = prod( rDeltaDeformationMatrix, UpdateMatrix);


               rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
               rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

               this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

               rPlasticMultiplier += DeltaGamma;
               double VolPlasticIncr = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  VolPlasticIncr += DeltaGamma * DeltaStressYieldCondition(i);
               rPlasticVolDef += VolPlasticIncr;

               double DevPlasticIncr = 0.0;
               for (unsigned int i = 0; i < 3; i++)
                  DevPlasticIncr += pow( DeltaGamma * DeltaStressYieldCondition(i) - VolPlasticIncr/3.0, 2.0);
               for (unsigned int i = 3; i < 6; i++)
                  DevPlasticIncr += 2.0 * pow( DeltaGamma *  DeltaStressYieldCondition(i) /2.0 , 2.0);
               DevPlasticIncr = sqrt(DevPlasticIncr);
               rPlasticDevDef += DevPlasticIncr;


               double hs = rhos * ( rPS) * (     VolPlasticIncr  + rChis*sqrt(2.0/3.0) * DevPlasticIncr );
               double ht = rhot * (-rPT) * (fabs(VolPlasticIncr) + rChit*sqrt(2.0/3.0) * DevPlasticIncr );

               rPS -= hs;
               rPT -= ht;

               rPCstar = rPS + (1.0 + k)*rPT;

               KRATOS_CATCH("")
            }

            //********************************************************************
            //********************************************************************
            // UpdateInternalVariables
            virtual void UpdateInternalVariables(ModelDataType& rValues, PlasticDataType& rVariables, const MatrixType& rStressMatrix) override
            {
               KRATOS_TRY

               mPreviousInternal.Variables[6] = mInternal.Variables[6];
               mInternal.Variables[6] = mInternal.Variables[6] + fabs( rVariables.Internal.Variables[1] - mInternal.Variables[1]);
               for (unsigned int i = 0; i < 6; i++) {
                  double & rCurrentPlasticVariable = rVariables.Internal.Variables[i]; 
                  double & rPreviousPlasticVariable    = mInternal.Variables[i];

                  mPreviousInternal.Variables[i] = rPreviousPlasticVariable;
                  rPreviousPlasticVariable = rCurrentPlasticVariable;
               }


               KRATOS_CATCH("")
            }

            //***************************************************************************************
            //***************************************************************************************
            // Correct Yield Surface Drift According to 
            virtual void ReturnStressToYieldSurface( ModelDataType & rValues, PlasticDataType & rVariables) override
            {
               KRATOS_TRY



               double Tolerance = 1e-7;

                  MatrixType StressMatrix;
                  this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);
               double YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);

               if ( fabs(YieldSurface) < Tolerance)
                  return;

               const ModelDataType & rModelData = rVariables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();
               double rhos = rMaterialProperties[RHOS];
               double rhot = rMaterialProperties[RHOT];
               double chis = rMaterialProperties[CHIS];
               double chit = rMaterialProperties[CHIT];
               double k =  rMaterialProperties[KSIM];
               // evaluate constitutive matrix and plastic flow
               double & rPlasticVolDef = rVariables.Internal.Variables[1]; 
               //double & rPlasticMultiplier = rVariables.Internal.Variables[0];
               double & rPlasticDevDef = rVariables.Internal.Variables[2];
               double & rPS     = rVariables.Internal.Variables[3];
               double & rPT     = rVariables.Internal.Variables[4];
               double & rPCstar = rVariables.Internal.Variables[5];

               for (unsigned int i = 0; i < 150; i++) {

                  Matrix ElasticMatrix(6,6);
                  noalias(ElasticMatrix) = ZeroMatrix(6,6);
                  this->mElasticityModel.CalculateConstitutiveTensor( rValues, ElasticMatrix);

                  VectorType DeltaStressYieldCondition = this->mYieldSurface.CalculateDeltaStressYieldCondition( rVariables, DeltaStressYieldCondition);
                  VectorType PlasticPotentialDerivative;
                  PlasticPotentialDerivative = DeltaStressYieldCondition; // LMV

                  MatrixType PlasticPotDerTensor;
                  PlasticPotDerTensor = ConstitutiveModelUtilities::StrainVectorToTensor( PlasticPotentialDerivative, PlasticPotDerTensor);
                  double H = this->mYieldSurface.GetHardeningRule().CalculateDeltaHardening( rVariables, H, PlasticPotDerTensor);

                  double DeltaGamma = YieldSurface;
                  DeltaGamma /= ( H + MathUtils<double>::Dot( DeltaStressYieldCondition, prod(ElasticMatrix, PlasticPotentialDerivative) ) );

                  MatrixType UpdateMatrix;
                  ConvertHenckyVectorToCauchyGreenTensor( -DeltaGamma * PlasticPotentialDerivative / 2.0, UpdateMatrix);

                  rValues.StrainMatrix = prod( UpdateMatrix, rValues.StrainMatrix);
                  rValues.StrainMatrix = prod( rValues.StrainMatrix, trans(UpdateMatrix));

                  this->mElasticityModel.CalculateStressTensor( rValues, StressMatrix);

                  double VolPlasticIncr = 0.0;
                  for (unsigned int i = 0; i < 3; i++)
                     VolPlasticIncr += DeltaGamma * DeltaStressYieldCondition(i);
                  rPlasticVolDef += VolPlasticIncr;

                  double DevPlasticIncr = 0.0;
                  for (unsigned int i = 0; i < 3; i++)
                     DevPlasticIncr += pow( DeltaGamma * DeltaStressYieldCondition(i) - VolPlasticIncr/3.0, 2.0);
                  for (unsigned int i = 3; i < 6; i++)
                     DevPlasticIncr += 2.0 * pow( DeltaGamma *  DeltaStressYieldCondition(i) /2.0 , 2.0);
                  DevPlasticIncr = DeltaGamma/fabs(DeltaGamma) * sqrt(DevPlasticIncr);
                  rPlasticDevDef += DevPlasticIncr;


                  double hs =  rhos * rPS * (VolPlasticIncr + chis * sqrt(2.0/3.0)* DevPlasticIncr);
                  double ht = rhot * rPT * ( -fabs(VolPlasticIncr)  + chit*sqrt(2.0/3.0)*DevPlasticIncr);
                  rPS -= hs;
                  rPT -= ht;

                  rPCstar = rPS + (1.0 + k)*rPT;


                  YieldSurface = this->mYieldSurface.CalculateYieldCondition( rVariables, YieldSurface);


                  if ( fabs( YieldSurface) < Tolerance) {
                     return;
                  }
               }


               KRATOS_CATCH("")
            }
         ///@}
         ///@name Protected  Access
         ///@{


         ///@}
         ///@name Protected Inquiry
         ///@{


         ///@}
         ///@name Protected LifeCycle
         ///@{


         ///@}

      private:
         ///@name Static Member Variables
         ///@{


         ///@}
         ///@name Member Variables
         ///@{


         ///@}
         ///@name Private Operators
         ///@{


         ///@}
         ///@name Private Operations
         ///@{


         ///@}
         ///@name Private  Access
         ///@{


         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Serialization
         ///@{    
         friend class Serializer;

         virtual void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Un accessible methods
         ///@{


         ///@}

   }; // Class NonlocalV2GensNovaModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@} 
   ///@name Input and output 
   ///@{


   ///@}

   ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_NONLOCAL_V2_GENS_NOVA_MODEL_H_INCLUDED  defined 


