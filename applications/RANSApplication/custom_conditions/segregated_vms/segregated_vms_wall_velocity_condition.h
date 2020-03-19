//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Suneth Warnakulasuriya
//

#ifndef KRATOS_SEGREGATED_VMS_WALL_VELOCITY_CONDITION_H
#define KRATOS_SEGREGATED_VMS_WALL_VELOCITY_CONDITION_H

// System includes

// External includes

// Project includes

// Application includes
#include "includes/condition.h"

namespace Kratos
{
///@addtogroup RANSApplication
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

/**
 * @brief
 *
 * @tparam TDim       Dimensionality of the condition (2D or 3D)
 * @tparam TNumNodes  Number of nodes in the condition
 */
template <unsigned int TDim, unsigned int TNumNodes = TDim>
class SegregatedVMSWallVelocityCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SegregatedVMSWallVelocityCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SegregatedVMSWallVelocityCondition);

    using BaseType = Condition;
    using NodeType = Node<3>;
    using PropertiesType = Properties;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit SegregatedVMSWallVelocityCondition(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    SegregatedVMSWallVelocityCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    SegregatedVMSWallVelocityCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    SegregatedVMSWallVelocityCondition(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    SegregatedVMSWallVelocityCondition(SegregatedVMSWallVelocityCondition const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor.
    ~SegregatedVMSWallVelocityCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    SegregatedVMSWallVelocityCondition& operator=(SegregatedVMSWallVelocityCondition const& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new SegregatedVMSWallVelocityCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SegregatedVMSWallVelocityCondition<TDim, TNumNodes>>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SegregatedVMSWallVelocityCondition<TDim, TNumNodes>>(
            NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition =
            Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override;

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& rCurrentProcessInfo) override;

    /// Calculate wall stress term for all nodes with SLIP set.
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo) override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& Values, int Step = 0) override;

    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& Values, int Step = 0) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Applies rans based wall law
     *
     * This method calculate left hand side matrix and right hand side vector for following equation
     *
     * \[
     *      u_\tau = max\left(C_\mu^0.25 \sqrt{k}, \frac{||\underline{u}||}{\frac{1}{\kappa}ln(y^+)+\beta}\right)
     * \]
     *
     * integration point value = \rho \frac{u_\tau^2}{||\underline{u}||}\underline{u}
     *
     * @param rLocalMatrix         Left hand side matrix
     * @param rLocalVector         Right hand side vector
     * @param rCurrentProcessInfo  Current process info from model part
     */
    void ApplyWallLaw(MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo);

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

    double mWallHeight;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SegregatedVMSWallVelocityCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                SegregatedVMSWallVelocityCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SegregatedVMSWallVelocityCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_SEGREGATED_VMS_WALL_VELOCITY_CONDITION_H
