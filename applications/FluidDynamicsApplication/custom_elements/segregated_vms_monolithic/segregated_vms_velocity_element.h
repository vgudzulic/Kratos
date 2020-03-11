//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_SEGREGATED_VMS_VELOCITY_ELEMENT_H_INCLUDED)
#define KRATOS_SEGREGATED_VMS_VELOCITY_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{
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

template <unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class SegregatedVMSVelocityElement : public Element
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of SegregatedVMSVelocityElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SegregatedVMSVelocityElement);

    /// definition of base type
    using BaseType = Element;

    /// definition of element type
    using ElementType = SegregatedVMSVelocityElement;

    /// definition of node type (default is: Node<3>)
    using NodeType = Node<3>;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    using PropertiesType = Properties;

    /// definition of the geometry type with given NodeType
    using GeometryType = Geometry<NodeType>;

    /// definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using VectorType = Vector;

    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using SizeType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using DofsArrayType = PointerVectorSet<Dof<double>, IndexedObject>;

    /// Type definition for integration methods
    using IntegrationMethod = GeometryData::IntegrationMethod;

    using GeometryDataType = GeometryData;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    static constexpr IndexType LocalSize = TDim * TNumNodes;
    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit SegregatedVMSVelocityElement(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    SegregatedVMSVelocityElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    SegregatedVMSVelocityElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    SegregatedVMSVelocityElement(IndexType NewId,
                                 GeometryType::Pointer pGeometry,
                                 PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.

    SegregatedVMSVelocityElement(SegregatedVMSVelocityElement const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor.

    ~SegregatedVMSVelocityElement() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * assignment operator: MANDATORY
     */

    /// Assignment operator.

    SegregatedVMSVelocityElement& operator=(SegregatedVMSVelocityElement const& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Informations
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SegregatedVMSVelocityElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SegregatedVMSVelocityElement>(NewId, pGeom, pProperties);
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        SegregatedVMSVelocityElement::Pointer p_new_elem =
            Kratos::make_intrusive<SegregatedVMSVelocityElement>(
                NewId, GetGeometry().Create(ThisNodes), pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * returns the used integration method. In the general case this is the
     * default integration method of the used geometry. I an other integration
     * method is used the method has to be overwritten within the element
     * @return default integration method of the used Geometry
     */
    IntegrationMethod GetIntegrationMethod() const override;

    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
    void GetValuesVector(Vector& Values, int Step = 0) override;

    /**
     * Getting method to obtain the time derivative of variable which defines the degrees of freedom
     */
    void GetFirstDerivativesVector(Vector& Values, int Step = 0) override;

    /**
     * Getting method to obtain the second time derivative of variable which defines the degrees of freedom
     */
    void GetSecondDerivativesVector(Vector& Values, int Step = 0) override;

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate Damp matrix and add velocity contribution to RHS
     * @param rDampingMatrix the velocity-proportional "damping" matrix
     * @param rRightHandSideVector the elemental right hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SegregatedVMSVelocityElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SegregatedVMSVelocityElement #" << Id();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
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
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
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

    void CalculateGeometryData(Vector& rGaussWeights,
                               Matrix& rNContainer,
                               GeometryType::ShapeFunctionsGradientsType& rDN_DX) const;

    double EvaluateInPoint(const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0) const;

    array_1d<double, 3> EvaluateInPoint(const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0) const;

    BoundedVector<double, TNumNodes> GetConvectionOperator(const array_1d<double, 3>& rVector,
                                                           const Matrix& rShapeDerivatives) const;

    array_1d<double, 3> GetVelocity(const Vector& rShapeFunction, const int Step = 0);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

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

}; // Class SegregatedVMSVelocityElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                SegregatedVMSVelocityElement<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SegregatedVMSVelocityElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}
} // namespace Kratos.
#endif // KRATOS_SEGREGATED_VMS_VELOCITY_ELEMENT_H_INCLUDED  defined
