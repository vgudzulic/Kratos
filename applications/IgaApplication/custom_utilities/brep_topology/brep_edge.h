//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//                   Thomas Oberbichler
//


#if !defined(KRATOS_BREP_EDGE_H_INCLUDED )
#define  KRATOS_BREP_EDGE_H_INCLUDED

// System includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

#include "brep_trimming_curve.h"

#include "../node_curve_geometry_3d.h"

// Kratos includes
#include "includes/node.h"
#include "includes/model_part.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /// Edge in global space.
    /** Detail class definition.
    */
    class BrepEdge : public IndexedObject, public Flags
    {
    public:

        /// Pointer definition of KratosNurbsBrepApplication
        KRATOS_CLASS_POINTER_DEFINITION(BrepEdge);

        /* Geometry Refinement Parameters are used to pass the full refinement information
        *  needed for face-refinement.
        */
        struct GeometryRefinementParametersCurve
        {
            /// knot variables
            Vector knot_insertions_u;
            int multiply_knots_u;
            double max_element_size_u;
            /// degree variables
            int order_elevation_p;
            int min_order_p;

            /* Constructor */
            GeometryRefinementParametersCurve() {
                knot_insertions_u = ZeroVector(0);
                multiply_knots_u = 0;
                max_element_size_u = 0.0;
                order_elevation_p = 0;
                min_order_p = 0;
            }
        };


        /* Used to separate the curve into trimmed ranges. */
        struct TrimmingRange
        {
            int trim_index;
            Vector range;

            TrimmingRange(const int& rTrimIndex, const Vector& rRange)
            {
                trim_index = rTrimIndex;
                range = rRange;
            }
        };

        /* Used to descibe the topology of edges. */
        struct EdgeTopology
        {
            int brep_id;
            int trim_index;
            bool relative_direction;

            EdgeTopology(const int& rBrepId, const int& rTrimIndex, const bool& rRelativeDirection)
            {
                brep_id = rBrepId;
                trim_index = rTrimIndex;
                relative_direction = rRelativeDirection;
            }
        };

        ///@name Life Cycle 
        ///@{ 

        bool IsCouplingEdge();

        ///Constructor
        BrepEdge(
            int rBrepId,
            std::vector<EdgeTopology>& rBrepEdgeTopologyVector,
            std::vector<TrimmingRange>& rTrimmingRangeVector,
            int& rDegree,
            Vector& rKnotVector,
            Vector& rActiveRange,
            std::vector<int>& rControlPointIds,
            ModelPart& rModelPart);

        /// Destructor.
        virtual ~BrepEdge() {};

        ///@} 
    protected:

    private:

        ///@name Member Variables
        ///@{ 
        // topology parameter
        std::vector<EdgeTopology>              m_brep_edge_topology_vector;
        std::vector<TrimmingRange>             m_trimming_range_vector;

        //anurbs variables
        std::shared_ptr<NodeCurveGeometry3D>   m_node_curve_geometry_3d;

        //3d curve parameters
        Vector                        m_active_range;
        std::vector<int>              m_control_point_ids;
        ModelPart&                    m_model_part;
        ///@}    
    }; // Class BrepEdge 

}  // namespace Kratos.

#endif // KRATOS_BREP_EDGE_H_INCLUDED defined