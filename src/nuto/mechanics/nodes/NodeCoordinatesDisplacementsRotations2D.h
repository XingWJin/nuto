// $Id: $
#ifndef NodeCoordinatesDisplacementsRotations2D_H
#define NodeCoordinatesDisplacementsRotations2D_H

#include "nuto/mechanics/nodes/NodeCoordinates2D.h"
#include "nuto/mechanics/nodes/NodeDisplacements2D.h"
#include "nuto/mechanics/nodes/NodeRotations2D.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having coordinates, displacements and rotations
class NodeCoordinatesDisplacementsRotations2D : public  NodeCoordinates2D,public NodeDisplacements2D, public NodeRotations2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeCoordinatesDisplacementsRotations2D() : NodeCoordinates2D (),
            NodeDisplacements2D(),
            NodeRotations2D()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements2D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeRotations2D);
    }
#endif  // ENABLE_SERIALIZATION
    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeCoordinates2D::SetGlobalDofs(rDOF);
        NodeDisplacements2D::SetGlobalDofs(rDOF);
        NodeRotations2D::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        NodeCoordinates2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::SetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        NodeCoordinates2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeDisplacements2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        NodeRotations2D::GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        NodeCoordinates2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeDisplacements2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
        NodeRotations2D::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeCoordinatesDisplacementsRotations2D");
    }
};
}

#endif //NodeCoordinatesDisplacementsRotations2D_H