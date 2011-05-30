// $Id$
#ifndef Voxel8N_H
#define Voxel8N_H

#include "nuto/mechanics/elements/Solid.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class ConstitutiveTangentLocal6x6;
class EngineeringStress3D;

//! @author Andrea Keßler, ISM
//! @date April 2010
//! @brief ... voxel element without nodes
class Voxel8N : public Solid
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    Voxel8N(NuTo::StructureBase* rStructure,int rElementID,int rNumLocalCoefficientMatrix0,
            NuTo::ElementData::eElementDataType rElementDataType , IpData::eIpDataType rIpDataType):
            	NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(), rIpDataType)
    {
        mVoxelID = rElementID;
        mVoxelLocation=0;
        mNodeIds=0;
        mNumLocalCoefficientMatrix0 =rNumLocalCoefficientMatrix0;
    }
    //! @brief destructor
    ~Voxel8N()
    {
        delete mVoxelLocation;
        mVoxelLocation=0;
        delete mNodeIds;
        mNodeIds=0;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        return NuTo::Element::VOXEL8N;
    }

    int GetVoxelID()const
    {
        return mVoxelID;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes() const
    {
        return 8;
    }

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr);

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctions() const
    {
        return 8;
    }

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctions(const double rLocalCoordinates[3],
                                 std::vector<double>& rShapeFunctions) const;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsLocal(
        const double rLocalCoordinates[3],
        std::vector<double>& rDerivativeShapeFunctions) const;

    //! @brief returns the number of local degrees of freedom
    //! @return number of local degrees of freedom
    inline int GetNumDofs() const
    {
        return 24;
    }

    //! @TODO change routine to find the nodes

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::GetNode] implement first.");
        /*
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<=7);
        return mNodes[rLocalNodeNumber];
         */
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber) const
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::GetNode] implement first.");
        /*
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<=7);
        return mNodes[rLocalNodeNumber];
        */
    }

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::SetNode] implement first.");
        /*
       assert(rLocalNodeNumber>=0 && rLocalNodeNumber<=7);
        mNodes[rLocalNodeNumber] = rNode;
        */
    }
    //! @brief ... reorder nodes is not possible in gris structure
    void ReorderNodes()
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::ReorderNodes] Grid nodes can not be reordered.");
    }


    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType();

    int GetNumLocalStiffnessMatrix();

    //! @brief Get voxel location for this element
    //! @return int number in x, y, z direction
    int* GetVoxelLocation();

    //! @brief Set voxel location for this element
    //! @param int* number in x, y, z direction
    void SetVoxelLocation(int * rVoxelLocation);

    //! @brief Get ids of all nodes of this element
    //! @brief sorted: bottom - top, each counterclockwise
    //! @return int* ids of all nodes
    int* GetNodeIds();

    //! @brief Set ids of all nodes for this element
    //! @brief sorted: bottom - top, each counterclockwise
    //! @param  int* ids of all nodes
    void SetNodeIds(int * rNodeIds);

 protected:
    //! @brief ... just for serialization
    Voxel8N(){};

    int mVoxelID;
    int* mVoxelLocation;
    int* mNodeIds;
    int mNumLocalCoefficientMatrix0;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Voxel8N)
#endif // ENABLE_SERIALIZATION
#endif //Voxel8N_H
