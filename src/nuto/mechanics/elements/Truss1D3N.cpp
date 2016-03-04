// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>

NuTo::Truss1D3N::Truss1D3N(NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType) :
        NuTo::Truss1D::Truss1D(rStructure, rElementDataType, GetStandardIntegrationType(),rIpDataType)
{
	if (rNodes.size()!=3)
    	throw MechanicsException("[NuTo::Truss3N::Truss3N] Exactly three nodes are required for this type of element.");
    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D3N::CalculateShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    ShapeFunctions1D::ShapeFunctionsTrussOrder2(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D3N::CalculateShapeFunctionsNonlocalTotalStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
   ShapeFunctions1D::ShapeFunctionsTrussOrder1(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss1D3N::CalculateShapeFunctionsNonlocalEqStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    ShapeFunctions1D::ShapeFunctionsTrussOrder1(rLocalCoordinates, rShapeFunctions);
}


//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D3N::CalculateDerivativeShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D3N::CalculateDerivativeShapeFunctionsNonlocalTotalStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss1D3N::CalculateDerivativeShapeFunctionsNonlocalEqStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss1D3N::GetNodeNonlocalTotalStrain(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[2];
    throw MechanicsException("[NuTo::Truss1D3N::GetNodeNonlocalTotalStrain] the interpolation order is only linear for the total strain.");
    return 0;
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Truss1D3N::GetNodeNonlocalTotalStrain(int rLocalNodeNumber)const
{
    assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[2];
    throw MechanicsException("[NuTo::Truss1D3N::GetNodeNonlocalTotalStrain] the interpolation order is only linear for the total strain.");
    return 0;
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss1D3N::GetNodeNonlocalEqStrain(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[2];
    throw MechanicsException("[NuTo::Truss1D3N::GetNodeNonlocalEqStrain] the interpolation order is only linear for the nonlocal eq strain.");
    return 0;
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Truss1D3N::GetNodeNonlocalEqStrain(int rLocalNodeNumber)const
{
    assert(rLocalNodeNumber==0 || rLocalNodeNumber==1 );
    if (rLocalNodeNumber==0)
        return mNodes[0];
    if (rLocalNodeNumber==1)
        return mNodes[2];
    throw MechanicsException("[NuTo::Truss1D3N::GetNodeNonlocalEqStrain] the interpolation order is only linear for the nonlocal eq strain.");
    return 0;
}

//! @brief returns the enum of the standard integration type for this element
NuTo::IntegrationType::eIntegrationType NuTo::Truss1D3N::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType1D2NGauss2Ip;
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Truss1D3N::ReorderNodes()
{
    std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[0];
    this->mNodes[0] = this->mNodes[2];
    this->mNodes[2] = tmp;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Truss1D3N::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<3; count++)
    {
        if (this->mNodes[count]==rOldPtr)
        {
            this->mNodes[count]=rNewPtr;
            break;
        }
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Truss1D3N::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D3N::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Truss1D3N::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Truss1D3N" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss1D);

    for(int i = 0; i < 3; i++)
    {
        std::uintptr_t& temp = reinterpret_cast<std::uintptr_t&>(mNodes[i]);
        ar & boost::serialization::make_nvp("mNode"+i, temp);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Truss1D3N" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Truss1D3N)
#endif // ENABLE_SERIALIZATION
