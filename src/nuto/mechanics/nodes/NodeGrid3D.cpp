// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeGrid3D.h"

//! @brief constructor
NuTo::NodeGrid3D::NodeGrid3D(int rNodeGridNum)  : NodeBase ()
{
	mNodeGridNum = rNodeGridNum;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeGrid3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeGrid3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeGrid3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeGrid3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeGrid3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeGrid3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
void NuTo::NodeGrid3D::serialize(Archive & ar, const unsigned int version)
{
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
	& BOOST_SERIALIZATION_NVP(mNodeGridNum);
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeGrid3D)
BOOST_CLASS_TRACKING(NuTo::NodeGrid3D, track_always)

#endif // ENABLE_SERIALIZATION


//! @brief returns the number of coordinates of the node
//! @return number of coordinates
int NuTo::NodeGrid3D::GetNumCoordinates()const
{
	return 3;
}

//! @brief set the grid node number of this node
//! @param rNodeGridNum  given grid node number
void NuTo::NodeGrid3D::SetNodeGridNum(int rNodeGridNum)
{
	mNodeGridNum = rNodeGridNum;
}

//! @brief get the grid node number of this node
//! @return NodeGridNum  grid node number
int NuTo::NodeGrid3D::GetNodeGridNum()const
{
	return mNodeGridNum;
}

//! @brief Get the Number of belonging elements
//! @return Number of belonging elements
int NuTo::NodeGrid3D::GetNumElems()
{
	return mNumElems;
}

//!@brief Set Number of belonging elements
//! @param Number of belonging elements
void NuTo::NodeGrid3D::SetNumElems(int rNumElems)
{
	mNumElems=rNumElems;
}


//! @brief Get Ids of the elements
//! @retutn int * Ids of the elements
int* NuTo::NodeGrid3D::GetElementIds()
{
	return mElementIds;
}

//! @brief Set Ids of the elements
//! @param int * Ids of the elements
void NuTo::NodeGrid3D::SetElementIds(int * rElementIds)
{
	mElementIds=new int[8];
	for (int count=0;count<8;++count)
		mElementIds[count]=rElementIds[count];

}

//! @brief writes the coordinates of a node to the prescribed pointer
//! @param rCoordinates coordinates
void NuTo::NodeGrid3D::GetCoordinates3D(double rCoordinates[3])const
{
	throw MechanicsException("[NuTo::NodeGrid3D::GetCoordinates3D] to be implemented.");
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeGrid3D::SetGlobalDofs(int& rDOF)
{
	//empty since coordinates are no DOFs
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid3D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid3D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid3D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid3D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid3D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	//empty since coordinates are no DOFs
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeGrid3D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	//empty since coordinates are no DOFs
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeGrid3D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	//empty since coordinates are no DOFs
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeGrid3D::GetNodeTypeStr()const
{
	return std::string("NodeGrid3D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeGrid3D::GetNodeType()const
{
    return Node::NodeGrid3D;
}
