// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeAccelerations2D.h"

// default constructor
NuTo::NodeAccelerations2D::NodeAccelerations2D() : NodeBase ()
{
	this->mAccelerations[0] = 0.0;
	this->mAccelerations[1] = 0.0;
}

// constructor
NuTo::NodeAccelerations2D::NodeAccelerations2D(const double rAccelerations[2]) : NodeBase ()
{
	this->mAccelerations[0] = rAccelerations[0];
	this->mAccelerations[1] = rAccelerations[1];
}

// returns the number of accelerations of the node
int NuTo::NodeAccelerations2D::GetNumAccelerations()const
{
	return 2;
}

// set the accelerations
void NuTo::NodeAccelerations2D::SetAccelerations2D(const double rAccelerations[2])
{
	this->mAccelerations[0] = rAccelerations[0];
	this->mAccelerations[1] = rAccelerations[1];
}

// writes the accelerations of a node to the prescribed pointer
void NuTo::NodeAccelerations2D::GetAccelerations2D(double rAccelerations[2])const
{
	rAccelerations[0] = this->mAccelerations[0];
	rAccelerations[1] = this->mAccelerations[1];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeAccelerations2D::GetNodeTypeStr()const
{
	return std::string("NodeAccelerations2D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeAccelerations2D::GetNodeType()const
{
    return NuTo::Node::NodeAccelerations2D;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeAccelerations2D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mAccelerations);
}
#endif // ENABLE_SERIALIZATION
