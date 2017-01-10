#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constraints/ConstraintLinearNodeRotations2D.h"
#include "math/SparseMatrixCSRGeneral.h"

NuTo::ConstraintLinearNodeRotations2D::ConstraintLinearNodeRotations2D(const NodeBase* rNode, double rValue) :
        ConstraintNode(rNode), ConstraintLinear()
{
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeRotations2D::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeRotations2D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeRotations2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    if (mNode->GetNum(Node::eDof::ROTATIONS)!=1)
        throw MechanicsException("[NuTo::ConstraintLinearNodeRotations2D::AddToConstraintMatrix] Node does not have a rotation component.");
    rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDof(Node::eDof::ROTATIONS, 0),1);

    curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeRotations2D::GetRHS(int& curConstraintEquation,Eigen::VectorXd& rRHS)const
{
    rRHS(curConstraintEquation) = mRHS;
    curConstraintEquation++;
}

NuTo::Node::eDof NuTo::ConstraintLinearNodeRotations2D::GetDofType() const
{
    return Node::eDof::ROTATIONS;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLinearNodeRotations2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeRotations2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeRotations2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeRotations2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeRotations2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeRotations2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearNodeRotations2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeRotations2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNode)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRHS);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeRotations2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeRotations2D)

void NuTo::ConstraintLinearNodeRotations2D::SetNodePtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mNodeMapCast)
{
    NuTo::ConstraintNode::SetNodePtrAfterSerialization(mNodeMapCast);
}
#endif // ENABLE_SERIALIZATION
