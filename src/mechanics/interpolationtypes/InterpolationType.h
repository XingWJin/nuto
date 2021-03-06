/*
 * InterpolationType.h
 *
 *  Created on: 31 Mar 2015
 *      Author: ttitsche
 */

#pragma once


#include <eigen3/Eigen/Core>
#include <boost/ptr_container/ptr_map.hpp>
#include <set>
#include <vector>

namespace NuTo
{
class InterpolationBase;
enum class eIntegrationType;

namespace Interpolation
{
enum class eShapeType;
enum class eTypeOrder;
} // namespace Interpolation

namespace Node
{
enum class eDof : unsigned char;
} // namespace Node

class InterpolationType
{
public:
    InterpolationType(NuTo::Interpolation::eShapeType rShapeType, int rDimension);

    virtual ~InterpolationType();

    const InterpolationBase& Get(const Node::eDof& rDofType) const;

    //! @brief adds a dof type and the corresponding interpolation of IGA type
    //! @param rDofType ... dof type
    void AddDofInterpolation(Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, const Eigen::VectorXi& rDegree,
                             const std::vector<Eigen::VectorXd>& rKnots, const Eigen::MatrixXd& rWeights);

    //! @brief adds a dof type and the corresponding interpolation order, calculate and store
    //! @param rDofType ... dof type
    //! @param rTypeOrder ... type and order of interpolation
    void AddDofInterpolation(Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    eIntegrationType GetStandardIntegrationType() const;

    //! @brief returns the shape type
    Interpolation::eShapeType GetShapeType() const;

    //! @brief clears the cached shape functions / N-matrices
    void ClearCache() const;

    //********************************************
    //             DOF METHODS
    //********************************************

    //! @brief defines whether or not the dof is active
    //! @param rIsActive ... true if active
    //! @param rDofType ... dof type
    void SetIsActive(bool rIsActive, Node::eDof rDofType);

    //! @brief returns whether or not the dof is active
    //! @param rDofType ... dof type
    bool IsActive(const Node::eDof& rDofType) const;

    //! @brief defines whether or not the dof is constitutive input
    //! @param rIsConstitutiveInput ... true if constitutive input
    //! @param rDofType ... dof type
    void SetIsConstitutiveInput(bool rIsConstitutiveInput, Node::eDof rDofType);

    //! @brief returns whether or not the dof is constitutive input
    //! @param rDofType ... dof type
    bool IsConstitutiveInput(const Node::eDof& rDofType) const;

    //! @brief returns true, if rDofType exists
    //! @param rDofType ... dof type
    bool IsDof(const Node::eDof& rDofType) const;

    //! @brief returns the number of dofs
    int GetNumDofs() const;

    //! @brief returns the number of active
    int GetNumActiveDofs() const;

    const std::set<Node::eDof>& GetActiveDofs() const;

    const std::set<Node::eDof>& GetDofs() const;

    std::set<Node::eDof> GetNodeDofs(int rNodeIndex) const;

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    const Eigen::VectorXd& GetNaturalNodeCoordinates(int rNodeIndex) const;

    int GetNumNodes() const;

    //! @brief returns the total number of nodes on a surface
    //! @param rSurface ... surface id
    int GetNumSurfaceNodes(int rSurface) const;

    //! @brief returns the node index of a specific DOF node on the surface
    //! @param rNodeIndex ... node index
    //! @param rSurface ... surface id
    int GetSurfaceNodeIndex(int rSurface, int rNodeIndex) const;

    //! @brief returns the node indices that span the surface
    //! 2 nodes for a 1D surface, 3 to 4 nodes on a 2D surface, (1 node for a 0D surface)
    //! calculates the shape functions of the parametrized surface and returns the index, where the value is 1
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... surface node indices
    Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const;

    //! @brief returns the number of surfaces
    int GetNumSurfaces() const;

    //! @brief returns the dof interpolation information as a string
    std::string Info() const;

    void PrintNodeIndices() const;
    void PrintNodeCoordinates() const;
    const Eigen::MatrixX2i& GetNodeRenumberingIndices() const;

private:
    //! @brief returns a nonconst reference to the object, stress that with the name
    //! @param rDofType ... dof type
    InterpolationBase& GetNonConst(Node::eDof rDofType);

    //! @brief returns whether or not the coordinate vectors rC1 and rC2 are equal
    //! @param rC1,rC2 ... coordinate vectors
    bool CoordinatesAreEqual(const Eigen::VectorXd& rC1, const Eigen::VectorXd& rC2) const;

    //! @brief Calculates index pairs that - if swapped - change the orientation of the element.
    //! It is implemented by reflecting each point at a plane at (0,0,0) with normal vector (1,-1,0) which is equal to
    //! swapping xi and eta coordinates
    //! @remark Different behavior for 1D: xi' = -xi. This could be done using polymorphism, but I think that bundling
    //! it here is sufficient.
    void UpdateNodeRenumberingIndices();

    Node::eDof GetDofWithHighestStandardIntegrationOrder() const;

    //! @brief map of single dof interpolations
    boost::ptr_map<Node::eDof, InterpolationBase> mInterpolations;

    //! @brief shape of the interpolation type
    const Interpolation::eShapeType mShapeType;

    //! @brief set of all dofs
    std::set<Node::eDof> mDofs;

    //! @brief set of active dofs
    std::set<Node::eDof> mActiveDofs;

    //! @brief number of dofs
    int mNumDofs;

    //! @brief number of active dofs
    int mNumActiveDofs;

    //! @brief contains a set of dofs for each local node
    std::vector<std::set<Node::eDof>> mNodeDofs;

    //! @brief contains local node coordinates
    std::vector<Eigen::VectorXd> mNodeCoordinates;

    //! @brief node renumbering indices that (if applied) change the orientation of the element
    Eigen::MatrixX2i mNodeRenumberingIndices;

    //! @brief vector (for each surface) of vectors (for each surface node) of surface node indices
    std::vector<std::vector<int>> mSurfaceNodeIndices;

    //! @brief dimension = Structure.GetDimension()
    const int mDimension;
};

} /* namespace NuTo */
