#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/elements/Element1DInXD.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/elements/EvaluateDataContinuum.h"


NuTo::Element1DInXD::Element1DInXD(const std::vector<NuTo::NodeBase*>& rNodes,
                                   const InterpolationType& rInterpolationType,
                                   const IntegrationTypeBase& integrationType, const DofStatus& dofStatus,
                                   int globalDimension)
    : NuTo::ContinuumElement<1>(rNodes, rInterpolationType, integrationType, dofStatus)
    , mGlobalDimension(globalDimension)
{
    mRotationMatrix = CalculateRotationMatrix();
}


Eigen::MatrixXd NuTo::Element1DInXD::CalculateRotationMatrix()
{
    const Eigen::VectorXd globalNodeCoordinates = ExtractGlobalNodeValues(0, Node::eDof::COORDINATES);

    assert((mGlobalDimension == 2 or mGlobalDimension == 3) and "Need 2d or 3d coordinates");

    // the rotation matrix consists of three basis vectors that define the new coordinate system.
    Eigen::Vector3d basisVector00 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector01 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector02 = Eigen::Vector3d::Zero();

    // basisVector00 is aligned with the truss. Note: Trusses have to be straight!
    basisVector00.block(0, 0, mGlobalDimension, 1) =
            globalNodeCoordinates.tail(mGlobalDimension) - globalNodeCoordinates.head(mGlobalDimension);

    // check if basisVector00 is linear independent to unitVectorZ and calculate the cross product to determine
    // basisVector01
    if (basisVector00(0, 0) > 1e-8 or basisVector00(1, 0) > 1e-8)
    {
        const Eigen::Vector3d unitVectorZ = Eigen::Vector3d::UnitZ();
        basisVector01 = basisVector00.cross(unitVectorZ);
    }
    else
    {
        const Eigen::Vector3d unitVectorX = Eigen::Vector3d::UnitX();
        basisVector01 = basisVector00.cross(unitVectorX);
    }

    // basisVector02 is othogonal to basisVector01 and basisVector00
    basisVector02 = basisVector00.cross(basisVector01);

    Eigen::Matrix3d rotationMatrix3d = Eigen::Matrix3d::Zero();

    basisVector00.normalize();
    basisVector01.normalize();
    basisVector02.normalize();

    rotationMatrix3d.col(0) = basisVector00;
    rotationMatrix3d.col(1) = basisVector01;
    rotationMatrix3d.col(2) = basisVector02;

    return rotationMatrix3d.block(0, 0, mGlobalDimension, mGlobalDimension);
}

Eigen::MatrixXd NuTo::Element1DInXD::CalculateTransformationMatrix(unsigned int rGlobalDimension,
                                                                   unsigned int rNumberOfNodes) const
{

    Eigen::MatrixXd transformationMatrix(rGlobalDimension * rNumberOfNodes, rGlobalDimension * rNumberOfNodes);
    transformationMatrix.setZero();

    for (unsigned int i = 0; i < rNumberOfNodes; ++i)
    {
        transformationMatrix.block(rGlobalDimension * i, rGlobalDimension * i, rGlobalDimension, rGlobalDimension) =
                mRotationMatrix;
    }


    return transformationMatrix;
}

Eigen::VectorXd NuTo::Element1DInXD::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    Eigen::VectorXd globalNodeValues = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    const unsigned int numNodes = GetInterpolationType().Get(rDofType).GetNumNodes();

    Eigen::VectorXd nodeValues(numNodes);

    for (unsigned iNode = 0; iNode < numNodes; ++iNode)
    {
        Eigen::VectorXd tmp =
                (mRotationMatrix.transpose() * globalNodeValues.segment(iNode * mGlobalDimension, mGlobalDimension));
        nodeValues[iNode] = tmp(0, 0);
    }

    return nodeValues;
}

const Eigen::VectorXd NuTo::Element1DInXD::ExtractGlobalNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{

    const InterpolationBase& interpolationTypeDof = GetInterpolationType().Get(rDofType);
    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode = GetNumDofsPerNode(rDofType);

    Eigen::VectorXd globalNodeValues(numDofsPerNode * numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase& node = *GetNode(iNode, rDofType);

        switch (rDofType)
        {
        case Node::eDof::COORDINATES:
            globalNodeValues.segment(iNode * numDofsPerNode, mGlobalDimension) = node.Get(Node::eDof::COORDINATES);
            break;
        case Node::eDof::DISPLACEMENTS:
            globalNodeValues.segment(iNode * numDofsPerNode, mGlobalDimension) =
                    node.Get(Node::eDof::DISPLACEMENTS, rTimeDerivative);
            break;
        default:
            throw Exception(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(rDofType));
        }
    }

    return globalNodeValues;
}

void NuTo::Element1DInXD::CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0,
                                                         EvaluateDataContinuum<1>& rData, int,
                                                         const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            auto& hessian0 = rHessian0(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            {
                const unsigned int numberOfNodes = mInterpolationType->Get(Node::eDof::DISPLACEMENTS).GetNumNodes();
                const unsigned int numberOfDofs = mGlobalDimension * numberOfNodes;

                Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(mGlobalDimension, numberOfNodes);

                const auto& tangentStressStrain = *static_cast<ConstitutiveMatrix<1, 1>*>(
                        constitutiveOutputMap.at(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)
                                .get());
                const Eigen::MatrixXd localHessian0 = rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose() *
                                                      tangentStressStrain * rData.mB.at(dofRow);
                Eigen::MatrixXd globalHessian0(numberOfDofs, numberOfDofs);
                globalHessian0.setZero();

                for (unsigned iRow = 0; iRow < localHessian0.rows(); ++iRow)
                {
                    for (unsigned iCol = 0; iCol < localHessian0.cols(); ++iCol)
                    {
                        globalHessian0(mGlobalDimension * iRow, mGlobalDimension * iCol) = localHessian0(iRow, iCol);
                    }
                }
                hessian0 += transformationMatrix * globalHessian0 * transformationMatrix.transpose();
            }
            break;
            default:
                throw Exception(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for "
                                                     "(" + Node::DofToString(dofRow) +
                                                             "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

void NuTo::Element1DInXD::CalculateElementOutputInternalGradient(
        BlockFullVector<double>& rInternalGradient, EvaluateDataContinuum<1>& rData, int, const ConstitutiveInputMap&,
        const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
        {
            const unsigned int numberOfNodes = mInterpolationType->Get(dofRow).GetNumNodes();

            Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(mGlobalDimension, numberOfNodes);

            const auto& engineeringStress = *static_cast<EngineeringStress<1>*>(
                    constitutiveOutputMap.at(Constitutive::eOutput::ENGINEERING_STRESS).get());
            const Eigen::VectorXd localInternalGradient =
                    rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose() * engineeringStress;
            Eigen::VectorXd globalInternalGradient(mGlobalDimension * numberOfNodes);
            globalInternalGradient.setZero();

            for (unsigned iNode = 0; iNode < numberOfNodes; ++iNode)
                globalInternalGradient[iNode * mGlobalDimension] = localInternalGradient[iNode];


            rInternalGradient[dofRow] += transformationMatrix * globalInternalGradient;
        }
        break;
        default:
            throw Exception(__PRETTY_FUNCTION__,
                            "Element output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");
        }
    }
}

int NuTo::Element1DInXD::GetNumDofsPerNode(Node::eDof rDofType) const
{
    switch (rDofType)
    {
    case NuTo::Node::eDof::COORDINATES:
        return mGlobalDimension;
    case NuTo::Node::eDof::DISPLACEMENTS:
        return mGlobalDimension;
    default:
        throw NuTo::Exception("[NuTo::Element1DInXD::GetNumDofsPerNode] dof type not found.");
    }
}

Eigen::VectorXd NuTo::Element1DInXD::InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates,
                                                          NuTo::Node::eDof rDofType) const
{
    return InterpolateDofGlobal(0, rNaturalCoordinates, rDofType);
}
Eigen::VectorXd NuTo::Element1DInXD::InterpolateDofGlobal(int rTimeDerivative,
                                                          const Eigen::VectorXd& rNaturalCoordinates,
                                                          NuTo::Node::eDof rDofType) const
{

    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    const Eigen::VectorXd nodalValues = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    const Eigen::MatrixXd matrixNLocal = interpolationType.MatrixN(rNaturalCoordinates);

    int numNodes = GetNumNodes(rDofType);
    int dimBlock = GetNumDofsPerNode(rDofType);


    Eigen::MatrixXd matrixNGlobal(dimBlock, numNodes * dimBlock);
    for (int iNode = 0, iBlock = 0; iNode < numNodes; ++iNode, iBlock += dimBlock)
    {
        matrixNGlobal.block(0, iBlock, dimBlock, dimBlock) =
                Eigen::MatrixXd::Identity(dimBlock, dimBlock) * matrixNLocal(iNode);
    }


    return matrixNGlobal * nodalValues;
}


void NuTo::Element1DInXD::CheckElement()
{
    unsigned int numIntegrationPoints = GetNumIntegrationPoints();

    // check number of integration points
    assert(numIntegrationPoints > 0);

    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::eDof::COORDINATES);

    double length = 0;
    for (unsigned int iIp = 0; iIp < numIntegrationPoints; ++iIp)
    {
        const auto coords = GetIntegrationType().GetLocalIntegrationPointCoordinates(iIp);
        const Eigen::MatrixXd& derivativeShapeFunctions =
                mInterpolationType->Get(Node::eDof::COORDINATES).DerivativeShapeFunctionsNatural(coords);
        Eigen::Matrix<double, 1, 1> detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates);
        assert(detJacobian(0, 0) > 0 and "Jacobian needs to be greater than 0");

        length += this->GetIntegrationPointWeight(iIp) * detJacobian(0, 0);
    }

    // check element length
    assert(length > 1e-14 and "element with zero length (check nodes)");
}
