/*
 * Element3D.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/Element3D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/math/FullMatrix.h"

NuTo::Element3D::Element3D(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType,
        IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIpDataType, rInterpolationType), mNodes(rNodes)
{
    mSection = 0;
    CheckElement();
}

NuTo::Error::eError NuTo::Element3D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    if (mStructure->GetHessianConstant(1) == false)
        throw MechanicsException("[NuTo::Element3D::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2) == false)
        throw MechanicsException("[NuTo::Element3D::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
        const SectionBase* section(GetSection());
        if (section == 0)
            throw MechanicsException("[NuTo::Element3D::Evaluate] no section allocated for element.");

        const std::set<Node::eAttributes>& dofs = mInterpolationType->GetDofs();
        const std::set<Node::eAttributes>& activeDofs = mInterpolationType->GetActiveDofs();

        int numActiveDofs = mInterpolationType->GetNumActiveDofs();

        // extract all node values and store them
        std::map<Node::eAttributes, Eigen::MatrixXd> nodalValues;
        for (auto dof : dofs)
        {
            ExtractNodeValues(nodalValues[dof], 0, dof);
        }

        //define inputs and outputs
        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*> constitutiveInputList;
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*> constitutiveOutputList;

        EngineeringStress3D engineeringStress3D;
        EngineeringStrain3D engineeringStrain3D;

        EngineeringStrain3D engineeringPlasticStrain3D;

        NuTo::ConstitutiveTangentLocal<6, 6> tangentStressStrain;

        DeformationGradient3D deformationGradient;

        //allocate damage
        Damage damage;

        //for the lumped mass calculation
        double total_mass = 0.;

        for (auto dof : dofs)
        {
            if (mInterpolationType->IsConstitutiveInput(dof) == false)
                continue;
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D] = &(deformationGradient);
            }
                break;
            default:
                throw MechanicsException(
                        "[NuTo::Element3D::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
            }
        }

        //define outputs
        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
                //on the global level
                if (mStructure->GetHessianConstant(0) == false)
                {
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &(engineeringStress3D);
                        }
                            break;
                        default:
                            throw MechanicsException(
                                    "[NuTo::Element3D::Evaluate] Constitutive output INTERNAL_GRADIENT for " + Node::AttributeToString(dof)
                                            + " not implemented.");

                        }
                    }
                }
                break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->GetFullMatrixDouble().setZero();
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D] = &tangentStressStrain;
                    }
                        break;
                    default:
                        throw MechanicsException(
                                "[NuTo::Element3D::Evaluate] Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof)
                                        + " not implemented.");

                    }
                }
            }
                break;
            case Element::HESSIAN_1_TIME_DERIVATIVE:
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                throw MechanicsException("[NuTo::Element3D::Evaluate] Constitutive output HESSIAN_1_TIME_DERIVATIVE not implemented.");
                break;
            case Element::HESSIAN_2_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                //there is only a constant mass part for the mechanics problem
            }
                break;
            case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                break;
            case Element::UPDATE_STATIC_DATA:
                constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
                break;
            case Element::UPDATE_TMP_STATIC_DATA:
                constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = 0;
                break;
            case Element::IP_DATA:
                switch (it->second->GetIpDataType())
                {
                case NuTo::IpData::ENGINEERING_STRAIN:
                    it->second->GetFullMatrixDouble().Resize(6, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D] = &(engineeringStrain3D);
                    break;
                case NuTo::IpData::ENGINEERING_STRESS:
                    it->second->GetFullMatrixDouble().Resize(6, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &(engineeringStress3D);
                    break;
                case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
                    it->second->GetFullMatrixDouble().Resize(6, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D] = &(engineeringPlasticStrain3D);
                    break;
                case NuTo::IpData::DAMAGE:
                    it->second->GetFullMatrixDouble().Resize(1, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::DAMAGE] = &(damage);
                    break;
                default:
                    throw MechanicsException("[NuTo::Element3D::Evaluate] this ip data type is not implemented.");
                }
                break;
            case Element::GLOBAL_ROW_DOF:
            {
                const Eigen::VectorXi& globalRowDofsEigen = CalculateGlobalRowDofs();
                std::vector<int> globalRowDofsStd(globalRowDofsEigen.data(), globalRowDofsEigen.data() + globalRowDofsEigen.rows());
                it->second->GetVectorInt() = globalRowDofsStd;
            }
                break;
            case Element::GLOBAL_COLUMN_DOF:
            {
                const Eigen::VectorXi& globalColumnDofsEigen = CalculateGlobalColumnDofs();
                std::vector<int> globalColumnDofsStd(globalColumnDofsEigen.data(), globalColumnDofsEigen.data() + globalColumnDofsEigen.rows());
                it->second->GetVectorInt() = globalColumnDofsStd;
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element3D::Evaluate] element output not implemented.");
            }
        }

        Eigen::Matrix3d invJacobian;
        double detJacobian;

        std::map<Node::eAttributes, Eigen::VectorXd> shapeFunctions;
        std::map<Node::eAttributes, Eigen::MatrixXd> derivativeShapeFunctions;

        // loop over the integration points
        for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
        {

            // calculate Jacobian
            const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural =
                    mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);

            this->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, nodalValues[Node::COORDINATES], invJacobian, detJacobian);

            // calculate shape functions and their derivatives

            for (auto dof : dofs)
            {
                if (dof == Node::COORDINATES)
                    continue;
                const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
                shapeFunctions[dof] = interpolationType.GetShapeFunctions(theIP);
                // this lazy product here is so much faster than any other implementation via a seperate method
                // possibly due to more efficient mallocs
                derivativeShapeFunctions[dof] = interpolationType.GetDerivativeShapeFunctionsNatural(theIP).lazyProduct(invJacobian);
            }

            // define constitutive inputs

            for (auto dof : dofs)
            {
                if (mInterpolationType->IsConstitutiveInput(dof) == false)
                    continue;
                switch (dof)
                {
                case Node::DISPLACEMENTS:
                {
//                    constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D] = &(deformationGradient);
                    deformationGradient = CalculateDeformationGradient(derivativeShapeFunctions.at(dof), nodalValues.at(dof));
                }
                    break;
                default:
                    throw MechanicsException(
                            "[NuTo::Element3D::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
                }
            }

            ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
            try
            {
                Error::eError error = constitutivePtr->Evaluate3D(this, theIP, constitutiveInputList, constitutiveOutputList);
                if (error != Error::SUCCESSFUL)
                    return error;
            } catch (NuTo::MechanicsException &e)
            {
                e.AddMessage("[NuTo::Element3D::Evaluate] error evaluating the constitutive model.");
                throw e;
            }

            //calculate output
            for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
            {
                switch (it->first)
                {
                case Element::INTERNAL_GRADIENT:
                {
                    //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
                    //on the global level
                    if (mStructure->GetHessianConstant(0) == false)
                    {
                        double factor(fabs(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
                        for (auto dof : activeDofs)
                        {
                            int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                            switch (dof)
                            {
                            case Node::DISPLACEMENTS:
                            {
                                AddDetJBtSigma(derivativeShapeFunctions.at(dof), engineeringStress3D, factor, startIndex, it->second->GetFullVectorDouble());
                            }
                                break;
                            default:
                                throw MechanicsException(
                                        "[NuTo::Element3D::Evaluate] Element output INTERNAL_GRADIENT for " + Node::AttributeToString(dof)
                                                + " not implemented.");

                            }
                        }
                    }
                }
                    break;
                case Element::HESSIAN_0_TIME_DERIVATIVE:
                {
                    //factor for the numerical integration
                    double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
                    if (tangentStressStrain.GetConstant() == false)
                        it->second->SetConstant(false);
                    if (tangentStressStrain.GetSymmetry() == false)
                        it->second->SetSymmetry(false);

                    for (auto dof : activeDofs)
                    {
                        int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            AddDetJBtCB(derivativeShapeFunctions.at(dof), tangentStressStrain, factor, startIndex, startIndex,
                                    it->second->GetFullMatrixDouble());
                        }
                            break;
                        default:
                            throw MechanicsException(
                                    "[NuTo::Element3D::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof)
                                            + " not implemented.");

                        }
                    }
                }
                    break;
                case Element::HESSIAN_1_TIME_DERIVATIVE:
                {
                    throw MechanicsException("[NuTo::Element3D::Evaluate] Element output HESSIAN_1_TIME_DERIVATIVE not implemented.");

                }
                    break;
                case Element::HESSIAN_2_TIME_DERIVATIVE:
                {
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {

                            double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))
                                            * constitutivePtr->GetDensity());
                            const Eigen::MatrixXd tmpMatrix = shapeFunctions.at(dof) * shapeFunctions.at(dof).transpose() * factor;
                            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result(it->second->GetFullMatrixDouble());
                            for (int count = 0; count < tmpMatrix.rows(); count++)
                            {
                                for (int count2 = 0; count2 < tmpMatrix.cols(); count2++)
                                {
                                    result(3*count,  3*count2)   += tmpMatrix(count,count2);
                                    result(3*count+1,3*count2+1) += tmpMatrix(count,count2);
                                    result(3*count+2,3*count2+2) += tmpMatrix(count,count2);
                                }
                            }
                        }
                            break;
                        default:
                            throw MechanicsException(
                                    "[NuTo::Element3D::Evaluate] Element output HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof)
                                            + " not implemented.");

                        }
                    }
                }
                    break;
                case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {

                            // calculate local mass matrix (the nonlocal terms are zero)
                            // don't forget to include determinant of the Jacobian and area
                            // detJ * area * density * HtH, :
                            double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))
                                            * constitutivePtr->GetDensity());
                            FullVector<double, Eigen::Dynamic>& result(it->second->GetFullVectorDouble());
                            total_mass += factor;
                            const Eigen::VectorXd& shapeFunction = shapeFunctions.at(dof);
                            //calculate for the translational dofs the diagonal entries
                            for (int count = 0; count < shapeFunction.rows(); count++)
                            {
                                result(3 * count) += shapeFunction[count] * shapeFunction[count] * factor;
                            }

                            if (theIP + 1 == GetNumIntegrationPoints())
                            {
                                //calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
                                double sum_diagonal(0);
                                for (int count = 0; count < shapeFunction.rows(); count++)
                                {
                                    sum_diagonal += result(3 * count);
                                }

                                //scale so that the sum of the diagonals represents the full mass
                                double scaleFactor = total_mass / sum_diagonal;
                                for (int count = 0; count < shapeFunction.rows(); count++)
                                {
                                    result(3*count) *= scaleFactor;
                                    result(3*count+1) = result(3*count);
                                    result(3*count+2) = result(3*count);
                                }
                            }
                        }
                            break;
                        default:
                            throw MechanicsException(
                                    "[NuTo::Element3D::Evaluate] Element output LUMPED_HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof)
                                            + " not implemented.");
                        }
                    }

                    break;
                case Element::UPDATE_STATIC_DATA:
                case Element::UPDATE_TMP_STATIC_DATA:
                    break;
                case Element::IP_DATA:
                    switch (it->second->GetIpDataType())
                    {
                    case NuTo::IpData::ENGINEERING_STRAIN:
                        //error = constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP * 6]), engineeringStrain3D.GetData(), 6 * sizeof(double));
                        break;
                    case NuTo::IpData::ENGINEERING_STRESS:
                        //error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP * 6]), engineeringStress3D.GetData(), 6 * sizeof(double));
                        break;
                    case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
                        //error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP * 6]), engineeringPlasticStrain3D.GetData(), 6 * sizeof(double));
                        break;
                    case NuTo::IpData::DAMAGE:
                        //error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]), damage.GetData(), sizeof(double));
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element3D::GetIpData] Ip data not implemented.");
                    }
                    break;
                case Element::GLOBAL_ROW_DOF:
                case Element::GLOBAL_COLUMN_DOF:
                    break;
                default:
                    throw MechanicsException("[NuTo::Element3D::Evaluate] element output not implemented.");
                }
            }

        }


    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::Element3D::Evaluate] Error evaluating element data of element" + ss.str() + ".");
        throw e;
    }
    return Error::SUCCESSFUL;
}

NuTo::Element::eElementType NuTo::Element3D::GetEnumType() const
{
    return NuTo::Element::ELEMENT3D;
}

int NuTo::Element3D::GetGlobalDimension() const
{
    return 3;
}

NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

const NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber, Node::eAttributes rDofType)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

const NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber, Node::eAttributes rDofType) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

void NuTo::Element3D::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    assert(rNode != nullptr);
    mNodes[rLocalNodeNumber] = rNode;
}

void NuTo::Element3D::ResizeNodes(int rNewNumNodes)
{
   if (rNewNumNodes == (int)mNodes.size())
       return;

   if (rNewNumNodes > (int)mNodes.size())
   {
       // just resize (enlarge)
       mNodes.resize(rNewNumNodes);
   }
   else
   {
       throw MechanicsException("[NuTo::Element3D::ResizeNodes] Resize that reduces the number of nodes is not implemented yet.");
   }
}

void NuTo::Element3D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (unsigned int count = 0; count < mNodes.size(); count++)
    {
        if (this->mNodes[count] == rOldPtr)
        {
            this->mNodes[count] = rNewPtr;
            break;
        }
    }
}

void NuTo::Element3D::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

const NuTo::SectionBase* NuTo::Element3D::GetSection() const
{
    return mSection;
}

NuTo::ConstitutiveStaticDataBase* NuTo::Element3D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain3D(this);
}

const Eigen::VectorXd NuTo::Element3D::GetIntegrationPointVolume() const
{

    Eigen::MatrixXd localNodeCoord = this->ExtractNodeValues(0, Node::COORDINATES);

    const InterpolationBase& interpolationType = mInterpolationType->Get(Node::COORDINATES);

    double detJac;
    Eigen::Matrix3d dummyJacobian;

    Eigen::VectorXd volume(GetNumIntegrationPoints());

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = interpolationType.GetDerivativeShapeFunctionsNatural(theIP);

        CalculateJacobian(derivativeShapeFunctionsNatural, localNodeCoord, dummyJacobian, detJac);

        volume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

const Eigen::MatrixXd NuTo::Element3D::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofs = interpolationTypeDof.GetNumDofs();
    int numDofsPerNode = numDofs / numNodes;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodalValues;
    nodalValues.resize(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            nodalValues.block<3, 1>(0, iNode) = node->GetCoordinates3D();
            break;

        case Node::DISPLACEMENTS:
            nodalValues.block<3, 1>(0, iNode) = node->GetDisplacements3D(rTimeDerivative);
            break;

        case Node::TEMPERATURES:
            nodalValues(0, iNode) = node->GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            nodalValues(0, iNode) = node->GetNonlocalEqStrain(rTimeDerivative);
            break;

        default:
            throw MechanicsException("[NuTo::Element3D::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) +".");
        }
    }
    return nodalValues;
}

void NuTo::Element3D::ExtractNodeValues(Eigen::MatrixXd& rNodeValues, int rTimeDerivative, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofs = interpolationTypeDof.GetNumDofs();
    int numDofsPerNode = numDofs / numNodes;

    rNodeValues.resize(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            rNodeValues.block<3, 1>(0, iNode) = node->GetCoordinates3D();
            break;

        case Node::DISPLACEMENTS:
            rNodeValues.block<3, 1>(0, iNode) = node->GetDisplacements3D(rTimeDerivative);
            break;

        case Node::TEMPERATURES:
            rNodeValues(0, iNode) = node->GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            rNodeValues(0, iNode) = node->GetNonlocalEqStrain(rTimeDerivative);
            break;

        default:
            throw MechanicsException("[NuTo::Element3D::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) +".");
        }
    }
}

const Eigen::VectorXi NuTo::Element3D::CalculateGlobalRowDofs() const
{
    int numActiveDos = mInterpolationType->GetNumActiveDofs();
    Eigen::VectorXi globalRowDofs(numActiveDos);

    for (auto dof : mInterpolationType->GetActiveDofs())
    {
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        int index = interpolationType.GetLocalStartIndex(); //

        for (int iNodeDof = 0; iNodeDof < interpolationType.GetNumNodes(); ++iNodeDof)
        {
            const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(0);
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(1);
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(2);
            }
                break;
            case Node::TEMPERATURES:
            {
                globalRowDofs[index++] = nodePtr->GetDofTemperature();
            }
                break;
            case Node::NONLOCALEQSTRAIN:
            {
                globalRowDofs[index++] = nodePtr->GetDofNonlocalEqStrain();
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element3D::CalculateGlobalRowDofs] Not implemented for " + Node::AttributeToString(dof) +".");

            }
        }
    }
    return globalRowDofs;
}

const Eigen::VectorXi NuTo::Element3D::CalculateGlobalColumnDofs() const
{
    int numNonlocalElements = GetNumNonlocalElements();
    if (numNonlocalElements == 0)
        return CalculateGlobalRowDofs();
    else
    {
        throw NuTo::MechanicsException("[NuTo::Element3D::CalculateGlobalColumnDofs] not implemented for nonlocal integral element formulation.");
    }
}

const NuTo::DeformationGradient3D NuTo::Element3D::CalculateDeformationGradient(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal,
        const Eigen::MatrixXd& rNodalDisplacements) const
{
    DeformationGradient3D deformationGradient;

    assert(rDerivativeShapeFunctionsLocal.rows() == rNodalDisplacements.cols());
    assert(rDerivativeShapeFunctionsLocal.cols() == rNodalDisplacements.rows());

    Eigen::Matrix3d deformationGradientEigen = rNodalDisplacements.lazyProduct(rDerivativeShapeFunctionsLocal);
    deformationGradientEigen.noalias() += Eigen::Matrix3d::Identity();
    Eigen::Map<Eigen::Matrix3d>(deformationGradient.mDeformationGradient) = deformationGradientEigen;

    return deformationGradient;
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element (that's the thermal solution)
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Element3D::AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal, const ConstitutiveTangentLocal<6, 6>& rConstitutiveTangent,
        double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{
    const double *C = rConstitutiveTangent.data();
    double x1,x2,y1,y2,z1,z2,x2x1,y2x1,z2x1,x2y1,y2y1,z2y1,x2z1,y2z1,z2z1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        assert((int)rDerivativeShapeFunctionsGlobal.size()>node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal(theNode1,0);
        y1 = rFactor * rDerivativeShapeFunctionsGlobal(theNode1,1);
        z1 = rFactor * rDerivativeShapeFunctionsGlobal(theNode1,2);
        node1mul3 +=rRow;
        node1mul3plus1 +=rRow;
        node1mul3plus2 +=rRow;
        for (int theNode2=0; theNode2<GetNumNodes(); theNode2++)
        {
            int node2mul3 = 3*theNode2;
            int node2mul3plus1 = node2mul3+1;
            int node2mul3plus2 = node2mul3plus1+1;
            node2mul3 +=rCol;
            node2mul3plus1 +=rCol;
            node2mul3plus2 +=rCol;

            assert((int)rDerivativeShapeFunctionsGlobal.size()>node2mul3plus2);
            x2 = rDerivativeShapeFunctionsGlobal(theNode2,0);
            y2 = rDerivativeShapeFunctionsGlobal(theNode2,1);
            z2 = rDerivativeShapeFunctionsGlobal(theNode2,2);

            x2x1 = x2*x1;
            y2x1 = y2*x1;
            z2x1 = z2*x1;
            x2y1 = x2*y1;
            y2y1 = y2*y1;
            z2y1 = z2*y1;
            x2z1 = x2*z1;
            y2z1 = y2*z1;
            z2z1 = z2*z1;
            assert(rCoefficientMatrix.GetNumRows()>node1mul3plus2 && rCoefficientMatrix.GetNumColumns()>node1mul3plus2);
            assert(rCoefficientMatrix.GetNumRows()>node2mul3plus2 && rCoefficientMatrix.GetNumColumns()>node2mul3plus2);

            rCoefficientMatrix(node1mul3,node2mul3)          +=x2x1*C[0] +x2y1*C[3] +x2z1*C[5] +y2x1*C[18]+y2y1*C[21]+y2z1*C[23]+z2x1*C[30]+z2y1*C[33]+z2z1*C[35];
            rCoefficientMatrix(node1mul3,node2mul3plus1)     +=y2x1*C[6] +y2y1*C[9] +y2z1*C[11]+x2x1*C[18]+x2y1*C[21]+x2z1*C[23]+z2x1*C[24]+z2y1*C[27]+z2z1*C[29];
            rCoefficientMatrix(node1mul3,node2mul3plus2)     +=z2x1*C[12]+z2y1*C[15]+z2z1*C[17]+y2x1*C[24]+y2y1*C[27]+y2z1*C[29]+x2x1*C[30]+x2y1*C[33]+x2z1*C[35];
            rCoefficientMatrix(node1mul3plus1,node2mul3)     +=x2y1*C[1] +x2x1*C[3] +x2z1*C[4] +y2y1*C[19]+y2x1*C[21]+y2z1*C[22]+z2y1*C[31]+z2x1*C[33]+z2z1*C[34];
            rCoefficientMatrix(node1mul3plus1,node2mul3plus1)+=y2y1*C[7] +y2x1*C[9] +y2z1*C[10]+x2y1*C[19]+x2x1*C[21]+x2z1*C[22]+z2y1*C[25]+z2x1*C[27]+z2z1*C[28];
            rCoefficientMatrix(node1mul3plus1,node2mul3plus2)+=z2y1*C[13]+z2x1*C[15]+z2z1*C[16]+y2y1*C[25]+y2x1*C[27]+y2z1*C[28]+x2y1*C[31]+x2x1*C[33]+x2z1*C[34];
            rCoefficientMatrix(node1mul3plus2,node2mul3)     +=x2z1*C[2] +x2y1*C[4] +x2x1*C[5] +y2z1*C[20]+y2y1*C[22]+y2x1*C[23]+z2z1*C[32]+z2y1*C[34]+z2x1*C[35];
            rCoefficientMatrix(node1mul3plus2,node2mul3plus1)+=y2z1*C[8] +y2y1*C[10]+y2x1*C[11]+x2z1*C[20]+x2y1*C[22]+x2x1*C[23]+z2z1*C[26]+z2y1*C[28]+z2x1*C[29];
            rCoefficientMatrix(node1mul3plus2,node2mul3plus2)+=z2z1*C[14]+z2y1*C[16]+z2x1*C[17]+y2z1*C[26]+y2y1*C[28]+y2x1*C[29]+x2z1*C[32]+x2y1*C[34]+x2x1*C[35];
        }
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Element3D::AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const EngineeringStress3D& rEngineeringStress, double rFactor,
        int rRow, FullVector<double, Eigen::Dynamic>& rResult) const
{
    const double *s = rEngineeringStress.GetData();

    double x1,y1,z1;
    for (int theNode1=0; theNode1<rDerivativeShapeFunctionsLocal.rows(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        assert(rResult.GetNumRows()>node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1,0);
        y1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1,1);
        z1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1,2);

        rResult(rRow + node1mul3)     +=x1*s[0]+y1*s[3]+z1*s[5];
        rResult(rRow + node1mul3plus1)+=y1*s[1]+x1*s[3]+z1*s[4];
        rResult(rRow + node1mul3plus2)+=z1*s[2]+y1*s[4]+x1*s[5];
    }
}

void NuTo::Element3D::ReorderNodes()
{
    // swap all nodes
//    unsigned int lastNode = mNodes.size() - 1;
//    for (unsigned int iNode = 0; iNode < lastNode; ++iNode, --lastNode)
//    {
////        std::cout << "Swapping node " << iNode << " and " << lastNode << "." << std::endl;
//        NodeBase* tmpNode = mNodes[iNode];
//        mNodes[iNode] = mNodes[lastNode];
//        mNodes[lastNode] = tmpNode;
//    }


}

void NuTo::Element3D::CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::MatrixXd& rNodeCoordinates,
        Eigen::Matrix3d& rInvJacobian, double& rDetJac) const
{
    /*       jacobian
     j0, j2,
     j1, j3 */

    assert(rDerivativeShapeFunctions.rows() == rNodeCoordinates.cols());
    assert(rDerivativeShapeFunctions.cols() == 3);
    assert(rNodeCoordinates.rows() == 3);

    Eigen::Matrix3d jacobian = rNodeCoordinates.lazyProduct(rDerivativeShapeFunctions);
    rDetJac = jacobian.determinant();

    if (rDetJac==0)
        throw MechanicsException("[NuTo::Element3D::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    rInvJacobian = jacobian.inverse();
}

void NuTo::Element3D::CheckElement()
{

    int numIntegrationPoints = GetNumIntegrationPoints();
    // check number of integration points
    if (numIntegrationPoints < 1)
    {
        throw MechanicsException("[NuTo::Element3D::CheckElement] invalid integration type.");
    }

    int theIP = 0;
    const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    Eigen::Matrix3d invJacobian;
    double detJacobian;

    CalculateJacobian(derivativeShapeFunctions, nodeCoordinates, invJacobian, detJacobian);
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    }

    double volume = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        CalculateJacobian(derivativeShapeFunctions, nodeCoordinates, invJacobian, detJacobian);
        if (detJacobian <= 0)
            throw MechanicsException("[NuTo::Element3D::CheckElement] Determinant of the Jacobian <= zero, no inversion possible.");
        volume += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Element3D::CheckElement] element with zero volume (check nodes).");
    }
}

