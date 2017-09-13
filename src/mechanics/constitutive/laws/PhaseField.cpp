#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"

#include "base/Logger.h"
#include "base/Exception.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#include "mechanics/nodes/NodeEnum.h"

#include <eigen3/Eigen/Core>

#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"




NuTo::PhaseField::PhaseField(const double youngsModulus, const double poissonsRatio, const double lengthScaleParameter,
                             const double fractureEnergy, const double artificialViscosity)
    : ConstitutiveBase()
    , mYoungsModulus(youngsModulus)
    , mPoissonsRatio(poissonsRatio)
    , mLengthScaleParameter(lengthScaleParameter)
    , mFractureEnergy(fractureEnergy)
    , mArtificialViscosity(artificialViscosity)
    , mLameLambda((youngsModulus * poissonsRatio) / ((1 + poissonsRatio) * (1 - 2 * poissonsRatio)))
    , mLameMu(youngsModulus / (2 * (1 + poissonsRatio)))
{
    assert(mYoungsModulus > 0.0);
    assert(mLengthScaleParameter > 0.0);
    assert(mFractureEnergy > 0.0);
    assert(mArtificialViscosity >= 0.0);
    assert(mPoissonsRatio >= 0.0 and mPoissonsRatio < 0.5);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::PhaseField::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLaw<PhaseField>>(*this, 0.);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
NuTo::ConstitutiveInputMap NuTo::PhaseField::GetConstitutiveInputs(const ConstitutiveOutputMap&) const
{
    ConstitutiveInputMap constitutiveInputMap;

    constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
    constitutiveInputMap[Constitutive::eInput::CRACK_PHASE_FIELD];

    return constitutiveInputMap;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double NuTo::PhaseField::Evaluate2DIsotropic(const double oldEnergyDensity,
                                             const ConstitutiveInputMap& rConstitutiveInput,
                                             const ConstitutiveOutputMap& rConstitutiveOutput)
{

    const auto& engineeringStrain =
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
    const auto& damage = *rConstitutiveInput.at(Constitutive::eInput::CRACK_PHASE_FIELD);

    // calculate coefficients
    double C11, C12, C33;
    std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mYoungsModulus, mPoissonsRatio);

    // calculate the effective stress
    Eigen::Vector3d effectiveStress;
    effectiveStress[0] = (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
    effectiveStress[1] = (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
    effectiveStress[2] = C33 * engineeringStrain[2];

    // calculate the elastic energy density
    const double elasticEnergyDensity = 0.5 * effectiveStress.dot(engineeringStrain);

    double currentEnergyDensity = std::max(elasticEnergyDensity, oldEnergyDensity);

    bool performUpdateAtEnd = false;
    constexpr double residualEnergyDensity = 1.e-12;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;
            engineeringStress[0] = factor * effectiveStress[0];
            engineeringStress[1] = factor * effectiveStress[1];
            engineeringStress[2] = factor * effectiveStress[2];

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

            // right coefficients are calculated above
            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;

            tangent(0, 0) = factor * C11;
            tangent(1, 0) = factor * C12;
            tangent(2, 0) = 0.;

            tangent(0, 1) = factor * C12;
            tangent(1, 1) = factor * C11;
            tangent(2, 1) = 0.;

            tangent(0, 2) = 0.;
            tangent(1, 2) = 0.;
            tangent(2, 2) = factor * C33;

            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);


            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;

            // plane strain
            engineeringStress3D[0] = factor * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
            engineeringStress3D[1] = factor * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
            engineeringStress3D[2] = factor * C12 * (engineeringStrain[0] + engineeringStrain[1]);
            engineeringStress3D[3] = 0.;
            engineeringStress3D[4] = 0.;
            engineeringStress3D[5] = factor * C33 * engineeringStrain[2];

            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() =
                    engineeringStrain.As3D(mPoissonsRatio, ePlaneState::PLANE_STRAIN);
            break;
        }

        case NuTo::Constitutive::eOutput::ELASTIC_ENERGY_DAMAGED_PART:
        {

            ConstitutiveIOBase& elasticEnergyLoadTerm = *itOutput.second;

            elasticEnergyLoadTerm[0] = currentEnergyDensity;

            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_PHASE_FIELD:
        {

            ConstitutiveIOBase& dStressDPhaseField = *itOutput.second;
            dStressDPhaseField.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);


            const double degradationFactor = -2 * (1 - damage[0]);


            dStressDPhaseField[0] = degradationFactor * effectiveStress[0];
            dStressDPhaseField[1] = degradationFactor * effectiveStress[1];
            dStressDPhaseField[2] = degradationFactor * effectiveStress[2];


            break;
        }

        case NuTo::Constitutive::eOutput::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN:
        {
            // the tangent  is equal to the damaged part of the engineering stress for loading and zero for unloading
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);


            if (elasticEnergyDensity == currentEnergyDensity)
            {
                // loading
                tangent[0] = effectiveStress[0];
                tangent[1] = effectiveStress[1];
                tangent[2] = effectiveStress[2];
            }
            else
            {
                // unloading
                tangent.SetZero();
            }

            break;
        }

        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        {
            throw Exception(__PRETTY_FUNCTION__,
                            "tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            continue;

        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            continue;

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    // return old/new history variables
    if (performUpdateAtEnd)
    {
        return currentEnergyDensity;
    }

    return oldEnergyDensity;
}

namespace NuTo
{
// template specialization needs an explicit namespace block...
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
void NuTo::PhaseField::Evaluate<1>(const ConstitutiveInputMap&, const ConstitutiveOutputMap&, Data&)
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
void NuTo::PhaseField::Evaluate<2>(const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput, Data& rStaticData)
{
    const auto& planeState =
            *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());
    if (planeState.GetPlaneState() != ePlaneState::PLANE_STRAIN)
        throw Exception(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found.");

    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::eInput::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw Exception(__PRETTY_FUNCTION__,
                        "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData =
            dynamic_cast<const ConstitutiveCalculateStaticData&>(*itCalculateStaticData->second);
    int index = calculateStaticData.GetIndexOfPreviousStaticData();

    double oldEnergyDensity = rStaticData.GetData(index);

    double energyDensity = 0.0;

    energyDensity = Evaluate2DIsotropic(oldEnergyDensity, rConstitutiveInput, rConstitutiveOutput);


    // update history variables
    rStaticData.SetData(energyDensity);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
void NuTo::PhaseField::Evaluate<3>(const ConstitutiveInputMap&, const ConstitutiveOutputMap&, Data&)
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}
}// namespace NuTo

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool NuTo::PhaseField::CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 0:
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::CRACKPHASEFIELD):
        case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
            return true;
        default:
            return false;
        }
        break;
    case 1:
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
            return true;
        default:
            return false;
        }
        break;
    default:
        return false;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double NuTo::PhaseField::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{

    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY:
        return mArtificialViscosity;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return mPoissonsRatio;
    case Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER:
        return mLengthScaleParameter;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return mYoungsModulus;

    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NuTo::PhaseField::SetParameterDouble(Constitutive::eConstitutiveParameter, double)
{
    throw Exception(__PRETTY_FUNCTION__, "Function must not be used. Use CTOR.");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
NuTo::Constitutive::eConstitutiveType NuTo::PhaseField::GetType() const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NuTo::PhaseField::CheckParameters() const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NuTo::PhaseField::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus         : " << mYoungsModulus << "\n";
    rLogger << "    Poisson's ratio         : " << mPoissonsRatio << "\n";
    rLogger << "    fracture energy         : " << mFractureEnergy << "\n";
    rLogger << "    artificial viscosity    : " << mArtificialViscosity << "\n";
    rLogger << "    length scale parameter  : " << mLengthScaleParameter << "\n";
}
