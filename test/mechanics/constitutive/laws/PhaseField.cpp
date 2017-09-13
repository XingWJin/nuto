#include "BoostUnitTest.h"

#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include <iostream>

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eInput;
using NuTo::Constitutive::eOutput;

double EvaluatePhaseFieldModel(const Eigen::Vector3d& rStrain, const double rCrackPhaseField,
                               NuTo::PhaseField* rPhaseField)
{
    constexpr int VoigtDim = NuTo::ConstitutiveIOBase::GetVoigtDim(2);

    NuTo::ConstitutiveInputMap myConstitutiveInputMap;
    myConstitutiveInputMap[eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<NuTo::ConstitutiveCalculateStaticData>(NuTo::eCalculateStaticData::EULER_BACKWARD);

    myConstitutiveInputMap[eInput::ENGINEERING_STRAIN] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eInput::ENGINEERING_STRAIN);
    auto& engineeringStrain =
            *static_cast<NuTo::EngineeringStrain<2>*>(myConstitutiveInputMap.at(eInput::ENGINEERING_STRAIN).get());
    engineeringStrain.AsVector() = rStrain;

    myConstitutiveInputMap[eInput::CRACK_PHASE_FIELD] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eInput::CRACK_PHASE_FIELD);
    auto& damage = *static_cast<NuTo::ConstitutiveScalar*>(myConstitutiveInputMap.at(eInput::CRACK_PHASE_FIELD).get());
    damage[0] = rCrackPhaseField;

    cout << "Inputs:" << endl;
    cout << "engineeringStrain.AsVector()" << endl << engineeringStrain.AsVector() << endl;
    cout << "damage.AsScalar()" << endl << damage.AsScalar() << endl;

    NuTo::ConstitutiveOutputMap myConstitutiveOutputMap;
    myConstitutiveOutputMap[eOutput::ENGINEERING_STRESS] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eOutput::ENGINEERING_STRESS);
    auto& engineeringStress =
            *static_cast<NuTo::EngineeringStress<2>*>(myConstitutiveOutputMap.at(eOutput::ENGINEERING_STRESS).get());
    engineeringStress.SetZero();

    myConstitutiveOutputMap[eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN);
    auto& tangentStressStrain = *static_cast<NuTo::ConstitutiveMatrix<VoigtDim, VoigtDim>*>(
            myConstitutiveOutputMap.at(eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
    tangentStressStrain.setZero();

    double oldEnergyDensity = 0.0;
    rPhaseField->Evaluate2DIsotropic(oldEnergyDensity, myConstitutiveInputMap,
                                                            myConstitutiveOutputMap);

    cout << "Outputs:" << endl;
    cout << "engineeringStress.AsVector()" << endl << engineeringStress.AsVector() << endl;
    cout << "tangentStressStrain.AsMatrix()" << endl << tangentStressStrain.AsMatrix() << endl;

    Eigen::Vector3d residual = tangentStressStrain.AsMatrix() * engineeringStrain - engineeringStress;

    return residual.norm();
}

BOOST_AUTO_TEST_CASE(check_d_stress_d_strain)
{
    constexpr double tolerance = 1.e-10;
    constexpr double youngsModulus = 2715;
    constexpr double poissonsRatio = 0.27;
    constexpr double lengthScale = 0.76;
    constexpr double fractureEnergy = 2.71;
    constexpr double artificialViscosity = 0.13;

    NuTo::PhaseField* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScale, fractureEnergy,
                                                        artificialViscosity);

    Eigen::Vector3d strain;
    strain[0] = 0.;
    strain[1] = 0.;
    strain[2] = 0.;

    double crackPhaseField = 0.;

    double residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 0.;
    strain[1] = 0.;
    strain[2] = 0.;

    crackPhaseField = 1.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.;
    strain[1] = 0.;
    strain[2] = 0.;

    crackPhaseField = 0.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = -1.;
    strain[1] = 0.;
    strain[2] = 0.;

    crackPhaseField = 0.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = -1.;
    strain[1] = 0.;
    strain[2] = 0.;

    crackPhaseField = 1.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.;
    strain[1] = 0.;
    strain[2] = 0.;

    crackPhaseField = 1.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.;
    strain[1] = 1.;
    strain[2] = 0.;

    crackPhaseField = 0.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.;
    strain[1] = 1.;
    strain[2] = 0.;

    crackPhaseField = 1.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 2.;
    strain[1] = 1.;
    strain[2] = 0.;

    crackPhaseField = 1.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 2.;
    strain[1] = -1.;
    strain[2] = 0.;

    crackPhaseField = 1.;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 0.;
    strain[1] = 0.;
    strain[2] = 2.;

    crackPhaseField = 1.0;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.;
    strain[1] = 1.;
    strain[2] = 2.;

    crackPhaseField = 1.0;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);

    cout << "**********************************************" << endl;

    strain[0] = 1.;
    strain[1] = -3.;
    strain[2] = 2.;

    crackPhaseField = 1.0;

    residualNorm = EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
    BOOST_CHECK_SMALL(residualNorm, tolerance);
}
