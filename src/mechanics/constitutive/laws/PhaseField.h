#pragma once

#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/staticData/DataContainer.h"

#include <functional>
namespace NuTo
{

class ConstitutiveIOBase;
class Logger;
class ElementBase;
class IPConstitutiveLawBase;

//! \class  PhaseField
//! \author Philip Huschke
//! \date   June 14, 2016
//! @brief  A phase-field model for brittle fracture

//! Recommended literature:
//!
//! Miehe et al. \n
//! "Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE
//! implementations"
//!
//! Ambati et al. \n
//! "A review on phase-field models of brittle fracture and a new fast hybrid formulation"
class PhaseField : public ConstitutiveBase
{
public:
    using StaticDataType = double;
    using Data = typename Constitutive::StaticData::DataContainer<double>;

    //! @brief      Constructor
    //! @param[in]  youngsModulus Young's Modulus
    //! @param[in]  poissonsRatio Poisson's Ratio
    //! @param[in]  lengthScaleParameter Parameter that corresponds to the band-width of the diffusive crack
    //! @param[in]  fractureEnergy Fracture energy/critical energy release rate
    //! @param[in]  artificialViscosity Parameter to improve robustness of the model (non-physical)
    PhaseField(double youngsModulus, double poissonsRatio, double lengthScaleParameter,
               double fractureEnergy, double artificialViscosity);

    //! @brief creates corresponding IPConstitutiveLaw
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override;

    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const override;

    //! @brief Evaluate the constitutive law
    //! @param[in] rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.)
    //! @param[out] rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.)
    //! @param[in] rStaticData static data
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput,
                  Data& rStaticData);

    //! @brief Determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow row dof
    //! @param rDofCol column dof
    //! @param rTimeDerivative time derivative
    bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

    //! @brief Gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief Sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief Gets the type of constitutive relationship
    //! @return Type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief Check parameters of the constitutive relationship
    void CheckParameters() const override;

    //! @brief Print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief Returns true, if a material model has tmp static data (which has to be updated before stress or stiffness
    //! are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
    {
        return false;
    }

    double Evaluate2DIsotropic(double oldEnergyDensity, const ConstitutiveInputMap& rConstitutiveInput,
                               const ConstitutiveOutputMap& rConstitutiveOutput);

protected:
    //! @brief Young's modulus \f$ E \f$
    const double mYoungsModulus;

    //! @brief Poisson's ratio \f$ \nu \f$
    const double mPoissonsRatio;

    //! @brief Length scale parameter \f$ l \f$
    const double mLengthScaleParameter;

    //! @brief Fracture energy \f$ G_f \f$
    const double mFractureEnergy;

    //! @brief Artificial viscosity to improve numerical robustness \f$ \eta \f$
    const double mArtificialViscosity;

    //! @brief First Lame parameter \f$ \lambda \f$
    const double mLameLambda;

    //! @brief Second Lame parameter \f$ \mu \f$
    const double mLameMu;

};
} // namespace NuTo
