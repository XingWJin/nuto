#pragma once

#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
namespace NuTo
{


/*! @brief This file contains a collection of static functions to
 *  calculate equivalent strains and their derivative with to
 *  strains.
 *
 *  A class is used (instead of a pure namespace) to hide some private
 *  functions.
 */

template <int TDim>
class EquivalentStrainModifiedMises
{
public:
    //! @param rStrain ... engineering strain
    //! @param rK ... k parameter, compressiveStrength / tensileStrength
    //! @param rNu ... poisson ratio
    //! @param rSectionType ... only needed for 2D: Plane strain/ plane stress
    EquivalentStrainModifiedMises(const EngineeringStrain<TDim>& rStrain, double rK, double rNu,
                                  ePlaneState planeState = ePlaneState::PLANE_STRESS);

    //! @brief calculates the modified mises equivalent strain

    //! @return modified mises equivalent strain
    double Get() const;

    ConstitutiveVector<ConstitutiveIOBase::GetVoigtDim(TDim)> GetDerivative() const;

private:
    double mK1;
    double mK2;
    EngineeringStrain<3> mStrain3D;
    double mI1;
    double mJ2;
    double mA;
    double mNu;
    ePlaneState mPlaneState = ePlaneState::PLANE_STRESS;
};

} // namespace NuTo
