#pragma once

#include "mechanics/cell/CellInterface.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/elements/ElementSimple.h"
#include "mechanics/nodes/DofContainer.h"
#include "mechanics/cell/Integrand.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
template <int TDim>
class Cell : public CellInterface
{
public:
    Cell(const ElementSimple& rCoordinateElement, DofContainer<ElementSimple*> rElements,
         const IntegrationTypeBase& rIntegrationType, const Integrand<TDim>& rIntegrand)
        : mCoordinateElement(rCoordinateElement)
        , mElements(rElements)
        , mIntegrationType(rIntegrationType)
        , mIntegrand()
    {
        for (int i = 0; i < rIntegrationType.GetNumIntegrationPoints(); i++)
            mIntegrand.push_back(rIntegrand.Clone());
    }

    //! @brief builds the internal gradien
    DofVector<double> Gradient() override
    {
        DofVector<double> gradient;
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            Jacobian<TDim> jacobian(mCoordinateElement.ExtractNodeValues(),
                                    mCoordinateElement.GetInterpolation().GetDerivativeShapeFunctions(ipCoords));
            CellIPData<TDim> cellipData(mElements, jacobian, ipCoords);
            gradient += mIntegrand[iIP].Gradient(cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return gradient;
    }

    DofVector<int> DofNumbering() override
    {
        return DofVector<int>();
    }

    //! @brief Extracts a vector (each IP) of vectors (several IPValues for the same integrion point) of IPValues
    std::vector<std::vector<IPValue>> IPValues() override
    {
        std::vector<std::vector<IPValue>> ipValues;

        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian<TDim> jacobian(mCoordinateElement.ExtractNodeValues(),
                                    mCoordinateElement.GetInterpolation().GetDerivativeShapeFunctions(ipCoords));
            CellIPData<TDim> cellipData(mElements, jacobian, ipCoords);
            ipValues.push_back(mIntegrand[iIP].IPValues(cellData, cellipData));
        }
        return ipValues;
    }

private:
    const ElementSimple& mCoordinateElement;
    DofContainer<ElementSimple*> mElements;
    const IntegrationTypeBase& mIntegrationType;
    boost::ptr_vector<Integrand<TDim>> mIntegrand;
};
} /* NuTo */
