#pragma once
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationQuadLinear : public InterpolationSimple
{
public:
    InterpolationQuadLinear(int rDofDimension)
        : mDofDimension(rDofDimension)
    {
    }

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationQuadLinear>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::ShapeFunctionsQuadOrder1(rNaturalIPCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(rNaturalIPCoords);
    }

    NaturalCoords GetLocalCoords(int rNodeId) const override
    {
        return ShapeFunctions2D::NodeCoordinatesQuadOrder1(rNodeId);
    }

    int GetNumNodes() const override
    {
        return 4;
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    int mDofDimension;
};
} /* NuTo */
