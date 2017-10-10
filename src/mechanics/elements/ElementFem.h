#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementInterface.h"
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/cell/Matrix.h"

namespace NuTo
{

class ElementFem : public ElementInterface
{
public:
    ElementFem(std::vector<NuTo::NodeSimple*> nodes, const InterpolationSimple& interpolation)
        : mNodes(nodes)
        , mInterpolation(interpolation)
    {
        assert(nodes.size() == interpolation.GetNumNodes());
        assert(nodes.front()->GetNumValues() == interpolation.GetDofDimension());
    }

    virtual NodeValues ExtractNodeValues() const override
    {
        const int dim = GetDofDimension();
        Eigen::VectorXd nodeValues(GetNumNodes() * dim);
        for (size_t i = 0; i < GetNumNodes(); ++i)
            nodeValues.segment(dim * i, dim) = mNodes[i]->GetValues();
        return nodeValues;
    }

    NMatrix GetNMatrix(NaturalCoords ipCoords) const override
    {
        return Matrix::N(Interpolation().GetShapeFunctions(ipCoords), Interpolation().GetNumNodes(),
                         Interpolation().GetDofDimension());
    }

    ShapeFunctions GetShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetShapeFunctions(ipCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(NaturalCoords ipCoords) const override
    {
        return Interpolation().GetDerivativeShapeFunctions(ipCoords);
    }

    int GetDofDimension() const override
    {
        return Interpolation().GetDofDimension();
    }

    int GetNumNodes() const override
    {
        return mNodes.size(); 
    }

    const InterpolationSimple& Interpolation() const
    {
        return mInterpolation;
    }

private:

    std::vector<NuTo::NodeSimple*> mNodes;
    std::reference_wrapper<const InterpolationSimple> mInterpolation;
};
} /* NuTo */
