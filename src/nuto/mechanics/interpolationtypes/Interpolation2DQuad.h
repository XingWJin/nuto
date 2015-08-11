/*
 * Interpolation2DQuad.h
 *
 *  Created on: 24 Apr 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION2DQUAD_H_
#define INTERPOLATION2DQUAD_H_

#include "nuto/mechanics/interpolationtypes/Interpolation2D.h"

namespace NuTo
{

/**
@brief 2D quadrilateral element with the following natural coordinate system and its surface parametrization
\f[\fbox{ \begin{tikzpicture}
  \draw[dotted, -latex] (0,0) -- (1.5,0) node[above]{$\xi$};
  \draw[dotted,-latex] (0,0) -- (0,1.5) node[above]{$\eta$};
  \draw[dashed] (-1,-1) node[below  left]{$(-1,-1)$}
              --( 1,-1) node[below right]{$( 1,-1)$}
              --( 1, 1) node[above right]{$( 1, 1)$}
              --(-1, 1) node[above  left]{$(-1, 1)$} -- cycle;
  \draw[|-latex] ( 0,-1) -- (.8,-1) node[midway, below] {$\alpha_0$};
  \draw[|-latex] ( 1, 0) -- ( 1,.8) node[midway, above right] {$\alpha_1$};
  \draw[|-latex] ( 0, 1) -- (-.8, 1) node[midway, below] {$\alpha_2$};
  \draw[|-latex] (-1, 0) -- (-1,-.8) node[midway, left] {$\alpha_3$};
  \node[rectangle, draw] at (1.5,2) {$\alpha_i = [-1,1]$};
\end{tikzpicture}   } \f]
**/
class Interpolation2DQuad: public Interpolation2D
{
public:
    Interpolation2DQuad(const StructureBase* rStructure, Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    IntegrationType::eIntegrationType GetStandardIntegrationType() const override;

    //! @brief returns the natural coordinates of the dof node
    //! @param rDofType ... dof type
    //! @param rNodeIndexDof ... node index of the dof type
    const Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    const Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns derivative shape functions in the local coordinate system
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... map of derivative shape functions in the natural coordinate system for all dofs
    const Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    const Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    //! @brief returns the derivative of the surface parametrization
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... derivative of the surface parametrization
    const Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    //! @brief returns the number of surfaces
    int GetNumSurfaces() const override
    {
        return 4;
    }

private:
    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */

#endif /* INTERPOLATION2DQUAD_H_ */
