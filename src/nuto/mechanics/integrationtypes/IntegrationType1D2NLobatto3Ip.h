// $Id$
#ifndef IntegrationType1D2NLobatto3Ip_H
#define IntegrationType1D2NLobatto3Ip_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#include "nuto/mechanics/integrationtypes/IntegrationType1D.h"

namespace NuTo
{
//! @author Jörg F. Unger, BAM
//! @date November 2013
//! @brief ... integration types in 1D with two nodes Lobatto integration and 3 integration points
class IntegrationType1D2NLobatto3Ip : public IntegrationType1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType1D2NLobatto3Ip();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType1D)
        & BOOST_SERIALIZATION_NVP(iPts)
        & BOOST_SERIALIZATION_NVP(weights);
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    void GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints()const;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum)const;

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    std::string GetStrIdentifier()const;

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    static std::string GetStrIdentifierStatic();

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP)const;
#endif // ENABLE_VISUALIZE
private:
    //! @brief ... integration points coordinates
    double iPts[3];

    //! @brief ... weights for the integration
    double weights[3];
};
} //namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType1D2NLobatto3Ip)
#endif

#endif //IntegrationType1D2NLobatto3Ip_H
