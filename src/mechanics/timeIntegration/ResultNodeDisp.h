// $Id: $

#pragma once


#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/timeIntegration/ResultNodeDof.h"


namespace NuTo
{

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class ResultNodeDisp : public ResultNodeDof
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ResultNodeDisp(const std::string& rIdent, int rNodeId);

    //! @brief calculate the relevant nodal dofs
    void CalculateValues(const StructureBase& rStructure, Eigen::Matrix<double, 1, Eigen::Dynamic>& rValues)const override;

    //! @brief number of data points per time step (e.g. number of displacement components of a node
    int GetNumData(const StructureBase& rStructure)const override;

    NuTo::eTimeIntegrationResultType GetResultType()const override;

    std::string GetTypeId() const
    {
    	return std::string("ResultNodeDisp");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const
    {

    }

protected:
};
}

