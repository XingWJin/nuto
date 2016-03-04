// $Id: LoadSurface3D.h 178 2009-12-11 20:53:12Z eckardt4 $
#ifndef LoadSurfaceBase2D_H
#define LoadSurfaceBase2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadBase.h"

namespace NuTo
{
class NodeBase;
class Element2D;
class StructureBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all surface loads in 3D
class LoadSurfaceBase2D : public LoadBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadSurfaceBase2D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId);

    //! @brief just for serialization
    LoadSurfaceBase2D(){ }

    //! @brief adds the load to global sub-vectors
    //! @param rLoadCase number of the current load case
    //! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const;

    //! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates ... global coordinates
    //! @param rNormal ... normal to the surface (pointing outwards)
    //! @param rLoadVector ... load vector
    virtual void CalculateSurfaceLoad(NuTo::FullVector<double,2>& rCoordinates,NuTo::FullVector<double,2>& rNormal,
    		NuTo::FullVector<double,2>& rLoadVector)const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief deserializes (load) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start load LoadSurfaceBase2D\n";
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadBase);

        std::vector<std::pair<std::uintptr_t, int> >  mElements2DAdress;
        ar & boost::serialization::make_nvp("mNodesAdress", mElements2DAdress);
        for(std::vector<std::pair<std::uintptr_t, int> >::const_iterator it = mElements2DAdress.begin(); it != mElements2DAdress.end(); it++)
        {
            const Element2D* tempElement2D = reinterpret_cast<const Element2D* >(it->first);
            std::pair<const Element2D*, int> tempPair(tempElement2D, it->second);
            mElements2D.push_back(tempPair);
        }

        ar & BOOST_SERIALIZATION_NVP(mIntegrationType2NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType3NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType4NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType5NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType3NPtrLobatto)
           & BOOST_SERIALIZATION_NVP(mIntegrationType4NPtrLobatto)
           & BOOST_SERIALIZATION_NVP(mIntegrationType5NPtrLobatto);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish load LoadSurfaceBase2D\n";
#endif
    }

    //! @brief serializes (saves) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start save LoadSurfaceBase2D\n";
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadBase);

        std::vector<std::pair<std::uintptr_t, int> >  mElements2DAdress;
        for(std::vector<std::pair<const Element2D*, int> >::const_iterator it = mElements2D.begin(); it != mElements2D.end(); it++)
        {
            std::uintptr_t tempAdressElement2D = reinterpret_cast<std::uintptr_t >(it->first);
            std::pair<std::uintptr_t, int> tempPair(tempAdressElement2D, it->second);
            mElements2DAdress.push_back(tempPair);
        }
        ar & boost::serialization::make_nvp("mNodesAdress", mElements2DAdress);

        ar & BOOST_SERIALIZATION_NVP(mIntegrationType2NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType3NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType4NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType5NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType3NPtrLobatto)
           & BOOST_SERIALIZATION_NVP(mIntegrationType4NPtrLobatto)
           & BOOST_SERIALIZATION_NVP(mIntegrationType5NPtrLobatto);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish load LoadSurfaceBase2D\n";
#endif
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Adress (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the adress in the map
    //! @param mNodeMapCast   std::map containing the old and new adresses
    virtual void SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast) override
    {
        for(std::vector<std::pair<const Element2D*, int> >::const_iterator it = mElements2D.begin(); it != mElements2D.end(); it++)
        {
            std::uintptr_t temp = reinterpret_cast<std::uintptr_t>(it->first);
            std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mElementMapCast.find(temp);
            if(itCast!=mElementMapCast.end())
            {
                Element2D** tempPtr = const_cast<Element2D**>(&(it->first));
                *tempPtr = reinterpret_cast<Element2D*>(itCast->second);
            }
            else
                throw MechanicsException("[NuTo::LoadSurfaceBase2D] The Element2D-Pointer could not be updated.");
        }
    }
#endif // ENABLE_SERIALIZATION

protected:
    std::vector<std::pair<const Element2D*, int> > mElements2D;
    IntegrationTypeBase* mIntegrationType2NPtr;
    IntegrationTypeBase* mIntegrationType3NPtr;
    IntegrationTypeBase* mIntegrationType4NPtr;
    IntegrationTypeBase* mIntegrationType5NPtr;
    IntegrationTypeBase* mIntegrationType3NPtrLobatto;
    IntegrationTypeBase* mIntegrationType4NPtrLobatto;
    IntegrationTypeBase* mIntegrationType5NPtrLobatto;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LoadSurfaceBase2D)
#endif

#endif //LoadSurfaceBase2D_H

