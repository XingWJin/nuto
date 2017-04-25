// $Id$ 

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "visualize/Component.h"
#include "visualize/VisualizeException.h"


NuTo::Visualize::Component::Component(NuTo::eVisualizeWhat visualizeComponent) :
        mComponent(visualizeComponent)
{
}

NuTo::eVisualizeWhat NuTo::Visualize::Component::GetComponentEnum() const
{
    return mComponent;
}

std::string NuTo::Visualize::Component::GetComponentName() const
{
    switch (mComponent)
    {
    case eVisualizeWhat::ACCELERATION:
        return "Accelerations";
    case eVisualizeWhat::ANGULAR_ACCELERATION:
        return "AngularAccelerations";
    case eVisualizeWhat::ANGULAR_VELOCITY:
        return "AngularVelocities";
    case eVisualizeWhat::BOND_STRESS:
        return "BondStress";
    case eVisualizeWhat::DAMAGE:
        return "Damage";
    case eVisualizeWhat::CRACK_PHASE_FIELD:
        return "CrackPhaseField";
    case eVisualizeWhat::DISPLACEMENTS:
        return "Displacements";
    case eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        return "EngineeringPlasticStrain";
    case eVisualizeWhat::ENGINEERING_STRAIN:
        return "EngineeringStrain";
    case eVisualizeWhat::ENGINEERING_STRESS:
        return "EngineeringStress";
    case eVisualizeWhat::HEAT_FLUX:
        return "HeatFlux";
    case eVisualizeWhat::LATTICE_STRAIN:
        return "LatticeStrain";
    case eVisualizeWhat::LATTICE_STRESS:
        return "LatticeStress";
    case eVisualizeWhat::LOCAL_EQ_STRAIN:
        return "LocalEqStrain";
    case eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        return "NonlocalEqStrain";
    case eVisualizeWhat::PARTICLE_RADIUS:
        return "ParticleRadius";
    case eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        return "PrincipalEngineeringStress";
    case eVisualizeWhat::RELATIVE_HUMIDITY:
        return "RelativeHumidity";
    case eVisualizeWhat::ROTATION:
        return "Rotations";
    case eVisualizeWhat::SHRINKAGE_STRAIN:
        return "ShrinkageStrains";
    case eVisualizeWhat::SLIP:
        return "Slip";
    case eVisualizeWhat::TEMPERATURE:
        return "Temperature";
    case eVisualizeWhat::THERMAL_STRAIN:
        return "ThermalStrain";
    case eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
        return "TotalInelasticEqStrain";
    case eVisualizeWhat::VELOCITY:
        return "Velocities";
    case eVisualizeWhat::WATER_VOLUME_FRACTION:
        return "WaterVolumeFraction";
    case eVisualizeWhat::ELECTRIC_POTENTIAL:
        return "ElectricPotential";
    case eVisualizeWhat::ELECTRIC_FIELD:
        return "ElectricField";
    case eVisualizeWhat::ELECTRIC_DISPLACEMENT:
        return "ElectricDisplacement";
    default:
        throw VisualizeException(__PRETTY_FUNCTION__, "Visualization component not implemented.");
    }
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Visualize::Component::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Visualize::Component::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Visualize::Component::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Visualize::Component::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Visualize::Component::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Visualize::Component::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Visualize::Component::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Visualize::Component" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Visualize::Component" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Visualize::Component)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Visualize::Component)
#endif // ENABLE_SERIALIZATION