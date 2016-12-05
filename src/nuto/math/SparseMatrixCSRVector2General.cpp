// $Id: SparseMatrixCSRGeneral.cpp 195 2009-12-16 09:13:29Z eckardt4 $

#include <string>

#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSRVector2General<double>::GetTypeId()const
{
    return std::string("SparseMatrixCSRVector2GeneralDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSRVector2General<int>::GetTypeId()const
{
    return std::string("SparseMatrixCSRVector2GeneralInt");
}

#ifdef ENABLE_SERIALIZATION
template<class T>
void SparseMatrixCSRVector2General<T>::Save ( const std::string &filename, std::string rType)const
{
    try
    {
        // transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int (*)(int))toupper);

        // open file
        std::ofstream ofs(filename.c_str(), std::ios_base::binary);
        if (!ofs.is_open())
        {
            throw MathException(__PRETTY_FUNCTION__, "Error opening file.");
        }

        // write data to file
        std::string typeIdString(this->GetTypeId());
        if (rType == "BINARY")
        {
            boost::archive::binary_oarchive oba(ofs, std::ios::binary);
            oba& boost::serialization::make_nvp("Object_type", typeIdString);
            oba& boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType == "XML")
        {
            boost::archive::xml_oarchive oxa(ofs, std::ios::binary);
            std::string tmpString(this->GetTypeId());
            oxa& boost::serialization::make_nvp("Object_type", typeIdString);
            oxa& boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType == "TEXT")
        {
            boost::archive::text_oarchive ota(ofs, std::ios::binary);
            ota& boost::serialization::make_nvp("Object_type", typeIdString);
            ota& boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MathException(__PRETTY_FUNCTION__, "File type not implemented.");
        }

        // close file
        ofs.close();
    }
    catch (boost::archive::archive_exception& e)
    {
        std::string s(__PRETTY_FUNCTION__ + "File save exception in boost - " + e.what());
        throw MathException(s);
    }
    catch (MathException& e)
    {
        throw;
    }
    catch (std::exception& e)
    {
        throw MathException(e.what());
    }
}

template<class T>
void SparseMatrixCSRVector2General<T>::Restore ( const std::string &filename,  std::string rType)
{
	try
	{
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

		// open file
		std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
		if(! ifs.is_open())
		{
			throw MathException("[NuTo::SparseMatrixCSRVector2General::Restore] Error opening file.");
		}

		std::string typeIdString;
		if (rType=="BINARY")
		{
			boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
			oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
			if ( typeIdString != this->GetTypeId() )
			{
				throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
			}
			oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		}
		else if (rType=="XML")
		{
			boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
			oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
			if ( typeIdString != this->GetTypeId() )
			{
				throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
			}
			oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		}
		else if (rType=="TEXT")
		{
			boost::archive::text_iarchive ota ( ifs, std::ios::binary );
			ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
			if ( typeIdString != this->GetTypeId() )
			{
				throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
			}
			ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		}
		else
		{
			throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Restore]File type not implemented" );
		}
		// close file
		ifs.close();
	}
	catch ( boost::archive::archive_exception& e )
	{
		std::string s ( std::string ( "[NuTo::SparseMatrixCSRVector2General::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
        throw;
	}
	catch ( std::exception &e )
	{
		throw MathException ( e.what() );
	}
}

template void SparseMatrixCSRVector2General<int>::Save (const std::string&, std::string) const;
template void SparseMatrixCSRVector2General<int>::Restore (const std::string&, std::string);
template void SparseMatrixCSRVector2General<double>::Save (const std::string&, std::string) const;
template void SparseMatrixCSRVector2General<double>::Restore (const std::string&, std::string);



#endif // ENABLE_SERIALIZATION

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SparseMatrixCSRVector2General<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SparseMatrixCSRVector2General<int>)
#endif  // ENABLE_SERIALIZATION
