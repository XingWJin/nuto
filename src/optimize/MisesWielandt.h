// $Id $
#pragma once


// parent
#include "optimize/Optimizer.h"

#include "optimize/OptimizeException.h"


namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... standard class for vonMises-Wielandt method
#ifdef ENABLE_MECHANICS
	class StructureGrid;
#endif // ENABLE_MECHANICS

class MisesWielandt : public virtual Optimizer
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    MisesWielandt(unsigned int rNumParameters) : Optimizer(rNumParameters,(unsigned int)0,(unsigned int) 0)
    {
        mAccuracyGradient = 1e-6;
        mMinDeltaObjBetweenRestarts = 1e-6;
        mMaxGradientCalls = INT_MAX,
        mMaxHessianCalls = INT_MAX,
        mMaxIterations = INT_MAX;
        mShowSteps = 100;
       	mUseDiagHessian =true;
       	mNumParameters=rNumParameters;
       	mObjectiveType=CONDITION_NUMBER_OF_MATRIX;
	}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
    	   & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
           & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
           & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
           & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
           & BOOST_SERIALIZATION_NVP(mMaxIterations)
           & BOOST_SERIALIZATION_NVP(mShowSteps)
           & BOOST_SERIALIZATION_NVP(mUseDiagHessian)
           & BOOST_SERIALIZATION_NVP(mNumParameters)
           & BOOST_SERIALIZATION_NVP(mObjectiveType);
    }
#endif // SWIG
#endif // ENABLE_SERIALIZATION

	enum eObjectiveType
	{
		MAX_EIGENVALUE_OF_PRECOND_MATRIX,
		MAX_EIGENVALUE_OF_MATRIX,
		SPECTRAL_RADIUS_OF_MATRIX,
		SPECTRAL_RADIUS_OF_PRECOND_MATRIX,
		CONDITION_NUMBER_OF_MATRIX,
		CONDITION_NUMBER_OF_PRECOND_MATRIX,
	 };

    int Optimize();


    void SetObjectiveType(std::string rObjectiveType)
    {
		std::string upperCaseObjectiveType;
		 std::transform(rObjectiveType.begin(), rObjectiveType.end(), std::back_inserter(upperCaseObjectiveType), (int(*)(int)) toupper);

		 if (upperCaseObjectiveType=="MAX_EIGENVALUE_OF_PRECOND_MATRIX")
			mObjectiveType=MAX_EIGENVALUE_OF_PRECOND_MATRIX;
		 else if(upperCaseObjectiveType=="MAX_EIGENVALUE_OF_MATRIX")
			mObjectiveType=MAX_EIGENVALUE_OF_MATRIX;
		 else if(upperCaseObjectiveType=="SPECTRAL_RADIUS_OF_MATRIX")
			 mObjectiveType=SPECTRAL_RADIUS_OF_MATRIX;
		 else if(upperCaseObjectiveType=="SPECTRAL_RADIUS_OF_PRECOND_MATRIX")
			 mObjectiveType=SPECTRAL_RADIUS_OF_PRECOND_MATRIX;
		 else if(upperCaseObjectiveType=="CONDITION_NUMBER_OF_MATRIX")
			 mObjectiveType=CONDITION_NUMBER_OF_MATRIX;
		 else if(upperCaseObjectiveType=="CONDITION_NUMBER_OF_PRECOND_MATRIX")
			 mObjectiveType=CONDITION_NUMBER_OF_PRECOND_MATRIX;
		else
		{
			throw OptimizeException("[NuTo::MisesWielandt::SetObjectiveType] ObjectiveType "+upperCaseObjectiveType +" does not exist.");
		}
    }

	inline void SetObjectiveType(eObjectiveType rObjectiveType)
    {
			 mObjectiveType=rObjectiveType;
    }
	inline eObjectiveType GetObjectiveType()
    {
    	return mObjectiveType;
    }

    inline void SetMaxGradientCalls(int rMaxGradientCalls)
    {
        mMaxGradientCalls = rMaxGradientCalls;
    }

    inline void SetMaxHessianCalls(int rMaxHessianCalls)
    {
        mMaxHessianCalls = rMaxHessianCalls;
    }

    inline void SetMaxIterations(int rMaxIterations)
    {
        mMaxIterations = rMaxIterations;
    }

    inline void SetAccuracyGradient(double rAccuracyGradient)
    {
        mAccuracyGradient = rAccuracyGradient;
    }

    inline void SetMinDeltaObjBetweenRestarts(double rMinDeltaObjBetweenRestarts)
    {
        mMinDeltaObjBetweenRestarts = rMinDeltaObjBetweenRestarts;
    }

    inline void SetShowSteps(int rShowSteps)
    {
        mShowSteps = rShowSteps;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const;


    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType);
#endif // ENABLE_SERIALIZATION

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
	virtual void Info()const;


protected:

	double mAccuracyGradient;
	double mMinDeltaObjBetweenRestarts;
	int    mMaxGradientCalls;
	int    mMaxHessianCalls;
	int    mMaxIterations;
	int    mShowSteps;
    bool   mUseDiagHessian;
    eObjectiveType 	mObjectiveType;
   	size_t mNumParameters;


};
} // namespace NuTo