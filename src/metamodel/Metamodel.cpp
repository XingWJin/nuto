// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include <time.h>
#include "metamodel/MinMaxTransformation.h"
#include "metamodel/ZeroMeanUnitVarianceTransformation.h"
#include "math/FullMatrix.h"
#include "metamodel/Metamodel.h"
// constructor
NuTo::Metamodel::Metamodel() : NuTo::NuToObject(), mRandomNumberGenerator(time (NULL)) // init random number generator with current time
{
}

void NuTo::Metamodel::AppendMinMaxTransformationInput(int rCoordinate, double rMin, double rMax)
{
    if (rCoordinate>=mSupportPoints.GetDimInput())
	{
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationInput - coordinate is out of range (larger than dimInput).");
	    
	}
	MinMaxTransformation *newTransformation = new MinMaxTransformation(rCoordinate, rMin, rMax);
	mSupportPoints.AppendTransformationInput(newTransformation);
}

void NuTo::Metamodel::AppendMinMaxTransformationInput(double rMin, double rMax)
{
    for (int count=0; count<mSupportPoints.GetDimInput(); count++)
		AppendMinMaxTransformationInput(count, rMin, rMax);
}

void NuTo::Metamodel::AppendMinMaxTransformationOutput(int rCoordinate, double rMin, double rMax)
{
    if (rCoordinate>=mSupportPoints.GetDimOutput())
	{
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimOutput).");
	    
	}
	MinMaxTransformation *newTransformation = new MinMaxTransformation(rCoordinate, rMin, rMax);
	mSupportPoints.AppendTransformationOutput(newTransformation);
}

void NuTo::Metamodel::AppendMinMaxTransformationOutput(double rMin, double rMax)
{
    for (int count=0; count<mSupportPoints.GetDimOutput(); count++)
		AppendMinMaxTransformationOutput(count, rMin, rMax);
}

// add zero mean, unit variance transformation to inputs
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationInput()
{
    for (int count = 0; count < this->mSupportPoints.GetDimInput(); count++)
    {
		this->AppendZeroMeanUnitVarianceTransformationInput(count);
    }
}
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationInput(int rCoordinate)
{
    if( (rCoordinate < 0) || (rCoordinate >= this->mSupportPoints.GetDimInput()) )
    {
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimInput).");
    }
	ZeroMeanUnitVarianceTransformation *newTransformation = new ZeroMeanUnitVarianceTransformation(rCoordinate);
	mSupportPoints.AppendTransformationInput(newTransformation);
}

// add zero mean, unit variance transformation to outputs
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationOutput()
{
    for (int count = 0; count < this->mSupportPoints.GetDimOutput(); count++)
    {
		this->AppendZeroMeanUnitVarianceTransformationOutput(count);
    }
}
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationOutput(int rCoordinate)
{
    if( (rCoordinate < 0) || (rCoordinate >= this->mSupportPoints.GetDimOutput()) )
    {
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimOutput).");
    }
	ZeroMeanUnitVarianceTransformation *newTransformation = new ZeroMeanUnitVarianceTransformation(rCoordinate);
	mSupportPoints.AppendTransformationOutput(newTransformation);
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::Metamodel::GetOriginalSupportPointsInput()const
{
    return mSupportPoints.GetOrigSupportPointsInput();
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::Metamodel::GetOriginalSupportPointsOutput()const
{
    return mSupportPoints.GetOrigSupportPointsOutput();
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::Metamodel::GetTransformedSupportPointsInput()const
{
	if (!mSupportPoints.IsTransformationBuild())
	    throw MetamodelException("Metamodel::GetTransformedSupportPoints - build the transformation first.");
    
	return mSupportPoints.GetTransformedSupportPointsInput();
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::Metamodel::GetTransformedSupportPointsOutput()const
{
	if (!mSupportPoints.IsTransformationBuild())
	    throw MetamodelException("Metamodel::GetTransformedSupportPoints - build the transformation first.");
    
	return mSupportPoints.GetTransformedSupportPointsOutput();
}

void NuTo::Metamodel::SetSupportPoints(int rDimInput, int rDimOutput, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rInputCoordinates, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rOutputCoordinates)
{
    if (rDimInput!=rInputCoordinates.GetNumRows())
	    throw MetamodelException("Metamodel::SetSupportPoints - dimension of input  must be equal to number of rows in the input matrix.");

    if (rDimOutput!=rOutputCoordinates.GetNumRows())
	    throw MetamodelException("Metamodel::SetSupportPoints - dimension of output  must be equal to number of rows in the output matrix.");

    if (rOutputCoordinates.GetNumColumns()!=rInputCoordinates.GetNumColumns())
	    throw MetamodelException("Metamodel::SetSupportPoints - number of samples (number of columns) in input and output matrix must be identical .");
		
	mSupportPoints.SetSupportPoints(rInputCoordinates,rOutputCoordinates);
}

void NuTo::Metamodel::BuildTransformation()
{
	if (!mSupportPoints.IsTransformationBuild())
		mSupportPoints.BuildTransformation();
}

void NuTo::Metamodel::InitRandomNumberGenerator(int rSeed)
{
    mRandomNumberGenerator.seed(rSeed);
}

double NuTo::Metamodel::RandomDouble()
{
    return std::uniform_real_distribution<double>(0.0, 1.0)(mRandomNumberGenerator); // interval [0, 1)  ( equals dSFMT close_open)
//    or
//    std::uniform_real_distribution<double> distr(0.0, 1.0); // interval [0, 1)  ( equals dSFMT close_open)
//    return distr(mRandomNumberGenerator);
}

void NuTo::Metamodel::Solve(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates)const
{
    //apply transformation of inputs
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rInputCoordinatesTransformed = rInputCoordinates;
    mSupportPoints.TransformForwardInput(rInputCoordinatesTransformed);

    //solve the submodule
    SolveTransformed(rInputCoordinatesTransformed, rOutputCoordinates);
    
    //apply transformation of outputs
    mSupportPoints.TransformForwardOutput(rOutputCoordinates);
}

void NuTo::Metamodel::SolveConfidenceInterval(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates,
                             NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinatesMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinatesMax)const
{
    //apply transformation of inputs
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rInputCoordinatesTransformed = rInputCoordinates;
    mSupportPoints.TransformForwardInput(rInputCoordinatesTransformed);
    
    //solve the submodule
    SolveConfidenceIntervalTransformed(rInputCoordinatesTransformed, rOutputCoordinates, rOutputCoordinatesMin, rOutputCoordinatesMax);
    
    //apply transformation of outputs
    mSupportPoints.TransformForwardOutput(rOutputCoordinates);
    mSupportPoints.TransformForwardOutput(rOutputCoordinatesMin);
    mSupportPoints.TransformForwardOutput(rOutputCoordinatesMax);
    
}
void NuTo::Metamodel::Build()
{
    if (mSupportPoints.GetDimOutput()<1)
        throw MetamodelException("NuTo::Metamodel::Build - number of outputs must be positive - set training data first.");
    
    if (mSupportPoints.GetDimInput()<1)
        throw MetamodelException("NuTo::Metamodel::Build - number of inputs must be positive - set training data first.");

    if (!mSupportPoints.IsTransformationBuild())
        mSupportPoints.BuildTransformation();
        
    BuildDerived();
}

void NuTo::Metamodel::Info()const
{
	mSupportPoints.Info();
}

// calculate the sample mean for each support point using original support point coordinates
void NuTo::Metamodel::GetOriginalSupportPointsMeanValue(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean) const
{
	try
	{
		this->mSupportPoints.GetMeanValueOriginalInput(rInputMean);
		this->mSupportPoints.GetMeanValueOriginalOutput(rOutputMean);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetOriginalSupportPointsMeanValue] error calculating sample means.");
		throw myException;
	}
}

// calculate the sample mean for each support point using transformed support point coordinates
void NuTo::Metamodel::GetTransformedSupportPointsMeanValue(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean) const
{
	try
	{
		this->mSupportPoints.GetMeanValueTransformedInput(rInputMean);
		this->mSupportPoints.GetMeanValueTransformedOutput(rOutputMean);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetTransformedSupportPointsMeanValue] error calculating sample means.");
		throw myException;
	}
}

// calculate the sample variance for each support point using original support point coordinates
void NuTo::Metamodel::GetOriginalSupportPointsVariance(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputVariance, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputVariance) const
{
	try
	{
		this->mSupportPoints.GetVarianceOriginalInput(rInputVariance);
		this->mSupportPoints.GetVarianceOriginalOutput(rOutputVariance);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetOriginalSupportPointsVariance] error calculating sample variance.");
		throw myException;
	}
}

// calculate the sample variance for each support point using transformed support point coordinates
void NuTo::Metamodel::GetTransformedSupportPointsVariance(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputVariance, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputVariance) const
{
	try
	{
		this->mSupportPoints.GetVarianceTransformedInput(rInputVariance);
		this->mSupportPoints.GetVarianceTransformedOutput(rOutputVariance);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetTransformedSupportPointsVariance] error calculating sample variance.");
		throw myException;
	}
}

// calculate the covariance matrix of the support points using original support point coordinates
void NuTo::Metamodel::GetOriginalSupportPointsCovarianceMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const
{
	try
	{
		this->mSupportPoints.GetCovarianceMatrixOriginal(rCovarianceMatrix);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetOriginalSupportPointsCovarianceMatrix] error calculating covariance matrix.");
		throw myException;
	}
}

// calculate the covariance matrix of the support points using transformed support point coordinates
void NuTo::Metamodel::GetTransformedSupportPointsCovarianceMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const
{
	try
	{
		this->mSupportPoints.GetCovarianceMatrixTransformed(rCovarianceMatrix);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetTransformedSupportPointsCovarianceMatrix] error calculating covariance matrix.");
		throw myException;
	}
}

// calculate Pearson's correlation matrix of the support points using original support point coordinates
void NuTo::Metamodel::GetOriginalSupportPointsPearsonCorrelationMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const
{
	try
	{
		this->mSupportPoints.GetPearsonCorrelationMatrixOriginal(rCorrelationMatrix);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetOriginalSupportPointsPearsonCorrelationMatrix] error calculating Pearson's correlation matrix.");
		throw myException;
	}
}

// calculate Pearson's correlation matrix of the support points using transformed support point coordinates
void NuTo::Metamodel::GetTransformedSupportPointsPearsonCorrelationMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const
{
	try
	{
		this->mSupportPoints.GetPearsonCorrelationMatrixTransformed(rCorrelationMatrix);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetTransformedSupportPointsPearsonCorrelationMatrix] error calculating Pearson's correlation matrix.");
		throw myException;
	}
}

// calculate confidence interval on the coefficients of Pearson's correlation matrix using original support point coordinates
void NuTo::Metamodel::GetOriginalSupportPointsPearsonCorrelationMatrixConfidenceInterval(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha) const
{
	try
	{
		this->mSupportPoints.GetPearsonCorrelationMatrixConfidenceIntervalsOriginal(rCorrelationMatrix, rMinCorrelationMatrix, rMaxCorrelationMatrix, rAlpha);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetOriginalSupportPointsPearsonCorrelationMatrixConfidenceInterval] error calculating confidence interval on the coefficients of Pearson's correlation matrix.");
		throw myException;
	}
}

// calculate confidence interval on the coefficients of Pearson's correlation matrix using transformed support point coordinates
void NuTo::Metamodel::GetTransformedSupportPointsPearsonCorrelationMatrixConfidenceInterval(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha) const
{
	try
	{
		this->mSupportPoints.GetPearsonCorrelationMatrixConfidenceIntervalsTransformed(rCorrelationMatrix, rMinCorrelationMatrix, rMaxCorrelationMatrix, rAlpha);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::Metamodel::GetTransformedSupportPointsPearsonCorrelationMatrixConfidenceInterval] error calculating confidence interval on the coefficients of Pearson's correlation matrix.");
		throw myException;
	}
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Metamodel::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Metamodel::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Metamodel" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mSupportPoints);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Metamodel" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Metamodel)
#endif  // ENABLE_SERIALIZATION