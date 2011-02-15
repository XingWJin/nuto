// $Id: SparseMatrixCSRVector2General.h 235 2010-04-22 09:25:38Z arnold2 $

#ifndef SPARSE_MATRIX_CSR_VECTOR2_GENERAL_H
#define SPARSE_MATRIX_CSR_VECTOR2_GENERAL_H

#include <algorithm>

#include "nuto/math/SparseMatrixCSRVector2General_Def.h"

#include "nuto/math/FullMatrix.h"
//#include "nuto/math/SparseMatrixCSRGeneral_Def.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/MathException.h"

//! @brief ... constructor
//! @param rNumRows_ ... number of rows
//! @param rNumColumns_ ... number of columns
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(int rNumRows_, int rNumColumns_) :
                             NuTo::SparseMatrixCSRVector2<T>(rNumRows_,rNumColumns_)
{
}

//! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
//! @param rFullMatrix ... input matrix (full storage)
//! @param rAbsoluteTolerance ... absolute tolerance
//! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(NuTo::FullMatrix<T>& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance):
                             NuTo::SparseMatrixCSRVector2<T>(rFullMatrix.GetNumRows(),rFullMatrix.GetNumColumns())
{
	double tolerance = rAbsoluteTolerance;
	if (rRelativeTolerance > 1e-14)
	{
		T maxValue(0), minValue(0);
		this->Max(maxValue);
		this->Min(minValue);
		if (fabs(maxValue)>fabs(minValue))
			tolerance += rRelativeTolerance * fabs(maxValue);
		else
			tolerance += rRelativeTolerance * fabs(minValue);

	}

	for (int row = 0; row < rFullMatrix.GetNumRows(); row++)
	{
		for (int col = 0; col < rFullMatrix.GetNumColumns(); col++)
		{
			if (rFullMatrix(row,col) > tolerance)
			{
				this->AddEntry(rFullMatrix(row,col),row,col);
			}
		}
	}
}

//! @brief ... create sparse matrix with vector of vector from standard CSR format
//! @param rCSRMatrix ... input matrix (full storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(const NuTo::SparseMatrixCSRGeneral<T>& rCSRMatrix) :
              NuTo::SparseMatrixCSRVector2<T>::SparseMatrixCSRVector2(rCSRMatrix.GetNumRows(),rCSRMatrix.GetNumColumns() )
{
    this->mNumColumns = rCSRMatrix.GetNumColumns();
    this->mColumns.resize(rCSRMatrix.GetNumRows());
    this->mValues.resize(rCSRMatrix.GetNumRows());
    std::vector<int>::const_iterator startIteratorColumns(rCSRMatrix.mColumns.begin());
    typename std::vector<T>::const_iterator startIteratorValues(rCSRMatrix.mValues.begin());
    int numEntries;
    for (unsigned int row = 0; row < this->mColumns.size(); row++)
    {
   	    numEntries = rCSRMatrix.mRowIndex[row+1]-rCSRMatrix.mRowIndex[row];
    	this->mColumns[row].resize(numEntries);
    	this->mValues[row].resize(numEntries);
    	std::copy(startIteratorColumns, startIteratorColumns+numEntries,this->mColumns[row].begin());
    	std::copy(startIteratorValues, startIteratorValues+numEntries,this->mValues[row].begin());
    	startIteratorColumns+=numEntries;
    	startIteratorValues+=numEntries;
     }
}

//! @brief ... returns whether the matrix is symmetric or unsymmetric
//! @return true if the matrix is symmetric and false if the matrix is unsymmetric
template<class T>
bool NuTo::SparseMatrixCSRVector2General<T>::IsSymmetric() const
{
	return false;
}

//! @brief ... add nonzero entry to matrix
//! @param rRow ... row of the nonzero entry (zero based indexing!!!)
//! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
//! @param rValue ... value of the nonzero entry
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::AddEntry(int rRow, int rColumn, T rValue)
{
	// check for overflow
	assert(rRow < INT_MAX);
	assert(rColumn < INT_MAX);
	assert(rRow >= 0);
	assert(rColumn >= 0);

	// check bounds
	if (rRow > (int)this->mValues.size() || rRow<0)
	{
		throw MathException("[SparseMatrixCSRVector2General::addEntry] row index is out of bounds.");
	}
	if (rColumn >= this->mNumColumns || rColumn<0)
	{
		throw MathException("[SparseMatrixCSRVector2General::addEntry] column index is out of bounds.");
	}
	if (this->mOneBasedIndexing)
	{
		rColumn++;
	}

/*
	// find position in matrix
	bool existingValue(false);
	int col_pos;

	for (col_pos = 0; (col_pos < (int)this->mColumns[rRow].size()); col_pos++)
	{
		if (this->mColumns[rRow][col_pos] < rColumn)
			continue;
		else
		{
			if (this->mColumns[rRow][col_pos] == rColumn)
			{
				existingValue = true;
			}
			break;
		}
	}
	// add value
	if (existingValue)
	{
		// add to existing value
		this->mValues[rRow][col_pos] += rValue;
	}
	else
	{
		// insert new value
		this->mColumns[rRow].insert(this->mColumns[rRow].begin() + col_pos, rColumn);
		this->mValues[rRow].insert(this->mValues[rRow].begin() + col_pos, rValue);
	}

*/
    typename std::vector<int>::iterator it = lower_bound(this->mColumns[rRow].begin(),this->mColumns[rRow].end(),rColumn);
	if (it==this->mColumns[rRow].end())
	{
		// insert new value
		unsigned int pos = it-this->mColumns[rRow].begin();
		this->mColumns[rRow].insert(it, rColumn);
		this->mValues[rRow].insert(this->mValues[rRow].begin() + pos, rValue);
	}
	else
	{
		if (*it==rColumn)
		{
			*(this->mValues[rRow].begin() + (it-this->mColumns[rRow].begin())) += rValue;
		}
		else
		{
			// insert new value
			int pos = it-this->mColumns[rRow].begin();
			this->mColumns[rRow].insert(it, rColumn);
			this->mValues[rRow].insert(this->mValues[rRow].begin() + pos, rValue);
		}
	}
}

//! @brief ... import matrix from slang object stored in  a text file
//! @param rFileName ... file name
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ImportFromSLangText(const char* rFileName)
{
	throw MathException("NuTo::SparseMatrixCSRVector2General::ImportFromSLangText] to be implemented.");
}

//! @brief ... write nonzero matrix entries into a full matrix
//! @param rFullMatrix ... the full matrix
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::WriteEntriesToFullMatrix(FullMatrix<T>& rFullMatrix) const
{
	if (this->mOneBasedIndexing)
	{
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				rFullMatrix(row, this->mColumns[row][col_count]-1) = this->mValues[row][col_count];
			}
		}
	}
	else
	{
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				rFullMatrix(row, this->mColumns[row][col_count]) = this->mValues[row][col_count];
			}
		}
	}
}

#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param rType ... type of file, either BINARY, XML or TEXT
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Save ( const std::string &filename, std::string rType)const
{
	try
	 {
		 //transform to uppercase
		 std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

		 // open file
		 std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
		 if(! ofs.is_open())
		 {
			 throw MathException("[NuTo::SparseMatrixCSRVector2General::Save] Error opening file.");
		 }

		 // write data to file
		 std::string typeIdString(this->GetTypeId());
		 if (rType=="BINARY")
		 {
			 boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
			 oba & boost::serialization::make_nvp ("Object_type", typeIdString );
			 oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		 }
		 else if (rType=="XML")
		 {
			 boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
			 std::string tmpString(this->GetTypeId());
			 oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
			 oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		 }
		 else if (rType=="TEXT")
		 {
			 boost::archive::text_oarchive ota ( ofs, std::ios::binary );
			 ota & boost::serialization::make_nvp("Object_type", typeIdString );
			 ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
		 }
		 else
		 {
			 throw MathException ( "[NuTo::SparseMatrixCSRVector2General::Save] File type not implemented." );
		 }

		 // close file
		 ofs.close();
	 }
	 catch ( boost::archive::archive_exception e )
	 {
		 std::string s ( std::string ( "[NuTo::SparseMatrixCSRVector2General::Save]File save exception in boost - " ) + std::string ( e.what() ) );
		 throw MathException ( s );
	 }
	 catch ( MathException &e )
	 {
		 throw e;
	 }
	 catch ( std::exception &e )
	 {
		 throw MathException ( e.what() );
	 }
}

//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Restore ( const std::string &filename,  std::string rType)
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
	catch ( boost::archive::archive_exception e )
	{
		std::string s ( std::string ( "[NuTo::SparseMatrixCSRVector2General::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MathException ( e.what() );
	}
}
#endif // ENABLE_SERIALIZATION


//! @brief ... add two matrices
//! @param rOther ... general sparse matrix stored in the CSRVector2 format
//! @return general sparse matrix stored in the CSR format
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::operator+ ( const NuTo::SparseMatrixCSRVector2General<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRVector2General::operator+] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator+] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRVector2General<T> result(this->GetNumRows(), this->GetNumColumns());
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		unsigned int thisPos(0), otherPos(0);
		std::vector<int>& thisColumnVec(this->mColumns[row]);
		const std::vector<int>& otherColumnVec(rOther.mColumns[row]);
		std::vector<T>& thisValueVec(this->mValues[row]);
		const std::vector<T>& otherValueVec(rOther.mValues[row]);
		while ((thisPos < thisColumnVec.size()) || (otherPos < otherColumnVec.size()))
		{
			int thisColumn;
			if (thisPos < thisColumnVec.size())
			{
				thisColumn = thisColumnVec[thisPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				thisColumn = this->mNumColumns;
			}

			int otherColumn;
			if (otherPos < otherColumnVec.size())
			{
				otherColumn = otherColumnVec[otherPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				otherColumn = this->mNumColumns;
			}

			int resultColumn;
			T resultValue;
			if (thisColumn < otherColumn)
			{
				resultColumn = thisColumn;
				resultValue = thisValueVec[thisPos];
				thisPos++;
			}
			else if (otherColumn < thisColumn)
			{
				resultColumn = otherColumn;
				resultValue = otherValueVec[otherPos];
				otherPos++;
			}
			else
			{
				resultColumn = thisColumn;
				resultValue = thisValueVec[thisPos] + otherValueVec[otherPos];
				thisPos++;
				otherPos++;
			}
			result.mColumns[row].push_back(resultColumn);
			result.mValues[row].push_back(resultValue);
		}
		assert(result.mColumns.size() == result.mValues.size());
	}
	return result;
}

//! @brief ... subtract two matrices
//! @param rOther ... general sparse matrix stored in the CSRVector2 format
//! @return general sparse matrix stored in the CSR format
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::operator- ( const NuTo::SparseMatrixCSRVector2General<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRVector2General::operator+] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator+] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRVector2General<T> result(this->GetNumRows(), this->GetNumColumns());
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		unsigned int thisPos(0), otherPos(0);
		std::vector<int>& thisColumnVec(this->mColumns[row]);
		const std::vector<int>& otherColumnVec(rOther.mColumns[row]);
		std::vector<T>& thisValueVec(this->mValues[row]);
		const std::vector<T>& otherValueVec(rOther.mValues[row]);
		while ((thisPos < thisColumnVec.size()) || (otherPos < otherColumnVec.size()))
		{
			int thisColumn;
			if (thisPos < thisColumnVec.size())
			{
				thisColumn = thisColumnVec[thisPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				thisColumn = this->mNumColumns;
			}

			int otherColumn;
			if (otherPos < otherColumnVec.size())
			{
				otherColumn = otherColumnVec[otherPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				otherColumn = this->mNumColumns;
			}

			int resultColumn;
			T resultValue;
			if (thisColumn < otherColumn)
			{
				resultColumn = thisColumn;
				resultValue = thisValueVec[thisPos];
				thisPos++;
			}
			else if (otherColumn < thisColumn)
			{
				resultColumn = otherColumn;
				resultValue = -otherValueVec[otherPos];
				otherPos++;
			}
			else
			{
				resultColumn = thisColumn;
				resultValue = thisValueVec[thisPos] - otherValueVec[otherPos];
				thisPos++;
				otherPos++;
			}
			result.mColumns[row].push_back(resultColumn);
			result.mValues[row].push_back(resultValue);
		}
		assert(result.mColumns.size() == result.mValues.size());
	}
	return result;
}

//! @brief ... subtract two matrices
//! @param rOther ... general sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::operator-=  ( const NuTo::SparseMatrixCSRVector2General<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRVector2General::operator-=] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator-=] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (unsigned int pos = 0; pos < rOther.mValues[row].size(); pos++)
		{
			this->AddEntry(row, rOther.mColumns[row][pos], - rOther.mValues[row][pos]);
		}
	}
	return *this;
}

//! @brief ... add two matrices
//! @param rOther ... general sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::operator+=  ( const NuTo::SparseMatrixCSRVector2General<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRVector2General::operator+=] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator+=] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (unsigned int pos = 0; pos < rOther.mValues[row].size(); pos++)
		{
			this->AddEntry(row, rOther.mColumns[row][pos], rOther.mValues[row][pos]);
		}
	}
	return *this;
}

//! @brief ... matrix - matrix multiplication
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return general sparse matrix stored in the CSR format
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::operator* ( const NuTo::SparseMatrixCSRVector2General<T> &rOther ) const
{
	if (this->GetNumColumns() != rOther.GetNumRows())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator*] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator*] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRVector2General<T> result(this->GetNumRows(), rOther.GetNumColumns());

	for (int thisRow = 0; thisRow < this->GetNumRows(); thisRow++)
	{
		const std::vector<T>& thisValueVec(this->mValues[thisRow]);
		const std::vector<int>& thisColumnVec(this->mColumns[thisRow]);
		for (unsigned int thisPos = 0; thisPos < thisValueVec.size(); thisPos++)
		{
			unsigned int thisColumn =thisColumnVec[thisPos];
			T thisValue = thisValueVec[thisPos];

			const std::vector<int>& otherColumnVec(rOther.mColumns[thisColumn]);
			const std::vector<T>& otherValueVec(rOther.mValues[thisColumn]);
			for (unsigned int otherPos = 0; otherPos < otherColumnVec.size(); otherPos++)
			{
				result.AddEntry(thisRow, otherColumnVec[otherPos], thisValue * otherValueVec[otherPos]);
			}
		}
	}
	return result;
}

//! @brief ... multiplies the matrix with an scalar value
//! @param rOther ... scalar value
//! @return ... the multiplied matrix (sparse csr storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::operator* ( const T &rOther ) const
{
	SparseMatrixCSRVector2General<T> result(*this);
	for (unsigned int thisRow=0; thisRow< this->mValues.size(); thisRow++)
	{
		std::vector<T>& resultValueVec(result.mValues[thisRow]);
		for (unsigned int thisColumn=0; thisColumn< this->mValues[thisRow].size(); thisColumn++)
		{
			resultValueVec[thisColumn]*=rOther;
		}
	}
	return result;
}

//! @brief ... multiply sparse matrix with a full matrix
//! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
//! @return ... full matrix
template<class T>
NuTo::FullMatrix<T> NuTo::SparseMatrixCSRVector2General<T>::operator* (const FullMatrix<T> &rMatrix) const
{
	if (this->GetNumColumns() != rMatrix.GetNumRows())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator*] invalid matrix dimensions.");
	}
	FullMatrix<T> result(this->GetNumRows(),rMatrix.GetNumColumns());
	if (this->HasOneBasedIndexing())
	{
		// loop over rows
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			// initialize result
			for (int matrixCol = 0; matrixCol < result.GetNumColumns(); matrixCol++)
			{
				result(row,matrixCol) = 0.0;
			}
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			// perform multiplication
			for (unsigned int pos = 0; pos < this->mColumns.size(); pos++)
			{
				int column = thisColumnVec[pos] - 1;
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	else
	{
		// loop over rows
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			// initialize result
			for (int matrixCol = 0; matrixCol < result.GetNumColumns(); matrixCol++)
			{
				result(row,matrixCol) = 0.0;
			}
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			// perform multiplication
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				int column = thisColumnVec[pos];
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	return result;
}

template<class T>
NuTo::FullMatrix<T> NuTo::SparseMatrixCSRVector2General<T>::TransMult(const NuTo::FullMatrix<T>& rMatrix) const
{
	if (this->GetNumRows() != rMatrix.GetNumRows())
	{
		throw MathException("[SparseMatrixCSRVector2General::TransMult] invalid matrix dimensions.");
	}
	FullMatrix<T> result(this->GetNumColumns(),rMatrix.GetNumColumns());
	if (this->HasOneBasedIndexing())
	{
		// loop over columns of transpose
		for (int column = 0; column < this->GetNumRows(); column++)
		{
			// perform multiplication
			const std::vector<T>& thisValueVec(this->mValues[column]);
			const std::vector<int>& thisColumnVec(this->mColumns[column]);
			for (unsigned int pos = 0; pos < thisValueVec.size(); pos++)
			{
				int row = thisColumnVec[pos] - 1;
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	else
	{
		// loop over columns of transpose
		for (int column = 0; column < this->GetNumRows(); column++)
		{
			// perform multiplication
			const std::vector<T>& thisValueVec(this->mValues[column]);
			const std::vector<int>& thisColumnVec(this->mColumns[column]);
			for (unsigned int pos = 0; pos < thisValueVec.size(); pos++)
			{
				int row = thisColumnVec[pos];
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	return result;
}

//! @brief ... calculate the transpose of the matrix (transpose row and columns)
//! @return ... transpose of this matrix (sparse csr storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::Transpose() const
{
	SparseMatrixCSRVector2General<T> transMatrix(this->GetNumColumns(), this->GetNumRows());
	if (this->mOneBasedIndexing)
	{
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				transMatrix.AddEntry(thisColumnVec[pos]-1, row, thisValueVec[pos]);
			}
		}
		transMatrix.SetOneBasedIndexing();
	}
	else
	{
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				transMatrix.AddEntry(thisColumnVec[pos], row, thisValueVec[pos]);
			}
		}
	}
	return transMatrix;
}

template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Gauss(NuTo::FullMatrix<T>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    throw MathException("SparseMatrixCSRVector2General::Gauss] : to be implemented");
}

//! @brief ... reorder columns of the matrix
//! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering)
{
	throw MathException("[SparseMatrixCSRVector2General::ReorderColumns] To be implemented.");
	/* copied from SparseMatrixCSR
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
		{
			this->mColumns[pos] = rMappingInitialToNewOrdering[this->mColumns[pos]];
		}

		// sort columns (simple bubble sort algorithm)
		int start = this->mRowIndex[row];
		int end = this->mRowIndex[row + 1] - 1;
		bool swapFlag;
		do
		{
			swapFlag = false;
			for (int pos = start; pos < end; pos++)
			{
				if (this->mColumns[pos] > this->mColumns[pos + 1])
				{
					swapFlag=true;
					int tmpInt = this->mColumns[pos];
					this->mColumns[pos] = this->mColumns[pos + 1];
					this->mColumns[pos + 1] = tmpInt;

					T tmpDouble = this->mValues[pos];
					this->mValues[pos] = this->mValues[pos + 1];
					this->mValues[pos + 1] = tmpDouble;
				}
			}
			end--;
		}
		while (swapFlag);
	}
*/
}

//! @brief ... Concatenate columns from another matrix
//! @param rOther ... other matrix with same number of rows
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ConcatenateColumns(const SparseMatrixCSRVector2General<T>& rOther)
{
/*    //! @brief value of nonzero matrix entries
    std::vector<std::vector<T> > mValues;
    //! @brief columns of nonzero matrix entries
    std::vector<std::vector<int> > mColumns;
    //! @brief ... number of columns
    int mNumColumns;
*/
    if (this->mOneBasedIndexing!=rOther.mOneBasedIndexing)
        throw MathException("[NuTo::SparseMatrixCSRVector2General<T>::ConcatenateColumns] index base (0 or 1) should be identical for both matrices");
    if (this->mValues.size()!=rOther.mValues.size())
        throw MathException("[NuTo::SparseMatrixCSRVector2General<T>::ConcatenateColumns] number of rows has to be identical for both matrices.");

    for (unsigned int theRow=0; theRow<this->mValues.size(); theRow++)
    {
        this->mValues[theRow].insert(this->mValues[theRow].end(), rOther.mValues[theRow].begin(), rOther.mValues[theRow].end());
        unsigned int oldSize(this->mColumns[theRow].size());
        this->mColumns[theRow].insert(this->mColumns[theRow].end(), rOther.mColumns[theRow].begin(), rOther.mColumns[theRow].end());
        //add to all columns of the newly added entries the columns of the original matrix
        for (unsigned int theCol=oldSize; theCol<this->mColumns[theRow].size(); theCol++)
        {
            this->mColumns[theRow][theCol]+=this->mNumColumns;
        }
    }
    this->mNumColumns+=rOther.mNumColumns;
}

//! @brief ... Concatenate rows from another matrix
//! @param rOther ... other matrix with same number of columns
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ConcatenateRows(const SparseMatrixCSRVector2General<T>& rOther)
{
    if (this->mOneBasedIndexing!=rOther.mOneBasedIndexing)
        throw MathException("[NuTo::SparseMatrixCSRVector2General<T>::ConcatenateRows] index base (0 or 1) should be identical for both matrices");
    if (this->mNumColumns!=rOther.mNumColumns)
        throw MathException("[NuTo::SparseMatrixCSRVector2General<T>::ConcatenateColumns] number of columns has to be identical for both matrices.");

    this->mValues.insert(this->mValues.end(), rOther.mValues.begin(), rOther.mValues.end());
    this->mColumns.insert(this->mColumns.end(), rOther.mColumns.begin(), rOther.mColumns.end());
}

#endif // SPARSE_MATRIX_CSR_VECTOR2_GENERAL_H
