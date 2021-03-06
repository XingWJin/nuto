#pragma once


#include "base/Exception.h"
#include <memory>
#include <eigen3/Eigen/Core>

namespace NuTo
{

template <int TRows, int TCols>
class ConstitutiveMatrix;
template <int TRows>
class ConstitutiveVector;
template <int TDim>
class EngineeringStrain;
namespace Constitutive
{
enum class eInput;
enum class eOutput;
std::string OutputToString(const eOutput);
}

class ConstitutiveIOBase
{

public:
    ConstitutiveIOBase() = default;
    ConstitutiveIOBase(const ConstitutiveIOBase& rOther) = default;

    virtual ~ConstitutiveIOBase() = default;

    //! Copy construct `this` into a `unique_ptr`.
    virtual std::unique_ptr<ConstitutiveIOBase> clone() = 0;

    //! Factory for polymorphic construction of constitutive outputs.
    //! @param outputType Determines which derived object is returned
    //! @return `unique_ptr` to the created object
    template <int TDim>
    static std::unique_ptr<ConstitutiveIOBase> makeConstitutiveIO(NuTo::Constitutive::eOutput outputType);
    //! Factory for polymorphic construction of constitutive inputs.
    //! @param inputType Determines which derived object is returned
    //! @return `unique_ptr` to the created object
    template <int TDim>
    static std::unique_ptr<ConstitutiveIOBase> makeConstitutiveIO(NuTo::Constitutive::eInput inputType);

    ConstitutiveIOBase& operator=(const ConstitutiveIOBase& rOther);

    //! @remark 1-->1,   2-->3,   3-->6
    static constexpr int GetVoigtDim(int D)
    {
        return D == 1 ? 1 : (D == 2 ? 3 : 6);
    }

    //! matrix access
    virtual double& operator()(int rRow, int rCol);
    virtual double operator()(int rRow, int rCol) const;

    //! vector access
    virtual double& operator[](int rRow);
    virtual double operator[](int rRow) const;

    virtual void SetZero()
    {
        throw Exception(__PRETTY_FUNCTION__, "not implemented for this constitutive input type");
    }
    virtual int GetNumRows() const
    {
        throw Exception(__PRETTY_FUNCTION__, "not implemented for this constitutive input type");
    }
    virtual int GetNumColumns() const
    {
        throw Exception(__PRETTY_FUNCTION__, "not implemented for this constitutive input type");
    }

    //! @brief copies itsself to an Eigen::MatrixXd
    //! @return matrix
    Eigen::MatrixXd CopyToEigenMatrix() const;

    /**************************************************************************
     *
     *  some "pretty asserts" to ensure type safety
     *
     ***************************************************************************/
    void AssertIsScalar(Constitutive::eOutput rOutputEnum, std::string rMethodName) const;
    // implementation in cpp file, since the dynamic_cast to ConstitutiveScalar
    // requires the full include instead of the forward declaration

    template <int TRows>
    void AssertIsVector(Constitutive::eOutput rOutputEnum, std::string rMethodName) const
    {
#ifndef NDEBUG
        AssertDimension<TRows, 1>(rOutputEnum, rMethodName);
        bool isNotVector = dynamic_cast<const ConstitutiveVector<TRows>*>(this) == nullptr;
        if (isNotVector)
            throw Exception(rMethodName, "Constitutive output " + Constitutive::OutputToString(rOutputEnum) +
                                                 " is not a ConstitutiveVector<>.");
#endif
    }

    template <int TRows, int TCols>
    void AssertIsMatrix(Constitutive::eOutput rOutputEnum, std::string rMethodName) const
    {
#ifndef NDEBUG
        AssertDimension<TRows, TCols>(rOutputEnum, rMethodName);
        bool isNotMatrix = dynamic_cast<const ConstitutiveMatrix<TRows, TCols>*>(this) == nullptr;
        if (isNotMatrix)
            throw Exception(rMethodName, "Constitutive output " + Constitutive::OutputToString(rOutputEnum) +
                                                 " is not a ConstitutiveMatrix<>.");
#endif
    }

    template <int TDim>
    EngineeringStrain<TDim>& AsEngineeringStrain();

    template <int TDim>
    const EngineeringStrain<TDim>& AsEngineeringStrain() const;

    virtual const EngineeringStrain<1>& AsEngineeringStrain1D() const
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "input/output is not engineering strain.");
    }
    virtual const EngineeringStrain<2>& AsEngineeringStrain2D() const
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "input/output is not engineering strain.");
    }
    virtual const EngineeringStrain<3>& AsEngineeringStrain3D() const
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "input/output is not engineering strain.");
    }

    virtual EngineeringStrain<1>& AsEngineeringStrain1D()
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "input/output is not engineering strain.");
    }
    virtual EngineeringStrain<2>& AsEngineeringStrain2D()
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "input/output is not engineering strain.");
    }
    virtual EngineeringStrain<3>& AsEngineeringStrain3D()
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "input/output is not engineering strain.");
    }

    void SetIsCalculated(bool rIsCalculated)
    {
        mIsCalculated = rIsCalculated;
    }

    bool GetIsCalculated() const
    {
        return mIsCalculated;
    }

private:
#ifndef NDEBUG
    template <int TRows, int TCols>
    void AssertDimension(Constitutive::eOutput rOutputEnum, const std::string& rMethodName) const
    {
        if (GetNumRows() != TRows || GetNumColumns() != TCols)
        {
            std::string exception;
            exception += "[" + rMethodName + "] \n";
            exception += "Dimension mismatch of constitutive output. \n";
            exception += "Dim(" + Constitutive::OutputToString(rOutputEnum) + ") = (";
            exception += std::to_string(GetNumRows()) + "x" + std::to_string(GetNumColumns()) + ") ";
            exception += "Expected: (" + std::to_string(TRows) + "x" + std::to_string(TCols) + ") \n";
            throw Exception(exception);
        }
    }
#endif

    //!@brief Is supposed to be set to <B>TRUE<B> by the constitutive law after
    //! calculation. Elements should check every output object before using its
    //! value.
    bool mIsCalculated = false;
};

} /* namespace NuTo */
