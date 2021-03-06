#pragma once

#include <set>
#include "mechanics/dofSubMatrixStorage/DofStatus.h"

namespace NuTo
{

namespace Node
{
enum class eDof : unsigned char;
} // namespace Node

//! @author Thomas Titscher, BAM
//! @date January 2016
//! @brief ... standard abstract class for all block storage classes
//! @remark ... all child classes are expected to hold the data of their subtypes, not references or pointers.
class BlockStorageBase
{

public:
    //! @brief ctor
    //! @param rDofStatus ... reference to DofStatus for automatic matrix resizing
    BlockStorageBase(const DofStatus& rDofStatus)
        : mDofStatus(rDofStatus)
    {
    }


    virtual ~BlockStorageBase();

    //! @brief gets the number of columns of the block storage
    //! @return number of columns
    int GetNumColumns() const;

    //! @brief gets the number of rows of the block storage
    //! @return number of rows
    int GetNumRows() const;

    //! @brief gets the number of columns of the block storage, only for active dof types
    //! @return number of columns
    int GetNumActiveColumns() const;

    //! @brief gets the number of rows of the block storage, only for active dof types
    //! @return number of rows
    int GetNumActiveRows() const;

#ifndef SWIG

    //! @brief gets the number of columns of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of columns
    virtual int GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const = 0;

    //! @brief gets the number of rows of the block storage for a specific set of dofs
    //! @param rDofTypes ... set of dofs
    //! @return number of rows
    virtual int GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const = 0;

#endif // SWIG

    //! @brief info method
    virtual void Info() const = 0;

    const NuTo::DofStatus& GetDofStatus() const
    {
        return mDofStatus;
    }


protected:
    const DofStatus& mDofStatus;
};

} /* namespace NuTo */
