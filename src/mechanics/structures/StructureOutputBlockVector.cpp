

#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"

NuTo::StructureOutputBlockVector::StructureOutputBlockVector(const DofStatus& rDofStatus, bool rAutomaticResize)
    : StructureOutputBase()
    , J(rDofStatus)
    , K(rDofStatus)
{
    if (rAutomaticResize)
        Resize(rDofStatus.GetNumActiveDofsMap(), rDofStatus.GetNumDependentDofsMap());
}

NuTo::StructureOutputBlockVector::~StructureOutputBlockVector()
{
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::
operator=(const NuTo::StructureOutputBlockVector& rOther)
{
    J = rOther.J;
    K = rOther.K;
    return *this;
}


void NuTo::StructureOutputBlockVector::AddElementVector(const NuTo::BlockFullVector<double>& rElementVector,
                                                        const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers)
{
    const auto& activeDofTypes = J.GetDofStatus().GetActiveDofTypes();
    const auto& numActiveDofTypeMap = J.GetDofStatus().GetNumActiveDofsMap();
    for (auto dofRow : activeDofTypes)
    {
        const auto& elementVector = rElementVector[dofRow];
        const auto& globalRowDofs = rGlobalRowDofNumbers[dofRow];

        int numActiveDofsRow = numActiveDofTypeMap.at(dofRow);

        assert(elementVector.rows() == globalRowDofs.rows());

        auto& activeRow = J[dofRow];
        auto& dependentRow = K[dofRow];

        for (int iRow = 0; iRow < globalRowDofs.rows(); ++iRow)
        {
            int globalRowDof = globalRowDofs[iRow];
            if (globalRowDof < numActiveDofsRow)
            {
                activeRow[globalRowDof] += elementVector[iRow];
            }
            else
            {
                dependentRow[globalRowDof - numActiveDofsRow] += elementVector[iRow];
            }
        }
    }
}

void NuTo::StructureOutputBlockVector::Resize(const std::map<Node::eDof, int>& rNumActiveDofsMap,
                                              const std::map<Node::eDof, int>& rNumDependentDofsMap)
{
    assert(rNumActiveDofsMap.size() == rNumDependentDofsMap.size());

    J.Resize(rNumActiveDofsMap);
    K.Resize(rNumDependentDofsMap);
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator+=(const StructureOutputBlockVector& rRhs)
{
    J += rRhs.J;
    K += rRhs.K;
    return *this;
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator-=(const StructureOutputBlockVector& rRhs)
{
    J -= rRhs.J;
    K -= rRhs.K;
    return *this;
}
NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator*=(double rRhs)
{
    J *= rRhs;
    K *= rRhs;
    return *this;
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator/=(double rRhs)
{
    J /= rRhs;
    K /= rRhs;
    return *this;
}


namespace NuTo
{
std::ostream& operator<<(std::ostream& rOut, const NuTo::StructureOutputBlockVector& rStructureOutputBlockVector)
{
    rOut << "Active Dofs" << std::endl;
    rOut << rStructureOutputBlockVector.J << std::endl;
    rOut << "Dependent Dofs" << std::endl;
    rOut << rStructureOutputBlockVector.K << std::endl;
    return rOut;
}
} // namespace NuTo

void NuTo::StructureOutputBlockVector::SetZero()
{
    J.SetZero();
    K.SetZero();
}
