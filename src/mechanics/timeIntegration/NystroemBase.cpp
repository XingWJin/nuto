#include "math/SparseDirectSolverMUMPS.h"

#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "mechanics/timeIntegration/NystroemBase.h"

#include "base/Timer.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NystroemBase::NystroemBase(StructureBase* rStructure)
    : TimeIntegrationBase(rStructure)
{
    mTimeStep = 0.;
    mUseDiagonalMassMatrix = true;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NystroemBase::Info() const
{
    TimeIntegrationBase::Info();
}


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
void NuTo::NystroemBase::Solve(double rTimeDelta)
{
    NuTo::Timer timer(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());
    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    if (mStructure->HasInteractingConstraints())
        throw Exception(__PRETTY_FUNCTION__, "not implemented for constrained systems including multiple dofs.");

    if (mTimeStep == 0.)
    {
        if (this->HasCriticalTimeStep())
        {
            mTimeStep = this->CalculateCriticalTimeStep();
        }
        else
        {
            throw Exception("[NuTo::NystroemBase::Solve] time step not set for unconditionally stable algorithm.");
        }
    }
    // allocate solver (usually not needed)
    NuTo::SparseDirectSolverMUMPS mySolver;

    std::cout << "modify computation of critical time step to include the dependence on the time integration scheme."
              << std::endl;
    // calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time
    // step

    std::cout << "time step " << mTimeStep << std::endl;
    std::cout << "number of time steps " << rTimeDelta / mTimeStep << std::endl;


    StructureOutputBlockVector outOfBalance(mStructure->GetDofStatus(), true);

    CalculateStaticAndTimeDependentExternalLoad();

    // store last converged displacements, velocities and accelerations
    auto dof_dt0 = mStructure->NodeExtractDofValues(0);
    auto dof_dt1 = mStructure->NodeExtractDofValues(1);

    StructureOutputBlockMatrix hessian2(mStructure->GetDofStatus(), true);

    if (mUseDiagonalMassMatrix)
    {
        hessian2 = mStructure->BuildGlobalHessian2Lumped();
        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= mStructure->GetDimension(); // since the mass is added to nodes in every direction
        std::cout << "the total mass is " << numericMass << std::endl;

        // invert the mass matrix
        hessian2.CwiseInvert();
    }
    else
    {
        // get full mass matrix
        //          StructureOutputBlockMatrix hessian2 = mStructure->BuildGlobalHessian2();
        std::cout << "ONLY LUMPED MASS MATRIX implemented" << std::endl;
    }

    double curTime = 0.;
    auto extLoad = CalculateCurrentExternalLoad(curTime);
    auto intForce = mStructure->BuildGlobalInternalGradient();

    std::vector<StructureOutputBlockVector> dof_dt2_tmp(this->GetNumStages(), mStructure->GetDofStatus());
    std::vector<double> stageDerivativeFactor(this->GetNumStages() - 1);

    while (curTime < rTimeDelta)
    {
        // calculate for delta_t = 0
        if (mStructure->GetVerboseLevel() > 5)
            std::cout << "curTime " << curTime << " (" << curTime / rTimeDelta
                      << ") max Disp = " << dof_dt0.J[Node::eDof::DISPLACEMENTS].maxCoeff() << std::endl;

        auto dof_dt0_new = dof_dt0 + dof_dt1 * mTimeStep;
        auto dof_dt1_new = dof_dt1;
        //          std::cout << "dof_dt0_new "<< dof_dt0_new << std::endl;
        //          std::cout << "dof_dt1_new "<< dof_dt1_new << std::endl;

        double prevTime(mTime);
        double prevCurTime(curTime);
        for (int countStage = 0; countStage < this->GetNumStages(); countStage++)
        {
            double deltaTimeStage = this->GetStageTimeFactor(countStage) * mTimeStep;
            this->GetStageDerivativeFactor(stageDerivativeFactor, countStage);
            //              std::cout << "countStage "<< countStage << std::endl;
            auto dof_dt0_tmp = dof_dt0 + dof_dt1 * deltaTimeStage;
            for (int countStage2 = 0; countStage2 < countStage; countStage2++)
            {
                if (stageDerivativeFactor[countStage2] != 0.)
                {
                    dof_dt0_tmp +=
                            dof_dt2_tmp[countStage2] * (stageDerivativeFactor[countStage2] * mTimeStep * mTimeStep);
                }
            }
            //              std::cout << "dof_dt0_tmp "<< dof_dt0_tmp << std::endl;

            if (this->HasTimeChanged(countStage) == true)
            {
                curTime = prevCurTime + deltaTimeStage;
                mTime = prevTime + deltaTimeStage;

                UpdateConstraints(curTime);

                // calculate external force
                extLoad = CalculateCurrentExternalLoad(curTime);
            }

            dof_dt0_tmp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_tmp.J);
            mStructure->NodeMergeDofValues(0, dof_dt0_tmp);
            mStructure->ElementTotalUpdateTmpStaticData();

            // calculate internal force (with update of history variables = true)
            intForce = mStructure->BuildGlobalInternalGradient();
            //              std::cout << "F-R " << (extLoad-intForce) << std::endl;

            if (mUseDiagonalMassMatrix)
            {
                // no system has to be solved
                dof_dt2_tmp[countStage] = hessian2 * (extLoad - intForce);
                //                  std::cout << "dof_dt2_tmp "<< dof_dt2_tmp[countStage] <<
                // std::endl;
            }
            else
            {
                // system of equations has to be solved
                //                  mySolver.Solve(fullMassMatrix, extLoad-intForce,
                // dof_dt2_tmp[countStage]);
                std::cout << "ONLY LUMPED MASS MATRIX implemented" << std::endl;
            }

            dof_dt0_new += dof_dt2_tmp[countStage] * (GetStageWeights1(countStage) * mTimeStep * mTimeStep);
            dof_dt1_new += dof_dt2_tmp[countStage] * (GetStageWeights2(countStage) * mTimeStep);
            //              std::cout << "dof_dt0_new "<< dof_dt0_new << std::endl;
            //              std::cout << "dof_dt1_new "<< dof_dt1_new << std::endl;
        }

        mTime = prevTime + mTimeStep;
        curTime = prevCurTime + mTimeStep;

        dof_dt0_new.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_new.J);
        mStructure->NodeMergeDofValues(0, dof_dt0_new);
        if (mStructure->GetNumTimeDerivatives() >= 1)
        {
            dof_dt1_new.K = mStructure->NodeCalculateDependentDofValues(dof_dt1_new.J);
            mStructure->NodeMergeDofValues(1, dof_dt1_new);
        }
        if (mStructure->GetNumTimeDerivatives() >= 2)
        {
            auto dof_dt2_new = (dof_dt1_new - dof_dt1) * (1. / mTimeStep);
            dof_dt2_new.K = mStructure->NodeCalculateDependentDofValues(dof_dt2_new.J);
            mStructure->NodeMergeDofValues(2, dof_dt2_new);
        }

        //          mStructure->ElementTotalUpdateTmpStaticData();
        mStructure->ElementTotalUpdateStaticData();

        dof_dt0 = dof_dt0_new;
        dof_dt1 = dof_dt1_new;

        //**********************************************
        // PostProcessing
        //**********************************************
        // postprocess data for plotting
        mTimeControl.SetCurrentTime(mTime);
        mPostProcessor->PostProcess(extLoad - intForce);
    }
}
