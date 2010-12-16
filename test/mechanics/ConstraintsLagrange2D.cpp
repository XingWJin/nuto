#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define MAXNUMNEWTONITERATIONS 20
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7

int main()
{
try
{
    //create structure
    NuTo::Structure myStructure(2);

    double lx(1.),ly(1.);

    //2x2 nodes 1x1 element grid
    //create nodes
    NuTo::FullMatrix<double> Coordinates(2,1);
    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 0.0;
    int node1 = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = lx;
    Coordinates(1,0) = 0.0;
    int node2 = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = ly;
    int node3 = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = lx;
    Coordinates(1,0) = ly;
    int node4 = myStructure.NodeCreate("displacements",Coordinates);

    //create elements
        NuTo::FullMatrix<int> Incidence(4,1);
    Incidence(0,0) = node1;
    Incidence(1,0) = node2;
    Incidence(2,0) = node4;
    Incidence(3,0) = node3;
    myStructure.ElementCreate("PLANE2D4N",Incidence);

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,1);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    double thickness(1.);
    myStructure.SectionSetThickness(mySection,thickness);

    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);

    //Create groups to apply the boundary conditions
    //left boundary
    int GrpNodesLeftBoundary = myStructure.GroupCreate("Nodes");
    int direction = 0; //either 0,1,2
    double min(0.);
    double max(0.);
    myStructure.GroupAddNodeCoordinateRange(GrpNodesLeftBoundary,direction,min,max);

    //right boundary
    int GrpNodesRightBoundary = myStructure.GroupCreate("Nodes");
    direction = 0; //either 0,1,2
    min=1.;
    max=1.;
    myStructure.GroupAddNodeCoordinateRange(GrpNodesRightBoundary,direction,min,max);

    //bottom boundary
    int GrpNodesBottomBoundary = myStructure.GroupCreate("Nodes");
    direction=1;
    min=0;
    max=0;
    myStructure.GroupAddNodeCoordinateRange(GrpNodesBottomBoundary,direction,min,max);

    //left lower node
    int GrpNodesBottomLeftNodeBoundary = myStructure.GroupIntersection(GrpNodesBottomBoundary,GrpNodesLeftBoundary);

    //fix bottom left corner node
    NuTo::FullMatrix<double> DirectionX(2,1);
    DirectionX.SetValue(0,0,1.0);
    DirectionX.SetValue(1,0,0.0);

    NuTo::FullMatrix<double>DirectionY(2,1);
    DirectionY.SetValue(0,0,0.0);
    DirectionY.SetValue(1,0,1.0);

    int constraintLHS = myStructure.ConstraintLagrangeSetDisplacementNodeGroup(GrpNodesLeftBoundary ,DirectionX, std::string("GREATER"),-0.5);
    int constraintRHS = myStructure.ConstraintLagrangeSetDisplacementNodeGroup(GrpNodesRightBoundary,DirectionX, std::string("EQUAL"),0.0);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(GrpNodesBottomLeftNodeBoundary,DirectionY, 0);

#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentEngineeringPlasticStrain();
    myStructure.AddVisualizationComponentPrincipalEngineeringStress();
    myStructure.ElementTotalUpdateTmpStaticData();
    myStructure.ExportVtkDataFile("ConstraintsLagrange2D.vtk");
#endif

    // init some result data
    NuTo::FullMatrix<double> PlotData(1,7);

    // start analysis
    double maxDisp(1);
    double deltaDispFactor(0.2);
    double maxDeltaDispFactor(0.2);
    double curDispFactor(0.2);

    //update conre mat
    myStructure.NodeBuildGlobalDofs();

    //update tmpstatic data with zero displacements
    myStructure.ElementTotalUpdateTmpStaticData();

    //init some auxiliary variables
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::FullMatrix<double> intForceVector;
    NuTo::FullMatrix<double> extForceVector;
    NuTo::FullMatrix<double> rhsVector;

    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
    if (PRINTRESULT)
        mySolver.SetShowTime(true);

    //calculate stiffness
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    //check stiffness
/*    {
        NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2_1, stiffnessMatrixCSRVector2_2;
        NuTo::FullMatrix<double> dispForceVector_1, intForceVector_1, intForceVector_2;
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2_1, dispForceVector_1);
        NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
        NuTo::FullMatrix<double> stiffnessMatrixCDFull(stiffnessMatrixCSRVector2);
        std::cout << "analytic solution " << std::endl;
        stiffnessMatrixFull.Info(10,13);
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector_1);
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        myStructure.NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        double delta=1e-5;
        for (int theDOF=0; theDOF<stiffnessMatrixFull.GetNumRows(); theDOF++)
        {
            displacementsActiveDOFsCheck(theDOF,0)=displacementsActiveDOFsCheck(theDOF,0)+delta;
            myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            myStructure.ElementTotalUpdateTmpStaticData();
            myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector_2);
            stiffnessMatrixCDFull.SetColumn(theDOF,(intForceVector_2-intForceVector_1)*(1/delta));
            displacementsActiveDOFsCheck(theDOF,0)=displacementsActiveDOFsCheck(theDOF,0)-delta;
        }
        std::cout << "central difference solution " << std::endl;
        stiffnessMatrixCDFull.Info(10,13);
        myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        myStructure.ElementTotalUpdateTmpStaticData();
    }
    exit(0);
*/
    //apply initial displacement
    double curDisp(maxDisp*curDispFactor);
    myStructure.ConstraintSetRHS(constraintRHS,-curDisp);

    //update conre mat
    myStructure.NodeBuildGlobalDofs();

    //update displacements of all nodes according to the new conre mat
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        myStructure.NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        myStructure.ElementTotalUpdateTmpStaticData();
    }

    //build global external load vector and RHS vector
    myStructure.BuildGlobalExternalLoadVector(extForceVector);
    //extForceVector.Info(10,13);
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
    //intForceVector.Info(10,13);
    rhsVector = extForceVector + dispForceVector - intForceVector;

    //calculate absolute tolerance for matrix entries to be not considered as zero
    double maxValue, minValue, ToleranceZeroStiffness;
    stiffnessMatrixCSRVector2.Max(maxValue);
    stiffnessMatrixCSRVector2.Min(minValue);
    //std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    myStructure.SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
    int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
    if (PRINTRESULT)
        std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;


    //repeat until max displacement is reached
    bool convergenceStatusLoadSteps(false);
    int loadstep(1);
    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    myStructure.NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
    while (!convergenceStatusLoadSteps)
    {

        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
        double normRHS(1.);
        double alpha(1.);
        int convergenceStatus(0);
        //0 - not converged, continue Newton iteration
        //1 - converged
        //2 - stop iteration, decrease load step
        while(convergenceStatus==0)
        {
            numNewtonIterations++;

            if (numNewtonIterations>MAXNUMNEWTONITERATIONS && alpha<0.5)
            {
                if (PRINTRESULT)
                {
                    std::cout << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << MAXNUMNEWTONITERATIONS << ")" << std::endl;
                }
                convergenceStatus = 2; //decrease load step
                break;
            }

            normRHS = rhsVector.Norm();

            // solve
            NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsDependentDOFs;
            if (PRINTRESULT)
            {
                NuTo::FullMatrix<double> stiffnessMatrixCSRVector2Full(stiffnessMatrixCSRVector2);
                stiffnessMatrixCSRVector2Full.Info(20,10);
            }

            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSR.SetOneBasedIndexing();
            mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);
            //deltaDisplacementsActiveDOFs.Trans().Info();

            // write displacements to node
            myStructure.NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

            //perform a linesearch
            alpha = 1.;
            do
            {
                //add new displacement state
                displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
                myStructure.NodeMergeActiveDofValues(displacementsActiveDOFs);
                myStructure.ElementTotalUpdateTmpStaticData();

                // calculate residual
                myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
                rhsVector = extForceVector - intForceVector;
                normResidual = rhsVector.Norm();
                if (PRINTRESULT)
                    std::cout << "alpha " << alpha << ", normResidual " << normResidual << ", normResidualInit "<< normRHS << ", normRHS*(1-0.5*alpha) " << normRHS*(1-0.5*alpha) << std::endl;
                alpha*=0.5;
            }
            while(alpha>1e-3 && normResidual>normRHS*(1-0.5*alpha) && normResidual>1e-5);
            if (normResidual>normRHS*(1-0.5*alpha) && normResidual>1e-5)
            {
                convergenceStatus=2;
                break;
            }

            maxResidual = rhsVector.Abs().Max();

            if (PRINTRESULT)
            {
                std::cout << std::endl << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<std::endl;
                //char cDummy[100]="";
                //std::cin.getline(cDummy, 100);
            }

            //check convergence
            if (normResidual<1e-5 || maxResidual<1e-5)
            {
                if (PRINTRESULT)
                {
                    std::cout << "Convergence after " << numNewtonIterations << " Newton iterations, curdispFactor " << curDispFactor << ", deltaDispFactor "<< deltaDispFactor << std::endl<< std::endl;
                }
                convergenceStatus=1;
                break;
            }

            //convergence status == 0 (continue Newton iteration)
            //build new stiffness matrix
            myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            if (PRINTRESULT)
            {
                int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
                std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;
            }
/*            //check stiffness
            {
                NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2_1, stiffnessMatrixCSRVector2_2;
                NuTo::FullMatrix<double> dispForceVector_1, intForceVector_1, intForceVector_2;
                myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2_1, dispForceVector_1);
                NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                NuTo::FullMatrix<double> stiffnessMatrixCDFull(stiffnessMatrixCSRVector2);
                stiffnessMatrixFull.Info(10,13);
                myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector_1);
                NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
                NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
                myStructure.NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
                double delta=1e-5;
                for (int theDOF=0; theDOF<stiffnessMatrixFull.GetNumRows(); theDOF++)
                {
                    displacementsActiveDOFsCheck(theDOF,0)=displacementsActiveDOFsCheck(theDOF,0)+delta;
                    myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                    myStructure.ElementTotalUpdateTmpStaticData();
                    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector_2);
                    stiffnessMatrixCDFull.SetColumn(theDOF,(intForceVector_2-intForceVector_1)*(1/delta));
                    displacementsActiveDOFsCheck(theDOF,0)=displacementsActiveDOFsCheck(theDOF,0)-delta;
                }
                stiffnessMatrixCDFull.Info(10,13);
                myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                myStructure.ElementTotalUpdateTmpStaticData();
            }
*/
        }

        if (deltaDispFactor<1e-7)
            throw NuTo::MechanicsException("Example ConcurrentMultiscale : No convergence, delta strain factor < 1e-7");

        if (convergenceStatus==1)
        {
            myStructure.ElementTotalUpdateStaticData();

            // visualize results
#ifdef ENABLE_VISUALIZE
            std::stringstream ss;
            ss << "ConstraintsLagrange2D" << loadstep << ".vtk";
            myStructure.ExportVtkDataFile(ss.str());
#endif
             //store result/plot data
            NuTo::FullMatrix<double> SinglePlotData(1,7);

            //displacements
            NuTo::FullMatrix<double> dispNode;
            myStructure.NodeGetDisplacements(node1,dispNode);
            SinglePlotData(0,0) = dispNode(0,0);
            myStructure.NodeGetDisplacements(node2,dispNode);
            SinglePlotData(0,1) = dispNode(0,0);

            //boundary force
            NuTo::FullMatrix<double> SupportingForce;
            myStructure.NodeGroupInternalForce(GrpNodesLeftBoundary,SupportingForce);
            SinglePlotData(0,2) = SupportingForce(0,0);
            myStructure.NodeGroupInternalForce(GrpNodesRightBoundary,SupportingForce);
            SinglePlotData(0,3) = SupportingForce(0,0);

            //lagrange multiplier
            NuTo::FullMatrix<double> lagrangeMultiplier;
            myStructure.ConstraintLagrangeGetMultiplier(constraintLHS,lagrangeMultiplier);
            SinglePlotData(0,4) = lagrangeMultiplier(0,0);
            myStructure.ConstraintLagrangeGetMultiplier(constraintRHS,lagrangeMultiplier);
            SinglePlotData(0,5) = lagrangeMultiplier(0,0);

            //number of Newton iterations
            SinglePlotData(0,6) = numNewtonIterations;

            PlotData.AppendRows(SinglePlotData);
            PlotData.WriteToFile("ConstraintsLagrange2DLoadDisp.txt"," ","#disp left and right, boundary force left and right, lagrange multiplier left and right","  ");
            if (PRINTRESULT)
            {
                std::cout << "disp left and right, boundary force left and right, lagrange multiplier left and right" << std::endl;
                SinglePlotData.Trans().Info();
            }

            myStructure.NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
            if (curDispFactor==1)
                convergenceStatusLoadSteps=true;
            else
            {
                //eventually increase load step
                if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3)
                {
                    deltaDispFactor*=1.5;
                    if (deltaDispFactor>maxDeltaDispFactor)
                        deltaDispFactor = maxDeltaDispFactor;
                }

                //increase displacement
                curDispFactor+=deltaDispFactor;
                if (curDispFactor>1)
                {
                    deltaDispFactor -= curDispFactor -1.;
                    curDispFactor=1;
                }

                curDisp = maxDisp*curDispFactor;

                //old stiffness matrix is used in first step of next load increment in order to prevent spurious problems at the boundary
                //std::cout << "press enter to next load increment, delta disp factor " << deltaDispFactor << " max delta disp factor " <<  maxDeltaDispFactor << std::endl << std::endl;
                //char cDummy[100]="";
                //std::cin.getline(cDummy, 100);
            }
            loadstep++;
        }
        else
        {
            assert(convergenceStatus==2);
            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
            //otherwise, the additional boundary displacements will result in an artificial localization in elements at the boundary
            curDispFactor-=deltaDispFactor;
            curDisp = maxDisp*curDispFactor;

            myStructure.ConstraintSetRHS(constraintRHS,-curDisp);

            // build global dof numbering
            myStructure.NodeBuildGlobalDofs();

            //set previous converged displacements
            myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
            myStructure.ElementTotalUpdateTmpStaticData();

            //decrease load step
            deltaDispFactor*=0.5;
            curDispFactor+=deltaDispFactor;
            curDisp = maxDisp*curDispFactor;

            //check for minimum delta (this mostly indicates an error in the software
            if (deltaDispFactor<MIN_DELTA_STRAIN_FACTOR)
            {
                throw NuTo::MechanicsException("Example ConstraintsLagrange1D : No convergence, delta strain factor < 1e-7");
            }

            std::cout << "press enter to reduce load increment" << std::endl;
            char cDummy[100]="";
            std::cin.getline(cDummy, 100);;
        }

        if (!convergenceStatusLoadSteps)
        {
            //update new displacement of RHS
            myStructure.ConstraintSetRHS(constraintRHS,-curDisp);

            // build global dof numbering and conre mat
            myStructure.NodeBuildGlobalDofs();

            //update stiffness in order to calculate new dispForceVector
            myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            if (PRINTRESULT)
                std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

            //update displacements of all nodes according to the new conre mat
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            myStructure.NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            myStructure.ElementTotalUpdateTmpStaticData();

            // calculate initial residual for next load step
            myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;

        }
    }
    NuTo::FullMatrix<double> PlotDataRef(6,7);
    PlotDataRef(0,0) = -0.0;
    PlotDataRef(1,0) = -0.2;
    PlotDataRef(2,0) = -0.4;
    PlotDataRef(3,0) = -0.5;
    PlotDataRef(4,0) = -0.5;
    PlotDataRef(5,0) = -0.5;

    PlotDataRef(0,1) = -0.0;
    PlotDataRef(1,1) = -0.2;
    PlotDataRef(2,1) = -0.4;
    PlotDataRef(3,1) = -0.6;
    PlotDataRef(4,1) = -0.8;
    PlotDataRef(5,1) = -1.0;


    PlotDataRef(0,2) = 0.0;
    PlotDataRef(1,2) = 0.0;
    PlotDataRef(2,2) = 0.0;
    PlotDataRef(3,2) = 0.1;
    PlotDataRef(4,2) = 0.3;
    PlotDataRef(5,2) = 0.5;

    PlotDataRef(0,3) = 0.0;
    PlotDataRef(1,3) = 0.0;
    PlotDataRef(2,3) = 0.0;
    PlotDataRef(3,3) = -0.1;
    PlotDataRef(4,3) = -0.3;
    PlotDataRef(5,3) = -0.5;

    PlotDataRef(0,4) = 0.0;
    PlotDataRef(1,4) = 0.0;
    PlotDataRef(2,4) = 0.0;
    PlotDataRef(3,4) = 0.05;
    PlotDataRef(4,4) = 0.15;
    PlotDataRef(5,4) = 0.25;

    PlotDataRef(0,5) = 0.0;
    PlotDataRef(1,5) = 0.0;
    PlotDataRef(2,5) = 0.0;
    PlotDataRef(3,5) = 0.05;
    PlotDataRef(4,5) = 0.15;
    PlotDataRef(5,5) = 0.25;

    PlotDataRef(0,6) = 0;
    PlotDataRef(1,6) = 1;
    PlotDataRef(2,6) = 1;
    PlotDataRef(3,6) = 2;
    PlotDataRef(4,6) = 1;
    PlotDataRef(5,6) = 1;

    if ((PlotDataRef-PlotData).Abs().Max()>1e-4)
    {
        std::cout<< "final results stored in load disp file as well" << std::endl;
        PlotData.Info();
        std::cout<< "reference results" << std::endl;
        PlotDataRef.Info();
        std::cout << "[ConstraintLagrange1D] result is not correct." << std::endl;
        return -1;
    }
    else
        std::cout<< "[ConstraintLagrange1D] nice, result is correct" << std::endl;
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
return 0;
}