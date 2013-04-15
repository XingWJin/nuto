import sys
import nuto
import os
import math

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = False

#system name and processor
system = sys.argv[1]+sys.argv[2]

#path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(sys.argv[3],"results",system,os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% start the real test file                                      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def SetParameters(rParameterTuple):
   global x,y
   x = rParameterTuple[0]
   y = rParameterTuple[1]
   if (printResult):
       print "[Python - Parameters routine]"
       print x ,y
       print 

def Objective():
   o = scaleFactorX/12*(x-2)*(x-2)*(x-2)*(x-2)+0.5*(y-3)*(y-3)
   #o = scaleFactorX*0.5*(x-2)*(x-2)+0.5*(y-3)*(y-3)
   if (printResult):
       print "[Python - Objective routine] (x,y)=(" + str(x) + "," + str(y) + ") : obj=" + str(o)
       print 
   return o
   #print rGradient.size()
def Gradient():
   rGradient = [scaleFactorX/3*(x-2)*(x-2)*(x-2),(y-3)]
   #rGradient = [scaleFactorX*(x-2),(y-3)]
   if (printResult):
       print "[Python - Gradient routine]"
       print rGradient
       print 
   return rGradient
   
def Hessian():
   rHessian = [scaleFactorX*(x-2)*(x-2),0,0,1]
   #rHessian = [scaleFactorX,0,0,1]
   if (printResult):
       print "[Python - Hessian routine]" 
       print rHessian
       print 
   return rHessian

x = 2.2
y = 15

scaleFactorX = 1e3

InitParameters = nuto.DoubleFullMatrix(2,1)
InitParameters.SetValue(0,0,x)
InitParameters.SetValue(1,0,y)

PythonCallback = nuto.CallbackHandlerPython()
PythonCallback.SetCallbackFunctions(SetParameters,Objective,Gradient,Hessian)

myOptimizer = nuto.ConjugateGradient(2)
#myOptimizer.Optimize()
#exit(0)
myOptimizer.SetCallback(PythonCallback)
myOptimizer.SetParameters(InitParameters)

myOptimizer.SetMaxFunctionCalls(1000000)
myOptimizer.SetMaxGradientCalls(1000)
myOptimizer.SetMaxHessianCalls(1000)
myOptimizer.SetMaxIterations(1000)
myOptimizer.SetMinObjective(0)
myOptimizer.SetMinDeltaObjBetweenRestarts(1e-6)
myOptimizer.SetAccuracyGradient(1e-12)

if (printResult):
    myOptimizer.SetVerboseLevel(10)
else:
    myOptimizer.SetVerboseLevel(0)

returnValue = myOptimizer.Optimize()
if (printResult):
    print "Final objective : ", myOptimizer.GetObjective()
if (createResult):
    f = open(pathToResultFiles+'Objective.txt','w')
    f.write('#Correct Objective\n')
    f.write(str(myOptimizer.GetObjective()))
    f.close()
else:
    f = open(pathToResultFiles+'Objective.txt','r')
    f.readline()
    objectiveExact = float(f.readline())
    f.close()
    if (math.fabs(myOptimizer.GetObjective()-objectiveExact)>1e-8):
        print '[' + system,sys.argv[0] + '] : objective is not correct.'
        error = True;

if (printResult):
    print "Final Set of Parameters"
    myOptimizer.GetParameters().Trans().Info()
    print ""
if createResult:
	myOptimizer.GetParameters().WriteToFile(pathToResultFiles+"Parameters.txt"," ","#Correct result","  ")
else:
    ParametersExact = nuto.DoubleFullMatrix(1,1)
    ParametersExact.ReadFromFile(pathToResultFiles+"Parameters.txt",1," ")
    if ((ParametersExact-myOptimizer.GetParameters()).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : Parameters is not correct.'
        error = True;


if (error):
    sys.exit(-1)
else:
    sys.exit(0)
