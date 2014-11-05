// $Id$
#ifndef TIMEINTEGRATIONENUM_H_
#define TIMEINTEGRATIONENUM_H_

namespace NuTo
{
namespace TimeIntegration
{
enum eConvergenceState
{
	CONTINUE_ITERATION=0,           //!< continue iteration to converge
	CONVERGED,                      //!< converged
    DECREASE_TIME_STEP              //!< no convergence, decrease load step
};

enum eResultType
{
	TIME=0,                         //!< time
	NODE_DISPLACEMENT,              //!< nodal diplacement
	NODE_ACCELERATION,              //!< nodal diplacement
	GROUP_NODE_FORCE
};

}
}
#endif /* TIMEINTEGRATIONENUM_H_ */