/*  File    : mixedmex.c
*  Abstract:
*
*      Example MEX-file for a single-input, single-output system
*      with mixed continuous and discrete states.
*
*      Syntax:  [sys, x0] = mixedmex(t,x,u,flag)
*
*  Copyright 1990-2013 The MathWorks, Inc.
*/

#define S_FUNCTION_NAME  testone
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "math.h"
#include<windows.h>

#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */
#define learnsp         0.25
#define SampleTime      0.02
#define SAMPLE_SIZE  150
#define  coef		0.02

double r_last = 0;
double x1 = 0;
double x2 = 0;
double x3 = 0;
double x2_l = 0;
double x3_l = 0;
double q1 = 0;
double q2 = 0;
double q3 = 0;
double y_k = 0;
int i = 0;
double Error2 = 5;
double p[6] = { 0 };
double weight_1[3][2] = { { 1,-1 },{ 1,-1 },{ 1,-1 } };
double weight_2[3] = { 1.38,0.0,0.78 };
double new_value_1[3][2] = { { 0,0 },{ 0,0 },{ 0,0 }, };
double new_value_2[3] = { 0 };

double symble_com(double m1, double m2, double m3, double m4)
{
	double t = (m1-m2)*(m3-m4);
	if (t < 0)
	{
		return -1;
	}
	else if (t>0)
		return 1;
	else
		return 0;
}
struct Pid_tuple
{
	double m_p;
	double m_i;
	double m_d;
};
struct PID_Unit
{
	struct Pid_tuple m_X;
	struct Pid_tuple m_Q;
	double m_r;
	double m_y;
	double m_u;
	double m_err;
};
struct PID_comput
{
	struct PID_Unit m_units[SAMPLE_SIZE];
	int m_cur;
	double m_err2;
};
int assign(struct PID_comput* self, int i, double *ptr)
{
	self->m_cur = i;
	self->m_units[self->m_cur].m_X.m_p = ptr[0];
	self->m_units[self->m_cur].m_X.m_i = ptr[1];
	self->m_units[self->m_cur].m_X.m_d = ptr[2];
	self->m_units[self->m_cur].m_Q.m_p = ptr[3];
	self->m_units[self->m_cur].m_Q.m_i = ptr[4];
	self->m_units[self->m_cur].m_Q.m_d = ptr[5];
	self->m_units[self->m_cur].m_r = ptr[6];
	self->m_units[self->m_cur].m_y = ptr[7];
	self->m_units[self->m_cur].m_u = ptr[8];
	self->m_units[self->m_cur].m_err = ptr[9];
	self->m_err2 += self->m_units[self->m_cur].m_err * self->m_units[self->m_cur].m_err;

	//printf("X1[%d]=%f", i, self->m_units[i].m_X.m_p);
	return 0;
}
double update_weight_2nd(double* self, double* new_value)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		self[i] += learnsp * new_value[i];
	}
	return *self;
}
double sum_once(double _1_, double _2_, double _3_, double _4_)
{
	return symble_com(_1_, _2_, _3_, _4_);
}
double sum_once_err(double err, double next_y, double curr_y, double curr_u, double prev_u)
{
	return err*sum_once(next_y, curr_y, curr_u, prev_u);
}
double sum(double* self, struct PID_Unit prev, struct PID_Unit curr, struct PID_Unit next)
{
	double value = sum_once_err(curr.m_err, next.m_y, curr.m_y, curr.m_u, prev.m_u);
	self[0] += value * curr.m_Q.m_p;
	self[1] += value * curr.m_Q.m_i;
	self[2] += value * curr.m_Q.m_d;
	return *self;
}
double avg_2nd(double* self, int count)		//average sum_gradient 
{
	int i;
	for (i = 0; i < 3; i++)
	{
		self[i] /= count;
	}
	return *self;
}
int update_weight_1st( )
{
	for (int i = 0; i < 3;i++)
	{
		weight_1[i][0] += learnsp * new_value_1[i][0];
		weight_1[i][1] += learnsp * new_value_1[i][1];
	}
	return 0;
}
int sum_1st( struct PID_Unit curr, struct PID_Unit prev, struct PID_Unit next)
{
	double value_p = sum_once_err(curr.m_err, next.m_y, curr.m_y, curr.m_u, prev.m_u) * symble_com(curr.m_Q.m_p, prev.m_Q.m_p, curr.m_X.m_p, prev.m_X.m_p);
	double value_i = sum_once_err(curr.m_err, next.m_y, curr.m_y, curr.m_u, prev.m_u) * symble_com(curr.m_Q.m_i, prev.m_Q.m_i, curr.m_X.m_i, prev.m_X.m_i);
	double value_d = sum_once_err(curr.m_err, next.m_y, curr.m_y, curr.m_u, prev.m_u) * symble_com(curr.m_Q.m_d, prev.m_Q.m_d, curr.m_X.m_d, prev.m_X.m_d);
	new_value_1[0][0] += value_p*curr.m_Q.m_p*curr.m_X.m_p*curr.m_r;
	new_value_1[0][1] += value_p*curr.m_Q.m_p*curr.m_X.m_p*curr.m_y;
	new_value_1[1][0] += value_p*curr.m_Q.m_i*curr.m_X.m_i*curr.m_r;
	new_value_1[1][1] += value_p*curr.m_Q.m_i*curr.m_X.m_i*curr.m_y;
	new_value_1[2][0] += value_p*curr.m_Q.m_d*curr.m_X.m_d*curr.m_r;
	new_value_1[2][1] += value_p*curr.m_Q.m_d*curr.m_X.m_d*curr.m_y;
	return 0;
}
int avg_1st( int count)
{
	for (int i = 0; i < 3;i++)
	{
		new_value_1[i][0] /= count;
		new_value_1[i][1] /= count;
	}
	return 0;
}



/*====================*
* S-function methods *
*====================*/

/* Function: mdlInitializeSizes ===============================================
* Abstract:
*    The sizes information is used by Simulink to determine the S-function
*    block's characteristics (number of inputs, outputs, states, etc.).
*/
static void mdlInitializeSizes(SimStruct *S)
{
	ssSetNumSFcnParams(S, 0);  /* Number of expected parameters */
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
		return; /* Parameter mismatch will be reported by Simulink */
	}

	ssSetNumContStates(S, 2);
	ssSetNumDiscStates(S, 2);

	if (!ssSetNumInputPorts(S, 1)) return;
	ssSetInputPortWidth(S, 0, 2);
	ssSetInputPortDirectFeedThrough(S, 0, 0);

	if (!ssSetNumOutputPorts(S, 1)) return;
	ssSetOutputPortWidth(S, 0, 1);

	ssSetNumSampleTimes(S, 2);
	ssSetNumRWork(S, 1);
	ssSetNumIWork(S, 0);
	ssSetNumPWork(S, 0);
	ssSetNumModes(S, 0);
	ssSetNumNonsampledZCs(S, 0);
	ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

	ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}



/* Function: mdlInitializeSampleTimes =========================================
* Abstract:
*    Two tasks: One continuous, one with discrete sample time of 1.0.
*
* Note: This S-function is block based one. To run this S-function in a
*       multi-tasking or auto mode, it will have a warning message suggesting
*       to convert this to PORT_BASED_SAMPLE_TIME S-function.
*       See sfun_multirate.c file as an example.
*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
	ssSetSampleTime(S, 1, SampleTime);

	ssSetOffsetTime(S, 0, 0.0);
	ssSetOffsetTime(S, 1, 0.0);
	ssSetModelReferenceSampleTimeDefaultInheritance(S);
}



#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
* Abstract:
*    Initialize both states to one
*/
static void mdlInitializeConditions(SimStruct *S)
{
	real_T *U_Con = ssGetContStates(S);				//define input constant
	real_T *U_Dis = ssGetRealDiscStates(S);			//define input discrete
	//real_T *lasty = ssGetRWork(S);

	/* initialize the states to 0.0 */
	U_Con[0] = 0.0;
	U_Con[1] = 0.0;
	U_Dis[0] = 0.0;
	U_Dis[1] = 0.0;

	/* Set initial output to initial value of delay */
	//lasty[0] = 0.0;
}


/* Function: mdlOutputs =======================================================
* Abstract:
*      y = x
*     If this is a sample hit for the discrete task, set the output to the
*     value of the discrete state, otherwise, maintain the current output.
*/

//double X1[150] = { 0 };
//double X2[150] = { 0 };
//double X3[150] = { 0 };


struct PID_comput computer;

static void mdlOutputs(SimStruct *S, int_T tid)
{

	real_T *y = ssGetOutputPortRealSignal(S, 0);
	real_T *U_Dis = ssGetRealDiscStates(S);
	//real_T *lasty = ssGetRWork(S);

    
	UNUSED_ARG(tid); /* not used in single tasking mode */
	p[6] = U_Dis[0];
	p[7] = U_Dis[1];
	x1 = weight_1[0][0] * U_Dis[0] + weight_1[0][1] * U_Dis[1];
    x2= weight_1[1][0] * U_Dis[0] + weight_1[1][1] * U_Dis[1];
    x3= weight_1[2][0] * U_Dis[0] + weight_1[2][1] * U_Dis[1];
	//x1 =  U_Dis[0] -  U_Dis[1];
	//x2 =  U_Dis[0] -  U_Dis[1];
	//x3 = U_Dis[0] - U_Dis[1];
    q1=x1;
	q2 = (x2_l + x2) * coef;
	q3 = (x3 - x3_l) / coef;
	y_k = weight_2[0] * q1 + weight_2[1] * q2 + weight_2[2] * q3;
	y[0] = y_k;
	p[8] = y_k;
	p[9]= U_Dis[0] - U_Dis[1];

  	//y[0] = InputDis[0];
	//     if (ssIsSampleHit(S, 1, tid)) {
	//         y[0]     = xD[0];
	//         lasty[0] = xD[0];
	//     } else {
	//         y[0] = lasty[0];
	//     }
}

#define MDL_UPDATE
/* Function: mdlUpdate ======================================================
* Abstract:
*      xD = xC
*/
static void mdlUpdate(SimStruct *S, int_T tid)
{
	real_T *U_Dis = ssGetRealDiscStates(S);
	real_T *U_Con = ssGetContStates(S);
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
	U_Con[0] = U(0);
	U_Con[1] = U(1);
	UNUSED_ARG(tid); /* not used in single tasking mode */

	if (ssIsSampleHit(S, 1, tid)) {
		U_Dis[0] = U_Con[0];
		U_Dis[1] = U_Con[1];
            x2_l=x2;
            x3_l=x3;
			//X1[i] =  x1;
			//X2[i] =  x2;
			//X3[i] =  x3;
			//printf("x1=%f\n",x1);
			 p[0] = x1;
			 p[1] = x2;
			 p[2] = x3;
			 p[3] = q1;
			 p[4] = q2;
			 p[5] = q3;
			 
			 if (U_Dis[0] != r_last)
			 {
				 i = 0;
				 computer.m_err2 = 0;
			 }
			 if (i < 149 && i >= 0)
			 {
				 assign(&computer, i, p);
				// printf("i=%d\n", i);
			 }
			 else if (i == 149)
			 {
				 assign(&computer, i, p);
				 printf("i=%d\n", i);
				// printf("X1[%d]=%f\n", i, computer.m_units[i].m_X.m_p);
				 if (computer.m_err2 < Error2)
				 {
					 for (int k = 1; k < 148; k++)
					 {
						 *new_value_2 = sum(new_value_2, computer.m_units[k - 1], computer.m_units[k], computer.m_units[k + 1]);
						  sum_1st( computer.m_units[k], computer.m_units[k - 1], computer.m_units[k + 1]);
					 }
					 *new_value_2 = avg_2nd(new_value_2, SAMPLE_SIZE);
					 *weight_2 = update_weight_2nd(weight_2, new_value_2);
					 avg_1st( SAMPLE_SIZE);
					 update_weight_1st();

					 printf("weight_2[3]={%f,%f,%f}\n", weight_2[0], weight_2[1], weight_2[2]);
					 printf("weight_1[3][2]={{%f,%f},{%f,%f},{%f,%f}}\n", weight_1[0][0], weight_1[0][1], weight_1[1][0], weight_1[1][1], weight_1[2][0], weight_1[2][1]);
				 }
				 Error2 = computer.m_err2;
				 computer.m_err2 = 0;
					 i = -1;
				 
			 }
			 else
				 i = -1;

			 i++;
			 r_last = U_Dis[0];
	}

}



#define MDL_DERIVATIVES
/* Function: mdlDerivatives =================================================
* Abstract:
*      xdot = u
*/
static void mdlDerivatives(SimStruct *S)
{
	//     real_T            *dx   = ssGetdX(S);
	//     InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
	// 
	//     dx[0] = U(0);
}



/* Function: mdlTerminate =====================================================
* Abstract:
*    No termination needed, but we are required to have this routine.
*/
static void mdlTerminate(SimStruct *S)
{
	UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
