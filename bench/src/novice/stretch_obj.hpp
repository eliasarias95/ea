/*************************************************************************
 * stretch_obj.hpp
 * Author:     Zijian Liu
 *************************************************************************/
#ifndef _STRETCH 
#define _STRETCH

#include <ucsl.h>
#define   T2D        1
#define   D2T        2
#define   ERROR      -1
#define   GOOD       0

class stretch_obj {

public:

/**************************************************************
Constructors
Parameters:
mode               T2D (1) or D2T (2)
stretch_vel        Define if you want to output stretched
                   velocity, 1 for yes and 0 for no
****************************************************************/
	stretch_obj(float in_dz_, int in_nz_, float mdl_dz_, int mdl_nz_, float out_dz_, int out_nz_, int nh_, int mode_, int stretch_vel_);	
	
    ~stretch_obj();	
	
/**************************************************************
function process
Parameters:
in                 Input data of single CDP, with the dimension as [pd->nh][pd->in_nz]
vel                Input velocity of single CDP, with the dimension as [pd->in_nz]
out                Output data of single CDP, with the dimension as [pd->nh][pd->out_nz]
out_nz             Sample total for output data
vel2               Output velocity of single CDP, with the dimension as [pd->out_nz]
                   This is only valid when stretch_vel=1
****************************************************************/
	int process(float **in_, float *vel_, float **out_, float *vel2_);
	
	
	/**************************************************************
	function process_mdl: stretch multiple models, from mdl_dtz to mdl_dzt
	****************************************************************/
	int process_mdls(float *mdl_, float *mdl2_, float *vel_);
	
	/**************************************************************
	function hor_stretch: stretch single horizon point
	****************************************************************/	
	float hor_stretch(float hor_in, float *vel_);
    
    int mode,stretch_vel;                      // T2D=1, D2T=2
    float in_dz;
    float out_dz;
    int in_nz;
    int nh;
    int out_nz;
    float mdl_dz;
    int mdl_nz;
    
	
private:
	
	int stretch_TT(float *vel_);
    void ccuint (float *X1, float *Y1, int N1, float *X2, float *Y2, int N2, int *IZ, float **ZZ, int INIT, int EXT);
	
	float *dat_ind, *in_ind, *resamp_ind;
    float **dat_fac, **mdl_fac;
    int *dat_tmp, *mdl_tmp;
    float *inv_vp, *inv_mdl;
    float *mdl_ind,*mdl_resamp;
    float *vel_resamp;


};

#endif 
