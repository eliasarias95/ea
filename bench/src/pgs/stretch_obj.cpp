/*************************************************************************
 * stretch_obj.cpp
 * Author:     Zijian Liu
 *************************************************************************/
 
#include <stretch_obj.hpp>
using namespace std;

/**************************************************************
Constructors
Parameters:
in_dz              Sample rate for input data (unit: meter/sec)
in_nz              Sample total for input data 
out_dz             Sample rate for output data (unit: meter/sec)
out_nz             Sample total for output data
nh                 Number of bins with in one CDP gather, 1 for stack input
mode               T2D (1) or D2T (2)
stretch_vel        Define if you want to output both stretched
                   data and velocity, 1 for yes and 0 for no
****************************************************************/
stretch_obj::stretch_obj(float in_dz_, int in_nz_, float mdl_dz_, int mdl_nz_, 
		            float out_dz_, int out_nz_, int nh_, int mode_, int stretch_vel_) {

   int i;
// passing parameters and check
   mode = mode_;
   in_nz = in_nz_;
   out_nz = out_nz_;
   in_dz = in_dz_;
   out_dz = out_dz_;
   stretch_vel = stretch_vel_;
   nh = nh_;
   mdl_dz = mdl_dz_;
   mdl_nz = mdl_nz_;
   
   if ((in_nz<0)||(out_nz<0)) {
      ucsl_error("STRETCH ERROR:Spcified sample total can not be negative!\n");
   }
   if ((in_dz<0.0)||(out_dz<0.0)) {
      ucsl_error("STRETCH ERROR:Spcified sample interval can not be negative!\n");
   }
   
   if (mode == T2D) {
	   if (in_dz>1.0) in_dz /= 1000.0;
	   if (mdl_dz>1.0) mdl_dz /= 1000.0;
   }
   if (mode == D2T) {
   	   if (out_dz>1.0) out_dz /= 1000.0;
      }
   
   
 // allocate related buffers
   dat_ind = in_ind = NULL;
   dat_fac = NULL;
   dat_fac = NULL;
   inv_vp = NULL;
   mdl_ind = NULL;
   mdl_resamp = NULL;
   vel_resamp = NULL;
   inv_mdl = NULL;
   
   dat_ind = (float*)mem_alloc(out_nz*sizeof(float));
   in_ind = (float*)mem_alloc(in_nz*sizeof(float));
   dat_fac = (float**)mem_alloc2(4, out_nz,sizeof(float));
   dat_tmp = (int*)mem_alloc(out_nz*sizeof(int));
   inv_vp = (float*)mem_alloc(in_nz*sizeof(float));
   mdl_resamp = (float*)mem_alloc(in_nz*sizeof(float));
   vel_resamp = (float*)mem_alloc(in_nz*sizeof(float));
   if ((mdl_dz!=in_dz)||(mdl_nz!=in_nz)) {
	   mdl_ind = (float*)mem_alloc(mdl_nz*sizeof(float));
       for (i=0;i<mdl_nz;i++) mdl_ind[i] = i*mdl_dz;
       mdl_tmp = (int*)mem_alloc(in_nz*sizeof(int));
       mdl_fac = (float**)mem_alloc2(4, in_nz,sizeof(float));
       resamp_ind = (float*)mem_alloc(in_nz*sizeof(float));
       for (i=0;i<in_nz;i++)  resamp_ind[i] = i*in_dz;
       inv_mdl = (float*)mem_alloc(mdl_nz*sizeof(float));
   }
   
   if ((dat_ind == NULL)||(in_ind == NULL))
	   ucsl_error("STRETCH ERROR:Insurficient Memory!\n");
   if (dat_fac == NULL)
	   ucsl_error("STRETCH ERROR:Insurficient Memory!\n");
   if ((dat_tmp == NULL)||(inv_vp == NULL)) 
	   ucsl_error("STRETCH ERROR:Insurficient Memory!\n");

   for (i=0;i<out_nz;i++)  dat_ind[i] = i*out_dz;  
   for (i=0;i<in_nz;i++)  in_ind[i] = -9999.0;
}

//destructors
stretch_obj::~stretch_obj(){
   if (dat_ind) free(dat_ind);
   if (in_ind) free(in_ind);
   if (dat_fac) free(dat_fac);
   if (dat_tmp) free(dat_tmp);
   if (inv_vp) free(inv_vp);
   if (mdl_ind) free(mdl_ind);
   if (mdl_resamp) free(mdl_resamp);
   if (vel_resamp) free(vel_resamp);
   if (inv_mdl) free(inv_mdl);
}

int stretch_obj::process(float **in_, float *vel_, float **out_, float *vel2_){
   int i,flag, iz;

   if ((mdl_nz!=in_nz)||(mdl_dz!=in_dz)){
	  for (i=0;i<in_nz;i++)  in_ind[i] = i*in_dz;
	  for (i=0;i<mdl_nz;i++) inv_mdl[i] = 1.0/vel_[i];
      ccuint(mdl_ind,inv_mdl,mdl_nz,resamp_ind,vel_resamp,in_nz,mdl_tmp,mdl_fac,1,1);
      for (i=0;i<in_nz;i++)  vel_resamp[i] = 1.0/vel_resamp[i];;
   } else {
	  memcpy(vel_resamp, vel_, in_nz*sizeof(float)); 
   }
   
   if (stretch_TT(vel_resamp)==ERROR) {
      ucsl_error("STRETCH ERROR:error when computing travel time map!\n");  
      goto err;
   }

   for (iz=0; iz<in_nz; iz++) {
      inv_vp[iz] = 1.0/vel_resamp[iz];
   }


   for (i=0; i<nh; i++) {
      ccuint(in_ind,in_[i],in_nz,dat_ind,out_[i],out_nz,dat_tmp,dat_fac,1,0);
   }
   if (stretch_vel) {
	   ccuint(in_ind,inv_vp,in_nz,dat_ind,vel2_,out_nz,dat_tmp,dat_fac,1,1);
	   for (iz=0; iz<out_nz; iz++) vel2_[iz] = 1.0/vel2_[iz];
   }

   return GOOD;
   
err:
   return ERROR;
}

int stretch_obj::process_mdls(float *mdl_, float *mdl2_, float *vel_){
   int i;
   
   if ((mdl_nz!=in_nz)||(mdl_dz!=in_dz)){
	  for (i=0;i<in_nz;i++)  in_ind[i] = i*in_dz;
      ccuint(mdl_ind,mdl_,mdl_nz,resamp_ind,mdl_resamp,in_nz,mdl_tmp,mdl_fac,1,1);
   } else {
	  memcpy(mdl_resamp, mdl_, in_nz*sizeof(float));	  
   }

   if (in_ind[0]<0.0) {
      if (stretch_TT(vel_resamp)==ERROR) {
         ucsl_error("STRETCH ERROR:error when computing travel time map!\n");  
         return ERROR;
      }
   }
   
   ccuint(in_ind,mdl_resamp,in_nz, dat_ind, mdl2_,out_nz,dat_tmp,dat_fac,1,1);

   return GOOD;
}

float stretch_obj::hor_stretch(float hor_in, float *vel_) {
   float fac;
   int il,ih;
   float hor_out;

   if (in_ind[0]<0.0) {
      if (stretch_TT(vel_)==ERROR) {
      ucsl_error("STRETCH ERROR:error when computing travel time map!\n");  
      return -999.0;
      }
   }   


   hor_in /= in_dz;
   fac = hor_in - floor(hor_in);
   il = int(floor(hor_in));
   if (il<0){
      ucsl_error("Index can not be minus!\n");  
	  return -999.0;
   }
   if (il<in_nz-1) ih = il+1;
   else ih=il;
   hor_out = in_ind[il]*(1-fac)+in_ind[ih]*fac;
  
   if (mode == D2T) hor_out *= 1000.0;
   return hor_out;
}



int stretch_obj::stretch_TT(float *vel) {
   int z,t;
   float sumt,sumz;
   
   if (mode == D2T) {
		sumt = 0.0;
		in_ind[0] = 0.0;
		for (z = 1; z < in_nz - 1; z++) {
			sumt = sumt + in_dz * (1.0 / vel[z-1] + 1.0 / vel[z]); //two way travel time
			in_ind[z] = sumt;
		}
		in_ind[in_nz - 1] = sumt + in_dz * (1.0 / vel[in_nz - 1]);
	} else { //T2D
		sumz = 0.0;
		in_ind[0] = 0.0;
		for (t = 1; t < in_nz - 1; t++) {
			sumz = sumz + 0.25 * in_dz * (vel[t-1] + vel[t]); //two-way trave time
			in_ind[t] = sumz;
		}
		in_ind[in_nz - 1] = sumz + in_dz * vel[in_nz - 1];
	}
   
   
   return GOOD;
}


void stretch_obj::ccuint (float *X1, float *Y1, int N1, float *X2, float *Y2, int N2, int *IZ, float **ZZ, int INIT, int EXT){
   float XX[64];
   float X1LO,X1HI;
   int J,I1,I2,II,I;
   int J1,J2,J3,J4;
   float DX1,DX2,DX3,DX4;
   float D12,D13,D14,D23,D34,D42;
   
   if ((N1<4)||(N2<1)) return;
   if (INIT) {
      X1LO = X1[0];
      X1HI = X1[N1-1];
      J  = 3;
      I2 = 0;
      while (I2<N2) {
         I1 = I2 + 1;
         I2 = I2 + 64;
         if (I2>N2) I2 = N2;
	 //DO FIX UP FOR OUT-OF-RANGE VALUES OF X2
	 II=0;
	 for (I=I1;I<=I2;I++) {
	    II=II+1;
	    if (X2[I-1]<X1LO) XX[II-1]=X1LO;
	    else if (X2[I-1] >= X1HI) XX[II-1] = X1HI;
	    else XX[II-1] = X2[I-1];
	 }
	//   CALCULATE IZ
         II = 0;
	 for (I=I1; I<=I2; I++) {
	    II = II+1;
	    while ((XX[II-1]>X1[J-1])&&(J<N1-1)) J=J+1;
	    IZ[I-1] = J-2;
	 }
	 // CALCULATE ZZ
	 II=0;
	 for (I=I1;I<=I2; I++) {
	    II=II+1; 
	    J1 = IZ[I-1];
            J2 = J1 + 1;
            J3 = J2 + 1;
            J4 = J3 + 1;

            DX1 = XX[II-1]-X1[J1-1];
            DX2 = XX[II-1]-X1[J2-1];
            DX3 = XX[II-1]-X1[J3-1];
            DX4 = XX[II-1]-X1[J4-1];
            D12 = X1[J1-1]-X1[J2-1];
            D13 = X1[J1-1]-X1[J3-1];
            D14 = X1[J1-1]-X1[J4-1];
            D23 = X1[J2-1]-X1[J3-1];
            D34 = X1[J3-1]-X1[J4-1];
            D42 = X1[J4-1]-X1[J2-1];

            ZZ[I-1][0]=DX2*DX3*DX4/(D12*D13*D14);
            ZZ[I-1][1]=DX1*DX3*DX4/(D12*D23*D42);
            ZZ[I-1][2]=DX1*DX2*DX4/(D13*D23*D34);
            ZZ[I-1][3]=DX1*DX2*DX3/(D14*D42*D34);
	 }
      }
   }   
      //---------------------
     //PERFORM INTERPOLATION
     //---------------------
   if (EXT==0) {
      for (I=1;I<=N2;I++) {
         J = IZ[I-1];
         if (X2[I-1]<=X1HI) {
            Y2[I-1] = ZZ[I-1][0]*Y1[J-1]+ZZ[I-1][1]*Y1[J]+ZZ[I-1][2]*Y1[J+1]+ZZ[I-1][3]*Y1[J+2];
         } else {
            Y2[I-1] = 0.0;
         }
      }
   }
   else {
      for (I=1;I<=N2;I++) {
         J = IZ[I-1];
         Y2[I-1] = ZZ[I-1][0]*Y1[J-1]+ZZ[I-1][1]*Y1[J]+ZZ[I-1][2]*Y1[J+1]+ZZ[I-1][3]*Y1[J+2];
      }
    }

   return;

}