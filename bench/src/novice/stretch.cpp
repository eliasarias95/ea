/*************************************************************************
 * Author:     Zijian Liu
 *************************************************************************/
#include <ucsl.h>
#include <stretch_obj.hpp>
#include <omp.h>

#define   ERROR      -1
#define   GOOD       0
#define   GDIM       4
#define   MDIM       3
#define   TBLOCK       1024

char *usage[] = {
"************************************************************",
" stretch  par=  [key=val]",
" Parameters:",
" gth_in=                input gathers/stack (required)",
" gth_out=               stretched gathers/stack (required)",
" vel_in=            Input velocity model (required)",
" vel_out=           stretched velocity model (optional)",
" mdl_in=            Input additional model (optional & repeatable)",
" mdl_out=           stretched additional model (optional & repeatable)",
" project=           optional project name",
" mode=              stretch mode=D2T/T2D",
" out_dzt=           sample rate for output data/model (required)",
" out_nzt=           Sample total for output data/model (optional, auto evaluate if missing input)",
" freecoords=        Input freecoords file (optional)",
" subline_seq=       min and max of subline range you want to process",
" xline_seq=         min and max of xline range you want to process",
"************************************************************",
NULL};

typedef struct {
    char hostname[128];
    int nprocs,iproc;
    stretch_obj **str;
    char *fn_in, *fn_out, *fn_vel, *fn_vel2;
    char *fn_hor;
    char **fn_mdl, **fn_mdl2;
    int nmdl;
    char *project;
    int mode;
    int out_nzt,in_nzt,vel_nzt;
    float out_dzt,in_dzt,vel_dzt;
    file_trace *fd_in,*fd_out,*fd_vel;
    file_trace **fd_mdl,**fd_mdl2;
    file_mdvd *fd_vel2;
    FILE *fp_hor;
    hdr_map *hmap_dat;
    hdr_map *hmap_vel,*hmap_vel2;
    hdr_map **hmap_mdl;
    axis **axis_in,**axis_out,**axis_vel,**axis_vel2;
    int nh;		    
    int mype;
    int npe;   
    
    int flag_hdr,flag_fc;
    int **mpi_range;
    int nhdr_dat,nhdr_vel,nhdr_vel2;
    int *nhdr_mdl;
    
    float ***traces_in, ***dtraces_in;
    float ***traces_out, ***dtraces_out;
    float **traces_vel, **dtraces_vel;
    float **traces_vel2, **dtraces_vel2;
    float ***traces_mdl, ***dtraces_mdl;
    float ***traces_mdl2, ***dtraces_mdl2;
    float *hor_map;   
    
    int submin,submax;
    int crsmin,crsmax;
    
} prog_data;

//*************************************************************************************
void initialize(prog_data *pd, parlist *par);
void process(prog_data *pd);
void clean_up(prog_data *pd);

void open_file(char *fn, file_trace **fd, hdr_map **hmap, axis ***ax, int dims);
int estimate_nzt(file_trace *fd, hdr_map *hmap, axis *axz, float dzt_out, int mode);
void stretch(stretch_obj *str, float **in, float *vel, float **out, float *vel2, 
                             float **mdl, float **mdl2, int nmdl, float *hor, int ind);
void check_cons(axis **ax_d, axis **ax_v, int submin, int submax, int crsmin, int crsmax);
int get_ntrace(axis **ax_v, int submin, int submax, int crsmin, int crsmax);
int local_index(int itr, axis **ax_out, axis **ax_in, int ndim);

//************************************************************************************

int main (int argc, char *argv[]) {
    parlist *par;
    prog_data *pd;
    int ipe;
    int dims;

    par = ucsl_initialize(usage, argv, argc, 1, 1);
    pd = (prog_data*)mem_alloc(sizeof(prog_data));
    initialize(pd, par);
    delete par;
    ucsl_abort_if_error();
    ucsl_barrier( MPI_COMM_WORLD);
    process(pd);
    ucsl_abort_if_error();
   
    for (ipe = 0; ipe < pd->npe; ipe++) {
		ucsl_barrier( MPI_COMM_WORLD);
		if (ipe == pd->mype) {
			ucsl_printf("****** PE %d - %s ******\n", pd->mype, pd->hostname);
			fflush( _UCSL_LOG);
		}
	}
   
    fflush(_UCSL_LOG);
    ucsl_barrier(MPI_COMM_WORLD);
    if(!pd->mype) ucsl_printf("Successful Task 1 of 1\n");
    clean_up(pd);

    ucsl_finalize();
    return 0;
}

/*************************************************************************
 * function initialize 
 *************************************************************************/
void initialize(prog_data *pd, parlist *par) {
	
	char *ctmp;
	int nmdl_tmp;
	int i,j,iax;
	int ncdp,nh,nhdr;
	int *node_traces;
	axis **axis_tmp;
	int wbz_exist;

//pointer initialize
    pd->str = NULL;
    pd->fn_in = NULL, pd->fn_out = NULL;
    pd->fn_vel = NULL, pd->fn_vel2 = NULL;
    pd->fn_hor = NULL;
    pd->fn_mdl = NULL, pd->fn_mdl2 = NULL;
    pd->fd_in = NULL,pd->fd_out = NULL;
    pd->fd_vel = NULL,pd->fd_vel2 = NULL;
    pd->fd_mdl = NULL,pd->fd_mdl2 = NULL;
    pd->fp_hor = NULL;
    pd->hmap_dat = NULL;
    pd->hmap_vel = NULL,pd->hmap_vel2 = NULL;
    pd->hmap_mdl = NULL;
    pd->axis_in = NULL,pd->axis_out = NULL;
    pd->axis_vel = NULL,pd->axis_vel2 = NULL;
    pd->traces_in = NULL, pd->dtraces_in = NULL;
    pd->traces_out = NULL, pd->dtraces_out = NULL;
    pd->traces_vel = NULL, pd->dtraces_vel = NULL;
    pd->traces_vel2 = NULL, pd->dtraces_vel2 = NULL;
    pd->traces_mdl = NULL, pd->dtraces_mdl = NULL;
    pd->traces_mdl2 = NULL, pd->dtraces_mdl2 = NULL;
    pd->hor_map = NULL;   
    pd->project = NULL;
	
	
//MPI and OpenMP initialize
    if (!par->get_status()) {
		ucsl_errorf("ERROR retrieving parameters!\n");
		ucsl_abort();
	}

	pd->mype = ucsl_mype();
	pd->npe = ucsl_npe();
	pd->nprocs = ucsl_get_n_cores();
	omp_set_num_threads(pd->nprocs);
	gethostname(pd->hostname, 128);

	if (!pd->mype) {
		ucsl_printf("************************************************************\n");
		ucsl_printf("stretch_ucsl\n");
		ucsl_printf("Version: %s\n", _UCSL_VERSION);
		ucsl_printf("************************************************************\n");
		ucsl_printf("Parameters:\n");
		par->print(_UCSL_LOG);
		ucsl_printf("************************************************************\n");
		ucsl_printf("There are %d nodes report to you.\n", pd->npe);
		ucsl_printf("Nodes:\n");
	}

	fflush( _UCSL_LOG);
	ucsl_barrier( MPI_COMM_WORLD);
	for (int ipe = 0; ipe < pd->npe; ipe++) {
		if (ipe == pd->mype) {
			ucsl_printf("%d\t%s\n", ipe, pd->hostname);
		}
		fflush(_UCSL_LOG);
		ucsl_barrier(MPI_COMM_WORLD);
	}

	ucsl_abort_if_error();		
	
//read-in parameters & open files
	par->get_string_req("gth_in", &pd->fn_in);
	par->get_string_req("gth_out", &pd->fn_out);
	par->get_string_req("vel_in", &pd->fn_vel);
	par->get_float_req("out_dzt", &pd->out_dzt);
	par->get_string("vel_out", &pd->fn_vel2);
	par->get_string("project", &pd->project);
	
	ctmp = NULL;
	par->get_string_req("mode", &ctmp);
	if(strcasecmp(ctmp,"D2T") == 0) {pd->mode = D2T;}
	else if (strcasecmp(ctmp,"T2D") == 0) {pd->mode = T2D;}
	else {ucsl_errorf("ERROR: Input mode must be either T2D or D2T!\n");}
	if (ctmp!=NULL) {free(ctmp); ctmp=NULL;}
	
	open_file(pd->fn_in, &(pd->fd_in), &(pd->hmap_dat), &(pd->axis_in), GDIM);
	open_file(pd->fn_vel, &(pd->fd_vel), &(pd->hmap_vel), &(pd->axis_vel), MDIM);
	pd->flag_hdr = 0;
	wbz_exist = 0;
	if (pd->mode == T2D) {
	   if (pd->hmap_vel->get_hdr_map_index("WBT")>0){
		   if ((pd->hmap_vel->get_hdr_map_index("WBZ"))<0){
			   pd->flag_hdr = 1;
			   if (!pd->mype) ucsl_printf("Program will stretch the header WBT to WBZ.\n");			
		   } else {
			   wbz_exist = 1;
		   }
	   }
	}	
		
	pd->nmdl = nmdl_tmp=0;
	pd->nmdl = get_trace_file_names(par, "mdl_in", &(pd->fn_mdl));
	nmdl_tmp = get_trace_file_names(par, "mdl_out", &(pd->fn_mdl2));
	if (pd->nmdl>0) {
		if (pd->nmdl!=nmdl_tmp) ucsl_errorf("ERROR: Totals are mismatch between input and output additoinal models!\n");
		pd->fd_mdl = (file_trace**)mem_alloc(pd->nmdl*sizeof(file_trace*));
		pd->fd_mdl2 = (file_trace**)mem_alloc(pd->nmdl*sizeof(file_trace*));
		pd->hmap_mdl = (hdr_map**)mem_alloc(pd->nmdl*sizeof(hdr_map*));
		
		for (i=0;i<pd->nmdl;i++){
     		open_file(pd->fn_mdl[i], &(pd->fd_mdl[i]), &(pd->hmap_mdl[i]), &(axis_tmp), MDIM);
			for (iax=0;iax<MDIM;iax++) {
				if (!(axis_tmp[iax]->is_equivalent(pd->axis_vel[iax]))) {
				   ucsl_errorf("ERROR: file %s has incompatible axes with input velocity.\n", pd->fn_mdl[i]);
				   ucsl_abort_if_error();
				}
			}
			for (iax=0; iax<MDIM; iax++) delete axis_tmp[iax];
		}		
	}
	
	pd->flag_fc = 0;	
	if ((pd->mode == T2D)&(!wbz_exist)&(!pd->flag_hdr)) {
	   par->get_string("freecoords", &pd->fn_hor);
	   if (pd->fn_hor) {
		   pd->fp_hor = fopen(pd->fn_hor, "r");
	       if (!pd->fp_hor) {		
   		      ucsl_errorf("ERROR: cannot open input freecoords file: %s!!!\n", &pd->fn_hor);	   
		   }
	       if (!pd->flag_hdr) pd->flag_fc = 1;
	       if (!pd->mype) ucsl_printf("Program will stretch the input horizon.\n");
	   }
	}
	
	
	pd->out_nzt = -1;
	par->get_int("out_nzt", &pd->out_nzt);
	if (pd->out_nzt<0){
		if (!pd->mype) ucsl_printf("Program will estimate the output sample total.\n");
		pd->out_nzt = estimate_nzt(pd->fd_vel,pd->hmap_vel,pd->axis_vel[0],pd->out_dzt,pd->mode);
        if (!pd->mype) ucsl_printf("Output sample total will be automatically derived as %d\n",pd->out_nzt);       		
	}
	
	pd->in_nzt = pd->axis_in[0]->n;
	pd->in_dzt = pd->axis_in[0]->d;
	pd->vel_nzt = pd->axis_vel[0]->n;
	pd->vel_dzt = pd->axis_vel[0]->d;
	pd->nh = pd->axis_in[1]->n;
	pd->nhdr_dat = pd->hmap_dat->get_n_hdr();
	pd->nhdr_vel = pd->hmap_vel->get_n_hdr();
	
	ucsl_abort_if_error();
	
//check sptial interval, sorting order and range	
    int count;
	count = par->get_n_val("subline_seq");
	if (count == 2) {
	   par->get_next_int("subline_seq", &pd->submin);
	   par->get_next_int("subline_seq", &pd->submax);
	   if (pd->submin>pd->submax)
	      ucsl_errorf("ERROR: error on input subline sequence!!!\n");		  
	} else if (count == 0) {
	   pd->submin = pd->axis_vel[MDIM-1]->o;
	   pd->submax = pd->axis_vel[MDIM-1]->e;
	} else {
	   ucsl_errorf("ERROR: illegal subline sequence input!!!\n");		  
	}
	
	count = par->get_n_val("xline_seq");
	if (count == 2) {
	   par->get_next_int("xline_seq", &pd->crsmin);
	   par->get_next_int("xline_seq", &pd->crsmax);
	   if (pd->crsmin>pd->crsmax)
	      ucsl_errorf("ERROR: error on input xline sequence!!!\n");		  
	} else if (count == 0) {
	   pd->crsmin = pd->axis_vel[MDIM-2]->o;
	   pd->crsmax = pd->axis_vel[MDIM-2]->e;
	} else {
	   ucsl_errorf("ERROR: illegal xline sequence input!!!\n");		  
	}
	
	check_cons(pd->axis_in, pd->axis_vel, pd->submin, pd->submax,
	             pd->crsmin, pd->crsmax);
	ucsl_abort_if_error();
	
//Build MPI database
	pd->mpi_range = (int**)mem_alloc2(2*sizeof(int), pd->npe, 1);
	node_traces = (int*)mem_alloc(pd->npe, sizeof(int));
    ncdp = get_ntrace(pd->axis_vel, pd->submin, pd->submax, pd->crsmin, pd->crsmax);
	if (ncdp<pd->npe) ucsl_errorf("ERROR: CDP total is smaller than nodes total, you shall not waste computational resouces!!!\n");
	for (i=0;i<pd->npe;i++) node_traces[i] = (int)(ncdp/pd->npe); 
	for (i=0;i<ncdp%pd->npe;i++) node_traces[i]++;	
	if (!pd->mype) { 
		ucsl_printf("************************************************************\n");
		ucsl_printf("There are %d cdps in total.\n", ncdp);
	}
	pd->mpi_range[0][0] = 0;
	pd->mpi_range[0][1] = node_traces[0]-1;
	for (i=1;i<pd->npe;i++) {
		pd->mpi_range[i][0] = pd->mpi_range[i-1][1]+1;
		pd->mpi_range[i][1] = pd->mpi_range[i][0]+node_traces[i]-1;		
	}
	if (!pd->mype) {
		for (i=0;i<pd->npe;i++) 
			ucsl_printf("Node %d will process cdps from %d to %d.\n", 
					 i,pd->mpi_range[i][0],pd->mpi_range[i][1]);
	}
	if (!pd->mype) ucsl_printf("************************************************************\n");
	
	free(node_traces);
	ucsl_abort_if_error();	

//Open output files
	
//output data
   pd->axis_out = (axis**)mem_alloc(GDIM*sizeof(axis*));
   pd->axis_out[0] = new axis(0.0, pd->out_dzt, pd->out_nzt);
   pd->axis_out[0]->set_label("samp");
   pd->axis_out[1] = new axis(pd->axis_in[1]);
   pd->axis_out[2] = new axis((float)pd->crsmin, (float)pd->crsmax, pd->axis_in[2]->d);
   pd->axis_out[2]->set_label("xline");
   pd->axis_out[3] = new axis((float)pd->submin, (float)pd->submax, pd->axis_in[3]->d);
   pd->axis_out[3]->set_label("subline");
   pd->fd_out = new file_trace(pd->fd_in, pd->fn_out, NULL, NULL, NULL, MPI_COMM_WORLD, 
			                 pd->axis_out, GDIM, _FILE_WRITE, 1);
							 
// vel2 
   pd->axis_vel2 = (axis**)mem_alloc(MDIM*sizeof(axis*));
   pd->axis_vel2[0] = new axis(0.0, pd->out_dzt, pd->out_nzt);
   pd->axis_vel2[1] = new axis(pd->axis_out[2]);
   pd->axis_vel2[2] = new axis(pd->axis_out[3]);
   if (pd->fn_vel2) {      
	  pd->hmap_vel2 = new hdr_map(pd->hmap_vel);
	  if (pd->flag_hdr) pd->hmap_vel2->add_hdr("WBZ", _HDR_FLOAT);
	  if (pd->flag_fc) {
	     if (pd->mode == T2D) pd->hmap_vel2->add_hdr("WBZ", _HDR_FLOAT);
	     if (pd->mode == D2T) pd->hmap_vel2->add_hdr("WBT", _HDR_FLOAT);
	  }
	  pd->nhdr_vel2 = pd->hmap_vel2->get_n_hdr();
	  pd->fd_vel2 = new file_mdvd(pd->fn_vel2, pd->project, "stretch", MPI_COMM_WORLD, pd->hmap_vel2, pd->axis_vel2, MDIM);
   }

// mdl
   if (pd->nmdl>0) {
      pd->nhdr_mdl = (int*)mem_alloc(pd->nmdl, sizeof(int));
	  for (i=0;i<pd->nmdl;i++) {
	     pd->nhdr_mdl[i] = pd->hmap_mdl[i]->get_n_hdr();
		 pd->fd_mdl2[i] = new file_trace(pd->fd_mdl[i], pd->fn_mdl2[i], NULL, NULL, NULL, MPI_COMM_WORLD, 
					            pd->axis_vel2, MDIM, _FILE_WRITE, 1);
		}
   }


//Allocate buffers
	pd->traces_in = (float***)mem_alloc3(pd->fd_in->get_bytes_trace(), pd->nh, TBLOCK, 1);
	pd->dtraces_in = pd->hmap_dat->get_data_ptr(pd->traces_in, pd->nh, TBLOCK);
	   
	pd->traces_out = (float***)mem_alloc3((pd->out_nzt+pd->nhdr_dat)*sizeof(float), pd->nh, TBLOCK, 1);
	pd->dtraces_out = pd->hmap_dat->get_data_ptr(pd->traces_out, pd->nh, TBLOCK);
	      
	pd->traces_vel = (float**)mem_alloc2(pd->fd_vel->get_bytes_trace(), TBLOCK, 1);
	pd->dtraces_vel = pd->hmap_vel->get_data_ptr(pd->traces_vel, TBLOCK);
	
	if (pd->fn_vel2) {
		pd->traces_vel2 = (float**)mem_alloc2((pd->out_nzt+pd->nhdr_vel2)*sizeof(float), TBLOCK, 1);
		pd->dtraces_vel2 = pd->hmap_vel2->get_data_ptr(pd->traces_vel2, TBLOCK);
	}
	
	if (pd->nmdl>0) {
		pd->dtraces_mdl = (float***)mem_alloc(pd->nmdl*sizeof(float**));
		pd->dtraces_mdl2 = (float***)mem_alloc(pd->nmdl*sizeof(float**));
		pd->traces_mdl = (float***)mem_alloc3(pd->fd_mdl[0]->get_bytes_trace(), TBLOCK, pd->nmdl, 1);
		pd->traces_mdl2 = (float***)mem_alloc3((pd->out_nzt+pd->nhdr_mdl[0])*sizeof(float), TBLOCK, pd->nmdl, 1);
		for (i=0;i<pd->nmdl;i++){
			pd->dtraces_mdl[i] = pd->hmap_mdl[i]->get_data_ptr(pd->traces_mdl[i], TBLOCK);
			pd->dtraces_mdl2[i] = pd->hmap_mdl[i]->get_data_ptr(pd->traces_mdl2[i], TBLOCK);
		}		
	}
	
	int nxl,nsl;
	nxl = pd->axis_out[GDIM-2]->n;
	nsl = pd->axis_out[GDIM-1]->n;
	if ((pd->flag_hdr)||(pd->flag_fc)) pd->hor_map = (float*)mem_alloc(nxl*nsl*sizeof(float));	
	pd->str = (stretch_obj**)mem_alloc(pd->nprocs*sizeof(stretch_obj*));
	int vel2_flag = 0;
	if (pd->fn_vel2) vel2_flag=1;
	for (i=0;i<pd->nprocs;i++) {
		pd->str[i] = new stretch_obj(pd->in_dzt, pd->in_nzt, pd->vel_dzt, pd->vel_nzt, 
	            pd->out_dzt, pd->out_nzt, pd->nh, pd->mode, vel2_flag);
	}
	
	ucsl_abort_if_error();
	
//Read freecoords file
	int ixl,isl,cc;
	int xl,sl;
        int xl_tmp,sl_tmp;
	float tz;
	
	ixl=isl=cc=0;
	if (pd->flag_fc) {
	   while (fscanf(pd->fp_hor, "%d %d %f\n", &sl, &xl, &tz) != EOF) {
			if ((pd->axis_out[GDIM - 2]->exist(xl)) && (pd->axis_out[GDIM - 1]->exist(sl))) {				
				if ((pd->mode == T2D) & (tz > 1.0)) tz *= 0.001;
				pd->hor_map[cc] = tz;
				cc++;
				ixl++;
				if (ixl == nxl) {
                    xl_tmp = pd->axis_out[GDIM - 2]->get_val(ixl-1);
                    sl_tmp = pd->axis_out[GDIM - 1]->get_val(isl);
					if ((xl_tmp != xl)||(sl_tmp != sl))
						ucsl_errorf("ERROR: Horizon SL/XL orders or spatial interval does not match with output data range!\n");
					ucsl_abort_if_error();
					ixl = 0;
					isl++;
				}
			}
		}
	    if (cc!=nxl*nsl) ucsl_errorf("ERROR: Horizon total does not match with total CDPs!\n");
		fclose(pd->fp_hor);
	}

	ucsl_abort_if_error();
}


void process(prog_data *pd) {
	int nprocs,iproc,mype,nprocs_tmp;
	int ntraces,begin,end;
	int nh,nmdl;
	int in_nzt,vel_nzt,out_nzt;
	float in_dzt,vel_dzt,out_dzt;
	
	int itr,idx_dat,idx_mdl;
	int nread,nread_tmp;
	int start;
	int i,j,ifile;

	float ***buf_in = NULL, ***buf_out = NULL, **buf_vel = NULL, ***buf_mdl = NULL;
	float **buf_vel2 = NULL,***buf_mdl2 = NULL;
	float *buf_hor = NULL;
	
	nprocs = pd->nprocs;
	mype = pd->mype;
	ntraces = pd->mpi_range[mype][1]-pd->mpi_range[mype][0]+1;
	begin = pd->mpi_range[mype][0];
	end = pd->mpi_range[mype][1];
	nh = pd->nh;
	nmdl = pd->nmdl;
	in_nzt = pd->in_nzt;
	out_nzt = pd->out_nzt;
	vel_nzt = pd->vel_nzt;
	in_dzt = pd->in_dzt;
	out_dzt = pd->out_dzt;
	vel_dzt = pd->vel_dzt;
	
	//TODO: shift pointers instead of memcpy? 
	buf_in = (float***) mem_alloc3(in_nzt, nh, nprocs, sizeof(float));
	buf_out = (float***) mem_alloc3(out_nzt, nh, nprocs, sizeof(float));
	buf_vel = (float**) mem_alloc2(vel_nzt, nprocs, sizeof(float));
	
	float **ptr_vel2;
    ptr_vel2 = (float**) mem_alloc(pd->nprocs*sizeof(float*));
    for (i = 0; i < pd->nprocs; i++) ptr_vel2[i] = NULL;		
	
	if (pd->fn_vel2) buf_vel2 = (float**) mem_alloc2(out_nzt, nprocs, sizeof(float));
	if (nmdl>0) {
		buf_mdl = (float***) mem_alloc3(vel_nzt, nmdl, nprocs, sizeof(float));
		buf_mdl2 = (float***) mem_alloc3(out_nzt, nmdl, nprocs, sizeof(float));
	}
		
	if ((pd->flag_hdr)||(pd->flag_fc)) buf_hor = (float*) mem_alloc(nprocs*sizeof(float));
	
		
	nread = TBLOCK;
	for (itr = begin; itr <= end; itr += TBLOCK) {
		if ((itr-begin) + nread > ntraces) nread = ntraces - (itr-begin);
		
		for (i=0;i<nread;i++) {
		   idx_dat = local_index(itr+i, pd->axis_out, pd->axis_in, GDIM);
		   idx_mdl = local_index(itr+i, pd->axis_vel2, pd->axis_vel, MDIM);
		   pd->fd_in->read_traces(pd->traces_in[i][0], nh, idx_dat * nh);
		   pd->fd_vel->read_traces(pd->traces_vel[i], 1, idx_mdl);
		   if (nmdl>0){			
		      for (j=0;j<nmdl;j++){
			     pd->fd_mdl[j]->read_traces(pd->traces_mdl[j][i], 1, idx_mdl);
			  }			
		   }
		}        
		
		nprocs_tmp = nprocs;
		nread_tmp = nread;
		start = 0;
		while (nread_tmp > 0) {
			if (nread_tmp < nprocs_tmp) nprocs_tmp = nread_tmp;	
			
			for (iproc = 0; iproc < nprocs_tmp; iproc++) {
				for (i = 0; i < nh; i++) {
					memcpy(buf_in[iproc][i], pd->dtraces_in[start + iproc][i], in_nzt * sizeof(float));					
				}
				
				memcpy(buf_vel[iproc], pd->dtraces_vel[start + iproc], vel_nzt*sizeof(float)); 
				if (pd->nmdl){
					for (i=0; i<nmdl; i++) {
						memcpy(buf_mdl[iproc][i], pd->dtraces_mdl[i][start + iproc], vel_nzt*sizeof(float));
					}
				}
				if (pd->flag_fc) buf_hor[iproc] = pd->hor_map[itr+start + iproc];
				if (pd->flag_hdr) {
					buf_hor[iproc] = pd->hmap_vel->get_hdr_float(pd->traces_vel[start + iproc], "WBT");
					if (buf_hor[iproc]>1.0) buf_hor[iproc] /= 1000.0;					
				}				
			}
			
			if (pd->fn_vel2) {
				for (i = 0; i < nprocs_tmp; i++) ptr_vel2[i] = buf_vel2[i];
			} else {				
				for (i = 0; i < nprocs_tmp; i++) ptr_vel2[i] = NULL;
			}
			
			#pragma omp parallel for default(shared) private (iproc)			
			for (iproc = 0; iproc < nprocs_tmp; iproc++) {				
				if (pd->nmdl>0) stretch(pd->str[iproc], buf_in[iproc], buf_vel[iproc], buf_out[iproc], ptr_vel2[iproc], buf_mdl[iproc], buf_mdl2[iproc], nmdl, buf_hor, iproc);
				else stretch(pd->str[iproc], buf_in[iproc], buf_vel[iproc], buf_out[iproc], ptr_vel2[iproc], NULL, NULL, nmdl, buf_hor, iproc);
			}

			for (iproc = 0; iproc < nprocs_tmp; iproc++) {
				for (i = 0; i < nh; i++) {
					memcpy(pd->dtraces_out[start + iproc][i], buf_out[iproc][i], out_nzt * sizeof(float));
					pd->hmap_dat->copy_hdrs(pd->traces_out[start + iproc][i], pd->traces_in[start + iproc][i]);
				}
				for (i=0; i<nmdl; i++) {
					memcpy(pd->dtraces_mdl2[i][start + iproc], buf_mdl2[iproc][i], out_nzt * sizeof(float));
					pd->hmap_mdl[i]->copy_hdrs(pd->traces_mdl2[i][start + iproc], pd->traces_mdl[i][start + iproc]);
				}
				if (pd->fn_vel2) {
					memcpy(pd->dtraces_vel2[start + iproc], buf_vel2[iproc], out_nzt * sizeof(float));
					pd->hmap_vel->copy_hdrs(pd->traces_vel2[start + iproc], pd->traces_vel[start + iproc]);
					if ((pd->flag_hdr)||(pd->flag_fc)) //TODO need verification
						pd->traces_vel2[start + iproc][pd->hmap_vel2->get_hdr_map_index("WBZ")] = buf_hor[iproc];
				}				
			}

			nread_tmp -= nprocs_tmp;
			start = start + nprocs_tmp;
		}

				
		pd->fd_out->write_traces(pd->traces_out[0][0], nread * nh, itr * nh);		
		if (pd->fn_vel2) pd->fd_vel2->write_traces(pd->traces_vel2[0], nread, itr);
		for (i=0;i<nmdl;i++) pd->fd_mdl2[i]->write_traces(pd->traces_mdl2[i][0], nread, itr);		
	}


//TODO: double check the mem_free, make sure it can free the multi-dimension buffer
	if (buf_in) mem_free3((void****)&buf_in);
	if (buf_out) mem_free3((void****)&buf_out);
	if (buf_vel) mem_free2((void***)&buf_vel);
	if (buf_vel2) mem_free2((void***)&buf_vel2);
	if (buf_mdl2) mem_free3((void****)&buf_mdl2);
	if (buf_mdl) mem_free3((void****)&buf_mdl);
	if (buf_hor) mem_free((void**)&buf_hor);
	if (ptr_vel2) free(ptr_vel2);
	
	ucsl_printf("Stretch process finishes sucessfully in node %d\n", mype);
}

void stretch(stretch_obj *str, float **in, float *vel, float **out, float *vel2, float **mdl, float **mdl2, int nmdl, float *hor, int ind){
	float val_in, val_out;
	int i;
	
	str->process(in, vel, out, vel2);
	if (hor) {
		val_in = hor[ind];
		val_out = str->hor_stretch(val_in, vel);
		hor[ind] = val_out;
	}
	if (nmdl>0){
		for (i=0; i<nmdl; i++) 
			str->process_mdls(mdl[i], mdl2[i], vel);
	}	
}


void open_file(char *fn, file_trace **fd, hdr_map **hmap, axis ***ax, int dims){
	int dims_tmp;
	int idim;
	
	fd[0] = new file_trace(fn, NULL, NULL, MPI_COMM_WORLD, _FILE_READ);
	hmap[0] = fd[0]->get_hdr_map();
	dims_tmp = fd[0]->get_dims();
	if (dims_tmp!=dims)  ucsl_errorf("ERROR: Dimensions of input file must be %d!\n", dims);
	
	ax[0] = (axis**)mem_alloc(dims*sizeof(axis*));
	for (idim=0; idim<dims; idim++) {
		ax[0][idim] = fd[0]->get_axis(idim);
	}	
}

int estimate_nzt(file_trace *fd, hdr_map *hmap, axis *axz, float dzt_out, int mode) {
	float *trace,*dtrace;
	int nhdr,nsamp;
	float vava;
	int iz;
	float dzt_in;
	int nsamp_out;
	float scale = 1.15;
	int res = 50;
	
	nhdr = hmap->get_n_hdr();
	nsamp = axz->n;
	dzt_in = axz->d;
	trace = (float*)mem_alloc(fd->get_bytes_trace());
	dtrace = &trace[nhdr];
	
	fd->read_traces(trace, 1, 0);
	
	vava = 0.0;
	for (iz=0;iz<nsamp;iz++) vava += dtrace[iz];
	vava /= (float)nsamp;
	
	if (mode == D2T) {
		if (dzt_out>1.0) dzt_out /= 1000.0;
		nsamp_out = 2.0*(int)((((nsamp*dzt_in)/vava)/dzt_out)*scale);
	}
	if (mode == T2D) {
		if (dzt_in>1.0) dzt_in /= 1000.0;
		nsamp_out = (int)(scale*vava*nsamp*dzt_in*0.5/dzt_out);
	}
	
	nsamp_out = ((int)(nsamp_out/res)+1)*res;
	
	free(trace);
	return nsamp_out;
}

void check_cons(axis **ax_d, axis **ax_v, int submin, int submax, int crsmin, int crsmax) {
   int iax;
   char *label;

   //check axis order, must be xline, subline
   label = ax_d[GDIM-2]->label;
   if (strcasecmp(label,"xline")&&strcasecmp(label,"XLINE"))
      ucsl_errorf("ERROR: The second axis of data must be xline!!!\n");
   label = ax_d[GDIM-1]->label;
   if (strcasecmp(label,"subline")&&strcasecmp(label,"SUBLINE"))
      ucsl_errorf("ERROR: The first axis of data must be subline!!!\n");
	  
   label = ax_v[MDIM-2]->label;  
   if (strcasecmp(label,"xline")&&strcasecmp(label,"XLINE"))
      ucsl_errorf("ERROR: The second axis of model(s) must be xline!!!\n");
   label = ax_v[MDIM-1]->label;
   if (strcasecmp(label,"subline")&&strcasecmp(label,"SUBLINE"))
      ucsl_errorf("ERROR: The first axis of model(s) must be subline!!!\n");
	  
   // check spatial interval
   if (ax_d[GDIM-2]->d != ax_v[MDIM-2]->d)
      ucsl_errorf("ERROR: The xline spatial interval are not the same between model(s) and data!!!\n");
   if (ax_d[GDIM-1]->d != ax_v[MDIM-1]->d)
      ucsl_errorf("ERROR: The subline spatial interval are not the same between model(s) and data!!!\n"); 
	  
   //check range
   if ((ax_d[GDIM-2]->o > crsmin) & (ax_d[GDIM-2]->e < crsmax))
      ucsl_errorf("ERROR: The xline range of data is smaller than what user specified!!!\n");
   if ((ax_v[MDIM-2]->o > crsmin) & (ax_v[MDIM-2]->e < crsmax))
      ucsl_errorf("ERROR: The xline range of model(s) is smaller than what user specified!!!\n");	  
   
   if ((ax_d[GDIM-1]->o > submin) & (ax_d[GDIM-1]->e < submax))
      ucsl_errorf("ERROR: The subline range of data is smaller than what user specified!!!\n");
   if ((ax_v[MDIM-1]->o > submin) & (ax_v[MDIM-1]->e < submax))
      ucsl_errorf("ERROR: The subline range of model(s) is smaller than what user specified!!!\n");	  
}

int get_ntrace(axis **ax_v, int submin, int submax, int crsmin, int crsmax){
   int nsub,ncrs;
   int ncdp;
   
   nsub = 1+(submax-submin)/(ax_v[MDIM-1]->d);
   ncrs = 1+(crsmax-crsmin)/(ax_v[MDIM-2]->d);
   
   ncdp = nsub*ncrs;
   
   return ncdp;
}

int local_index(int itr, axis **ax_out, axis **ax_in, int ndim) {
   int nsub_out, ncrs_out;
   int nsub_in, ncrs_in;
   int idx_sub_out, idx_crs_out;
   int idx_sub_in, idx_crs_in;
   float sub,crs;
   int itr_in;
   
   nsub_out = ax_out[ndim-1]->n;
   ncrs_out = ax_out[ndim-2]->n;
   nsub_in = ax_in[ndim-1]->n;
   ncrs_in = ax_in[ndim-2]->n;
   
   idx_sub_out = itr/ncrs_out;
   idx_crs_out = itr%ncrs_out;
   
   sub = ax_out[ndim-1]->get_val(idx_sub_out);
   crs = ax_out[ndim-2]->get_val(idx_crs_out);
   
   idx_sub_in = ax_in[ndim-1]->get_index(sub);
   idx_crs_in = ax_in[ndim-2]->get_index(crs);
   
   itr_in = idx_sub_in*ncrs_in+idx_crs_in;

   return itr_in;
}

void clean_up(prog_data *pd) {
	int i;
	
	delete pd->fd_in;
	delete pd->fd_vel;
	if (pd->fn_vel2) delete pd->fd_vel2;
	if (pd->nmdl>0) {
		for (i=0;i<pd->nmdl;i++){
			delete pd->fd_mdl[i];
			delete pd->fd_mdl2[i];
		}			
	}
	ucsl_printf("Clean up sucessfull\n");	
}