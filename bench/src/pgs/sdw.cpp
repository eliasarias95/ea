#include <ucsl.h>
#include <sdw_obj.h>
#include <cart_volume_io.h>
#include <ctime>
#ifdef _USE_OMP
#include <omp.h>
#endif
#include <mpi.h>

char *usage[] = {
  "************************************************************",
  " ufilter par=[key=val]",
  " Parameters:",
  " in=",
  " out=",
  " q:=",
  " type=(alpha-trim,median,boxcar,gaussian,kuwahara,std_dev,set_mask,curv)",
  " h1=2",
  " h2=2",
  " h3=2",
  " slowness=0",
  " niter=1",
  " nthread=number of cores",
  " alpha=(if type == alpha-trim)",
  " mask_val=1",
  " mask_option=eq (eq,lt,le,gt,ge)",
  "************************************************************",
  NULL};

typedef struct {
  char *fn_in;
  char *fn_out;
  char *fn_qc;
  char *mask_option;
  file_trace *fd_in;
  axis ** axes_in;
  char *filter_type;
  float mask_val;
  float binil, binxl;
  int pad3;
  int h1;
  int h2;
  int h3;
  int niter; 
  int naxes;
  int slowness;
  float alpha;
  float threshold;
  int mype;
  int nthread;
  int verbose;
} prog_data;

int init_cart_volume(cart_volume ** vol_out, axis **axes_in, 
    size_t const& naxes, int h1, int h2, int h3, MPI_Comm comm);
void do_filter3d(prog_data *pd, parlist *par); 
void do_filter2d(prog_data *pd, parlist *par); 

int main (int argc, char *argv[]) {
  double global_time_start, global_time_end;
  prog_data *pd;
  parlist *par;

  par = ucsl_initialize(usage, argv, argc, 1);
  if(!par->get_status()) return _ERROR_PAR;
  global_time_start = MPI_Wtime();

  pd = (prog_data*)mem_alloc(sizeof(prog_data));
  pd->mype = ucsl_mype();

  par->get_string_req("in", &(pd->fn_in));
  par->get_string_req("out", &(pd->fn_out));
  pd->fn_qc = NULL;
  par->get_string("qc", &(pd->fn_qc));
  par->get_int_def("h1", &(pd->h1), 2);
  par->get_int_def("h2", &(pd->h2), 2);
  par->get_int_def("h3", &(pd->h3), 2);
  //add checks for positive h values and also check pad against size of
  //total model
  par->get_int_def("pad_subline", &(pd->pad3),pd->h3);
  //par->get_string_req("type", &(pd->filter_type));
  par->get_int_def("nthread", &(pd->nthread), ucsl_get_n_cores());
  par->get_int_def("verbose", &(pd->verbose), 1); 
  par->get_string_def("mask_option", &(pd->mask_option), "eq");
  par->get_float_def("mask_val", &(pd->mask_val), 1.f);
  if(pd->verbose && !pd->mype) {
    ucsl_print("**********************************************************\n");
    ucsl_printf("ufilter version %s\n", _UCSL_VERSION);
    par->print(_UCSL_LOG);
    ucsl_printf("using %d threads\n", pd->nthread);
    ucsl_print("**********************************************************\n");
  }
#ifdef _USE_OMP
  omp_set_num_threads(pd->nthread);
#endif 
  ucsl_abort_if_error();
  pd->fd_in = new file_trace(pd->fn_in, _FILE_READ | _FILE_WRITE);
  ucsl_abort_if_error();
  pd->naxes = pd->fd_in->get_axes(&(pd->axes_in));

  if(pd->naxes == 3 && pd->axes_in[2]->n > 1) do_filter3d(pd, par);
  else if((pd->naxes == 2 || (pd->naxes == 3 && pd->axes_in[2]->n == 1)) && 
      !pd->mype) do_filter2d(pd, par);

  delete pd->fd_in;
  delete par;

  global_time_end = MPI_Wtime();
  double const global_elapsed_time = global_time_end - global_time_start;
  ucsl_printf("(PE: %d) Elapsed time: %.3f\n", pd->mype, global_elapsed_time);

  ucsl_finalize();

  return EXIT_SUCCESS;
}

void do_filter2d(prog_data *pd, parlist *par) {
  file_trace *fd_out;
  hdr_map *hdr;
  float **tbuf, **dbuf;
  float **data, **sdata, **slopes;
  axis *ax1, *axs1, *axs2;
  int i, j;
  int res;
  float binxl;

  int bytes_trace = pd->fd_in->get_bytes_trace(); // 96
  int n1 = pd->axes_in[0]->n; // 21
  int n2 = pd->axes_in[1]->n; // 15
  float d1 = pd->axes_in[0]->d;
  float d2 = pd->axes_in[1]->d;
  float o1 = pd->axes_in[0]->o;
  float o2 = pd->axes_in[1]->o;

  axs1 = new axis(0.0f,1.0f,n1); // axis for shifts in 1st dimension
  axs2 = new axis(0.0f,1.0f,n2); // axis for shifts in 2nd dimension
  ax1  = new axis(o1,d1,n1,0,0); // axis for data in 1st dimension

  tbuf   = (float**)mem_alloc2(bytes_trace,n2,1);
  data   = (float**)mem_alloc2(n1,n2,sizeof(float));
  sdata  = (float**)mem_alloc2(n1,n2,sizeof(float));
  slopes = (float**)mem_alloc2(n1,n2,sizeof(float));

  hdr = pd->fd_in->get_hdr_map();
  if(hdr) dbuf = hdr->get_data_ptr(tbuf,n2);
  else dbuf = tbuf;

  pd->fd_in->read_traces(tbuf[0],n2,0);

  for(j=0; j<n2; ++j) {
    for(i=0; i<n1; ++i) {
      data[j][i] = dbuf[j][i];
    }
  }

  memcpy(sdata[0],data[0],n1*sizeof(float));
  memcpy(sdata[n2-1],data[n2-2],n1*sizeof(float));
  for (int i2=1; i2<n2-1; ++i2) {
    memcpy(sdata[i2],data[i2-1],n1*sizeof(float));
  }

  //Smooth dynamic warping routine
  int k = 10;
  double pmax = 2.0;
  double r1 = 0.1, r2 = 0.2;
  // use axes for shifts when creating sdw_obj
  sdw_obj *sdw = new sdw_obj(k,-pmax,pmax,axs1,axs2);
  sdw->setStrainLimits(-r1,r1,-r2,r2);
  sdw->setSmoothness((double)pd->h1,(double)pd->h2);
  // use axes for data when finding shifts
  sdw->findShifts(ax1,data,ax1,sdata,slopes);

  float sum = 0.0f;
  for(j=0; j<n2; ++j) {
    for(i=0; i<n1; ++i) {
      sum += slopes[j][i];
      dbuf[j][i] = slopes[j][i];
    }
  }
  std::cout << "mean slope= " << sum/(8.0f*20.0f) << "\n";
  fd_out = new file_trace(pd->fd_in,pd->fn_out,NULL);
  fd_out->write_traces(tbuf[0],n2,0);
  delete fd_out;

  mem_free2((void***)&tbuf);
}


void do_filter3d(prog_data *pd, parlist *par) {
  /*filter_base *filter;
  cart_volume * vol(0);
  char * hdr_buffer(0);
  int res;
  float binxl, binil;

  res = init_cart_volume(&vol, pd->axes_in, pd->naxes, pd->h1, pd->h2, pd->h3, MPI_COMM_WORLD);
  if(res) ucsl_errorf("(me: %d) Failed to initialize cart_volume\n", pd->mype);
  ucsl_abort_if_error();

  res = fill_cart_volume(vol, &hdr_buffer, pd->fd_in, MPI_COMM_WORLD);
  if(res) ucsl_errorf("(me: %d) Failed to read input\n", pd->mype);
  ucsl_abort_if_error();

  if(pd->slowness) vol->invert();
  ucsl_abort_if_error();

  if(strcmp(pd->filter_type, "boxcar") == 0)
    filter = new filter_boxcar(pd->h1, pd->h2, pd->h3);
  else if(strcmp(pd->filter_type, "gaussian") == 0)
    filter = new filter_gaussian(pd->h1, pd->h2, pd->h3);
  else if(strcmp(pd->filter_type, "kuwahara") == 0) 
    filter = new filter_kuwahara(pd->h1, pd->h2, pd->h3, pd->nthread);
  else if(strcmp(pd->filter_type, "median") == 0)
    filter = new filter_median(pd->h1, pd->h2, pd->h3, pd->nthread);
  else if(strcmp(pd->filter_type, "curv") == 0) {
    par->get_float_req("binxl", &binxl);
    par->get_float_req("binil", &binil);
    filter = new filter_curv(vol->ax1->d, vol->ax2->d*binxl, vol->ax3->d*binil, 1.f);
  }
  else if(strcmp(pd->filter_type, "alpha-trim") == 0) {
    if(pd->alpha == ALPHA_UNDEFINED) 
      ucsl_errorf("ERROR: alpha must be defined for alpha-trim filter\n");
    filter = new filter_alphatrim(pd->h1, pd->h2, pd->h3, pd->alpha, pd->nthread);
  }
  else if(strcmp(pd->filter_type, "set_mask") == 0)
    filter = new filter_set_mask(pd->mask_option, pd->mask_val);
  else ucsl_errorf("type = '%s' is not a supported filter.\n", pd->filter_type);
  ucsl_abort_if_error();

  if(!pd->mype) filter->print(_UCSL_LOG);
  filter_loop loop;
  res = loop.run(filter, pd->niter, vol);
  if(res) ucsl_errorf("(PE: %d) Failed to apply filter\n", pd->mype);
  ucsl_abort_if_error();
  delete filter;

  file_trace * fd_out = new file_trace(pd->fd_in, pd->fn_out, MPI_COMM_WORLD);

  if(pd->slowness) vol->invert();
  ucsl_abort_if_error();

  res = dump_cart_volume(fd_out, vol, hdr_buffer, MPI_COMM_WORLD);
  if(res) ucsl_errorf("(PE: %d) Failed to write filtered output!\n", pd->mype);
  delete fd_out;

  delete vol;
  if(hdr_buffer) mem_free((void **)&hdr_buffer);
  */
}

int init_cart_volume(cart_volume ** vol_out, axis **axes_in, size_t const& naxes, 
    int h1, int h2, int h3, MPI_Comm comm) {
  axis **axes_filter = new axis*[naxes];
  int filter_hl_list[3];
  filter_hl_list[0] = h1;
  filter_hl_list[1] = h2;
  filter_hl_list[2] = h3;

  for(size_t i(0); i < naxes; ++i)
  {
    axes_filter[i] = new axis(axes_in[i]->o,
        axes_in[i]->d,
        axes_in[i]->n,
        filter_hl_list[i],
        0);

    axes_filter[i]->set_label(axes_in[i]->label);
  }

  int const mype = ucsl_mype();
  int const npe = ucsl_npe();

  int const estimated_distribution =
    cart_volume::compute_n(axes_filter[2]->n, mype, npe);

  if(estimated_distribution <= filter_hl_list[2])
  {
    ucsl_errorf("The ghost region is too small.\n");
    return 1;
  }

  (*vol_out) = new cart_volume(axes_filter[0], axes_filter[1],
      axes_filter[2], comm);
  return 0;
}
