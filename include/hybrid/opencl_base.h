#ifndef INCLUDE_OPENCL_BASE_
#define INCLUDE_OPENCL_BASE_
#include "common.h"
#ifdef USE_OPENCL

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <stdarg.h>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

typedef struct OpenCLProfile
{
  //Profiling
  cl_event event;
  double total_time;
  double kern_time;
} OpenCLProfile;

typedef struct OpenCL
{
  cl_uint device; //Selected device
  cl_uint pf_id;  //Selected platform

  cl_uint n_platforms;
  cl_uint n_devices;
  cl_platform_id platforms[32];
  cl_device_id devices[32];

  cl_program prog;
  cl_context context;
  cl_command_queue queue;

  int nvidia;
  int problem_size;
  int problem_sections;
  int offset;

  size_t global_size;
  size_t maxwgsize;
  size_t nthreads;
  size_t ngroups;
  size_t nblocks;

  //Build log
  char log[65536];

  //Default profiler
  OpenCLProfile profile;

  //Constant data
  char* constants;
} OpenCL;

void clCheck(int error, const char *msg);
int get_thread_count(int n, int max);
int get_global_count(int n, int threads);
char* cl_readfile(const char* filename);

OpenCL* cl_new(unsigned int N, unsigned int sections, unsigned int threadmax, unsigned int device_id);
void cl_size(OpenCL* cl, unsigned int N, unsigned int sections, unsigned int threadmax, unsigned int device_id);
void cl_build(OpenCL* cl, const char* srcptr[], const char* options, unsigned int count);
void cl_checklocalmem(OpenCL* cl, size_t required);
cl_mem cl_new_buffer(OpenCL* cl, cl_mem_flags flags, size_t size);
cl_kernel cl_get_kernel(OpenCL* cl, const char* name, int count, ...);
void cl_profile(OpenCL* cl, OpenCLProfile* profile);
void cl_delete(OpenCL* cl);

void cl_print_platform_info(OpenCL* cl, int i);
void cl_print_device_info(OpenCL* cl, int i);

void cl_define_integer_constant(OpenCL* cl, const char* name, int value);
void cl_define_real_constant(OpenCL* cl, const char* name, Real value);
void cl_define_constants(OpenCL* cl, const char* format, ...);
void cl_append_constant(OpenCL* cl, const char* str);

//Macros to make calls for enqueueing buffers and running kernels a bit cleaner
#define cl_set_arg(cl, ker, pos, arg) \
  clCheck(clSetKernelArg(ker, pos, sizeof(arg), &arg), "clSetKernelArg");

#define cl_enqueue_kernel(cl, ker, eve) \
  clCheck(clEnqueueNDRangeKernel(cl->queue, ker, 1, NULL, &cl->global_size, &cl->nthreads, 0, NULL, eve), "clEnqueueNDRangeKernel")

#define cl_enqueue_write(cl, buf, siz, ptr) \
  clCheck(clEnqueueWriteBuffer(cl->queue, buf, CL_TRUE, 0, siz, ptr, 0, NULL, NULL), "clEnqueueWriteBuffer")

#define cl_enqueue_read(cl, buf, siz, ptr) \
  clCheck(clEnqueueReadBuffer(cl->queue, buf, CL_FALSE, 0, siz, ptr, 0, NULL, NULL), "clEnqueueReadBuffer")

#define cl_run(cl) \
  clCheck(clFinish(cl->queue), "clFinish")

#define cl_run_profile(cl) \
  clCheck(clFinish(cl->queue), "clFinish"); cl_profile(cl, cl->pr)

#endif //USE_OPENCL
#endif // INCLUDE_OPENCL_BASE_
