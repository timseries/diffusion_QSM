#include "hybrid/opencl_base.h"

void clCheck(int error, const char *msg)
{
   if (error != CL_SUCCESS)
   {
      fprintf(stderr, "%s, error code %d\n", msg, error);
      abort();
   }
}

//Calculate number of threads to use based on total problem size and max thread size
int get_thread_count(int n, int max)
{
  int i;
  for (i=max; i>0; i--)
  {
    int rem = n % i;
    if (rem == 0)
      return i;
  }
  return 0;
}

//Calculates global work size given problem size and thread count
int get_global_count(int n, int threads)
{
  return threads * ceil((Real)n/(Real)threads);
}

char* cl_readfile(const char* filename)
{
  //Read kernel
  FILE* fp;
  size_t srcsize;
  char* kernel_src = NULL;
  fp = fopen(filename, "r");
  if (!fp) return NULL;
  fseek(fp, 0, SEEK_END);
  srcsize = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  kernel_src = (char*)malloc(srcsize+1);
  kernel_src[srcsize] = 0;
  if (srcsize != fread(kernel_src, sizeof(char), srcsize, fp)) 
  { 
     free(kernel_src);
     return NULL;
  }
  fclose(fp);
  return kernel_src;
}

void error_callback(const char *errinfo, const void *private_info, size_t cb, void *user_data)
{
  printf("OpenCL Error: %s\n", errinfo);
}

OpenCL* cl_new(unsigned int N, unsigned int sections, unsigned int threadmax, unsigned int device_id)
{
  //device, kernel, context, queue, size, sections, global_size, maxwgsize, nthreads, ngroups, in, out, event, kern_time, wall_time, clk
  if (threadmax == 0) threadmax = 64; //Default to 64 threads per workgroup
  if (threadmax > N/sections) threadmax = N/sections;
  OpenCL* cl = (OpenCL*)malloc(sizeof(OpenCL));
  memset(cl, 0, sizeof(OpenCL));
  cl->device = device_id;
  cl->problem_size = N;

  char buffer[2048];
  clGetPlatformIDs(32, cl->platforms, &cl->n_platforms);
  cl->pf_id = 0;  //Selected platform, using first available for now
  clGetPlatformInfo(cl->platforms[cl->pf_id], CL_PLATFORM_NAME, 2048, buffer, NULL);
  if (device_id==0 && cl->n_platforms > 0)
    fprintf(stderr, "%d OpenCL platform(s) found, using [%d] %s ===\n", cl->n_platforms, cl->pf_id, buffer);
  cl->nvidia = 0;
  if (strstr(buffer, "NVIDIA")) cl->nvidia = 1; 
  //cl_print_platform_info(cl, cl->pf_id);
  if (cl->n_platforms == 0) 
  {
    fprintf(stderr, "[PROC %d] No OpenCL platform(s) found!\n", device_id);
    return NULL;   //No platforms found
  }

  /*if (device_id==0)  //DEBUG: Show all platform specs on root proc
  {
    int ii;
    for (ii=0; ii<cl->n_platforms; ii++)
      cl_print_platform_info(cl, ii);
  }*/

  cl->n_devices = 0;
  clGetDeviceIDs(cl->platforms[cl->pf_id], CL_DEVICE_TYPE_ALL, 32, cl->devices, &cl->n_devices);

  //Ensure selected device within range
  if (cl->device >= cl->n_devices)
    cl->device %= cl->n_devices;

  clGetDeviceInfo(cl->devices[cl->device], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
  fprintf(stderr, "[PROC %d of %d] === %d OpenCL device(s) found on platform %d: (Device %d selected: %s)\n", device_id, sections, cl->n_devices, cl->pf_id, cl->device, buffer);
  /*if (device_id==0)  //DEBUG: Show all selected device specs on root proc
    cl_print_device_info(cl, cl->device);
    */

  //Save some useful device info
  clGetDeviceInfo(cl->devices[cl->device], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl->maxwgsize), &cl->maxwgsize, NULL);

  cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cl->platforms[cl->pf_id], 0};
  // Note that nVidia's OpenCL requires the platform property

#if 0
    /* Calculate number of threads to use - Note, if done here can't check maxwgsize */
    if (threadmax == 0) threadmax = 60;
    cl->nthreads = get_thread_count(N, threadmax); //cl_ maxwgsize);
    cl->ngroups = N / cl->nthreads;
    printf("Calculated thread count %d for total work size %d, (%d groups) max work group size %d\n", 
          cl->nthreads, cl->global_size, cl->ngroups, cl->maxwgsize);
#endif

printf("N %d sections %d threadmax %d device %d\n", N, sections, threadmax, device_id);
  cl_size(cl, N, sections, threadmax, device_id);

  cl_int error;
  cl->context = clCreateContext(properties, 1, &cl->devices[cl->device], &error_callback, NULL, &error);
  clCheck(error, "clCreateContext");
  cl->queue = clCreateCommandQueue(cl->context, cl->devices[cl->device], CL_QUEUE_PROFILING_ENABLE, &error);
  clCheck(error, "clCreateCommandQueue");

  return cl;
}

void cl_size(OpenCL* cl, unsigned int N, unsigned int sections, unsigned int threadmax, unsigned int device_id)
{
  //Calculates total work size to issue based on provided thread count and number of section divisions
  bool print = true;
  if (!sections) 
  {
    sections = cl->problem_sections;
    print = false;
  }
  cl->global_size = get_global_count(N/sections, threadmax);
  cl->offset = cl->global_size * device_id;
  if (sections > 1)
  {
    if (device_id+1 == sections)
    {
      //Adjust global count for the last device
      int remain = N-(sections-1)*cl->global_size;
      if (print) fprintf(stderr, "Last proc, original global_size %lu, remaining items %d\n", cl->global_size, remain);
      cl->global_size = get_global_count(remain, threadmax);
    }
  }

  cl->problem_sections = sections;
  cl->nthreads = threadmax;
  cl->ngroups = cl->global_size / cl->nthreads;
  cl->nblocks = ceil((N/(Real)cl->nthreads));
  if (print) fprintf(stderr, "[%d] Calculated global work size %lu for thread count %lu problem size %d, max work group size %lu, total blocks %lu\n", device_id, cl->global_size, cl->nthreads, N/sections, cl->maxwgsize, cl->nblocks);
}

// Build source
void cl_build(OpenCL* cl, const char* srcptr[], const char* options, unsigned int count)
{
  /*
  int k;
  for (k=0; k<count; k++)
  {
    if (srcptr[k]) printf("---------------------\n%s", srcptr[k]);
  }
  getchar();*/

  // Submit the source code of the kernel to OpenCL
  cl_int error;
  cl->prog=clCreateProgramWithSource(cl->context, count, srcptr, NULL, &error);
  clCheck(error, "Unable to create Program Object");

  //Set options and defines
  char* f_options = (char*)malloc(strlen(options) + 512);
  //User options + standard defines
  const char* opt_format = "%s -DN=%d -Dnthreads=%d -Dnblocks=%d -cl-fast-relaxed-math";
  //const char* opt_format = "%s -DN=%d -Dnthreads=%d -cl-fast-relaxed-math";
  sprintf(f_options, opt_format, options, cl->problem_size, cl->nthreads, cl->nblocks);
#ifdef SINGLE_PRECISION
  //Pass on the single precision switch
  strcat(f_options, " -DSINGLE_PRECISION");
#endif
  if (cl->nvidia)
    strcat(f_options, " -cl-nv-verbose");   //-cl-nv-maxrregcount=32
  //else
  //  strcat(f_options, " -g");

  // Compile (after this we could extract the compiled version)
  fprintf(stderr, "Compiling: %s\n", f_options);
  error = clBuildProgram(cl->prog, 1, &cl->devices[cl->device], f_options, NULL, NULL);
  free(f_options);
  //Always get the build log for compile info from cl-nv-verbose
  //if (error != CL_SUCCESS)
  {
    size_t length;
    clCheck(clGetProgramBuildInfo(cl->prog, cl->devices[cl->device], CL_PROGRAM_BUILD_STATUS, sizeof(cl->log), cl->log, &length), "Failed to get build log");
    int status = *(int*)cl->log;
    if (status) fprintf(stderr, "\nclGetProgramBuildInfo returned cl_build_status: %d\n", status);
    clCheck(clGetProgramBuildInfo(cl->prog, cl->devices[cl->device], CL_PROGRAM_BUILD_LOG, sizeof(cl->log), cl->log, &length), "Error getting program build info. Buffer may be too small.");
    if (length < sizeof(cl->log))
      cl->log[length]='\0';
  }
  if (error != CL_SUCCESS)
  {
    fprintf(stderr, "%s\n", cl->log);
    clCheck(error, "clBuildProgram");
    exit(error*100-2);
  }

  //Free memory
  int i;
  for (i=0; i<count; i++)
  {
    if (srcptr[i]) free((char*)srcptr[i]);
    //If constants passed as a source string, clear pointer
    if (srcptr[i] == cl->constants) cl->constants = NULL;
  }
}

void cl_checklocalmem(OpenCL* cl, size_t required)
{
  //Local memory check: 
  cl_ulong local_mem_size;
  clGetDeviceInfo(cl->devices[cl->device], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem_size), &local_mem_size, NULL);
  fprintf(stderr, "%lu bytes of local memory required, %lu bytes available\n", required, local_mem_size);
  if (local_mem_size > 0 && local_mem_size < required)
  {
    fprintf(stderr, "Insufficient local memory! aborting\n");
    exit(1);
  }
}

cl_mem cl_new_buffer(OpenCL* cl, cl_mem_flags flags, size_t size)
{
  cl_int error;
  cl_mem buffer = clCreateBuffer(cl->context, flags, size, NULL, &error);
  clCheck(error, "clCreateBuffer");
  return buffer;
}

cl_kernel cl_get_kernel(OpenCL* cl, const char* name, int count, ...)
{
  int i;
  va_list ap;
  cl_int error;
  cl_kernel kernel = clCreateKernel(cl->prog, name, &error);
  clCheck(error, "clCreateKernel");

  // Set the kernel arguments
  va_start(ap, count); 
  for (i=0; i < count; i++)
  {
    cl_mem buffer = va_arg(ap, cl_mem);
    clCheck(clSetKernelArg(kernel, i, sizeof(buffer), &buffer), "clSetKernelArg");
  }
  va_end(ap);

  //Optional info
  size_t wgsize;
  size_t wgmultiple = 0;
  clCheck(clGetKernelWorkGroupInfo(kernel, cl->devices[cl->device], CL_KERNEL_WORK_GROUP_SIZE, 
          sizeof(wgsize), &wgsize, NULL), "clGetKernelWorkGroupInfo");

#ifdef CL_VERSION_1_1
  //OpenCL 1.1 required? This gets "warp" size? work groups should be a multiple of this
  clCheck(clGetKernelWorkGroupInfo(kernel, cl->devices[cl->device], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
          sizeof(wgmultiple), &wgmultiple, NULL), "clGetKernelWorkGroupInfo");
#endif //CL_VERSION_1_1

  //fprintf(stderr, "Kernel %s, Preferred work group size multiple %d, work group size %d\n", name, wgmultiple, wgsize);

  return kernel;
}

void cl_profile(OpenCL* cl, OpenCLProfile* profile)
{
  //Get profiling info
  if (!profile) profile = &cl->profile;
  cl_ulong start;
  cl_ulong end;
  clCheck(clGetEventProfilingInfo(profile->event, CL_PROFILING_COMMAND_START, sizeof(start), &start, NULL), "clGetEventProfilingInfo(COMMAND_END)");
  clCheck(clGetEventProfilingInfo(profile->event, CL_PROFILING_COMMAND_END, sizeof(end), &end, NULL), "clGetEventProfilingInfo(COMMAND_END)");
  profile->kern_time += (double)(end - start) / 1e9; // Convert nanoseconds to secs
  profile->total_time += profile->kern_time;
}

void cl_delete(OpenCL* cl)
{
  if (cl->constants) free(cl->constants);
  clReleaseProgram(cl->prog);
  clReleaseCommandQueue(cl->queue);
  clReleaseContext(cl->context);
  free(cl);
}

void cl_print_platform_info(OpenCL* cl, int i)
{
  char buffer[2048];
  printf("  -- %d --\n", i);
  clGetPlatformInfo(cl->platforms[i], CL_PLATFORM_PROFILE, sizeof(buffer), buffer, NULL);
  printf("  PROFILE = %s\n", buffer);
  clGetPlatformInfo(cl->platforms[i], CL_PLATFORM_VERSION, sizeof(buffer), buffer, NULL);
  printf("  VERSION = %s\n", buffer);
  clGetPlatformInfo(cl->platforms[i], CL_PLATFORM_NAME, sizeof(buffer), buffer, NULL);
  printf("  NAME = %s\n", buffer);
  clGetPlatformInfo(cl->platforms[i], CL_PLATFORM_VENDOR, sizeof(buffer), buffer, NULL);
  printf("  VENDOR = %s\n", buffer);
  clGetPlatformInfo(cl->platforms[i], CL_PLATFORM_EXTENSIONS, sizeof(buffer), buffer, NULL);
  printf("  EXTENSIONS = %s\n", buffer);
}

void cl_print_device_info(OpenCL* cl, int i)
{
  cl_uint buf_uint;
  cl_ulong buf_ulong;
  size_t buf_size;
  char buffer[2048];
  printf("  -- %d --\n", i);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
  printf("  DEVICE_NAME = %s\n", buffer);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL);
  printf("  DEVICE_VENDOR = %s\n", buffer);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL);
  printf("  DEVICE_VERSION = %s\n", buffer);
  clGetDeviceInfo(cl->devices[i], CL_DRIVER_VERSION, sizeof(buffer), buffer, NULL);
  printf("  DRIVER_VERSION = %s\n", buffer);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_EXTENSIONS, sizeof(buffer), buffer, NULL);
  printf("  DEVICE_EXTENSIONS = %s\n", buffer);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL);
  printf("  DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(buf_size), &buf_size, NULL);
  printf("  DEVICE_MAX_WORK_GROUP_SIZE = %lu\n", buf_size);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL);
  printf("  DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
  printf("  DEVICE_GLOBAL_MEM_SIZE = %llu\n", (unsigned long long)buf_ulong);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
  printf("  DEVICE_LOCAL_MEM_SIZE = %llu\n", (unsigned long long)buf_ulong);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
  printf("  DEVICE_MAX_MEM_ALLOC_SIZE = %llu\n", (unsigned long long)buf_ulong);
  clGetDeviceInfo(cl->devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
}

//Utility functions to build constants into kernel, 
//Place calls to these before compiling kernel
void cl_define_integer_constant(OpenCL* cl, const char* name, int value)
{
  char temp[256];
  sprintf(temp, "__constant int %s = %d;\n", name, value);
  cl_append_constant(cl, temp);
}

void cl_define_real_constant(OpenCL* cl, const char* name, Real value)
{
  char temp[256];
#ifdef SINGLE_PRECISION
  sprintf(temp, "__constant Real %s = %.8ff;\n", name, value);
#else
  sprintf(temp, "__constant Real %s = %.16lf;\n", name, value);
#endif
  cl_append_constant(cl, temp);
}

void cl_define_constants(OpenCL* cl, const char* format, ...)
{
  char temp[512];
  va_list ap;                    // Pointer to arguments list
  if (format == NULL) return;    // No format string
  va_start(ap, format);          // Parse format string for variables
  vsprintf(temp, format, ap);    // Convert symbols
  va_end(ap);
  cl_append_constant(cl, temp);
}

void cl_append_constant(OpenCL* cl, const char* str)
{
  //Reallocate ensuring enough memory for old content + new
  size_t size = (cl->constants ? strlen(cl->constants)+1 : 1);
  cl->constants = (char*)realloc(cl->constants, strlen(str) + size * sizeof(char));
  if (size==1)
    strcpy(cl->constants, str);
  else
    strcat(cl->constants, str);
}

