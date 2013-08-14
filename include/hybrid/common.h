//Shared definitions between OpenCL and C code
#ifndef INCLUDE_COMMON_H_
#define INCLUDE_COMMON_H_

//Read a value of any type as an int
#define ASINT(A) *(int*)(&A)

//Number of atom types (changes to kernels etc will be required to support > 2)
#define Comp 2

//Define data types for switching precision on both CPU and GPU
//For single precision mode, #define SINGLE_PRECISION
#ifdef SINGLE_PRECISION
  #define Real float
  #define Int unsigned int
  #define MPI_Real MPI_FLOAT
  #define RealFormat "%f"
  #ifdef __OPENCL_VERSION__
    #define Real3 float3
    #define Real4 float4
    #define Int4 uint4
    #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
  #endif
#else
  #define Real double
  #define Int ulong
  #define MPI_Real MPI_DOUBLE
  #define RealFormat "%lf"
  #ifdef __OPENCL_VERSION__
    #define Real3 double3
    #define Real4 double4
    #define Int4 ulong4
    //Enable FP64 (double precision) support
    #if defined(cl_khr_fp64)  // Khronos extension available?
      #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #elif defined(cl_amd_fp64)  // AMD extension available?
      #pragma OPENCL EXTENSION cl_amd_fp64 : enable
    #endif
    //Enable int64 atomics (not supported by amd)
    #pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable
  #endif
#endif


#ifdef USE_OPENMP
  #define OPENMP 1
#else
  #define OPENMP 0
#endif

//Macros
//TODO(timseries): consider inlining these or including these or something else, as per google coding standards on preprocessor macros.

#define printall printf
#define printroot printf

//Constants
#define SIZEOF_CHAR 1
#define LOADDATAVERSION 1

// The models typedef used for computing the kernel basis functions.
// MODEL_SPHERICAL: use the spherical model only.
// MODLE_MIXED: use the spherical+cylindrical model.
typedef enum {
  MODEL_SPHERICAL, MODEL_MIXED
} models;

typedef struct OutputCL
{
  Real rms_x;
  Real rms_diff_x;
} OutputCL;

#endif //INCLUDE_COMMON_H_
