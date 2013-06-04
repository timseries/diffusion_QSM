#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "mpi.h"

#define LOADDATAVERSION 1

#define usedtype float
#define MPI_USEDTYPE MPI_FLOAT

#define printroot if (rank == 0) printf
#define printall printf

using namespace std;

//****************************************************************************************************************
// ENUM DEFINITIONS
//****************************************************************************************************************

//----------------------------------------------------------------------------------------------------------------
// models - type of dipole model to use
//   - MODEL_SPHERICAL = only spherical model
//   - MODEL_MIXED = spherical and cylindrical models
typedef enum {
	MODEL_SPHERICAL, MODEL_MIXED
} models;

//****************************************************************************************************************
// CLASS DEFINITIONS
//****************************************************************************************************************

//----------------------------------------------------------------------------------------------------------------
// DataSpec - specifications of the data, such as size, B0 value, TE value
class DataSpec {
public:
	int		size[3];		// data size
	int		N;
	int		yoffset;
	int		zoffset;
	double	B0;
	double	bhat[3];
	double	caxis[3];
	int		nBG;
	int		nFG;
	
	unsigned	start;
	unsigned	end;
	
	DataSpec() {}
	~DataSpec() {}
};

class ModelMap {
public:
	int			*mask;
	usedtype	*x;
	usedtype	*y;
	usedtype	*z;
	int			ncyls;
	ModelMap() {
		mask = NULL;
		x = NULL;
		y = NULL;
		z = NULL;
		ncyls = 0;
	}
	~ModelMap() {}
	bool Process(DataSpec &dspec);
	
	void close() {
		if (mask != NULL)	free(mask);
		if (x != NULL)		free(x);
		if (y != NULL)		free(y);
		if (z != NULL)		free(z);		
	}
};

class Kernel {
public:
	models			model;
	usedtype		*skernel, *ckernel;
	usedtype		*kernel;
	usedtype		threshold;
	usedtype		B0;
	
	// constants and arrays for mixed model on-the-fly calculations
	usedtype		CYL2alpha;
	usedtype		CYL3alpha;
	usedtype		CYL4alpha3;
	usedtype		CYLa;
	usedtype		*ctr;
	usedtype		*sin2beta;
	usedtype		*gx;
	usedtype		*gy;
	usedtype		*gz;
	
	void CreateSphericalKernel(DataSpec &dspec);
	void CreateCylindricalKernel(DataSpec &dspec);
	void InitMixedModel(DataSpec &dspec);
	
	int				yoffset;
	int				zoffset;
	
	ModelMap		modelmap;
	int				nnz;
	int				size;
	int				halfsize;
	int				N;
	
	Kernel() {
		kernel = NULL; 
		nnz = 0;
		skernel = NULL;
		ckernel = NULL;
		ctr = NULL;
		sin2beta = NULL;
		gx = NULL;
		gy = NULL;
		gz = NULL;
	}
	~Kernel() {}
	
	bool Create(models &model, DataSpec &dspec, usedtype &threshold);
	usedtype Get(int x, int y, int z, int o);
	usedtype GetCyl(int mix, int x, int y, int z);
	void close(void) {	
		modelmap.close();
		
		if (kernel != NULL) free(kernel);
		if (skernel != NULL) free(skernel);
		if (ckernel != NULL) free(ckernel);
		if (ctr != NULL) free(ctr);
		if (sin2beta != NULL) free(sin2beta);
		if (gx != NULL) free(gx);
		if (gy != NULL) free(gy);
		if (gz != NULL) free(gz);
	}
	//void output(void);
};

/*
class MemLog {
	char				TimePoints[32 * 50];
	PetscLogDouble		mem[50];
	int					maxlogs;
	int					n;
public:
	MemLog() { 
		maxlogs = 50;
		n = 0;
	}
	~MemLog() {}
	
	void AddLog(const char *tp) {
		PetscLogDouble	x;
		memcpy(&TimePoints[n*32], tp, 32);
		TimePoints[n*32 + 31] = '\0';
		PetscMallocGetCurrentUsage(&x);
		mem[n] = x;
		n++;
	}
	
	void PrintLog() {
		if (rank == 0) {
			printf("\nMemory Log (in MB)\n");
			printf("===========================================================================\n");
			printf("%-32s   Memory (MB)", "TimePoint");
			printf("\n---------------------------------------------------------------------------\n");
			for (int m = 0; m < n; m++) {
				printf("\n%-32s   %5.1f", &TimePoints[m*32], mem[m]/1024/1024);
			}
			printf("\n\n");
		}
	}
};
*/

class ArgHandler {
	int argc;
	char** args;
public:
	ArgHandler(){}
	~ArgHandler(){}
	
	void Init(int argc, char** args) {
		this->argc = argc;
		this->args = args;
	}
	
	bool GetArg(const char* option, char *&str) {
		for (int i = 0; i < argc; i++) {
			if (strcmp(args[i], option) == 0) {
				if (i == argc-1) {
					return false;
				}
				str = (char*) calloc(strlen(args[i+1])+1, sizeof(char));
				strcpy(str, args[i+1]);
				return true;
			}
		}
		return false;
	}
	
	bool GetArg(const char* option, usedtype &value) {
		for (int i = 0; i < argc; i++) {
			if (strcmp(args[i], option) == 0) {
				if (i == argc-1) {
					return false;
				}
				value = atof(args[i+1]);
				return true;
			}
		}
		return false;
	}
	
	bool GetArg(const char* option, int &value) {
		for (int i = 0; i < argc; i++) {
			if (strcmp(args[i], option) == 0) {
				if (i == argc-1) {
					return false;
				}
				value = atoi(args[i+1]);
				return true;
			}
		}
		return false;
	}
};

class Output {
	MPI_File	binfile;
	MPI_File	matfile;
	FILE		*errfile;
	bool		initialised;
	char		*tmpstr;
	
public:
	char		*outdir;
	
	Output(){initialised = false;};
	~Output(){};
	
	void Init();
	void LocalArray(int onproc, usedtype* array, int ndims, int* dims, const char* arrayname);
	void LocalArray(int onproc, bool* array, int ndims, int* dims, const char* arrayname);
	void DistrArray(usedtype* array, int localsize, int ndims, int* dims, const char* arrayname);
	void Close();
	
};

int				size;
int				rank;

//----------------------------------------------------------------------------------------------------------------
// Load functions
bool loadDeltaB(DataSpec &dspec, usedtype *&DeltaB);
bool loadMask(DataSpec &dspec, bool *&Mask);
void checkEndianness(MPI_File &fptr, bool &flgByteSwap);
void checkVersion(MPI_File &fptr, bool flgByteSwap);
void byteswap(char* buf, int buflength, int dataTypeBytes);

//----------------------------------------------------------------------------------------------------------------
// Main solve function
bool FullPass(models &model, DataSpec &dspec, usedtype *&DeltaB, bool *&mask, Kernel &kernel, usedtype *&chi);

Output			myout;
ArgHandler		arghandler;

int	main(int argc, char** args) {	
	
	MPI_Init(&argc, &args);
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	printroot("\n------------------------------------------\n");
	printroot("MPI Environment\n");
	printroot("Number of processes: %d\n", size);
	
	models					model;
	
	usedtype				*deltab = NULL;
	usedtype				*chi = NULL;
	bool					*mask = NULL;
	
	Kernel					kernel;
	DataSpec				dspec;
	
	char					*filepath = NULL;
	int						tmpInt;
	
	usedtype				threshold;
	
	bool					flg;	
	char					*str = NULL;
	
	MPI_File				fptr;
	
	printroot("\n------------------------------------------\n");
	printroot("INITIALISING\n");
		
	// Initialise argument handler
	arghandler.Init(argc, args);
#ifdef DEBUG
	printroot("   Argument handler initialised\n");
#endif

	// Initialise output object
	myout.Init();
#ifdef DEBUG
	printroot("   Output object initialised\n");
#endif
	
	// GET THRESHOLD
	flg = arghandler.GetArg("-threshold", threshold);
	if (!flg) threshold = 1e-6;
	printroot("   Threshold = %0.3e\n", threshold);
	
	// GET MODEL
	flg = arghandler.GetArg("-model", str);
	if (!flg) {
		model = MODEL_MIXED;
	}
	else {
		switch (str[0]) {
			case 's':
				model = MODEL_SPHERICAL;
				break;
			case 'm':
				model = MODEL_MIXED;
				break;
			default:
				printf("Unrecognised model");
				return 1;
				break;
		}
	}
	free(str); str = NULL;
	
	printroot("   Model: %s\n", model == MODEL_SPHERICAL ? "spherical" : "mixed");

	dspec.caxis[0] = 1;
	dspec.caxis[1] = 0;
	dspec.caxis[2] = 0;	
	
	// LOAD DATA
	if (!loadDeltaB(dspec, deltab)) goto exitnow;
	myout.DistrArray(deltab, dspec.end - dspec.start, 3, dspec.size, "deltab");
	
	// LOAD MASK
	if (!loadMask(dspec, mask)) goto exitnow;
	myout.LocalArray(0, mask, 3, dspec.size, "mask");
	
	// LOAD MODELMASK
	if (!kernel.modelmap.Process(dspec)) goto exitnow;
	
	//============================================================================================================================== 
	// CREATE CHI VECTOR
	chi = (usedtype*) calloc(dspec.N, sizeof(usedtype));
	
	flg = arghandler.GetArg("-chi", filepath);
	if (flg) {
		if (rank == 0) {
			MPI_File_open(MPI_COMM_SELF, filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
			MPI_File_read(fptr, chi, dspec.N, MPI_USEDTYPE, MPI_STATUS_IGNORE);
			MPI_File_close(&fptr);
		}
		MPI_Bcast(chi, dspec.N, MPI_USEDTYPE, 0, MPI_COMM_WORLD);
		free(filepath);
		filepath = NULL;
	}
	else {
		memset(chi, 0, dspec.N*sizeof(usedtype));
	}

	// IF NUMBER OF BACKGROUND VOXELS IS BELOW <THRESHOLD> ADD EXTRA ZERO PADDING
		
	//============================================================================================================================== 
	// CREATE KERNEL
	if (!kernel.Create(model, dspec, threshold)) goto exitnow;		
		
	//==============================================================================================================================
	// FULL DATA PROCESSING
	
	if (!FullPass(model, dspec, deltab, mask, kernel, chi)) goto exitnow;
	
	//==============================================================================================================================
	// OUTPUT CHI
	myout.LocalArray(0, chi, 3, dspec.size, "chi");
	
	if (rank == 0) {
		filepath = (char*) calloc(strlen(myout.outdir) + 9, sizeof(char));
		sprintf(filepath, "%s/chi.bin", myout.outdir);
		MPI_File_open(MPI_COMM_SELF, filepath, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
		MPI_File_write(fptr, chi, dspec.N, MPI_USEDTYPE, MPI_STATUS_IGNORE);
		MPI_File_close(&fptr);
		free(filepath);
		filepath = NULL;
	}
	
	
	// DESTROY EVERYTHING
exitnow:
	printroot("\n------------------------------------------\n");
	printroot("CLEANING UP\n");
	printroot("Freeing memory ...\n");
	if (deltab != NULL)			free(deltab);
	if (chi != NULL)			free(chi);
	if (mask != NULL)			free(mask);
	kernel.close();
	printroot("Finished\n", rank); fflush(stdout);
	
	myout.Close();
	
	MPI_Finalize();
}

//==============================================================================================================================
// OUTPUT CLASS

void Output::Init() {
	bool flg;
	
	flg = arghandler.GetArg("-out", outdir);
	if (!flg) {
		outdir = (char*) calloc(7, sizeof(char));
		sprintf(outdir, "output");
	}

	tmpstr = (char*) calloc(1024, sizeof(char));
	
	// redirect stderr to error files
	sprintf(tmpstr, "%s/err%03d.txt", outdir, rank);
	if ((errfile = freopen(tmpstr, "w", stderr)) == NULL) {
		printall("[%d] Could not redirect stderr\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
#ifdef DEBUG
	printroot("   Output::Init : redirected stderr\n");
#endif
	
	// Open binary file
	if (rank == 0) {
		sprintf(tmpstr, "%s/out.bin", outdir);
		MPI_File_open(MPI_COMM_SELF, tmpstr, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &binfile);
#ifdef DEBUG
		printroot("   Output::Init : opened binary file\n");
#endif
		MPI_File_set_size(binfile, 0);
	
	// Open matlab file on process 0
		sprintf(tmpstr, "%s/out.m", outdir);
		MPI_File_open(MPI_COMM_SELF, tmpstr, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &matfile);
		MPI_File_set_size(matfile, 0);
	
#ifdef DEBUG
		printroot("   Output::Init : opened matlab file\n");
#endif
		
		int endchkval = 1;
		// Write precision type to mat file
		if (sizeof(usedtype) == sizeof(float)) {
			sprintf(tmpstr, "TmpChiMapOut.precision = 'single';\n");
		}
		else {
			sprintf(tmpstr, "TmpChiMapOut.precision = 'double';\n");
		}
		MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
		
		// Write endian check value to binary file
		MPI_File_write(binfile, &endchkval, 1, MPI_INT, MPI_STATUS_IGNORE);
		
		// Write file open command to matlab
		sprintf(tmpstr, "TmpChiMapOut.fid = fopen('out.bin', 'r', 'b');\n");
		sprintf(tmpstr, "%sif fread(TmpChiMapOut.fid, 1, 'int32') ~= 1\n", tmpstr);
		sprintf(tmpstr, "%s\tfclose(TmpChiMapOut.fid);\n", tmpstr);
		sprintf(tmpstr, "%s\tTmpChiMapOut.fid = fopen('out.bin', 'r', 'l');\n", tmpstr);
		sprintf(tmpstr, "%s\tfread(TmpChiMapOut.fid, 1, 'int32');\n", tmpstr);
		sprintf(tmpstr, "%send\n", tmpstr);
		MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	
	initialised = true;
}

void Output::LocalArray(int onproc, usedtype* array, int ndims, int* dims, const char* arrayname) {
	int			ndims0;
	int			*dims0;
	int			arraynamelength;
	char		*arrayname0;
	int			n;
	usedtype	*array0;
	
	if (rank != onproc && rank != 0) return;
	
	// if data is not on process 0, send data to process 0
	if (onproc != 0) {
		
		// send data
		if (rank == onproc) {
			n = 1;
			for (int i = 0; i < ndims; i++)
				n *= dims[i];
			
			arraynamelength = strlen(arrayname);
			
			MPI_Send(&ndims,			1,					MPI_INT,		0, 0, MPI_COMM_WORLD);
			MPI_Send(dims,				ndims,				MPI_INT,		0, 1, MPI_COMM_WORLD);
			MPI_Send(&arraynamelength,	1,					MPI_INT,		0, 2, MPI_COMM_WORLD);
			MPI_Send((char*)arrayname,	arraynamelength,	MPI_CHAR,		0, 3, MPI_COMM_WORLD);
			MPI_Send(array,				n,					MPI_USEDTYPE,	0, 4, MPI_COMM_WORLD);
		}
		
		// receive data
		if (rank == 0) {
			MPI_Recv(&ndims0,			1,					MPI_INT,		onproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			dims0 = (int*) calloc(ndims0, sizeof(int));
			MPI_Recv(dims0,				ndims0,				MPI_INT,		onproc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&arraynamelength,	1,					MPI_INT,		onproc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			arrayname0 = (char*) calloc(arraynamelength, sizeof(char));
			MPI_Recv(arrayname0,		arraynamelength,	MPI_CHAR,		onproc, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			n = 1;
			for (int i = 0; i < ndims0; i++)
				n *= dims0[i];
			array0 = (usedtype*) calloc(n, sizeof(usedtype));
			MPI_Recv(array0,			n,					MPI_USEDTYPE,	onproc, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	// otherwise, set values and pointers as though process 0 had sent data to itself
	else {
		ndims0 = ndims;
		dims0 = dims;
		arrayname0 = (char*) arrayname;
		array0 = array;
		n = 1;
		for (int i = 0; i < ndims0; i++)
			n *= dims0[i];		
	}

	// On process 0, write the data to file
	if (rank == 0) {
		
		// write array to binary file
		MPI_File_write(binfile, array0, n, MPI_USEDTYPE, MPI_STATUS_IGNORE);
		
		// write read-in matlab code to matlab file
		if (ndims > 1) {
			sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname0, n, dims0[0]);
			for (int i = 1; i < ndims0; i++) {
				sprintf(tmpstr, "%s %d", tmpstr, dims0[i]);
			}
			sprintf(tmpstr, "%s]);\n", tmpstr);
		}
		else {
			sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname0, n);
		}
		
		MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
		
		// if data was sent from another process, free the memory that was allocated
		if (onproc != 0) {
			free(dims0);
			free(arrayname0);
			free(array0);
		}
		
	}
	
	// THE FOLLOWING IS THE TRUE MPI IMPLEMENTATION. DUE TO A BUG ON AVOCA, THIS HAS BEEN REPLACED BY THE ABOVE
	/* 
	if (rank != onproc) return;
	
	int n = 1;
	for (int i = 0; i < ndims; i++)
		n *= dims[i];
	
	MPI_File_write_shared(binfile, array, n, MPI_USEDTYPE, MPI_STATUS_IGNORE);
	
	if (ndims > 1) {
		sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname, n, dims[0]);
		for (int i = 1; i < ndims; i++) {
			sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
		}
		sprintf(tmpstr, "%s]);\n", tmpstr);
	}
	else {
		sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname, n);
	}
	
	MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
	 */
}

void Output::LocalArray(int onproc, bool* array, int ndims, int* dims, const char* arrayname) {
	int			ndims0;
	int			*dims0;
	int			arraynamelength;
	char		*arrayname0;
	int			n;
	bool		*array0;
	
	if (rank != onproc && rank != 0) return;
	
	// if data is not on process 0, send data to process 0
	if (onproc != 0) {
		
		// send data
		if (rank == onproc) {
			n = 1;
			for (int i = 0; i < ndims; i++)
				n *= dims[i];
			
			arraynamelength = strlen(arrayname);
			
			MPI_Send(&ndims,			1,					MPI_INT,		0, 0, MPI_COMM_WORLD);
			MPI_Send(dims,				ndims,				MPI_INT,		0, 1, MPI_COMM_WORLD);
			MPI_Send(&arraynamelength,	1,					MPI_INT,		0, 2, MPI_COMM_WORLD);
			MPI_Send((char*)arrayname,	arraynamelength,	MPI_CHAR,		0, 3, MPI_COMM_WORLD);
			MPI_Send(array,				n,					MPI_CHAR,		0, 4, MPI_COMM_WORLD);
		}
		
		// receive data
		if (rank == 0) {
			MPI_Recv(&ndims0,			1,					MPI_INT,		onproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			dims0 = (int*) calloc(ndims0, sizeof(int));
			MPI_Recv(dims0,				ndims0,				MPI_INT,		onproc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&arraynamelength,	1,					MPI_INT,		onproc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			arrayname0 = (char*) calloc(arraynamelength, sizeof(char));
			MPI_Recv(arrayname0,		arraynamelength,	MPI_CHAR,		onproc, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			n = 1;
			for (int i = 0; i < ndims0; i++)
				n *= dims0[i];
			array0 = (bool*) calloc(n, sizeof(bool));
			MPI_Recv(array0,			n,					MPI_CHAR,		onproc, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	
	// otherwise, set values and pointers as though process 0 had sent data to itself
	else {
		ndims0 = ndims;
		dims0 = dims;
		arrayname0 = (char*) arrayname;
		array0 = array;
		n = 1;
		for (int i = 0; i < ndims0; i++)
			n *= dims0[i];		
	}
	
	// On process 0, write the data to file
	if (rank == 0) {
		
		// write array to binary file
		MPI_File_write(binfile, array0, n, MPI_CHAR, MPI_STATUS_IGNORE);
		
		// write read-in matlab code to matlab file
		if (ndims > 1) {
			sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, 'uint8'), [%d", arrayname0, n, dims0[0]);
			for (int i = 1; i < ndims0; i++) {
				sprintf(tmpstr, "%s %d", tmpstr, dims0[i]);
			}
			sprintf(tmpstr, "%s]);\n", tmpstr);
		}
		else {
			sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, 'uint8');\n", arrayname0, n);
		}
		
		MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
		
		// if data was sent from another process, free the memory that was allocated
		if (onproc != 0) {
			free(dims0);
			free(arrayname0);
			free(array0);
		}
		
	}
	
	// THE FOLLOWING IS THE TRUE MPI IMPLEMENTATION. DUE TO A BUG ON AVOCA, THIS HAS BEEN REPLACED BY THE ABOVE
	/*
	if (rank != onproc) return;
	int n = 1;
	for (int i = 0; i < ndims; i++)
		n *= dims[i];
	
	MPI_File_write_shared(binfile, array, n, MPI_CHAR, MPI_STATUS_IGNORE);
	
	if (ndims > 1) {
		sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, 'uint8'), [%d", arrayname, n, dims[0]);
		for (int i = 1; i < ndims; i++) {
			sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
		}
		sprintf(tmpstr, "%s]);\n", tmpstr);
	}
	else {
		sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, 'uint8');\n", arrayname, n);
	}
	
	MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
	 */
}

void Output::DistrArray(usedtype* array, int localsize, int ndims, int* dims, const char* arrayname) {
	int			n;
	usedtype	*array0 = NULL;
	MPI_Request req[2];
	
	// Send data to process 0
	MPI_Isend(&localsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req[0]);
	MPI_Isend(array, localsize, MPI_USEDTYPE, 0, 1, MPI_COMM_WORLD, &req[1]);
		
	// On process 0, receive and write data to file.
	if (rank == 0) {
		
		n = 0; // current length of array0
		
		// Receive array length and data from each process in order, and write to binary file
		for (int p = 0; p < size; p++) {
			MPI_Recv(&localsize, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			// realloc memory if the localsize has increased
			if (n < localsize) {
				array0 = (usedtype*) realloc(array0, localsize*sizeof(usedtype));
				n = localsize;
			}
			MPI_Recv(array0, localsize, MPI_USEDTYPE, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_File_write(binfile, array0, localsize, MPI_USEDTYPE, MPI_STATUS_IGNORE);
		}
		
		
		// Write read-in matlab code to the matlab file
		int n = 1;
		for (int i = 0; i < ndims; i++)
			n *= dims[i];
		
		if (ndims > 1) {
			sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname, n, dims[0]);
			for (int i = 1; i < ndims; i++) {
				sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
			}
			sprintf(tmpstr, "%s]);\n", tmpstr);
		}
		else {
			sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname, n);
		}
		
		MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	
	// Free array
	if (array0 != NULL)	free(array0);
	
	// Wait for requests to complete
	MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
	
	// THE FOLLOWING IS THE TRUE MPI IMPLEMENTATION. DUE TO A BUG ON AVOCA, THIS HAS BEEN REPLACED BY THE ABOVE
	/*
	int n = 1;
	for (int i = 0; i < ndims; i++)
		n *= dims[i];
	
	MPI_File_write_ordered(binfile, array, localsize, MPI_USEDTYPE, MPI_STATUS_IGNORE);
	
	if (rank == 0) {
		if (ndims > 1) {
			sprintf(tmpstr, "%s = reshape(fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision), [%d", arrayname, n, dims[0]);
			for (int i = 1; i < ndims; i++) {
				sprintf(tmpstr, "%s %d", tmpstr, dims[i]);
			}
			sprintf(tmpstr, "%s]);\n", tmpstr);
		}
		else {
			sprintf(tmpstr, "%s = fread(TmpChiMapOut.fid, %d, TmpChiMapOut.precision);\n", arrayname, n);
		}
		
		MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	 */
}	

void Output::Close() {
	
	MPI_Offset offset, disp;
	
	if (rank == 0) {
		sprintf(tmpstr, "fclose(TmpChiMapOut.fid);\n");
		sprintf(tmpstr, "%sclear TmpChiMapOut\n", tmpstr);
		//MPI_File_write_shared(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write(matfile, tmpstr, strlen(tmpstr), MPI_CHAR, MPI_STATUS_IGNORE);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Close binary file
//	MPI_File_get_position_shared(binfile, &offset);
//	MPI_File_get_byte_offset(binfile, offset, &disp);
//	MPI_File_set_size(binfile, disp); 
	if (rank == 0)
		MPI_File_close(&binfile);
	
	// Close matlab file
//	MPI_File_get_position_shared(matfile, &offset);
//	MPI_File_get_byte_offset(matfile, offset, &disp);
//	MPI_File_set_size(matfile, disp); 
	if (rank == 0)
		MPI_File_close(&matfile);
	
	// Close error log file and remove if empty
	if (ftell(errfile) == 0) {
		fclose(errfile);
		sprintf(tmpstr, "%s/err%03d.txt", outdir, rank);
		remove(tmpstr);
	}
	else {
		fclose(errfile);
	}
	
	// Close files
	initialised = false;
	free(tmpstr);
}

//==============================================================================================================================
// PROCESS MODEL ARRAY

bool ModelMap::Process(DataSpec &dspec) {
	int			err;
	char		*filepath;
	MPI_File	fptr;
	double		buf[5];
	bool		flgByteSwap;
	MPI_Offset	offset, disp;
	
	bool		flg;
	usedtype	threshold;
	double		*mx, *my, *mz;
	
	printroot( "Loading model map ...\n");
	flg = arghandler.GetArg("-modelmap", filepath);
	if (!flg) {
		printroot( "Model map file not specified\n");
		return false;
	}
	printroot( "   file: %s\n", filepath);
	
	// open file for reading
	err = MPI_File_open(MPI_COMM_SELF, filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
	if (err) {printroot( "Could not open model mask file"); return false;}
		
	// check endianness
	checkEndianness(fptr, flgByteSwap);
	
	// check version
	checkVersion(fptr, flgByteSwap);
	
	// read in header
	MPI_File_read(fptr, buf, 5, MPI_DOUBLE, MPI_STATUS_IGNORE);
		
	// Check parameters match
	if (flgByteSwap) byteswap((char*)buf, 5, sizeof(double));
		
	if (dspec.size[0] != buf[1]) {printroot( "Size of first dimension of ModelMask does not match DeltaB"); return false; }
	if (dspec.size[1] != buf[2]) {printroot( "Size of second dimension of ModelMask does not match DeltaB"); return false; }
	if (dspec.size[2] != buf[3]) {printroot( "Size of third dimension of ModelMask does not match DeltaB"); return false; }
	if (buf[4] != 3) {printroot( "Size of fourth dimension of ModelMask is not 3"); return false; }
	if (buf[0] != buf[1] * buf[2] * buf[3] * buf[4]) {printroot( "Number of elements ModelMask does not match array sizes"); return false; }
		
	// Load model data
	mx = (double*) calloc(dspec.N, sizeof(double));
	my = (double*) calloc(dspec.N, sizeof(double));
	mz = (double*) calloc(dspec.N, sizeof(double));
	
	if (mx == NULL || my == NULL || mz == NULL) {
		printroot("Not enough memory available.");
		return false;
	}
	
	MPI_File_read(fptr, mx, dspec.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(fptr, my, dspec.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read(fptr, mz, dspec.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
	
	if (flgByteSwap) byteswap((char*)mx, dspec.N, sizeof(double));
	if (flgByteSwap) byteswap((char*)my, dspec.N, sizeof(double));
	if (flgByteSwap) byteswap((char*)mz, dspec.N, sizeof(double));
	
	MPI_File_close(&fptr);
		
	// GET THRESHOLD
	flg = arghandler.GetArg("-mt", threshold);
	if (!flg) threshold = 0.2;
	
	printroot( "   threshold = %0.3f\n", threshold);

	// CREATE MASK
	mask = (int*) calloc(dspec.N, sizeof(int));
	
	int ix = 0;
	for (int i = 0; i < dspec.N; i++) {
		mask[i] = (sqrt(mx[i]*mx[i] + my[i]*my[i] + mz[i]*mz[i]) > threshold) ? ix++ : -1;
	}
	
	// GATHER X Y AND Z ARRAYS
	ncyls = ix;
	printroot( "   number cylinders = %d\n", ncyls);
	x = (usedtype*) calloc(ncyls, sizeof(usedtype));
	y = (usedtype*) calloc(ncyls, sizeof(usedtype));
	z = (usedtype*) calloc(ncyls, sizeof(usedtype));
	for (int i = 0; i < dspec.N; i++) {
		if (mask[i] >= 0) {
			x[mask[i]] = mx[i];
			y[mask[i]] = my[i];
			z[mask[i]] = mz[i];
		}
	}
	printroot( "   ix = %d\n", ix);
	
	free(mx);
	free(my);
	free(mz);
	
	return true;
}

//==============================================================================================================================
// KERNEL FUNCTIONS

bool Kernel::Create(models &model, DataSpec &dspec, usedtype &threshold) {
	printroot("Creating kernel\n");
	
	this->model = model;
	this->threshold = threshold;
	B0 = dspec.B0;
	
	if (model == MODEL_SPHERICAL) {
		size = ceil(pow(1.0/threshold, 1.0/3.0));
	}
//	else if (model == MODEL_CYLINDRICAL) {
//		size = ceil(pow(1.0*threshold, -1.0/2.0));
//	}
	else { //model == MODEL_MIXED
		size = ceil(pow(1.0/threshold, 1.0/3.0)); // Create a spherical kernel. Cylinders will be calculated on the fly.
	}

	//if (iseven(size)) size++;
	if (size/2 == size/2.0) size++;
	
	long maxdatasize = dspec.size[0] > dspec.size[1] & dspec.size[0] > dspec.size[2] ? dspec.size[0] : dspec.size[1] > dspec.size[2] ? dspec.size[1] : dspec.size[2];
	if (size > maxdatasize * 2 + 1) size = maxdatasize * 2 + 1;
	
	printroot("   threshold = %0.3e\n", threshold);
	printroot("   kernel size = %d\n", size);
	halfsize = size/2;
	
	N = size*size*size;
	
	nnz = N;
	
	if (model == MODEL_SPHERICAL) {
		skernel = (usedtype*) calloc(nnz, sizeof(usedtype));
	}
//	else if (model == MODEL_CYLINDRICAL) {
//		ckernel = (usedtype*) calloc(nnz, sizeof(usedtype));
//	}
	else {
		skernel = (usedtype*) calloc(nnz, sizeof(usedtype));
		ctr = (usedtype*) calloc(modelmap.ncyls, sizeof(usedtype));
		sin2beta = (usedtype*) calloc(modelmap.ncyls, sizeof(usedtype));
		gx = (usedtype*) calloc(modelmap.ncyls, sizeof(usedtype));
		gy = (usedtype*) calloc(modelmap.ncyls, sizeof(usedtype));
		gz = (usedtype*) calloc(modelmap.ncyls, sizeof(usedtype));
	}
	
	yoffset = size;
	zoffset = size*size;
	
	// Create kernel
	if (model == MODEL_SPHERICAL)
		CreateSphericalKernel(dspec);
//	else if (model == MODEL_CYLINDRICAL)
//		CreateCylindricalKernel(dspec);
	else {
		CreateSphericalKernel(dspec);
		InitMixedModel(dspec);
	}
	
	printroot("   non-zeros = %d\n", nnz);
	
	return true;
}

void Kernel::CreateSphericalKernel(DataSpec &dspec) {
	usedtype			SPHa = dspec.B0 / (4 * M_PI);
	usedtype			minkernelvalue, maxkernelvalue = 0;
	
	for (int x = -halfsize; x <= halfsize; x++) {
		for (int y = -halfsize; y <= halfsize; y++) {
			for (int z = -halfsize; z <= halfsize; z++) {
				int p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
				if (x == 0 && y == 0 && z == 0) {
					skernel[p] = 0;
					nnz--;
				}
				else {
					usedtype rmag	= sqrt(x*x + y*y + z*z);
					usedtype rx = x/rmag;
					usedtype ry = y/rmag;
					usedtype rz = z/rmag;
					
					usedtype br = dspec.bhat[0]*rx + dspec.bhat[1]*ry + dspec.bhat[2]*rz;
					skernel[p] = SPHa * ( 3*br*br - 1) / (rmag*rmag*rmag);
					
					if (fabs(skernel[p]) > maxkernelvalue)
						maxkernelvalue = fabs(skernel[p]);
				}
				
			}			
		}		
	}
	printroot( "   max value : %0.20e\n", maxkernelvalue);
	minkernelvalue = threshold*maxkernelvalue;	
	printroot( "   threshold value : %0.20e\n", minkernelvalue);
	for (int i = 0; i < N; i++) {
		if (fabs(skernel[i]) < minkernelvalue) {
			skernel[i] = 0;
			nnz--;
		}
	}	
}

void Kernel::CreateCylindricalKernel(DataSpec &dspec) {
	// returns the number of non-zeros
	
	usedtype	a = 0.47;
	usedtype	CYL2alpha = 2*a;
	usedtype	CYL3alpha = 3*a;
	usedtype	CYL4alpha3 = 1/(4*a*a*a);
	usedtype	CYLa = dspec.B0 / 2 / M_PI; 
	
	usedtype	cmag = sqrt(dspec.caxis[0]*dspec.caxis[0] + dspec.caxis[1]*dspec.caxis[1] + dspec.caxis[2]*dspec.caxis[2]);
	usedtype	cx = dspec.caxis[0]/cmag;
	usedtype	cy = dspec.caxis[1]/cmag;
	usedtype	cz = dspec.caxis[2]/cmag;
	
	usedtype	cosbeta = cx*dspec.bhat[0] + cy*dspec.bhat[1] + cz*dspec.bhat[2];
	
	usedtype	CYLsin2beta = 1 - cosbeta*cosbeta;
	
	usedtype	gx = dspec.bhat[0] - cosbeta * cx;
	usedtype	gy = dspec.bhat[1] - cosbeta * cy;
	usedtype	gz = dspec.bhat[2] - cosbeta * cz;
	
	usedtype	minkernelvalue, maxkernelvalue = 0;
	
	for (int x = -halfsize; x <= halfsize; x++) {
		for (int y = -halfsize; y <= halfsize; y++) {
			for (int z = -halfsize; z <= halfsize; z++) {
				int p = (x+halfsize) + (y+halfsize)*yoffset + (z+halfsize)*zoffset;
				if (x == 0 && y == 0 && z == 0) {
					ckernel[p] = dspec.B0/6 * (3 * cosbeta*cosbeta - 1);
				}
				else {
					usedtype rc = x*cx + y*cy + z*cz;
					
					if (fabs(rc) <= CYL2alpha) {
						usedtype r2 = 1/(x*x + y*y + z*z - rc*rc);
						usedtype rg = x*gx + y*gy + z*gz;
						
						usedtype h = CYL2alpha-fabs(rc);
						usedtype q = h*h * (CYL3alpha-h) * CYL4alpha3;
						
						ckernel[p] = CYLa * ( 2*rg*rg*r2 - CYLsin2beta) * r2 * q;
					}
					else {
						ckernel[p] = 0.0;
					}
					
					if (fabs(ckernel[p]) > maxkernelvalue)
						maxkernelvalue = fabs(ckernel[p]);
				}
			}			
		}		
	}
	
	minkernelvalue = threshold*dspec.B0;	
	for (int i = 0; i < N; i++) {
		if (fabs(ckernel[i]) < minkernelvalue) {
			ckernel[i] = 0;
			nnz--;
		}
	}
}

void Kernel::InitMixedModel(DataSpec &dspec) {
	usedtype	a = 0.475;
	CYL2alpha = 2*a;
	CYL3alpha = 3*a;
	CYL4alpha3 = 1/(4*a*a*a);
	CYLa = dspec.B0 / 2 / M_PI; 
	
	for (int i = 0; i < modelmap.ncyls; i++) {
		usedtype	cmag = sqrt(modelmap.x[i]*modelmap.x[i] + modelmap.y[i]*modelmap.y[i] + modelmap.z[i]*modelmap.z[i]);
		usedtype	cx = modelmap.x[i]/cmag;
		usedtype	cy = modelmap.y[i]/cmag;
		usedtype	cz = modelmap.z[i]/cmag;
		
		modelmap.x[i] = cx;
		modelmap.y[i] = cy;
		modelmap.z[i] = cz;
		
		usedtype	cosbeta = cx*dspec.bhat[0] + cy*dspec.bhat[1] + cz*dspec.bhat[2];
		
		sin2beta[i] = 1 - cosbeta*cosbeta;
		
		ctr[i] = dspec.B0/6 * (3 * cosbeta*cosbeta - 1);
		gx[i] = dspec.bhat[0] - cosbeta * cx;
		gy[i] = dspec.bhat[1] - cosbeta * cy;
		gz[i] = dspec.bhat[2] - cosbeta * cz;
	}
}

/*
void Kernel::output(void) {
	ExportToMFile(0, "loadkernel", skernel, N, "kernel", "kernel = reshape(kernel, [%d %d %d]);", size, size, size);
} 
 */

usedtype Kernel::Get(int x, int y, int z, int o) {
	int ix = x + y*yoffset + z*zoffset;
	int mix;
	if (model == MODEL_MIXED &&  (mix = modelmap.mask[o]) != -1) {
		x -= halfsize;
		y -= halfsize;
		z -= halfsize;
		if (x == 0 && y == 0 && z == 0) {
			return ctr[mix];
		}
		else {
			usedtype rc = x*modelmap.x[mix] + y*modelmap.y[mix] + z*modelmap.z[mix];			
			if (fabs(rc) <= CYL2alpha) {
				usedtype r2 = 1/(x*x + y*y + z*z - rc*rc);
				usedtype rg = x*gx[mix] + y*gy[mix] + z*gz[mix];
				
				usedtype h = CYL2alpha-fabs(rc);
				usedtype q = h*h * (CYL3alpha-h) * CYL4alpha3;
				
				usedtype retval = CYLa * ( 2*rg*rg*r2 - sin2beta[mix]) * r2 * q;
				if (fabs(retval) < threshold*B0)
					return 0.0;
				else
					return retval;
			}
			else {
				return 0.0;
			}		
		}
	}
//	else if (model == MODEL_CYLINDRICAL) {
//		return ckernel[ix];		
//	}
	else {
		return skernel[ix];
	}
}

usedtype Kernel::GetCyl(int mix, int x, int y, int z) {
	usedtype rc = x*modelmap.x[mix] + y*modelmap.y[mix] + z*modelmap.z[mix];			
	if (fabs(rc) <= CYL2alpha) {
		usedtype r2 = 1/(x*x + y*y + z*z - rc*rc);
		usedtype rg = x*gx[mix] + y*gy[mix] + z*gz[mix];
		
		usedtype h = CYL2alpha-fabs(rc);
		usedtype q = h*h * (CYL3alpha-h) * CYL4alpha3;
		
		usedtype retval = CYLa * ( 2*rg*rg*r2 - sin2beta[mix]) * r2 * q;
		if (fabs(retval) < threshold*B0)
			return 0.0;
		else
			return retval;
	}
	else {
		return 0.0;
	}			
}

//==============================================================================================================================
// FULL DATA PASS 

bool FullPass(models &model, DataSpec &dspec, usedtype *&DeltaB, bool *&mask, Kernel &kernel, usedtype *&chi) {
	// Landweber iteration
	// x1_{n+1} = x_{n} - alpha * A' * (A * x_{n} - b) - beta * D' * (D * x_{n})
	// where b     = deltab, 
	//       x     = chi, 
	//       alpha = 0.5 * tau, 
	//       beta  = 0.5 * tau,
	//       tau   = 2/norm(A), and 
	//       D     = represents a 3D laplacian filter
	
	
	usedtype				tau = 0.15, alpha = 0.75, beta = 0.25;
	
	bool					PreCalcCylinders = true;
	usedtype				*cylColumns;
	
	usedtype				*Ax_b;	
	usedtype				*AtAx_b;
	usedtype				*Dx;
	usedtype				*DtDx;
	usedtype				*x;
	usedtype				*new_x;
	
	usedtype				D3_96 = 3.0/96.0;
	usedtype				D10_96 = 10.0/96.0;
	
	int						*FGindices; // indices corresponding to foreground elements
	
	int						dStart, dEnd;
	int						dN;	
	int						recvcounts, displs;
	
	//stopping criteria
	usedtype				relative_threshold = 1e-6;
	usedtype				absolute_threshold = 1e-16;
	int						max_iters = 1000;
	double					rms_x, rms_diff_x; //, old_rms_diff_x;
	int						iteration = 0;
	
	double					mem;
	
	int						o, ox, oy, oz;
	int						p, px, py, pz;
	int						rx, ry, rz;
	
	int						roffset;
	
	MPI_File				fptr;
	char					*fname;
	
	time_t					tsecs, tmins;
	time_t					tStart, tEnd, tsave;
	time_t					tIterStart1, tIterStart2, tIterEnd1, tIterEnd2;
	
	bool					flg;
	
	printroot("Full pass ...\n");
	
	fname = (char*) calloc(256, sizeof(char));

	flg = arghandler.GetArg("-rth", relative_threshold);
	flg = arghandler.GetArg("-ath", absolute_threshold);
	flg = arghandler.GetArg("-maxiters", max_iters);
	
	printroot("   Relative threshold = %0.3e\n", relative_threshold);
	printroot("   Absolute threshold = %0.3e\n", absolute_threshold);
	printroot("   Maximum iterations = %d\n", max_iters);
	
	flg = arghandler.GetArg("-tau", tau);
	flg = arghandler.GetArg("-alpha", alpha);
	if (flg)
		beta = 1 - alpha;
	
	printroot("   tau = %0.3f\n", tau);
	printroot("   alpha = %0.3f\n", alpha);
	printroot("   beta = %0.3f\n", beta);
	
	//==================================================================================================================
	// Get start and end of local portion of Deltab array
	dStart = dspec.start;
	dEnd = dspec.end;
	dN = dEnd - dStart;
	
	//==================================================================================================================
	// Create x array
	printroot("   Creating x array ...\n");
	x = (usedtype*) calloc(dspec.nFG, sizeof(usedtype));
	
 	//==================================================================================================================
	// Create foreground indices array
	printroot("   Creating foreground indices array ...\n");
	FGindices = (int*) calloc(dspec.N, sizeof(int));
	o = 0;
	for (p = 0; p < dspec.N; p++) {
		if (mask[p]) {
			FGindices[o] = p;
			o++;
		}
	}
	
 	//==================================================================================================================
	// Initialise x
	
	flg = arghandler.GetArg("-x", fname);
	
	if (flg) {
		int err = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
		if (err) {
			printroot("Cannot find file %s\n", fname);
			return false;
		}
		
		MPI_File_read(fptr, &iteration, 1, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_read(fptr, x, dspec.nFG, MPI_USEDTYPE, MPI_STATUS_IGNORE);
		
		MPI_File_close(&fptr);
		
		iteration++;
	}
	else {
		o = 0;
		for (p = 0; p < dspec.N; p++) {
			if (mask[p]) {
				x[o] = chi[p];
				o++;
			}
		}
	}
	
	//==================================================================================================================
	// Create Ax_b array
	printroot("   Creating Ax_b array ...\n");
	Ax_b = (usedtype*) calloc(dEnd - dStart, sizeof(usedtype));
	
	//==================================================================================================================
	// Create AtAx_b array
	printroot("   Creating AtAx_b array ...\n");
	AtAx_b = (usedtype*) calloc(dspec.nFG, sizeof(usedtype));
	
	//==================================================================================================================
	// Create Dx array
	printroot("   Creating Dx array ...\n");
	Dx = (usedtype*) calloc(dEnd - dStart, sizeof(usedtype));
	
	//==================================================================================================================
	// Create DtDx array
	printroot("   Creating DtDx array ...\n");
	DtDx = (usedtype*) calloc(dspec.nFG, sizeof(usedtype));
	
	//==================================================================================================================
	// Create new_x array
	printroot("   Creating new_x array ...\n");
	new_x = (usedtype*) calloc(dspec.nFG, sizeof(usedtype));
	memset(new_x, 0, dspec.nFG * sizeof(usedtype));
	
	//==================================================================================================================
	// Create cylColumns array
	double M = (double)kernel.modelmap.ncyls * dN * sizeof(usedtype);
	cylColumns = (usedtype*) malloc(M);
	if (cylColumns == NULL) {
		printroot("   Not enough memory to pre-calculate cylinder kernels.\n");
		printroot("   Cylinder kernels will be calculated on-the-fly, which will take longer to process.\n");
		PreCalcCylinders = false;
	}
	else {
		
		printroot("   Creating cylColumns array ...\n");
		
		// Initialise cylColumns array
		for (p = dStart; p < dEnd; p++) {
			
			pz = p / dspec.zoffset;
			py = (p - pz * dspec.zoffset) / dspec.yoffset;
			px = p - py * dspec.yoffset - pz * dspec.zoffset;
			
			for (o = 0; o < dspec.nFG; o++) {
				int mix = kernel.modelmap.mask[FGindices[o]];
				if (mix != -1) {
					oz = FGindices[o] / dspec.zoffset;
					oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
					ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;
					
					rx = px - ox + kernel.halfsize;
					ry = py - oy + kernel.halfsize;
					rz = pz - oz + kernel.halfsize;					
					
					cylColumns[mix * dN + p - dStart] = kernel.Get(rx, ry, rz, FGindices[o]);
				}
			}
		}
		
		PreCalcCylinders = true;
	}
	
	//==================================================================================================================
	// Start iterations
	
	printroot("   Starting iterations ...\n");
	//old_rms_diff_x = -1;
	
	tStart = time(NULL);
	tsave = tStart + 3600;
	int dn = -1; //kernel.halfsize - 1;
	int dz = 0; //kernel.halfsize;
	int dp = 1; //kernel.halfsize + 1;
	
	do {
		
		tIterStart1 = time(NULL);
		printroot("      Iteration %d: ", iteration);
		
		// Ax_b = A * x - b
		
		memset(Ax_b, 0, (dEnd-dStart) * sizeof(usedtype));
		memset(Dx, 0, (dEnd-dStart) * sizeof(usedtype));
		for (o = 0; o < dspec.nFG; o++) {
			
			if (x[o]) {
				
				oz = FGindices[o] / dspec.zoffset;
				oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
				ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;
				
				pz = dStart / dspec.zoffset;
				py = (dStart - pz * dspec.zoffset) / dspec.yoffset;
				px = dStart - py * dspec.yoffset - pz * dspec.zoffset - 1;			
				
				ry = py - oy;
				rz = pz - oz;
				
				roffset = (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset;
				
				for (p = dStart; p < dEnd; p++) {
					
					px++;
					if (px == dspec.size[0]) {
						px = 0;
						py++;
						roffset += kernel.yoffset;
						if (py == dspec.size[1]) {
							py = 0;
							pz++;
							rz = pz - oz;
							roffset = kernel.halfsize + kernel.yoffset*kernel.halfsize + (rz+kernel.halfsize)*kernel.zoffset;
						}
						ry = py - oy;
					}				
					rx = px - ox;
					
					if (rx >= -kernel.halfsize && rx <= kernel.halfsize && ry >= -kernel.halfsize && ry <= kernel.halfsize && rz >= -kernel.halfsize && rz <= kernel.halfsize) {
						// Linear system
						int mix = kernel.modelmap.mask[FGindices[o]];
						if (mix == -1) { // spherical kernel
							Ax_b[p - dStart] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * x[o];
						}
						else if (PreCalcCylinders) {
							Ax_b[p - dStart] += cylColumns[mix*dN + p-dStart] * x[o];
						}
						else if (rx == 0 && ry == 0 && rz == 0) {
							Ax_b[p - dStart] += kernel.ctr[mix] * x[o];
						}
						else {
							Ax_b[p - dStart] += kernel.GetCyl(mix, rx, ry, rz) * x[o];
						}
						
						// Laplacian
						if (rx == dn){
							if (ry == dn) {
								if (rz == dz)
									Dx[p-dStart] += D3_96 * x[o];
							}
							else if (ry == dz) {
								if (rz == dn)
									Dx[p-dStart] += D3_96 * x[o];
								else if (rz == dz)
									Dx[p-dStart] += D10_96 * x[o];
								else if (rz == dp)
									Dx[p-dStart] += D3_96 * x[o];
							}
							else if (ry == dp) {
								if (rz == dz)
									Dx[p-dStart] += D3_96 * x[o];
							}
						}
						else if (rx == dz) {
							if (ry == dn) {
								if (rz == dn)
									Dx[p-dStart] += D3_96 * x[o];
								else if (rz == dz)
									Dx[p-dStart] += D10_96 * x[o];
								else if (rz == dp)
									Dx[p-dStart] += D3_96 * x[o];
							}
							else if (ry == dz) {
								if (rz == dn)
									Dx[p-dStart] += D10_96 * x[o];
								else if (rz == dz)
									Dx[p-dStart] -= x[o];
								else if (rz == dp)
									Dx[p-dStart] += D10_96 * x[o];
							}
							else if (ry == dp) {
								if (rz == dn)
									Dx[p-dStart] += D3_96 * x[o];
								else if (rz == dz)
									Dx[p-dStart] += D10_96 * x[o];
								else if (rz == dp)
									Dx[p-dStart] += D3_96 * x[o];
							}
							
						}
						else if (rx == dp) {
							if (ry == dn) {
								if (rz == dz)
									Dx[p-dStart] += D3_96 * x[o];
							}
							else if (ry == dz) {
								if (rz == dn)
									Dx[p-dStart] += D3_96 * x[o];
								else if (rz == dz)
									Dx[p-dStart] += D10_96 * x[o];
								else if (rz == dp)
									Dx[p-dStart] += D3_96 * x[o];
							}
							else if (ry == dp) {
								if (rz == dz)
									Dx[p-dStart] += D3_96 * x[o];
							}						
						}
					}				
				}
			}
		}
		
		for (p = 0; p < dEnd - dStart; p++) {
			Ax_b[p] -= DeltaB[p];
		}
					
//		MPI_File fptr;
//		sprintf(fname, "%s/Dx.bin", outdir);
//		MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
//		MPI_File_write_ordered(fptr, Dx, dN, MPIU_SCALAR, MPI_STATUS_IGNORE);
//		MPI_File_close(&fptr);
//		
//		sprintf(fname, "%s/Ax_b.bin", outdir);
//		MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
//		MPI_File_write_ordered(fptr, Ax_b, dN, MPIU_SCALAR, MPI_STATUS_IGNORE);
//		MPI_File_close(&fptr);
		
		tIterEnd1 = time(NULL);
		tsecs = tIterEnd1 - tIterStart1;
		tmins = floor(tsecs/60);
		tsecs -= tmins*60;		
		printroot("(%ldmin %lds)", tmins, tsecs);
		tIterStart2 = time(NULL);
		
		// AtAx_b = A' * Ax_b
		// DtDx = D' * Dx
		memset(AtAx_b, 0, dspec.nFG*sizeof(usedtype));
		memset(DtDx, 0, dspec.nFG*sizeof(usedtype));
		for (o = 0; o < dspec.nFG; o++) {
			
			oz = FGindices[o] / dspec.zoffset;
			oy = (FGindices[o] - oz * dspec.zoffset) / dspec.yoffset;
			ox = FGindices[o] - oy * dspec.yoffset - oz * dspec.zoffset;
			
			pz = dStart / dspec.zoffset;
			py = (dStart - pz * dspec.zoffset) / dspec.yoffset;
			px = dStart - py * dspec.yoffset - pz * dspec.zoffset - 1;			
			
			ry = py - oy;
			rz = pz - oz;
			
			roffset = (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset;
			
			for (p = dStart; p < dEnd; p++) {
				
				px++;
				if (px == dspec.size[0]) {
					px = 0;
					py++;
					if (py == dspec.size[1]) {
						py = 0;
						pz++;
						rz = pz - oz;
					}
					ry = py - oy;
				}				
				rx = px - ox;
		 
				if (rx >= -kernel.halfsize && rx <= kernel.halfsize && ry >= -kernel.halfsize && ry <= kernel.halfsize && rz >= -kernel.halfsize && rz <= kernel.halfsize) {
		 
					// Linear system
					int mix = kernel.modelmap.mask[FGindices[o]];
					if (mix == -1) { // spherical kernel
						AtAx_b[o] += kernel.skernel[rx+kernel.halfsize + (ry+kernel.halfsize)*kernel.yoffset + (rz+kernel.halfsize)*kernel.zoffset] * Ax_b[p - dStart];
					}
					else if (PreCalcCylinders) {
						AtAx_b[o] += cylColumns[mix*dN + p-dStart] * Ax_b[p - dStart];
					}
					else if (rx == 0 && ry == 0 && rz == 0) {
						AtAx_b[o] += kernel.ctr[mix] * Ax_b[p - dStart];
					}
					else {
						AtAx_b[o] += kernel.GetCyl(mix, rx, ry, rz) * Ax_b[p - dStart];
					}
	
					// Laplacian
					if (rx == dn){
						if (ry == dn) {
							if (rz == dz)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
						else if (ry == dz) {
							if (rz == dn)
								DtDx[o] += D3_96 * Dx[p-dStart];
							else if (rz == dz)
								DtDx[o] += D10_96 * Dx[p-dStart];
							else if (rz == dp)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
						else if (ry == dp) {
							if (rz == dz)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
					}
					else if (rx == dz) {
						if (ry == dn) {
							if (rz == dn)
								DtDx[o] += D3_96 * Dx[p-dStart];
							else if (rz == dz)
								DtDx[o] += D10_96 * Dx[p-dStart];
							else if (rz == dp)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
						else if (ry == dz) {
							if (rz == dn)
								DtDx[o] += D10_96 * Dx[p-dStart];
							else if (rz == dz)
								DtDx[o] -= Dx[p-dStart];
							else if (rz == dp)
								DtDx[o] += D10_96 * Dx[p-dStart];
						}
						else if (ry == dp) {
							if (rz == dn)
								DtDx[o] += D3_96 * Dx[p-dStart];
							else if (rz == dz)
								DtDx[o] += D10_96 * Dx[p-dStart];
							else if (rz == dp)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
						
					}
					else if (rx == dp) {
						if (ry == dn) {
							if (rz == dz)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
						else if (ry == dz) {
							if (rz == dn)
								DtDx[o] += D3_96 * Dx[p-dStart];
							else if (rz == dz)
								DtDx[o] += D10_96 * Dx[p-dStart];
							else if (rz == dp)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}
						else if (ry == dp) {
							if (rz == dz)
								DtDx[o] += D3_96 * Dx[p-dStart];
						}						
				}
				}
			}
		}
		
		tIterEnd2 = time(NULL);
		tsecs = tIterEnd2 - tIterStart2;
		tmins = floor(tsecs/60);
		tsecs -= tmins*60;
		printroot(" (%ldmin %lds)", tmins, tsecs);
		
		// Reduce AtAx_b
		
		//printroot("      reducing AtAx_b ...\n");
		MPI_Allreduce(MPI_IN_PLACE, AtAx_b, dspec.nFG, MPI_USEDTYPE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, DtDx, dspec.nFG, MPI_USEDTYPE, MPI_SUM, MPI_COMM_WORLD);
		
//		if (rank == 0) {
//			sprintf(fname, "%s/DtDx.bin", outdir);
//			MPI_File_open(PETSC_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
//			MPI_File_write_ordered(fptr, DtDx, dspec.nFG, MPIU_SCALAR, MPI_STATUS_IGNORE);
//			MPI_File_close(&fptr);
//
//			sprintf(fname, "%s/AtAx_b.bin", outdir);
//			MPI_File_open(PETSC_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
//			MPI_File_write_ordered(fptr, AtAx_b, dspec.nFG, MPIU_SCALAR, MPI_STATUS_IGNORE);
//			MPI_File_close(&fptr);
//		}
		
		// Calculate new x and rms values
		rms_x = 0;
		rms_diff_x = 0;
		for (o = 0; o < dspec.nFG; o++) {
			new_x[o] = x[o] - alpha * tau * AtAx_b[o] - beta * tau * DtDx[o];
			rms_x += new_x[o] * new_x[o];
			rms_diff_x += (new_x[o] - x[o]) * (new_x[o] - x[o]);
		}
		
		rms_x = sqrt(rms_x / dspec.nFG);
		rms_diff_x = sqrt(rms_diff_x / dspec.nFG);
		
		tIterEnd2 = time(NULL);
		tsecs = tIterEnd2 - tIterStart1;
		tmins = floor(tsecs/60);
		tsecs -= tmins*60;
		
		printroot(" rms_diff_x = %0.3e rms_x = %0.3e rms_diff_x / rms_x  = %0.3e (%ldmin %lds)\n", 
					 rms_diff_x, rms_x, rms_diff_x / rms_x, tmins, tsecs);
		
		//old_rms_diff_x = rms_diff_x;
		
		// Copy results to x
		memcpy(x, new_x, dspec.nFG*sizeof(usedtype));
		
		tEnd = time(NULL);
		if (rank == 0 && (tEnd > tsave || iteration < 2 || iteration == max_iters - 1)) {
			MPI_File fptr;
			sprintf(fname, "%s/x_iter%06d.bin", myout.outdir, iteration);
			MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fptr);
			MPI_File_set_size(fptr, 0);
			MPI_File_write(fptr, &iteration, 1, MPI_INT, MPI_STATUS_IGNORE);
			MPI_File_write(fptr, x, dspec.nFG, MPI_USEDTYPE, MPI_STATUS_IGNORE);
			MPI_File_close(&fptr);
			
			tsave = tEnd + 3600;
		}
		
		iteration++;
		
	} while (rms_diff_x / rms_x > relative_threshold && rms_diff_x > absolute_threshold && iteration < max_iters);
	
	tEnd = time(NULL);
	tsecs = tEnd - tStart;
	tmins = floor(tsecs/60);
	
	if (rms_diff_x / rms_x <= relative_threshold)
		printroot("   Reached relative threshold in %ldmin %lds\n", tmins, tsecs - tmins*60);
	else if (rms_diff_x <= absolute_threshold)
		printroot("   Reached absolute threshold in %ldmin %lds\n", tmins, tsecs - tmins*60);
	else
		printroot("   Reached maximum iterations in %ldmin %lds\n", tmins, tsecs - tmins*60);
	
	printroot("   Completed solve.\n ");
	
	// Create chi vector
	printroot("Compiling chi vector ...\n");
	memset(chi, 0, dspec.N*sizeof(usedtype));
	
	for (o = 0; o < dspec.nFG; o++) {
		chi[FGindices[o]] = x[o];
	}
	
	// Free memory
	printroot("   Finished. Cleaning up ...\n");
	free(x);
	free(FGindices);
	free(Ax_b);
	free(AtAx_b);
	free(new_x);
	free(cylColumns);
	free(Dx);
	free(DtDx);
	free(fname);
	
	return true;
}

//==============================================================================================================================
// LOAD FUNCTIONS

bool loadDeltaB(DataSpec &dspec, usedtype *&DeltaB) {
	// reads in DeltaB to a vector
	
	bool		flg;
	char		*filepath;
	int			err;		
	MPI_File	fptr;		
	MPI_Status	status;
	MPI_Offset	offset, disp;
	
	bool		flgByteSwap;
	double		buf[11];
	double		bmag;
	
	double		*castbuf;
	
	int			ElsPerProc;
	
	printroot("Loading DeltaB ...\n");
	flg = arghandler.GetArg("-DeltaB", filepath);
	if (!flg) {
		printroot("Data file not specified\n");
		return false;
	}
	printroot("   file: %s\n", filepath);
	
	// open file for reading
	err = MPI_File_open(MPI_COMM_SELF, filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
	if (err) {
		printroot("Could not open file");
		return false;
	}
	
	// check endianness
	checkEndianness(fptr, flgByteSwap);
	
	// check version
	checkVersion(fptr, flgByteSwap);
	
	// read in header
	MPI_File_read(fptr, buf, 11, MPI_DOUBLE, &status);
	
	// byte swap header if required
	if (flgByteSwap) byteswap((char*)buf,11,sizeof(double));
	
	// process header into relevant variables
	dspec.N = (int) buf[0];
	dspec.size[0] = (int) buf[1];
	dspec.size[1] = (int) buf[2];
	dspec.size[2] = (int) buf[3];
	dspec.B0 = buf[4];
	memcpy(dspec.bhat, &buf[5], 3*sizeof(double));
	
	// normalise bhat
	bmag = sqrt(dspec.bhat[0]*dspec.bhat[0] + dspec.bhat[1]*dspec.bhat[1] + dspec.bhat[2]*dspec.bhat[2]);
	dspec.bhat[0] /=bmag;
	dspec.bhat[1] /=bmag;
	dspec.bhat[2] /=bmag;
	
	// Create vector and read in DeltaB
	ElsPerProc = dspec.N / size;
	dspec.start = 0;
	for (int p = 0; p < rank; p++) {
		dspec.start += ElsPerProc + ((p < dspec.N - ElsPerProc * size) ? 1 : 0);
	}
	dspec.end = dspec.start + ElsPerProc + ((rank < dspec.N - ElsPerProc * size) ? 1 : 0);
	
	//printall("[%d] start = %d end = %d\n", rank, dspec.start, dspec.end); 
	
	DeltaB = (usedtype*) calloc(dspec.end - dspec.start, sizeof(usedtype));
	
	MPI_File_seek(fptr, dspec.start*sizeof(double), MPI_SEEK_CUR);
	
	if (sizeof(usedtype) == sizeof(double)) {
		MPI_File_read(fptr, DeltaB, dspec.end-dspec.start, MPI_DOUBLE, &status);
		if (flgByteSwap) byteswap((char*)DeltaB, dspec.end-dspec.start, sizeof(double));
	}
	else {
		castbuf = (double*) calloc(dspec.end-dspec.start, sizeof(double));
		MPI_File_read(fptr, castbuf, dspec.end-dspec.start, MPI_DOUBLE, &status);
		if (flgByteSwap) byteswap((char*)castbuf, dspec.end-dspec.start, sizeof(double));
		for (int i = 0; i < dspec.end - dspec.start; i++)
			DeltaB[i] = (float) castbuf[i];
		free(castbuf);
	}
	
	// Close file pointer
	MPI_File_close(&fptr);
	
	// set offsets
	dspec.yoffset = dspec.size[0];
	dspec.zoffset = dspec.size[0] * dspec.size[1];
	
	printroot("   size = %d %d %d\n", dspec.size[0], dspec.size[1], dspec.size[2]);
	printroot("   number of elements = %d\n", dspec.N);
	
	free(filepath);
	
	return true;
}

bool loadMask(DataSpec &dspec, bool *&mask) {
	// Collective function loads mask onto each process
	
	int			err;
	char		*filepath;
	bool		flg;
	MPI_File	fptr;
	double		buf[4];
	bool		flgByteSwap;
	MPI_Offset	offset;
	
	printroot("Loading mask ...\n");
	
	flg = arghandler.GetArg("-mask", filepath);
	if (!flg) {
		printroot("Mask file not specified\n");
		return false;
	}
	printroot("   file: %s\n", filepath);
	
	mask = (bool*) calloc(dspec.N, sizeof(bool));
	
	// open file for reading
	if (rank == 0) {
		err = MPI_File_open(MPI_COMM_SELF, filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fptr);
		if (err) {
			printroot("Could not open mask file"); 
			return false;	
		}
		
		// check endianness
		checkEndianness(fptr, flgByteSwap);
		
		// check version
		checkVersion(fptr, flgByteSwap);
		
		// read in header
		MPI_File_read(fptr, buf, 4, MPI_DOUBLE, MPI_STATUS_IGNORE);
		
		// byte swap header if required
		if (flgByteSwap) byteswap((char*)buf,4,sizeof(double));
		
		// Check parameters match
		if (dspec.N != buf[0]) {printroot("Number of elements in Mask does not match DeltaB"); return false;}
		if (dspec.size[0] != buf[1]) {printroot("Size of first dimension of Mask does not match DeltaB"); return false;}
		if (dspec.size[1] != buf[2]) {printroot("Size of second dimension of Mask does not match DeltaB"); return false;}
		if (dspec.size[2] != buf[3]) {printroot("Size of third dimension of Mask does not match DeltaB"); return false;}
		
		MPI_File_read(fptr, mask, dspec.N, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fptr);
		
	}
	MPI_Bcast(mask, dspec.N, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	dspec.nBG = 0; dspec.nFG = 0;
	for (int i=0; i<dspec.N; i++) {
		if (mask[i])
			dspec.nFG++;
		else
			dspec.nBG++;
	}
	printroot("   number of foreground voxels = %d\n", dspec.nFG);
	printroot("   number of background voxels = %d\n", dspec.nBG);
	printroot("   background percentage = %0.1f%%\n", 100.0*dspec.nBG/dspec.N);
	
	free(filepath);
	
	return true;
}

void checkEndianness(MPI_File &fptr, bool &flgByteSwap) {
	// This function checks the first flag in the data file for whether an endian change has occurred.
	// It compares the endianness against the local process, not the server process.
	// Therefore subsequent reads of data from the file should either checked and converted before sending
	// to other processes, or sent and then checked and converted.
	
	// INPUTS:	MPI_File	&fptr
	// OUTPUTS:	bool		&flgByteSwap
	
	MPI_Status	status;
	unsigned	endianCheckValue;
	
	MPI_File_read(fptr, &endianCheckValue, 1, MPI_UNSIGNED, &status);
	
	if (endianCheckValue == 0x12345678) //Hex value 0x12345678
		flgByteSwap = false;
	else if (endianCheckValue == 0x78563412) //Hex value 0x78563412
		flgByteSwap = true;
	else {
		printroot("Could not recognise endian check flag\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	
}

void checkVersion(MPI_File &fptr, bool flgByteSwap) {
	MPI_Status	status;
	unsigned	versionCheckValue;
	bool		check;
	
	MPI_File_read(fptr, &versionCheckValue, 1, MPI_UNSIGNED, &status);
	
	if (flgByteSwap) byteswap((char*)&versionCheckValue, 1, sizeof(unsigned));
	
	if (versionCheckValue != LOADDATAVERSION) {
		printroot("Data file version is not compatible\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
}

void byteswap(char* buf, int buflength, int dataTypeBytes) {
	// buf			 : array
	// buflength	 : number of elements in array
	// dataTypeBytes : number of bytes in datatype
	char swappee[dataTypeBytes];
	char swapped[dataTypeBytes];
	
	for (int i=0; i<buflength; i++) {
		memcpy(swappee, &buf[i*dataTypeBytes] , dataTypeBytes);
		for (int j=0; j<dataTypeBytes/2; j++) {			
			swapped[j] = swappee[dataTypeBytes-1-j];
			swapped[dataTypeBytes-1-j] = swappee[j];
		}
		memcpy(&buf[i*dataTypeBytes], swapped, dataTypeBytes);
	}	
}

