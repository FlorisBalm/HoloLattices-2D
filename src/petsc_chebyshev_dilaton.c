//
// Created by balm on 29/07/2019.
//
#define N_DOF 9

#define MAX_DERIV 2

#define NEAREST_NEIGHBOURS 2
#define NEAREST_NEIGHBOURS_FULL_STENCIL ((2*NEAREST_NEIGHBOURS)+1)


// These may not all be required as they are included in eachother, but it doesn't hurt.
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscviewerhdf5.h>
#include <petscsys.h>

// Use some of these libraries. Check, they may be obsolete.
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Include the headers for all the bulk and boundary conditions
#include "equations/output/HeaderEOMs_I.h"
#include "equations/output/HeaderEOMs_2_1.h"
#include "equations/output/HeaderEOMs_2_0.h"

#include "coefs/output_inlined/HeaderCoefs_I.h"
#include "coefs/output_inlined/HeaderCoefs_2_1.h"
#include "coefs/output_inlined/HeaderCoefs_2_0.h"

extern PetscErrorCode GetDerivativeCoefficients(
        PetscReal deriv_centre,
        PetscReal* points,
        PetscInt len_points,
        PetscReal ret[MAX_DERIV+1][NEAREST_NEIGHBOURS_FULL_STENCIL]
        );

inline PetscReal ChebZ(PetscReal dz){
    return (1.0 - PetscCosScalar(PETSC_PI*dz))/2.0;
}

//Magic from Fornberg paper
PetscErrorCode GetDerivativeCoefficients(
        PetscReal deriv_centre,
        PetscReal* points,
        PetscInt len_points,
        PetscReal ret[MAX_DERIV+1][NEAREST_NEIGHBOURS_FULL_STENCIL]
        ){
    
//    for (int i = 0; i < MAX_DERIV+1; ++i){
//        for (int j = 0; j < NEAREST_NEIGHBOURS_FULL_STENCIL; ++j){
//            ret[i][j] = 0.0;
//        }
//    }
    PetscInt n = len_points - 1;

    PetscReal c1,c2,c3,c4,c5;
    PetscInt mn;

    c1 = 1.0;
    c4 = points[0] - deriv_centre;
    ret[0][0] = 1.0;

    PetscInt max_deriv = MAX_DERIV;
    PetscInt i,j,k;
    for (i = 1 ; i < n + 1; ++i){
        mn = i < max_deriv ? i : max_deriv;
        c2 = 1.0;
        c5 = c4;
        c4 = points[i] - deriv_centre;
        for (j = 0; j < i; ++j){
            c3 = points[i] - points[j];
            c2 = c2 * c3;
            if (j+1 == i){
                for (k = mn; k > 0; --k){
                    ret[k][i] = c1*(k*ret[k-1][i-1] - c5*ret[k][i-1])/c2;
                }
                ret[0][i] = -c1*c5*ret[0][i-1]/c2;

            }
            for (k = mn; k > 0 ; --k){
                ret[k][j] = (c4*ret[k][j] - k*ret[k-1][j])/c3;

            }
            ret[0][j] = (c4*ret[0][j])/c3;
        }
        c1 = c2;
    }
    return 0;
}

// Representation of all the parameters that are given to the code
typedef struct {
    PetscScalar mu, Q, ax, ay, Gx, Gy, nperiodsx, nperiodsy, phasex, phasey, B,c1;
} Parameters;

// This represents the global info of the simulation. I don't always use this. 
typedef struct {
    DMDAStencilType stencil;
    PetscInt ni, nj, nk;
    PetscInt dof;
    PetscScalar dx, dy, dz;
    PetscScalar x_scale, y_scale;
} GridInfo;

// To allow access to the Distributed Memory object, so we can see what parts are local to this process.
typedef struct {
    Parameters *params;
    GridInfo *grid;
    DM da;                /* distributed array data structure */
} AppCtx;


// At each point, we'll have N_DOF number of fields to solve for.

typedef struct {
    PetscScalar f[N_DOF];
} Field;

extern PetscErrorCode LoadParametersFromFile(AppCtx *user, Parameters* x, const char* fileName);

extern PetscErrorCode FormRNGuess(AppCtx *user, Vec X);
extern PetscErrorCode LoadInitialGuessFromFile(AppCtx* user, Vec X, const char*);

extern PetscErrorCode FormJacobianLocal(DMDALocalInfo *, Field***, Mat, Mat, void*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo *, Field***, Field ***f, void *);

int main(int argc, char **argv) {

    // Declaring basic setup
    SNES snes;
    Vec x;

    Mat J ,j;

    AppCtx *user;
    MPI_Comm comm;

    PetscErrorCode ierr;

    
    // This is the basic petsc loop for error codes. If you want, you can do each call as ierr=....; CHKERRQ(ierr); macro and have 
    // have some debugging that way
    ierr = PetscInitialize(&argc, &argv, (char *) 0, NULL);
    if (ierr) return ierr;


    comm = PETSC_COMM_WORLD;

    DM da;
    SNESCreate(comm, &snes);

    // Default parameters
    GridInfo grid;
    grid.stencil = DMDA_STENCIL_BOX;
    grid.ni = grid.nj = 20;
    grid.nk= 40;
    grid.dof= N_DOF;
    grid.dx = 1.0/(double)grid.ni;
    grid.dy = 1.0/(double)grid.nj;
    grid.dz = 1.0/(double)(grid.nk-1);

    // Number of nearest neighbours on each side.
    PetscInt stencilWidth = 2;
    
    // This sets up the distributed memory distributed array, and decides which BC are used.
    // Note that this can all be changed through the commandline - as we are using "DMSetFromOptions" below.
    DMDACreate3d(
            comm, 
            DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE,
            DMDA_STENCIL_BOX, // Need box for mixed derivatives. 
            grid.ni,grid.nj,grid.nk,
            PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
            N_DOF, stencilWidth,
            NULL, NULL, NULL, &da);
    DMSetFromOptions(da);

    DMSetUp(da);
    SNESSetDM(snes,da);
    SNESSetFromOptions(snes);


    DMCreateGlobalVector(da, &x);

    // These names are not used, but they could be changed or used in output.
//    for (int i = 0 ; i < N_DOF; ++i){
//        DMDASetFieldName(da, i,  sprintf("f%2d",i));
//    }
    // PetscNew makes objects, and also takes care of MPI details
    PetscNew(&user);

    Parameters params;

    user->params = &params;
    user->da = da;
    user->grid = &grid;

    // This allows to circulary access the parameters. 
    DMSetApplicationContext(da, user);

    // Allocate storage for the matrices. They preallocate 
    // according to the finite stencilwidth and type above
    DMCreateMatrix(da, &J);
    DMCreateMatrix(da, &j);


    // Check whether we want to load from file
    PetscBool loadFromFile = PETSC_FALSE;
    char fileName[PETSC_MAX_PATH_LEN];
    PetscOptionsGetString(NULL,NULL, "-load_from_file", fileName, PETSC_MAX_PATH_LEN,&loadFromFile);

    PetscBool useRNExample= PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-guess_rn", &useRNExample, NULL);

    PetscBool outputToFileProvided = PETSC_FALSE;
    char outputFileName[PETSC_MAX_PATH_LEN];
    PetscOptionsGetString(NULL,NULL, "-output_to_file", outputFileName, PETSC_MAX_PATH_LEN,&outputToFileProvided);


    if (loadFromFile){
        LoadParametersFromFile(user, &params,fileName);
        LoadInitialGuessFromFile(user,x, fileName);
    }else{
        //Default Params
        params.ax = params.ay = 0.00;
        params.Gx = params.Gy = 2.0 * PETSC_PI;
        params.Q = 1;
        params.mu = PetscSqrtScalar(3*params.Q*(1+params.Q));
        params.nperiodsx=params.nperiodsy = 1;
        params.phasex=params.phasey=0;
        params.B = 0;
        params.c1 = 0;

        FormRNGuess(user, x);
    }

    // Potentially override parameters specified in file
    PetscOptionsGetScalar(NULL,NULL, "-param_ax", &params.ax,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_ay", &params.ay,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_Gx",&params.Gx,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_Gy",&params.Gy,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_mu",&params.mu,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_Q",&params.Q,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_phasex",&params.phasex,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_phasey",&params.phasey,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_nperiodsx",&params.nperiodsx,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_nperiodsy",&params.nperiodsy,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_B",&params.B,NULL);
    PetscOptionsGetScalar(NULL,NULL, "-param_c1",&params.c1,NULL);


    // Check options. Due to MPI, this prints once
    PetscPrintf(PETSC_COMM_WORLD,"Running using params: mu %f Q %f ax %f ay %f Gx %f Gy %f npx %f npy %f fx %f fy %f B %f c1 %f \n ", params.mu, params.Q, params.ax, params.ay, params.Gx, params.Gy, params.nperiodsx, params.nperiodsy, params.phasex, params.phasey, params.B, params.c1);


    // This sets the jacobian and nonlinear function, the local form is used here 
    // because this is easier to deal with in the function itself
    DMDASNESSetFunctionLocal(da, INSERT_VALUES,
                             (PetscErrorCode(*)(DMDALocalInfo *, void *, void *, void *)) FormFunctionLocal,
                             (void *) user);

    DMDASNESSetJacobianLocal(da, (PetscErrorCode(*)(DMDALocalInfo*, void*, Mat, Mat, void*))FormJacobianLocal,
                             (void*)user);


    // The actual main body of the computation itself
    PetscInt its;
    ierr=SNESSolve(snes, NULL, x);CHKERRQ(ierr);
    SNESGetIterationNumber(snes, &its);
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes, &reason);
    PetscInt exitcode = 0;
    if ( reason < 0){
        exitcode=reason;
        PetscPrintf(PETSC_COMM_WORLD, "Exiting due to error in SNES solution\n");
        goto FINALIZE;
    }
    KSP kspSoln;
    SNESGetKSP(snes,&kspSoln);
    PetscInt kspIts;
    KSPGetIterationNumber(kspSoln, &kspIts);
    KSPType kspType;
    KSPGetType(kspSoln, &kspType);


    if (kspIts < 2 && strcmp(kspType, "preonly")!=0){
        exitcode=-1;
        PetscPrintf(PETSC_COMM_WORLD, "Exiting due to not enough KSP iterations, has to be wrong if PC type is not PREONLY\n");
        goto FINALIZE;
    }




    // From here it is mainly IO:
    // ObjectSetName detemrines what name something will
    // get in the final HDF file

    PetscObjectSetName((PetscObject)x, "result");

    // This needs to be set to a specific size so the output works. This is done with sizeof so
    // that this will work regardless of the exact type of scalar used in the computation.
    // Typically, PetscScalar will be double, but it can be other types.
    Vec paramVec;
    VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE, sizeof(params)/sizeof(params.mu),  &paramVec); 

    for(unsigned long i  =0; i < sizeof(params)/sizeof(params.mu); ++i){
        VecSetValue(paramVec, (PetscInt)i, ((PetscScalar*)&params)[i], INSERT_VALUES);
    } 

    // These steps are needed to actually make sure all processes have the correct data
    VecAssemblyBegin(paramVec);
    VecAssemblyEnd(paramVec);




    // I have opted to use HDF5, this is easy for structured data, can be compressed inside the file
    // itself and it is easy to manipulate using many languages (mathematica,python etc.).
    // The only difficulty is that it does not have a standard way of incorporating complex numbers,
    // so standards may differ, but these are often easily resolved.
    PetscViewer h5view;
    PetscObjectSetName((PetscObject)paramVec, "parameters");

    // Make up a timestamped name in case no valid name has been provided.
    if (outputToFileProvided){
        // No need to do anything here
    } else {
        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0){
            time_t rawtime;
            struct tm *info;
            time(&rawtime);
            info = localtime(&rawtime);
            strftime(outputFileName, PETSC_MAX_PATH_LEN, "data/Output-%F-%H%M%S.h5", info);
            MPI_Bcast(outputFileName, (PetscInt)(strlen(outputFileName)+1), MPI_BYTE, 0, PETSC_COMM_WORLD);
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }

    // Writing to hdf5 works by opening a "viewer" that's the HDF5 file, and then 
    // "viewing" an object through there. 
    // Loading works similarly, but using the VecLoad command
    PetscViewerHDF5Open(PETSC_COMM_WORLD, outputFileName , FILE_MODE_WRITE,&h5view);

    VecView(x, h5view);
    VecView(paramVec, h5view);

    // Default finalisation
//FINALIZE:

    PetscViewerDestroy(&h5view);

FINALIZE:

    VecDestroy(&x);
    VecDestroy(&paramVec);
    SNESDestroy(&snes);
    MatDestroy(&J);
    MatDestroy(&j);
    DMDestroy(&user->da);

    PetscFinalize();
    return exitcode;
}

    
// Basic, load "result" from file
PetscErrorCode LoadInitialGuessFromFile(AppCtx *user, Vec x, const char* fileName){

    PetscViewer h5Viewer;

    PetscViewerHDF5Open(PETSC_COMM_WORLD, fileName, FILE_MODE_READ, &h5Viewer);
    PetscObjectSetName((PetscObject)x, "result");
    VecLoad(x, h5Viewer);
    PetscViewerDestroy(&h5Viewer);

    return 0;
}

// Loading parameters is a bit trickier, as there we need to read into a vector first
// before we can access it and put it into the struct that the rest of the program uses.
PetscErrorCode LoadParametersFromFile(AppCtx *user, Parameters* x, const char* fileName){

    PetscViewer h5Viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, fileName, FILE_MODE_READ, &h5Viewer);

    Vec temp;
    // Create this as a local object on all processors!
    VecCreate(PETSC_COMM_SELF, &temp);
    VecSetType(temp, VECSEQ);
    //By setting the name, you specify which object of the HDF5 to load
    PetscObjectSetName((PetscObject)temp, "parameters");

    VecLoad(temp, h5Viewer);
    PetscScalar* reinterp = (PetscScalar*)x;


    PetscInt indices[12] = {0};
    for (int i = 0; i < 12;++i){
        indices[i] = i;
    }

    VecGetValues(temp, 12, indices, reinterp);

    // Cleanup
    VecDestroy(&temp);
    PetscViewerDestroy(&h5Viewer);

    return 0;
}

PetscErrorCode FormRNGuess(AppCtx *user, Vec X){

    //GridInfo* grid = user->grid;
    PetscInt xs, ys,zs, xm,ym,zm, i,j,k;

    // The corners are the rectangular bits of the domain that each local
    // process has access to. Do not try to write to positions outside of this domain,
    // this will at best not do anything, at worst it will segfault
    // The second triple of arguments are WIDTHS, NOT endpoints
    DMDAGetCorners(user->da, &xs, &ys, &zs, &xm,&ym,&zm);

    // Read the field into a local array
    // NOTE: this gives you read access beyond your local bounds,
    // meaning that you can access the "ghost points" that are needed
    // to do finite difference differentation

    Field ***x;
    DMDAVecGetArray(user->da, X, &x);

    // This is by default the best looping order; z->y->x->field, this is the PETSC indexing
    // order and gives also the best performance due to cache locality

    for(k=zs ; k < zs+zm; ++k){
        for(j = ys; j < ys + ym; ++j){
            for(i = xs; i < xs+xm; ++i){
                    // Standard RN ansatz
                    x[k][j][i].f[ 0] = 1.0; // psi
                    x[k][j][i].f[ 1] = 1.0; // phi
                    x[k][j][i].f[ 2] = 1.0; // qtt
                    x[k][j][i].f[ 3] = 1.0; // qxx

                    x[k][j][i].f[ 4] = 1.0; // qyy

                    x[k][j][i].f[ 5] = 1.0; // qzz
                    x[k][j][i].f[ 6] = 0.0; // qxy
                    x[k][j][i].f[ 7] = 0.0; // qxz
                    x[k][j][i].f[ 8] = 0.0; // qyz
            }
        }
    }

    // Restore the array we got to make sure it gets entered into the vector
    DMDAVecRestoreArray(user->da, X, &x);

    return 0;
}

// In the function signature, x is the current guess for the field, and f is the result
// of evaluating the equations of motion at each grid point, so we need to write to f

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, Field ***x, Field ***f, void *ptr) {

    AppCtx *user = (AppCtx *) ptr;
    Parameters *param = user->params;
    //GridInfo *grid = user->grid;

    PetscScalar mu, Q, Gx, Gy, ax, ay, nperiodsx, nperiodsy, phasex, phasey,B,c1;

    mu = param->mu;
    Q = param->Q;
    Gx = param->Gx;
    Gy = param->Gy;
    ax = param->ax;
    ay = param->ay;
    nperiodsx = param->nperiodsx;
    nperiodsy = param->nperiodsy;
    phasex = param->phasex;
    phasey = param->phasey;
    B = param->B;
    c1 = param->c1;

    PetscScalar x_scale,y_scale;

    x_scale = 2*PETSC_PI*nperiodsx/(Gx);
    y_scale = 2*PETSC_PI*nperiodsy/(Gy);



    PetscInt i, j, k, mx, my, mz, klim;
    PetscInt is, js, ks, ie, je, ke;

    PetscReal x_coord, y_coord, z_coord;

    mx = info->mx;
    my = info->my;
    mz = info->mz;

    PetscReal dx,dy,dz;
    dx =x_scale* 1.0/((PetscScalar)mx);
    dy =y_scale* 1.0/((PetscScalar)my);
    dz = 1.0/((PetscScalar)mz-1.0);

    klim = mz - 1;

    is = info->xs;
    js = info->ys;
    ks = info->zs;
    ie = is + info->xm;
    je = js + info->ym;
    ke = ks + info->zm;

    PetscInt nDOF = N_DOF;

    PetscScalar derivs[N_DOF][10]= {{0}};

    for (k = ks; k < ke; ++k) {

        PetscReal derivativeCoefsZ [MAX_DERIV+1][NEAREST_NEIGHBOURS_FULL_STENCIL] = {0};
        z_coord = ChebZ(((double)k)*dz);

        if (k == 0){
            PetscReal points[NEAREST_NEIGHBOURS+1] = {0};

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+1; ++idx){
                points[idx] = ChebZ(idx*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+1, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    x_coord = (double)i * dx;



                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] =0. +1.*x[k][j][i].f[n];
                        derivs[n][1] =0. +(0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] =0. +(0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] =0. +1.*derivativeCoefsZ[1][0]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[2 + k][j][i].f[n];
                        derivs[n][4] =0. -(0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] =0. -(0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] =0. +1.*derivativeCoefsZ[2][0]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[2 + k][j][i].f[n];
                        derivs[n][7] =0. +(0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] =0. +(0.08333333333333333*derivativeCoefsZ[1][0]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] =0. +(0.08333333333333333*derivativeCoefsZ[1][0]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][2 + j][i].f[n])/dy;
                    }


                    f[k][j][i].f[ 0] = eom_2_0_0(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 1] = eom_2_0_1(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 2] = eom_2_0_2(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 3] = eom_2_0_3(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 4] = eom_2_0_4(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 5] = eom_2_0_5(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 6] = eom_2_0_6(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 7] = eom_2_0_7(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 8] = eom_2_0_8(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 9] = eom_2_0_9(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[10] = eom_2_0_10(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[11] = eom_2_0_11(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[12] = eom_2_0_12(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[13] = eom_2_0_13(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[14] = eom_2_0_14(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                }
            }
        } else if (k==1){  // Z = 0 + dz
            PetscReal points[NEAREST_NEIGHBOURS+2] = {0};

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+2; ++idx){
                points[idx] = ChebZ(idx*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+2, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    x_coord = (double)i * dx;

                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][3]*x[2 + k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][3]*x[2 + k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][2 + j][i].f[n])/dy;
                    }

                    f[k][j][i].f[ 0] = eom_0(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 1] = eom_1(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 2] = eom_2(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 3] = eom_3(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 4] = eom_4(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 5] = eom_5(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 6] = eom_6(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 7] = eom_7(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 8] = eom_8(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 9] = eom_9(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 10] = eom_10(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 11] = eom_11(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 12] = eom_12(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 13] = eom_13(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 14] = eom_14(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                }
            }
        } else if (k == klim-1){ // z = 1-dz
            PetscReal points[NEAREST_NEIGHBOURS+2] = {0};

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+2; ++idx){
                points[idx] = ChebZ((klim - (NEAREST_NEIGHBOURS+1) + idx)*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+2, derivativeCoefsZ);


            for (j = js; j < je; ++j) {
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    x_coord = (double)i * dx;
                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][3]*x[1 + k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][3]*x[1 + k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][2 + j][i].f[n])/dy;
                    }

                    f[k][j][i].f[ 0] = eom_0(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 1] = eom_1(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 2] = eom_2(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 3] = eom_3(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 4] = eom_4(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 5] = eom_5(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 6] = eom_6(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 7] = eom_7(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 8] = eom_8(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 9] = eom_9(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 10] = eom_10(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 11] = eom_11(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 12] = eom_12(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 13] = eom_13(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 14] = eom_14(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                }
            }
        } else if (k == klim ) { // z = 1
            PetscReal points[NEAREST_NEIGHBOURS+1] = {0};

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+1; ++idx){
                points[idx] = ChebZ((klim - (NEAREST_NEIGHBOURS) + idx)*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+1, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    x_coord = (double)i * dx;

                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][2 + j][i].f[n])/dy;
                    }

                    f[k][j][i].f[ 0] = eom_2_1_0(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 1] = eom_2_1_1(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 2] = eom_2_1_2(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 3] = eom_2_1_3(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 4] = eom_2_1_4(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 5] = eom_2_1_5(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 6] = eom_2_1_6(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 7] = eom_2_1_7(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 8] = eom_2_1_8(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 9] = eom_2_1_9(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 10] = eom_2_1_10(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 11] = eom_2_1_11(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 12] = eom_2_1_12(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 13] = eom_2_1_13(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 14] = eom_2_1_14(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                }
            }
        } else {
            PetscReal points[NEAREST_NEIGHBOURS_FULL_STENCIL] = {0};

            for (int idx = 0; idx < NEAREST_NEIGHBOURS_FULL_STENCIL; ++idx){
                points[idx] = ChebZ(((k - NEAREST_NEIGHBOURS) + idx)*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS_FULL_STENCIL, derivativeCoefsZ);
            for (j = js; j < je; ++j) {
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    x_coord = (double)i * dx;

                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][3]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][4]*x[2 + k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][3]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][4]*x[2 + k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][2 + j][i].f[n])/dy;                    
                    }
                    
                    f[k][j][i].f[ 0] = eom_0(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 1] = eom_1(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 2] = eom_2(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 3] = eom_3(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 4] = eom_4(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 5] = eom_5(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 6] = eom_6(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 7] = eom_7(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                    f[k][j][i].f[ 8] = eom_8(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 9] = eom_9(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 10] = eom_10(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 11] = eom_11(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 12] = eom_12(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 13] = eom_13(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
//                    f[k][j][i].f[ 14] = eom_14(derivs, mu,Q, Gx, Gy, ax, ay, x_coord, y_coord, z_coord, nperiodsx,nperiodsy,phasex, phasey, B, c1);
                }


            }
        }
        memset(derivativeCoefsZ, 0, sizeof derivativeCoefsZ);
    }
    // No "restorearray or anything needed here
    return 0;
}

// Now the resulting matrices J and jac are where we need to write to. J is the Jacobian, jac would be
// the preconditioner matrix. We are not using that, instead we are letting Petsc figure that out.
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, Field ***x,Mat J, Mat jac, void *ptr) {


    AppCtx *user = (AppCtx *) ptr;
    Parameters *param = user->params;
    //GridInfo *grid = user->grid;

    PetscScalar mu, Q, Gx, Gy, ax, ay, nperiodsx, nperiodsy, phasex, phasey,B, c1;

    // This is needed to be able to write to the matrix
    MatSetOption(jac, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    MatSetOption(jac, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

    mu = param->mu;
    Q = param->Q;
    Gx = param->Gx;
    Gy = param->Gy;
    ax = param->ax;
    ay = param->ay;
    nperiodsx = param->nperiodsx;
    nperiodsy = param->nperiodsy;
    phasex = param->phasex;
    phasey = param->phasey;
    B = param->B;
    c1 = param->c1;



    PetscInt i, j, k, mx, my, mz, klim;
    PetscInt is, js, ks, ie, je, ke;

    double x_coord, y_coord, z_coord;
    PetscInt nDOF= N_DOF;

    mx = info->mx;
    my = info->my;
    mz = info->mz;


    PetscScalar x_scale,y_scale;

    x_scale = 2*PETSC_PI*nperiodsx/(Gx);
    y_scale = 2*PETSC_PI*nperiodsy/(Gy);

    PetscScalar dx,dy,dz;
    dx = x_scale*1.0/((PetscScalar)mx);
    dy = y_scale*1.0/((PetscScalar)my);
    dz = 1.0/((PetscScalar)mz-1.0);

    klim = mz - 1;

    is = info->xs;
    js = info->ys;
    ks = info->zs;
    ie = is + info->xm;
    je = js + info->ym;
    ke = ks + info->zm;

    MatStencil col[61], row;
    PetscReal val[61];

    klim = mz - 1;

    PetscScalar derivs[N_DOF][10]= {{0}};
    PetscScalar points[NEAREST_NEIGHBOURS_FULL_STENCIL] = {0};

    PetscScalar cVal[10] = {0};

    // Here we do much the same thing as in the equations of motion: We compute the derivatives, then we 
    // use those as well as the current guess function to evaluate the jacobian.
    PetscReal derivativeCoefsZ [MAX_DERIV+1][NEAREST_NEIGHBOURS_FULL_STENCIL] = {{0}};
    for (k = ks; k < ke; ++k) {

        row.k = k;
        z_coord = ChebZ(k*dz);
        if (k == 0){

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+1; ++idx){
                points[idx] = ChebZ(idx*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+1, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                row.j = j;
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    row.i = i;
                    x_coord = (double)i * dx;



                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. +1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. +(0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. +(0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. +1.*derivativeCoefsZ[1][0]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[2 + k][j][i].f[n];
                        derivs[n][4] = 0. -(0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. -(0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. +1.*derivativeCoefsZ[2][0]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[2 + k][j][i].f[n];
                        derivs[n][7] = 0. +(0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. +(0.08333333333333333*derivativeCoefsZ[1][0]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. +(0.08333333333333333*derivativeCoefsZ[1][0]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[2 + k][2 + j][i].f[n])/dy;
                    }


                    #include "HeadersSecondOrder/ExpressionCoefficients_2_0_0.h"
                }
            }
        } else if (k==1){  // Z = 0 + dz

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+2; ++idx){
                points[idx] = ChebZ(idx*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+2, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                row.j = j;
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    row.i = i;
                    x_coord = (double)i * dx;

                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][3]*x[2 + k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][3]*x[2 + k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][3]*x[2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][3]*x[2 + k][2 + j][i].f[n])/dy;
                    }

                    #include "HeadersSecondOrder/ExpressionCoefficients_2_0_1.h"
                }
            }
        } else if (k == klim-1){ // z = 1-dz

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+2; ++idx){
                points[idx] = ChebZ((k - (NEAREST_NEIGHBOURS) + idx)*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+2, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                row.j = j;
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    row.i = i;
                    x_coord = (double)i * dx;
                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][3]*x[1 + k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][3]*x[1 + k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][2 + j][i].f[n])/dy;
                    }

                    #include "HeadersSecondOrder/ExpressionCoefficients_2_1_1.h"
                }
            }
        } else if (k == klim ) { // z = 1

            for (int idx = 0; idx < NEAREST_NEIGHBOURS+1; ++idx){
                points[idx] = ChebZ((k- (NEAREST_NEIGHBOURS) + idx)*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS+1, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                row.j = j;
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    row.i = i;
                    x_coord = (double)i * dx;
                    
                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = 0. + 1.*x[k][j][i].f[n];
                        derivs[n][1] = 0. + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = 0. + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = 0. + 1.*derivativeCoefsZ[1][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[k][j][i].f[n];
                        derivs[n][4] = 0. - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = 0. - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = 0. + 1.*derivativeCoefsZ[2][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[k][j][i].f[n];
                        derivs[n][7] = 0. + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][2 + i].f[n])/dx;
                        derivs[n][9] = 0. + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][2 + j][i].f[n])/dy;
                    }

                    #include "HeadersSecondOrder/ExpressionCoefficients_2_1_0.h"
                }
            }
        } else {

            for (int idx = 0; idx < NEAREST_NEIGHBOURS_FULL_STENCIL; ++idx){
                points[idx] = ChebZ(((k - NEAREST_NEIGHBOURS) + idx)*dz);
            }

            GetDerivativeCoefficients(z_coord, points, NEAREST_NEIGHBOURS_FULL_STENCIL, derivativeCoefsZ);

            for (j = js; j < je; ++j) {
                row.j = j;
                y_coord = (double)j * dy;
                for (i = is; i < ie; ++i) {
                    row.i = i;
                    x_coord = (double)i * dx;


                    for(int n = 0; n < nDOF; ++n){
                        derivs[n][0] = + 1.*x[k][j][i].f[n];
                        derivs[n][1] = + (0.08333333333333333*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*x[k][j][2 + i].f[n])/dx;
                        derivs[n][2] = + (0.08333333333333333*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*x[k][2 + j][i].f[n])/dy;
                        derivs[n][3] = + 1.*derivativeCoefsZ[1][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[1][3]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[1][4]*x[2 + k][j][i].f[n];
                        derivs[n][4] = - (0.08333333333333333*x[k][j][-2 + i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][-1 + i].f[n])/(dx*dx) - (2.5*x[k][j][i].f[n])/(dx*dx) + (1.3333333333333333*x[k][j][1 + i].f[n])/(dx*dx) - (0.08333333333333333*x[k][j][2 + i].f[n])/(dx*dx);
                        derivs[n][5] = - (0.08333333333333333*x[k][-2 + j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][-1 + j][i].f[n])/(dy*dy) - (2.5*x[k][j][i].f[n])/(dy*dy) + (1.3333333333333333*x[k][1 + j][i].f[n])/(dy*dy) - (0.08333333333333333*x[k][2 + j][i].f[n])/(dy*dy);
                        derivs[n][6] = + 1.*derivativeCoefsZ[2][0]*x[-2 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][1]*x[-1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][2]*x[k][j][i].f[n] + 1.*derivativeCoefsZ[2][3]*x[1 + k][j][i].f[n] + 1.*derivativeCoefsZ[2][4]*x[2 + k][j][i].f[n];
                        derivs[n][7] = + (0.006944444444444444*x[k][-2 + j][-2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-2 + j][-1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-2 + j][1 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][-2 + j][2 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][-1 + j][-2 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][-1 + j][-1 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][-1 + j][1 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][-1 + j][2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][1 + j][-2 + i].f[n])/(dx*dy) - (0.4444444444444444*x[k][1 + j][-1 + i].f[n])/(dx*dy) + (0.4444444444444444*x[k][1 + j][1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][1 + j][2 + i].f[n])/(dx*dy) - (0.006944444444444444*x[k][2 + j][-2 + i].f[n])/(dx*dy) + (0.05555555555555555*x[k][2 + j][-1 + i].f[n])/(dx*dy) - (0.05555555555555555*x[k][2 + j][1 + i].f[n])/(dx*dy) + (0.006944444444444444*x[k][2 + j][2 + i].f[n])/(dx*dy);
                        derivs[n][8] = + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][j][2 + i].f[n])/dx + (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][j][-2 + i].f[n])/dx - (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][j][-1 + i].f[n])/dx + (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][j][1 + i].f[n])/dx - (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][j][2 + i].f[n])/dx;
                        derivs[n][9] = + (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][0]*x[-2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][0]*x[-2 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][1]*x[-1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][1]*x[-1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][2]*x[k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][2]*x[k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][3]*x[1 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][3]*x[1 + k][2 + j][i].f[n])/dy + (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][-2 + j][i].f[n])/dy - (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][-1 + j][i].f[n])/dy + (0.6666666666666666*derivativeCoefsZ[1][4]*x[2 + k][1 + j][i].f[n])/dy - (0.08333333333333333*derivativeCoefsZ[1][4]*x[2 + k][2 + j][i].f[n])/dy;                    
                    }


                    #include "HeadersSecondOrder/ExpressionsCoefficientsInternal.h"
                }


            }
        }

        memset(derivativeCoefsZ, 0.0, sizeof derivativeCoefsZ);
        memset(points, 0.0, sizeof points);
    }
    MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
    if(jac != J){
        MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    }
    return 0;

}
