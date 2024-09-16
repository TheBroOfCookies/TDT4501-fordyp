#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <omp.h>

#include <mpi.h>

#include "main.h"
#include "kernels.h"

#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

/* Function prototypes */
void options ( int argc, char **argv );
void config_print ( double elapsed_time );
void time_step ( int_t iteration );
void border_exchange_all ( void );
void border_exchange ( real_t *startS, real_t *startR, real_t *endS, real_t *endR );
void insert_source ( int_t ts, source_t type );

/* Weight coefficients, inner diff. loop */
#define HALF 8
const real_t W[HALF] = {1.2627, -0.1312, 0.0412, -0.0170, 0.0076, -0.0034, 0.0014, -0.0005};                 

/* Simulation parameters */
int_t
    Nx=0, Ny=0, Nz=0, Nt=0, st=0;

int_t
    tNx=0, tNy=0, tNz=0;

int
    rank, size,
    rcv0_rank, rcv0_min,
    rcv1_rank, rcv1_min;

source_t
    source_type = STRESS;

/* Discretization, grid */
const real_t
    kDx = 10.0, kDy = 10.0, kDz = 10.0, kDt = 0.001;

/* Discretization, model parameters */
const real_t
    kRho = 1.0e3, kVp = 2.2e3, kVs = 1.0e3;

/* Discretization, source */
const real_t
    kF0 = 5.0,
    kT0 = 0.3;
direction_t source_direction = XDIR;    // Used in insert_source (only)

int_t source_x, source_y, source_z; // Depend on input domain size
real_t *source; // Depends on interation count

/* Problem wrapper structures */
model_t *model __attribute__((aligned(64))) = NULL;
mesh_t *mesh = NULL;
receiver_t *recv = NULL;


/* Main setup & loop */
int
main ( int argc, char **argv )
{
    /* Command line arguments: resolution x,y,z, timesteps, sourcetype
     * sourcetype is magic constant determining force initialization
     */
    

    // MPI setup

    MPI_Init (&argc, &argv);

    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );

    int_t BroadcastBuffer[5];
    if(rank == 0){
        if (ceil(log2(size)) != floor(log2(size))) {
            fprintf( stderr, "\nInvalid number of mpi ranks. Please choose a power of 2\n\n" );
            exit(1);
        }
        options ( argc, argv );
        BroadcastBuffer[0] = tNx;
        BroadcastBuffer[1] = tNy;
        BroadcastBuffer[2] = tNz;
        BroadcastBuffer[3] = Nt;
        BroadcastBuffer[4] = st;
        printf("Halo size = %d\n",HALO);
        
    }

    MPI_Bcast( BroadcastBuffer, 5, MPI_INT64_T, 0, MPI_COMM_WORLD);
    tNx = BroadcastBuffer[0];
    tNy = BroadcastBuffer[1];
    tNz = BroadcastBuffer[2];
    Nt = BroadcastBuffer[3];

    switch (BroadcastBuffer[4]) {
        case 1:  source_type=STRESS; break;
        case 2: source_type=F_MONOP; break;
        case 3:   source_type=F_DIP; break;
    }

    Nx = tNx;   //simple x decomposition
    Ny = tNy;
    Nz = tNz/size;

    
    

    

    /* Source coordinates: initialized here b/c not constant with
     * different domain sizes
     */
    source_x = HNx/2, source_y = HNy/2, source_z = HNz/2;
    source = malloc ( Nt * sizeof(real_t) );
    for ( int_t t=0; t<Nt; t++ )
    {
        real_t arg = pow ( M_PI * kF0 * (kDt * t - kT0), 2.0 );
        source[t] = 1.0e4 * (2.0 * arg - 1.0) * exp(-arg);
    }

#ifdef SAVE_RECEIVERS
    recv = malloc(sizeof(receiver_t));
    int_t nrecvs;
    switch ( tNz )
    {
        case 64:   nrecvs = 8;  break;
        case 128:  nrecvs = 16; break;
        case 256:  nrecvs = 24; break;
        case 512:  nrecvs = 32; break;
        default:
            nrecvs = 0;
            fprintf ( stderr,
                "Warning: receiver setup not supported for Nz=%ld\n", Nz
            );
            break;
    }
    /* START Parallel section for saving recievers*/
    if(size > 1) {

        if (tNz != 64) {
            fprintf( stderr, "\nTotal number of gridpoints in Z-direction must be 64 to support SAVE_RECIEVERS and MPI parallelizaiton.\n" );
            exit(1);
        }
        
        for ( int_t r=0; r<size; r++ ) {
            int maxz = (r+1)*Nz-1;
            int minz = r*Nz;
            if(maxz >= 42 && minz <= 42) {
                
                rcv1_rank = r;
                rcv1_min = minz;

            }
            if(maxz >= 12 && minz <= 12) {
                
                rcv0_rank = r;
                rcv0_min = minz;
            }
        }

        nrecvs = nrecvs/2;
        if (rank == rcv0_rank || rank == rcv1_rank){
            printf("Rank %d deisgnated reciever\n", rank);
            receiver_init ( recv, nrecvs, false, true, true, true );
            printf("Init compelete Rank %d\n", rank);
            receiver_setup ( recv );
            printf("Setup compelete Rank %d\n", rank);
        } /* END Parallel section for saving recievers*/
    } else {
        receiver_init ( recv, nrecvs, false, true, true, true );
        receiver_setup ( recv );
    }
#endif

    model = malloc (sizeof(model_t));
    mesh = malloc (sizeof(mesh_t));
    model_init ( model );
    mesh_init ( mesh );
    model_set_uniform ( model );
    

    struct timeval t_start, t_end;
    gettimeofday ( &t_start, NULL );
    for ( int_t t=0; t<Nt; t++ ) {
        if(size > 1) border_exchange_all();
        time_step ( t );
    }
    gettimeofday ( &t_end, NULL );

#ifdef SAVE_RECEIVERS
    if(size > 1){
        if (rcv1_rank == rank || rcv0_rank == rank){
            receiver_write_MPI ( recv );
        }
    } else {
        receiver_write ( recv );
    }
#endif

    if (rank==0){config_print (WALLTIME(t_end) - WALLTIME(t_start) );}

    mesh_destroy ( mesh );
    model_destroy ( model );
    free ( mesh );
    free ( model );
    free ( source );
#ifdef SAVE_RECEIVERS
    receiver_destroy ( recv );
    free ( recv );
#endif
    MPI_Finalize ();
    exit ( EXIT_SUCCESS );
}


void
insert_source ( int_t ts, source_t type )
{
    // Convenience aliases
    int_t z=source_z, y=source_y, x=source_x;
    real_t s = source[ts];
    // Determine source type, act accordingly
    switch ( type )
    {
        case STRESS:
            SXX(z,y,x) += s * kDt; //printf("SXX %f, rank %d\n", SXX(z,y,x), rank);
            SYY(z,y,x) += s * kDt;
            SZZ(z,y,x) += s * kDt;
            break;
        case F_MONOP:   // Monopole force
            switch ( source_direction )
            {
                case XDIR: VX(z,y,x) = s * kDt * (1.0 / RHO(z,y,x)); break;
                case YDIR: VY(z,y,x) = s * kDt * (1.0 / RHO(z,y,x)); break;
                case ZDIR: VZ(z,y,x) = s * kDt * (1.0 / RHO(z,y,x)); break;
            }
            break;
        case F_DIP:   // Dipole force
            VX(z,y,x+1) += s * kDt * (1.0 / RHO(z,y,x)) / (2.0*kDx);
            VX(z,y,x-1) -= s * kDt * (1.0 / RHO(z,y,x)) / (2.0*kDx);
            VY(z,y+1,x) += s * kDt * (1.0 / RHO(z,y,x)) / (2.0*kDy);
            VY(z,y-1,x) -= s * kDt * (1.0 / RHO(z,y,x)) / (2.0*kDy);
            VZ(z+1,y,x) += s * kDt * (1.0 / RHO(z,y,x)) / (2.0*kDz);
            VZ(z-1,y,x) -= s * kDt * (1.0 / RHO(z,y,x)) / (2.0*kDz);
            break;
    }
}


void
time_step ( int_t ts )
{
    //if ( ts % 100 == 0 )
    //    printf ( "Iter %ld\n", ts );

    /* Insert source: stress, monopole force, or dipole force */
    insert_source ( ts, source_type );    

    /* Collect samples, if enabled */
#ifdef SAVE_RECEIVERS
    receiver_save ( ts, recv );
#endif

    /* Call pattern: to, from, scale
     * scale is 1.0 / kD(dim)
     * scale param is really given by name of differentiator
     */

    /* Time step: */
    /* Vx:
     *    dx forward  (del1,sxx)
     *    dz backward (del2,sxz)
     *    dy backward (del3,sxy)
     *    compute vx
     */
    FDD(DEL1,SXX,XF);
    FDD(DEL2,SXZ,ZB);
    FDD(DEL3,SXY,YB);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz; k++)
        for (int_t j = 0; j < HNy; j++)
            for (int_t i = 0; i < HNx-1; i++)
                VX(k,j,i) += kDt *
                    ( 2.0 / (RHO(k,j,i) + RHO(k,j,i+1)) ) *
                    ( DEL1(k,j,i) + DEL2(k,j,i) + DEL3(k,j,i) );


    /* Vy:
     *    dy forward  (del1,syy)
     *    dz backward (del2,syz)
     *    dx backward (del3,sxy)
     *    compute vy
     */
    FDD(DEL1,SYY,YF);
    FDD(DEL2,SYZ,ZB);
    FDD(DEL3,SXY,XB);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz; k++)
        for (int_t j = 0; j < HNy-1; j++)
            for (int_t i = 0; i < HNx; i++)
                VY(k,j,i) += kDt *
                    ( 2.0 / (RHO(k,j,i) + RHO(k,j+1,i)) ) *
                    ( DEL1(k,j,i) + DEL2(k,j,i) + DEL3(k,j,i) );

    /* Vz:
     *    dz forward  (del1,szz)
     *    dx backward (del2,sxz)
     *    dy backward (del3,syz)
     *    compute vz
     */
    FDD(DEL1,SZZ,ZF);
    FDD(DEL2,SXZ,XB);
    FDD(DEL3,SYZ,YB);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz-1; k++)
        for (int_t j = 0; j < HNy; j++)
            for (int_t i = 0; i < HNx; i++)
                VZ(k,j,i) += kDt *
                    ( 2.0 / (RHO(k,j,i) + RHO(k+1,j,i)) ) *
                    ( DEL1(k,j,i) + DEL2(k,j,i) + DEL3(k,j,i) );

    /* Sxx, Syy, Szz: (stress?)
     *    dz backward (del1,vz)
     *    dx backward (del2,vx)
     *    dy backward (del3,vy)
     *    compute sxx_syy_szz
     */
    FDD(DEL1,VZ,ZB);
    FDD(DEL2,VX,XB);
    FDD(DEL3,VY,YB);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz-1; k++)
        for (int_t j = 0; j < HNy; j++)
            for (int_t i = 0; i < HNx; i++)
            {
                SXX(k,j,i) += kDt * (
                      ( LAMBDA(k,j,i) + 2.0*MU(k,j,i) ) * DEL2(k,j,i)
                    + ( LAMBDA(k,j,i) * (DEL1(k,j,i)+DEL3(k,j,i)) )
                );
                SYY(k,j,i) += kDt * (
                      ( LAMBDA(k,j,i) + 2.0*MU(k,j,i) ) * DEL3(k,j,i)
                    + ( LAMBDA(k,j,i) * (DEL1(k,j,i)+DEL2(k,j,i)) )
                );
                SZZ(k,j,i) += kDt * (
                      ( LAMBDA(k,j,i) + 2.0*MU(k,j,i) ) * DEL1(k,j,i)
                    + ( LAMBDA(k,j,i) * (DEL2(k,j,i)+DEL3(k,j,i)) )
                );
            }

    /* Sxy:
     *    dy forward (del1,vx)
     *    dx forward (del2,vy)
     *    compute sxy
     */
    FDD(DEL1,VX,YF);
    FDD(DEL2,VY,XF);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz; k++)
        for (int_t j = 0; j < HNy-1; j++)
            for (int_t i = 0; i < HNx-1; i++)
                SXY(k,j,i) += kDt * (
                    (MU(k,j,i) + MU(k,j,i+1) + MU(k,j+1,i) + MU(k,j+1,i+1))
                    * 0.25 * (DEL1(k,j,i)+DEL2(k,j,i))
                );


    /* Syz:
     *    dz forward (del1,vy)
     *    dy forward (del2,vz)
     *    compute syz
     */
    FDD(DEL1,VY,ZF);
    FDD(DEL2,VZ,YF);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz-1; k++)
        for (int_t j = 0; j < HNy-1; j++)
            for (int_t i = 0; i < HNx; i++)
                SYZ(k,j,i) += kDt * (
                    (MU(k,j,i) + MU(k,j+1,i) + MU(k+1,j,i) + MU(k+1,j+1,i))
                    * 0.25 * (DEL1(k,j,i)+DEL2(k,j,i))
                );

    /* Sxz:
     *    dx forward (del1,vz)
     *    dz forward (del2,vx)
     *    comute sxz
     */
    FDD(DEL1,VZ,XF);
    FDD(DEL2,VX,ZF);
    #pragma omp parallel for
    for (int_t k = 0; k < HNz-1; k++)
        for (int_t j = 0; j < HNy; j++)
            for (int_t i = 0; i < HNx-1; i++)
                SXZ(k,j,i) += kDt * (
                    (MU(k,j,i) + MU(k,j,i+1) + MU(k+1,j,i) + MU(k+1,j,i+1))
                    * 0.25 * (DEL1(k,j,i)+DEL2(k,j,i))
                );
}
void
border_exchange_all ( void )
{
    border_exchange(&VX(HNz-2,0,0), &VX(HNz-1,0,0), &VX(1,0,0), &VX(0,0,0));
    border_exchange(&VY(HNz-2,0,0), &VY(HNz-1,0,0), &VY(1,0,0), &VY(0,0,0));
    border_exchange(&VZ(HNz-2,0,0), &VZ(HNz-1,0,0), &VZ(1,0,0), &VZ(0,0,0));

    border_exchange(&SXX(HNz-2,0,0), &SXX(HNz-1,0,0), &SXX(1,0,0), &SXX(0,0,0));
    border_exchange(&SYY(HNz-2,0,0), &SYY(HNz-1,0,0), &SYY(1,0,0), &SYY(0,0,0));
    border_exchange(&SZZ(HNz-2,0,0), &SZZ(HNz-1,0,0), &SZZ(1,0,0), &SZZ(0,0,0));

    border_exchange(&SXY(HNz-2,0,0), &SXY(HNz-1,0,0), &SXY(1,0,0), &SXY(0,0,0));
    border_exchange(&SYZ(HNz-2,0,0), &SYZ(HNz-1,0,0), &SYZ(1,0,0), &SYZ(0,0,0));
    border_exchange(&SXZ(HNz-2,0,0), &SXZ(HNz-1,0,0), &SXZ(1,0,0), &SXZ(0,0,0));

}
void
border_exchange ( real_t *upS, real_t *upR, real_t *downS, real_t *downR )
{
    // Exchange border for each valid array
    int count = HNx*HNy;
    if(rank != size - 1){  //Rank 0
        //printf("Rank %d sending...\n", rank);
        MPI_Send(upS, count, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
        //printf("Rank %d receiving...\n", rank);
        MPI_Recv(upR, count, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Sendrecv(&VX(HNz-2,0,0), count, MPI_FLOAT, rank + 1, 0, 
        //             &VX(HNz-1,0,0), count, MPI_FLOAT, rank + 1, 0, 
        //             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } 
    if (rank != 0) { //Rank size-1
        //printf("Rank %d receiving...\n", rank);
        MPI_Recv(downR, count, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("Rank %d sending...\n", rank);
        MPI_Send(downS, count, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
        //MPI_Sendrecv(&VX(1,0,0), count, MPI_FLOAT, rank - 1, 0, 
        //             &VX(0,0,0), count, MPI_FLOAT, rank - 1, 0, 
        //             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    //printf("Rank %d Exchange complete\n", rank);
}



void
model_init ( model_t *model )
{
    /* Initialize everything to NULL */
    model->input =
    model->vp = model->vs = 
    model->rho = model->lambda = model->mu = 
    model->l = model->m = NULL;

    /* Allocate the elements in use */
    CLEAR(model->input);
    CLEAR_HALO(model->rho);
    CLEAR_HALO(model->lambda);
    CLEAR_HALO(model->mu);
}


void
model_set_uniform ( model_t *model )
{ 
    real_t
        mu     = kRho * pow(kVs,2.0),
        lambda = ( pow(kVp,2.0) * kRho ) - ( 2.0 * mu );
    for ( int_t z=0; z<HNz; z++ )
        for ( int_t y=0; y<HNy; y++ )
            for ( int_t x=0; x<HNx; x++ )
                RHO(z,y,x) = kRho, LAMBDA(z,y,x) = lambda, MU(z,y,x) = mu;
}


void
model_destroy ( model_t *model )
{
    free ( model->input );
    free ( model->vp );
    free ( model->vs );
    free ( model->rho );
    free ( model->lambda );
    free ( model->mu );
    free ( model->l );
    free ( model->m );
}


void
mesh_init ( mesh_t *mesh )
{
    mesh->sxx  = mesh->syy  = mesh->szz =
    mesh->sxy  = mesh->syz  = mesh->sxz =
    mesh->vx   = mesh->vy   = mesh->vz =
    mesh->del1 = mesh->del2 = mesh->del3 = NULL;
    CLEAR_HALO(mesh->sxx);
    CLEAR_HALO(mesh->syy);
    CLEAR_HALO(mesh->szz);
    CLEAR_HALO(mesh->sxy);
    CLEAR_HALO(mesh->syz);
    CLEAR_HALO(mesh->sxz);
    CLEAR_HALO(mesh->vx);
    CLEAR_HALO(mesh->vy);
    CLEAR_HALO(mesh->vz);
    CLEAR_HALO(mesh->del1);
    CLEAR_HALO(mesh->del2);
    CLEAR_HALO(mesh->del3);
}


void
mesh_destroy ( mesh_t *mesh )
{
    free ( mesh->sxx );
    free ( mesh->syy );
    free ( mesh->szz );
    free ( mesh->sxy );
    free ( mesh->syz );
    free ( mesh->sxz );
    free ( mesh->vx );
    free ( mesh->vy );
    free ( mesh->vz );
    free ( mesh->del1 );
    free ( mesh->del2 );
    free ( mesh->del3 );
}


void
receiver_init ( receiver_t *recv, int_t n, bool p, bool vx, bool vy, bool vz )
{
    recv->p = recv->vx = recv->vy = recv->vz = NULL;
    recv->n = n;

    // In case Nz didn't match any of our receiver setups, leave struct blank
    if ( n == 0 )
        return;

    if ( p )
    {
        recv->p = malloc ( Nt*n*sizeof(real_t) );
        memset ( recv->p, 0, Nt*n*sizeof(real_t) ); 
    }
    if ( vx )
    {
        recv->vx = malloc ( Nt*n*sizeof(real_t) );
        memset ( recv->vx, 0, Nt*n*sizeof(real_t) ); 
    }
    if ( vy )
    {
        recv->vy = malloc ( Nt*n*sizeof(real_t) );
        memset ( recv->vy, 0, Nt*n*sizeof(real_t) ); 
    }
    if ( vz )
    {
        recv->vz = malloc ( Nt*n*sizeof(real_t) );
        memset ( recv->vz, 0, Nt*n*sizeof(real_t) ); 
    }
    // Coordinate arrays
    recv->x = malloc ( n*sizeof(int_t) );
    recv->y = malloc ( n*sizeof(int_t) );
    recv->z = malloc ( n*sizeof(int_t) );
}


void
receiver_setup ( receiver_t *recv )
{
    // In case setup doesn't have a config for Nz
    if ( recv->n == 0 )
        return;
    // Convenience aliases
    int_t x=source_x, y=source_y, z=source_z;
    if(size > 1) {
        if (rcv1_rank == rank){
            recv->x[0] = x+20, recv->y[0] = y+20, recv->z[0] = 42 - rcv1_min;
            recv->x[1] = x-20, recv->y[1] = y+20, recv->z[1] = 42 - rcv1_min;
            recv->x[2] = x-20, recv->y[2] = y-20, recv->z[2] = 42 - rcv1_min;
            recv->x[3] = x+20, recv->y[3] = y-20, recv->z[3] = 42 - rcv1_min;
        }
        if (rcv0_rank == rank){
            recv->x[0] = x+20, recv->y[4] = y+20, recv->z[4] = 12 - rcv0_min;   //actual 4
            recv->x[1] = x-20, recv->y[5] = y+20, recv->z[5] = 12 - rcv0_min;   //actual 5
            recv->x[2] = x-20, recv->y[6] = y-20, recv->z[6] = 12 - rcv0_min;   //actual 6
            recv->x[3] = x+20, recv->y[7] = y-20, recv->z[7] = 12 - rcv0_min;   //actual 7
        }
        return;
    }

    /* Incremental configuration: break statements omitted on purpose
     * "// Falls through"-comments silence gcc warning when
     * -Wimplicit-fallthrough=2
     */
    switch ( Nz )
    {
        case 1024:
            recv->x[32] = x+400, recv->y[32] = y+400, recv->z[32] = z+400;
            recv->x[33] = x-400, recv->y[33] = y+400, recv->z[33] = z+400;
            recv->x[34] = x-400, recv->y[34] = y-400, recv->z[34] = z+400;
            recv->x[35] = x+400, recv->y[35] = y-400, recv->z[35] = z+400;
            recv->x[36] = x+400, recv->y[36] = y+400, recv->z[36] = z-400;
            recv->x[37] = x-400, recv->y[37] = y+400, recv->z[37] = z-400;
            recv->x[38] = x-400, recv->y[38] = y-400, recv->z[38] = z-400;
            recv->x[39] = x+400, recv->y[39] = y-400, recv->z[39] = z-400;
            // Falls through
        case 512:
            recv->x[24] = x+200, recv->y[24] = y+200, recv->z[24] = z+200;
            recv->x[25] = x-200, recv->y[25] = y+200, recv->z[25] = z+200;
            recv->x[26] = x-200, recv->y[26] = y-200, recv->z[26] = z+200;
            recv->x[27] = x+200, recv->y[27] = y-200, recv->z[27] = z+200;
            recv->x[28] = x+200, recv->y[28] = y+200, recv->z[28] = z-200;
            recv->x[29] = x-200, recv->y[29] = y+200, recv->z[29] = z-200;
            recv->x[30] = x-200, recv->y[30] = y-200, recv->z[30] = z-200;
            recv->x[31] = x+200, recv->y[31] = y-200, recv->z[31] = z-200;
            // Falls through
        case 256:
            recv->x[16] = x+100, recv->y[16] = y+100, recv->z[16] = z+100;
            recv->x[17] = x-100, recv->y[17] = y+100, recv->z[17] = z+100;
            recv->x[18] = x-100, recv->y[18] = y-100, recv->z[18] = z+100;
            recv->x[19] = x+100, recv->y[19] = y-100, recv->z[19] = z+100;
            recv->x[20] = x+100, recv->y[20] = y+100, recv->z[20] = z-100;
            recv->x[21] = x-100, recv->y[21] = y+100, recv->z[21] = z-100;
            recv->x[22] = x-100, recv->y[22] = y-100, recv->z[22] = z-100;
            recv->x[23] = x+100, recv->y[23] = y-100, recv->z[23] = z-100;
            // Falls through
        case 128:
            recv->x[8]  = x+40, recv->y[8]  = y+40, recv->z[8]  = z+40;
            recv->x[9]  = x-40, recv->y[9]  = y+40, recv->z[9]  = z+40;
            recv->x[10] = x-40, recv->y[10] = y-40, recv->z[10] = z+40;
            recv->x[11] = x+40, recv->y[11] = y-40, recv->z[11] = z+40;
            recv->x[12] = x+40, recv->y[12] = y+40, recv->z[12] = z-40;
            recv->x[13] = x-40, recv->y[13] = y+40, recv->z[13] = z-40;
            recv->x[14] = x-40, recv->y[14] = y-40, recv->z[14] = z-40;
            recv->x[15] = x+40, recv->y[15] = y-40, recv->z[15] = z-40;
            // Falls through
        case 64:
            recv->x[0] = x+20, recv->y[0] = y+20, recv->z[0] = z+20;
            recv->x[1] = x-20, recv->y[1] = y+20, recv->z[1] = z+20;
            recv->x[2] = x-20, recv->y[2] = y-20, recv->z[2] = z+20;
            recv->x[3] = x+20, recv->y[3] = y-20, recv->z[3] = z+20;
            recv->x[4] = x+20, recv->y[4] = y+20, recv->z[4] = z-20;
            recv->x[5] = x-20, recv->y[5] = y+20, recv->z[5] = z-20;
            recv->x[6] = x-20, recv->y[6] = y-20, recv->z[6] = z-20;
            recv->x[7] = x+20, recv->y[7] = y-20, recv->z[7] = z-20;
            break;
    }
}


void
receiver_save ( int_t ts, receiver_t *recv )
{
    for ( int_t r=0; r<recv->n; r++ )
    {
        // Convenience aliases
        int_t k = recv->z[r], j = recv->y[r], i = recv->x[r];
        // Save pressure? 
        if ( recv->p )
            recv->p[r*Nt+ts] = (1.0/3.0) * (SXX(k,j,i)+SYY(k,j,i)+SZZ(k,j,i));
        if ( recv->vx )
            recv->vx[r*Nt+ts] = VX(k,j,i);
        if ( recv->vy )
            recv->vy[r*Nt+ts] = VY(k,j,i);
        if ( recv->vz )
            recv->vz[r*Nt+ts] = VZ(k,j,i);
    }
}

void
receiver_write_MPI ( receiver_t *recv )
{
    
    FILE *out = fopen ( "receivers0.csv", "w" );
    if (rank == rcv1_rank) {
        fclose ( out );
        out = fopen ( "receivers1.csv", "w" );
    }
    fprintf ( out, "%ld\n%ld\n", recv->n, Nt );
    for ( int_t r=0; r<recv->n; r++ )
    {
        if ( recv->p )
        {
            fprintf ( out, "P\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->p[r*Nt+t] );
        }
        if ( recv->vx )
        {
            fprintf ( out, "Vx\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->vx[r*Nt+t] );
        }
        if ( recv->vy )
        {
            fprintf ( out, "Vy\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->vy[r*Nt+t] );
        }
        if ( recv->vz )
        {
            fprintf ( out, "Vz\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->vz[r*Nt+t] );
        }
        fprintf ( out, "\n" );
    }
    fclose ( out );
    return;
}


void
receiver_write ( receiver_t *recv )
{
    FILE *out = fopen ( "receivers.csv", "w" );
    fprintf ( out, "%ld\n%ld\n", recv->n, Nt );
    for ( int_t r=0; r<recv->n; r++ )
    {
        if ( recv->p )
        {
            fprintf ( out, "P\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->p[r*Nt+t] );
        }
        if ( recv->vx )
        {
            fprintf ( out, "Vx\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->vx[r*Nt+t] );
        }
        if ( recv->vy )
        {
            fprintf ( out, "Vy\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->vy[r*Nt+t] );
        }
        if ( recv->vz )
        {
            fprintf ( out, "Vz\n" );
            for ( int_t t=0; t<Nt; t++ )
                fprintf ( out, "%e\n", recv->vz[r*Nt+t] );
        }
        fprintf ( out, "\n" );
    }
    fclose ( out );
}


void
receiver_destroy ( receiver_t *recv )
{
    // In case no receiver was allocated, unsupported choice of Nz
    if ( recv->n == 0 )
        return;
    free ( recv->p );
    free ( recv->vx );
    free ( recv->vy );
    free ( recv->vz );
    free ( recv->x );
    free ( recv->y );
    free ( recv->z );
}


void
config_print ( double elapsed_time )
{
    switch ( source_type )
    {
        case STRESS:  printf ( "Source type:\tStress monopole\n" ); break;
        case F_MONOP: printf ( "Source type:\tForce monopole\n" );  break;
        case F_DIP:   printf ( "Source type:\tForce dipole\n" );    break;
    }
    printf ( "Grid size\n" );
    printf ( "  Total\t\t%ld x %ld x %ld\n", tNx, tNy, tNz );
    printf ( "  Per rank\t%ld x %ld x %ld\n", Nx, Ny, Nz );
    printf ( "Ranks\t\t%d\n", size );
    printf ( "Iterations\t%ld\n", Nt );
    printf ( "Compute time\t%.4lf\n", elapsed_time );
    printf ( "Total effective MLUPS\t%.5lf\n",
        Nt*tNx*tNy*tNz*1.0e-6 / elapsed_time
    );
    printf ( "       Per rank MULPS\t%.5lf\n",
        Nt*tNx*tNy*tNz*1.0e-6 / (elapsed_time*size)
    );
}


void
options ( int argc, char **argv )
{
    if (argc!=6)
    {
        fprintf ( stderr, "Arguments: Nx Ny Nz Nt source_type\n" );
        exit ( EXIT_FAILURE );
    }
    else
    {
        tNx = strtol ( argv[1], NULL, 10 );
        tNy = strtol ( argv[2], NULL, 10 );
        tNz = strtol ( argv[3], NULL, 10 );
        Nt = strtol ( argv[4], NULL, 10 );
        switch ( strtol ( argv[5], NULL, 10 ) )
        {
            case STRESS:  source_type=STRESS; break;
            case F_MONOP: source_type=F_MONOP; break;
            case F_DIP:   source_type=F_DIP; break;
            default:
                fprintf ( stderr, "Error: unsupported source type\n" );
                exit ( EXIT_FAILURE );
                break;
        }
    }
}
