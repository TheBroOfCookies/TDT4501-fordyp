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
#define MPI_RANK_ROOT  ( rank == 0 )
#define MPI_RANK_TEST  ( rank == 37 )
#define ARRAY_INDEX(z,y,x,arr) (arr[(z)*HNy*HNx+(y)*HNx+(x)])


/* Function prototypes */
void options ( int argc, char **argv );
void config_print ( double elapsed_time );
void time_step ( int_t iteration );
void insert_source ( int_t ts, source_t type );

void mpi_types_create ( void );
void mpi_types_destroy ( void );
int_t point_to_rank(int_t x, int_t y, int_t z);
void allocate_receivers ( int_t *receivers_to_rank, int_t *recv_points_per_rank, int_t nrecvs );
int_t is_valid_coord (int coords[3]);
void find_neighbours ( void );
void fake_source ( void );
void border_exchange_all ( void );
void border_exchange_3d ( real_t *array );

void mpi_test_save ( real_t *array );

/* Weight coefficients, inner diff. loop */
#define HALF 8
const real_t W[HALF] = {1.2627, -0.1312, 0.0412, -0.0170, 0.0076, -0.0034, 0.0014, -0.0005};                 

/* Simulation parameters */
int_t
    Nx=0, Ny=0, Nz=0, Nt=0, st=0,
    tNx=0, tNy=0, tNz=0;

int
    rank, size, source_rank,
    n_dims = 3, cart_dims[3], cart_coords[3],
    t1_neighbours[6], t2_neighbours[12], t3_neighbours[8];
    //t1 flat neighbours, t2 edge neighbours, t3 corner neighbours, see readme for more info

int_t
    rcv0_min, rcv1_min,
    max_x, min_x,
    max_y, min_y,
    max_z, min_z;

MPI_Comm
    cart_com;

MPI_Datatype
    z_grid, 
    z_face_grid, y_face_grid, x_face_grid,
    z_edge_grid, y_edge_grid, x_edge_grid,
    corner_grid;

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
int_t *recv_to_rank = NULL;
int_t *points_per_rank = NULL;


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
    if (MPI_RANK_ROOT) {
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

    MPI_Dims_create ( size, n_dims, cart_dims );
    MPI_Cart_create ( MPI_COMM_WORLD, n_dims, cart_dims, (int[3]){0,0,0}, 0, &cart_com );

    MPI_Comm_rank ( cart_com, &rank );
    MPI_Cart_coords ( cart_com, rank, n_dims, cart_coords );

    if (cart_dims[0] != 0) { Nz = tNz/cart_dims[0]; } 
    else { Nz = tNz; }

    if (cart_dims[1] != 0) { Ny = tNy/cart_dims[1]; } 
    else { Ny = tNy; }

    if (cart_dims[2] != 0) { Nx = tNx/cart_dims[2]; } 
    else { Nx = tNx; }
    
    if (MPI_RANK_ROOT) {
        printf("In total: tNx %ld tNy %ld tNz %ld\n", tNx, tNy, tNz);
        printf("Per rank:  Nx %ld  Ny %ld  Nz %ld\n",Nx, Ny, Nz);
    }


    //max and min global coordiantes for this rank
    min_x = cart_coords[2]*Nx;
    max_x = (cart_coords[2]+1)*Nx-1;  
    min_y = cart_coords[1]*Ny;
    max_y = (cart_coords[1]+1)*Ny-1;  
    min_z = cart_coords[0]*Nz;
    max_z = (cart_coords[0]+1)*Nz-1;  

    if (MPI_RANK_ROOT) printf("Rank\tmin_x\tmax_x\tmin_y\tmax_y\tmin_z\tmax_z\tCarthesian cooridnates\n");
    else sleep(1);
    printf("%d\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\t%d\t%d\n", rank, min_x, max_x, min_y, max_y, min_z, max_z, cart_coords[0], cart_coords[1], cart_coords[2]);
    sleep(1);
    MPI_Barrier(cart_com);
    find_neighbours( );
    
    /* Source coordinates: initialized here b/c not constant with
     * different domain sizes
     */
    source_x = TNx/2, source_y = TNy/2, source_z = TNz/2;

    source = malloc ( Nt * sizeof(real_t) );
    for ( int_t t=0; t<Nt; t++ )
    {
        real_t arg = pow ( M_PI * kF0 * (kDt * t - kT0), 2.0 );
        source[t] = 1.0e4 * (2.0 * arg - 1.0) * exp(-arg);
    }

    source_rank = point_to_rank(source_x-HALO, source_y-HALO, source_z-HALO);

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
    if (size > 1) {
        if (tNz != 64) { //only 64x64x64 currently supported for saving recievers with MPI ranks
            fprintf( stderr, "\nTotal number of gridpoints in Z-direction must be 64 to support SAVE_RECIEVERS and MPI parallelizaiton.\n" );
            exit(1);
        }
        
    }
    
    recv_to_rank = malloc(sizeof(int_t)*nrecvs);
    points_per_rank = calloc(size, sizeof(int_t));
    allocate_receivers( recv_to_rank, points_per_rank, nrecvs);

    receiver_init ( recv, points_per_rank[rank], false, true, true, true );
    receiver_setup ( recv, nrecvs );
#endif

    model = malloc (sizeof(model_t));
    mesh = malloc (sizeof(mesh_t));
    model_init ( model );
    mesh_init ( mesh );
    model_set_uniform ( model );

    mpi_types_create ( );
    
    struct timeval t_start, t_end;
    gettimeofday ( &t_start, NULL );
    /* for ( int_t t=0; t<Nt; t++ ) {
        if (size > 1) border_exchange_all( );
        time_step ( t );
    } */

    fake_source( );
    #ifdef SAVE_RECEIVERS
    receiver_save ( 0, recv );
    #endif
    border_exchange_3d( mesh->vx );
    MPI_Barrier( cart_com );
    mpi_test_save ( mesh->vx );

    gettimeofday ( &t_end, NULL );


    double elapsed_time = WALLTIME(t_end) - WALLTIME(t_start);
#ifdef SAVE_RECEIVERS
    if (size > 1){
        receiver_write_MPI ( recv, elapsed_time );
    } else {
        receiver_write ( recv, elapsed_time );
    }
#endif
    if (rank==0) config_print (elapsed_time );

    mesh_destroy ( mesh );
    model_destroy ( model );
    free ( mesh );
    free ( model );
    free ( source );
#ifdef SAVE_RECEIVERS
    receiver_destroy ( recv );
    free ( recv );
#endif
    mpi_types_destroy ( );
    MPI_Finalize ( );
    exit ( EXIT_SUCCESS );
}


void 
mpi_types_create ( void ) 
{
    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { HALO, Ny, Nx }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &z_face_grid );
    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { Nz, HALO, Nx }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &y_face_grid );
    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { Nz, Ny, HALO }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &x_face_grid );

    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { Nz, HALO, HALO }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &z_edge_grid );
    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { HALO, Ny, HALO }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &y_edge_grid );
    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { HALO, HALO, Nx }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &x_edge_grid );

    MPI_Type_create_subarray ( n_dims,  (int[3]) { HNz, HNy, HNx }, (int[3]) { HALO, HALO, HALO }, (int[3]) {0,0,0}, MPI_ORDER_C, MPI_FLOAT, &corner_grid );

    MPI_Type_commit ( &z_face_grid );
    MPI_Type_commit ( &y_face_grid );
    MPI_Type_commit ( &x_face_grid );

    MPI_Type_commit ( &z_edge_grid );
    MPI_Type_commit ( &y_edge_grid );
    MPI_Type_commit ( &x_edge_grid );

    MPI_Type_commit ( &corner_grid );
}

int_t 
point_to_rank ( int_t x, int_t y, int_t z ) 
{
    for ( int r=0; r<size; r++ ) {  //Determining which rank contains the given point
        int r_coords[3];
        MPI_Cart_coords ( cart_com, r, n_dims, r_coords );
        int_t r_min_x = r_coords[2]*Nx;
        int_t r_max_x = (r_coords[2]+1)*Nx-1;
        int_t r_min_y = r_coords[1]*Ny;
        int_t r_max_y = (r_coords[1]+1)*Ny-1;
        int_t r_min_z = r_coords[0]*Nz;
        int_t r_max_z = (r_coords[0]+1)*Nz-1;
        if (r_max_z >= z && r_min_z <= z) {
            if (r_max_y >= y && r_min_y <= y) {
                if (r_max_x >= x && r_min_x <= x) {
                    //if (MPI_RANK_ROOT) printf("Rank: %d, contains x=%ld, y=%ld, z=%ld\n", r, x, y, z);
                    return r;
                }
            }
        }
    }
}

void
allocate_receivers ( int_t *receivers_to_rank, int_t *recv_points_per_rank, int_t nrecvs)
{
    int_t z=source_z, y=source_y, x=source_x;
    for ( int r=0; r<size; r++ ) {  //Determining which ranks will have the recievers for 64x64x64
        int r_coords[3];
        MPI_Cart_coords ( cart_com, r, n_dims, r_coords );
        int_t r_min_x = r_coords[2]*Nx;
        int_t r_max_x = (r_coords[2]+1)*Nx-1;
        int_t r_min_y = r_coords[1]*Ny;
        int_t r_max_y = (r_coords[1]+1)*Ny-1;
        int_t r_min_z = r_coords[0]*Nz;
        int_t r_max_z = (r_coords[0]+1)*Nz-1;
        if (r_max_z >= z+20 && r_min_z <= z+20) {
            if (r_max_y >= y+20 && r_min_y <= y+20) {
                if (r_max_x >= x+20 && r_min_x <= x+20) {
                    receivers_to_rank[0] = r; //recv->x[0] = x+20, recv->y[0] = y+20, recv->z[0] = z+20;
                    recv_points_per_rank[r]++;
                } if (r_max_x >= x-20 && r_min_x <= x-20) {
                    receivers_to_rank[1] = r; //recv->x[1] = x-20, recv->y[1] = y+20, recv->z[1] = z+20;
                    recv_points_per_rank[r]++;
                }
            } if (r_max_y >= y-20 && r_min_y <= y-20) {
                if (r_max_x >= x-20 && r_min_x <= x-20) {
                    receivers_to_rank[2] = r; //recv->x[2] = x-20, recv->y[2] = y-20, recv->z[2] = z+20;
                    recv_points_per_rank[r]++;
                } if (r_max_x >= x+20 && r_min_x <= x+20) {
                    receivers_to_rank[3] = r; //recv->x[3] = x+20, recv->y[3] = y-20, recv->z[3] = z+20;
                    recv_points_per_rank[r]++;
                }
            }
        } if (r_max_z >= z-20 && r_min_z <= z-20) {
            if (r_max_y >= y+20 && r_min_y <= y+20) {
                if (r_max_x >= x+20 && r_min_x <= x+20) {
                    receivers_to_rank[4] = r; //recv->x[4] = x+20, recv->y[4] = y+20, recv->z[4] = z-20;
                    recv_points_per_rank[r]++;
                } if (r_max_x >= x-20 && r_min_x <= x-20) {
                    receivers_to_rank[5] = r; //recv->x[5] = x-20, recv->y[5] = y+20, recv->z[5] = z-20;
                    recv_points_per_rank[r]++;
                }
            } if (r_max_y >= y-20 && r_min_y <= y-20) {
                if (r_max_x >= x-20 && r_min_x <= x-20) {
                    receivers_to_rank[6] = r; //recv->x[6] = x-20, recv->y[6] = y-20, recv->z[6] = z-20;
                    recv_points_per_rank[r]++;
                } if (r_max_x >= x+20 && r_min_x <= x+20) {
                    receivers_to_rank[7] = r; //recv->x[7] = x+20, recv->y[7] = y-20, recv->z[7] = z-20;
                    recv_points_per_rank[r]++;
                }
            }
        }
    }
    //print result of receiver allocation
    if ( MPI_RANK_ROOT ) {
        printf("Receiver point -> rank\n");
        for ( int_t i = 0; i < nrecvs; i++ ) {
            printf("%ld->%ld\t", i, recv_to_rank[i]);
        }
        printf("\n");
        printf("Rank -> number of receiver points\n");
        for ( int_t i = 0; i < size; i++ ) {
            if ( points_per_rank[i] > 0 ) printf("%ld->%ld\t", i, points_per_rank[i]);
        }
        printf("\n");
    }
}

int_t
is_valid_coord ( int coord[3] ) 
{
    if (coord[0] < cart_dims[0] && coord[0] >= 0) {
        if (coord[1] < cart_dims[1] && coord[1] >= 0) {
            if (coord[2] < cart_dims[2] && coord[2] >= 0) {
                return 1;
            }   
        }
    }
    return 0;
}

void
find_neighbours ( ) 
{
    //carthesian offsets from the current rank to all "face" neighbours
    //        offsets for     z  y  x
    int_t t1_offsets[6*3] = { 1, 0, 0,   -1, 0, 0,   
                              0, 1, 0,    0,-1, 0,
                              0, 0, 1,    0, 0,-1 };
    for ( int_t k=0; k<6; k++ ) {
        int coords[3] = {cart_coords[0]+t1_offsets[3*k], cart_coords[1]+t1_offsets[3*k+1], cart_coords[2]+t1_offsets[3*k+2]};
        if ( is_valid_coord(coords) ) { MPI_Cart_rank(cart_com, coords, &t1_neighbours[k]); }
        else { t1_neighbours[k] = MPI_PROC_NULL; }
    }

     //carthesian offsets from the current rank to all "edge" neighbours
    int_t t2_offsets[12*3] = { 1, 1, 0,   1,-1, 0,   1, 0, 1,   1, 0,-1,
                              -1, 1, 0,  -1,-1, 0,  -1, 0, 1,  -1, 0,-1,  
                               0, 1, 1,   0,-1, 1,   0, 1,-1,   0,-1,-1 };
    for ( int_t j=0; j<12; j++ ) {
        int coords[3] = {cart_coords[0]+t2_offsets[3*j], cart_coords[1]+t2_offsets[3*j+1], cart_coords[2]+t2_offsets[3*j+2]};
        if ( is_valid_coord(coords) ) { MPI_Cart_rank(cart_com, coords, &t2_neighbours[j]); }
        else { t2_neighbours[j] = MPI_PROC_NULL; }
    }

    //carthesian offsets from the current rank to all "corner" neighbours
    int_t t3_offsets[8*3] = { 1, 1, 1,   1, 1,-1,   1,-1, 1,   1,-1,-1,
                             -1, 1, 1,  -1, 1,-1,  -1,-1, 1,  -1,-1,-1 };
    for ( int_t i=0; i<8; i++ ) {
        int coords[3] = {cart_coords[0]+t3_offsets[3*i], cart_coords[1]+t3_offsets[3*i+1], cart_coords[2]+t3_offsets[3*i+2]};
        if ( is_valid_coord(coords) ) { MPI_Cart_rank(cart_com, coords, &t3_neighbours[i]); }
        else { t3_neighbours[i] = MPI_PROC_NULL; }
    }

    if (MPI_RANK_TEST) {
        printf("---------------------------\n");
        printf("Neighbours of rank %d\n", rank);

        printf("T1: \n");
        for ( int_t k=0; k<6; k++ ) {
            if ( k%2 == 0 && k != 0 ) printf("\n");
            if (t1_neighbours[k] != -2){
                printf("%d\t", t1_neighbours[k]);
            } else {
                printf("-\t");
            } 
        }

        printf("\nT2: \n");
        for ( int_t j=0; j<12; j++ ) {
            if ( j%4 == 0 && j != 0 ) printf("\n");
            if (t2_neighbours[j] != -2) {
                printf("%d\t", t2_neighbours[j]);
            } else {
                printf("-\t");
            } 
        }

        printf("\nT3: \n");
        for ( int_t i=0; i<8; i++ ) {
            if ( i%4 == 0 && i != 0 ) printf("\n");
            if (t3_neighbours[i] != -2) {
                printf("%d\t", t3_neighbours[i]);;
            } else {
                printf("-\t");
            } 
        }
        printf("\n");
        printf("---------------------------\n");
    }

}

void
fake_source ( )
{
    for ( int_t z=0; z<HNz; z++ ) {
        for ( int_t y=0; y<HNy; y++ ) {
            for ( int_t x=0; x<HNx; x++ ) {
                if(x < HALO || x >= Nx+HALO || y < HALO || y >= Ny+HALO || z < HALO || z > Nz+HALO) {
                    VX(z,y,x) = -1;
                    VY(z,y,x) = -1;
                    VZ(z,y,x) = -1;
                } else {
                    VX(z,y,x) = rank;
                    VY(z,y,x) = rank;
                    VZ(z,y,x) = rank;
                }
            }
        }
    }
}

void
insert_source ( int_t ts, source_t type )
{
    // Convenience aliases
    int_t z=source_z, y=source_y, x=source_x;
    real_t s = source[ts];
    // Determine source type, act accordingly
    // Determine correct rank for force application
    if (rank != source_rank) return;
    //adjust global coord to local coord
    x = x - (min_x);
    y = y - (min_y);
    z = z - (min_z);
    switch ( type )
    {
        case STRESS:
            SXX(z,y,x) += s * kDt;
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
    border_exchange_3d( mesh->vx );
    border_exchange_3d( mesh->vy );
    border_exchange_3d( mesh->vz );

    border_exchange_3d( mesh->sxx );
    border_exchange_3d( mesh->syy );
    border_exchange_3d( mesh->szz );

    border_exchange_3d( mesh->sxy );
    border_exchange_3d( mesh->syz );
    border_exchange_3d( mesh->sxz );
}
/*
int_t t1_neighbour_order[6*3] = { 1, 0, 0,   -1, 0, 0,   
                                  0, 1, 0,    0,-1, 0,
                                  0, 0, 1,    0, 0,-1 };

int_t t2_neighbour_order[12*3] = { 1, 1, 0,   1,-1, 0,   1, 0, 1,   1, 0,-1,
                                  -1, 1, 0,  -1,-1, 0,  -1, 0, 1,  -1, 0,-1,  
                                   0, 1, 1,   0,-1, 1,   0, 1,-1,   0,-1,-1 };

int_t t3_neighbour_order[8*3] = { 1, 1, 1,   1, 1,-1,   1,-1, 1,   1,-1,-1,
                                 -1, 1, 1,  -1, 1,-1,  -1,-1, 1,  -1,-1,-1 };

translation from t*_neighbour_order to t*_snd/t*_rcv
t1
snd 1 = N*  0 = HALO   -1 = HALO
rcv 1 = N*+HALO  0 = HALO  -1 = 0

t2
snd 1 = N*  0 = HALO  -1 = HALO
rcv 1 = N*+HALO  0 = HALO  -1 = 0

t3
snd  1 = N*   -1 = HALO
rcv  1 = N*+HALO   -1 = 0 
*/

void
border_exchange_3d ( real_t *array )
{
    MPI_Datatype t1_types[6] = { z_face_grid, z_face_grid, 
                                 y_face_grid, y_face_grid, 
                                 x_face_grid, x_face_grid };
    MPI_Datatype t2_types[12] = { x_edge_grid, x_edge_grid, y_edge_grid, y_edge_grid, 
                                  x_edge_grid, x_edge_grid, y_edge_grid, y_edge_grid,
                                  z_edge_grid, z_edge_grid, z_edge_grid, z_edge_grid };

    int_t t1_snd[6*3] = { Nz, HALO, HALO,        HALO, HALO, HALO,   
                          HALO, Ny, HALO,        HALO, HALO, HALO,
                          HALO, HALO, Nx,        HALO, HALO, HALO };
    int_t t1_rcv[6*3] = { Nz+HALO, HALO, HALO,   0, HALO, HALO,   
                          HALO, Ny+HALO, HALO,   HALO, 0, HALO,
                          HALO, HALO, Nx+HALO,   HALO, HALO, 0 };

    int_t t2_snd[12*3] = { Nz, Ny, HALO,             Nz, HALO, HALO,     Nz, HALO, Nx,             Nz, HALO, HALO,
                           HALO, Ny, HALO,           HALO, HALO, HALO,   HALO, HALO, Nx,           HALO, HALO, HALO,  
                           HALO, Ny, Nx,             HALO, HALO, Nx,     HALO, Ny, HALO,           HALO, HALO, HALO, };
    int_t t2_rcv[12*3] = { Nz+HALO, Ny+HALO, HALO,   Nz+HALO, 0, HALO,   Nz+HALO, HALO, Nx+HALO,   Nz+HALO, HALO, 0,
                           0, Ny+HALO, HALO,         0, 0, HALO,         0, HALO, Nx+HALO,         0, HALO, 0,  
                           HALO, Ny+HALO, Nx+HALO,   HALO, 0, Nx+HALO,   HALO, Ny+HALO, 0,         HALO, 0, 0 };

    int_t t3_snd[8*3] = { Nz, Ny, Nx,                   Nz, Ny, HALO,           Nz, HALO, Nx,           Nx, HALO, HALO,
                          HALO, Ny, Nx,                 HALO, Ny, HALO,         HALO, HALO, Nx,         HALO, HALO, HALO };
    int_t t3_rcv[8*3] = { Nz+HALO, Ny+HALO, Nx+HALO,    Nz+HALO, Ny+HALO, 0,    Nz+HALO, 0, Nx+HALO,    Nz+HALO, 0, 0,
                          0, Ny+HALO, Nx+HALO,          0, Ny+HALO, 0,          0, 0, Nx+HALO,          0, 0, 0 };
    
    for ( int_t k=0; k<6; k++ ) {
        real_t *snd = &ARRAY_INDEX(t1_snd[3*k], t1_snd[3*k+1], t1_snd[3*k+2], array);
        real_t *rcv = &ARRAY_INDEX(t1_rcv[3*k], t1_rcv[3*k+1], t1_rcv[3*k+2], array);
        MPI_Sendrecv(snd, 1, t1_types[k], t1_neighbours[k], 0, rcv, 1, t1_types[k], t1_neighbours[k], 0, cart_com, MPI_STATUS_IGNORE);
    }

    for ( int_t j=0; j<12; j++ ) {
        real_t *snd = &ARRAY_INDEX(t2_snd[3*j], t2_snd[3*j+1], t2_snd[3*j+2], array);
        real_t *rcv = &ARRAY_INDEX(t2_rcv[3*j], t2_rcv[3*j+1], t2_rcv[3*j+2], array);
        if (t2_neighbours[j] > -1) MPI_Sendrecv(snd, 1, t2_types[j], t2_neighbours[j], 0, rcv, 1, t2_types[j], t2_neighbours[j], 0, cart_com, MPI_STATUS_IGNORE);
    }

    for ( int_t i=0; i<8; i++ ) {
        real_t *snd = &ARRAY_INDEX(t3_snd[3*i], t3_snd[3*i+1], t3_snd[3*i+2], array);
        real_t *rcv = &ARRAY_INDEX(t3_rcv[3*i], t3_rcv[3*i+1], t3_rcv[3*i+2], array);
        if (t3_neighbours[i] > -1) MPI_Sendrecv(snd, 1, corner_grid, t3_neighbours[i], 0, rcv, 1, corner_grid, t3_neighbours[i], 0, cart_com, MPI_STATUS_IGNORE);
    }
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
receiver_setup ( receiver_t *recv, int_t nrecvs )
{
    // In case setup doesn't have a config for Nz
    if ( recv->n == 0 )
        return;
    // Convenience aliases
    int_t x=source_x, y=source_y, z=source_z;
    if (size > 1) {
        int_t p_num = 0;
        for (int_t p = 0; p < nrecvs; p++){
            if (recv_to_rank[p] == rank){
                switch (p) {
                    case 0: recv->x[p_num] = x+20 - min_x, recv->y[p_num] = y+20 - min_y, recv->z[p_num] = z+20 - min_z; p_num++; break;
                    case 1: recv->x[p_num] = x-20 - min_x, recv->y[p_num] = y+20 - min_y, recv->z[p_num] = z+20 - min_z; p_num++; break;
                    case 2: recv->x[p_num] = x-20 - min_x, recv->y[p_num] = y-20 - min_y, recv->z[p_num] = z+20 - min_z; p_num++; break;
                    case 3: recv->x[p_num] = x+20 - min_x, recv->y[p_num] = y-20 - min_y, recv->z[p_num] = z+20 - min_z; p_num++; break;
                    case 4: recv->x[p_num] = x+20 - min_x, recv->y[p_num] = y+20 - min_y, recv->z[p_num] = z-20 - min_z; p_num++; break;
                    case 5: recv->x[p_num] = x-20 - min_x, recv->y[p_num] = y+20 - min_y, recv->z[p_num] = z-20 - min_z; p_num++; break;
                    case 6: recv->x[p_num] = x-20 - min_x, recv->y[p_num] = y-20 - min_y, recv->z[p_num] = z-20 - min_z; p_num++; break;
                    case 7: recv->x[p_num] = x+20 - min_x, recv->y[p_num] = y-20 - min_y, recv->z[p_num] = z-20 - min_z; p_num++; break;
                }
            }
        }
        return;
    }

    /* Incremental configuration: break statements omitted on purpose
     * "// Falls through"-comments silence gcc warning when
     * -Wimplicit-fallthrough=2
     */
    switch ( Nz )
    {
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
receiver_write_MPI ( receiver_t *recv, double elapsed_time )
{
    for ( int r = 0; r < size; r++ ){
        if (rank != r) {
            MPI_Barrier(cart_com);
        } else {
            char filename[50];
            sprintf(filename, "receivers/receivers_mpi%d.csv", r);
            FILE *out = fopen ( filename, "w" );
            fprintf( out, "%d,%lf\n", size, elapsed_time);
            fprintf ( out, "%ld,%d,%ld\t", recv->n, HALO, Nt );
            fprintf( out, "%ld,%ld,%ld\n", tNx, tNy, tNz);
            //fprintf( out, "%d,%d,%d\n", rank, rank, rank);
            
            for ( int_t r=0; r<recv->n; r++ )
            {
                fprintf ( out, "Vx\n" );
                for ( int_t t=0; t<Nt; t++ )
                    fprintf ( out, "%e\n", recv->vx[r*Nt+t] );
                fprintf ( out, "Vy\n" );
                for ( int_t t=0; t<Nt; t++ )
                    fprintf ( out, "%e\n", recv->vy[r*Nt+t] );
                fprintf ( out, "Vz\n" );
                for ( int_t t=0; t<Nt; t++ )
                    fprintf ( out, "%e\n", recv->vz[r*Nt+t] );
                fprintf ( out, "\n" );
            }
            fclose ( out );
            MPI_Barrier(cart_com);
        }
    }
    return;
}


void
receiver_write ( receiver_t *recv, double elapsed_time )
{
    FILE *out = fopen ( "receivers/receivers.csv", "w" );
    fprintf( out, "1,%lf\n", elapsed_time);
    fprintf ( out, "%ld,%d,%ld\t", recv->n, HALO, Nt );
    fprintf( out, "%ld,%ld,%ld\n", tNx, tNy, tNz);
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

    free ( recv_to_rank );
    free ( points_per_rank );
}

void 
mpi_types_destroy ( void )
{
    MPI_Type_free ( &z_face_grid );
    MPI_Type_free ( &y_face_grid );
    MPI_Type_free ( &x_face_grid );

    MPI_Type_free ( &z_edge_grid );
    MPI_Type_free ( &y_edge_grid );
    MPI_Type_free ( &x_edge_grid );

    MPI_Type_free ( &corner_grid );
}

void mpi_test_save (real_t *array) {
    if (MPI_RANK_TEST) {
        for ( int_t k=0; k<HNz; k++ ) {
            char filename[50];
            sprintf(filename, "receivers/3d_%ldborder_exchange.csv", k);
            FILE *out = fopen ( filename, "w" );
            fprintf( out, "Rank: %d, size: %d, slice: %ld\n", rank, size, k);
            fprintf ( out, "%ld,%d,%ld\t", recv->n, HALO, Nt );
            fprintf( out, "%ld,%ld,%ld\n", tNx, tNy, tNz);
            for ( int_t j=1; j<=HNy; j++ ) {
                    fprintf ( out, "\n" );
                    for ( int_t i=0; i<HNx; i++ ){
                        int num = (int)ARRAY_INDEX(k, HNy-j, i, array);  //mirrored y coordinate
                        if (num != -1){
                            fprintf ( out, "%d\t", num );
                        } else {
                            fprintf ( out, "-\t");
                        }
                    }
            }
            fprintf ( out, "\n" );
            fclose ( out );
        }
    }
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
        Nt*tNx*tNy*tNz*1.0e-6 / elapsed_time);
    printf ( "       Per rank MLUPS\t%.5lf\n",
        Nt*tNx*tNy*tNz*1.0e-6 / (elapsed_time*size));
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
        Nt = strtol  ( argv[4], NULL, 10 );
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
