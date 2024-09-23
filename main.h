#include <stdint.h>
#include <stdbool.h>

#define HALO 1

#define HNx (Nx+2*HALO)
#define HNy (Ny+2*HALO)
#define HNz (Nz+2*HALO)

// Model indexing macros
#define INPUT(z,y,x) (model->rho[(z)*Ny*Nx+(y)*Nx+(x)])
#define RHO(z,y,x) (model->rho[(z)*HNy*HNx+(y)*HNx+(x)])
#define LAMBDA(z,y,x) (model->lambda[(z)*HNy*HNx+(y)*HNx+(x)])
#define MU(z,y,x) (model->mu[(z)*HNy*HNx+(y)*HNx+(x)])

// Mesh indexing macros
#define SXX(z,y,x) (mesh->sxx[(z)*HNy*HNx+(y)*HNx+(x)])
#define SYY(z,y,x) (mesh->syy[(z)*HNy*HNx+(y)*HNx+(x)])
#define SZZ(z,y,x) (mesh->szz[(z)*HNy*HNx+(y)*HNx+(x)])
#define SXY(z,y,x) (mesh->sxy[(z)*HNy*HNx+(y)*HNx+(x)])
#define SYZ(z,y,x) (mesh->syz[(z)*HNy*HNx+(y)*HNx+(x)])
#define SXZ(z,y,x) (mesh->sxz[(z)*HNy*HNx+(y)*HNx+(x)])
#define VX(z,y,x) (mesh->vx[(z)*HNy*HNx+(y)*HNx+(x)])
#define VY(z,y,x) (mesh->vy[(z)*HNy*HNx+(y)*HNx+(x)])
#define VZ(z,y,x) (mesh->vz[(z)*HNy*HNx+(y)*HNx+(x)])
#define DEL1(z,y,x) (mesh->del1[(z)*HNy*HNx+(y)*HNx+(x)])
#define DEL2(z,y,x) (mesh->del2[(z)*HNy*HNx+(y)*HNx+(x)])
#define DEL3(z,y,x) (mesh->del3[(z)*HNy*HNx+(y)*HNx+(x)])

#define MEMSET_CLEAR(d,s) memset(d,0,s*sizeof(real_t))
#define MEMSET_WIPE(d) memset(&d(0,0,0),0,Nz*Ny*Nx*sizeof(real_t))

#define FTOUCH_CLEAR(d,s) do {              \
    _Pragma("omp parallel for")             \
    for ( int_t i=0; i<s; i++ )             \
                d[i] = 0.0;                 \
} while ( false )
#define FTOUCH_WIPE(d) do {                 \
    _Pragma("omp parallel for")             \
    for ( int_t k=0; k<Nz; k++ )            \
        for ( int_t j=0; j<Ny; j++ )        \
            for ( int_t i=0; i<Nx; i++ )    \
                    d(k,j,i) = 0.0;         \
} while ( false )

#ifdef USE_MEMSET
#define CLEAR_METHOD(d,s) MEMSET_CLEAR(d,s)
#define WIPE_METHOD(d) MEMSET_WIPE(d)
#else
#define CLEAR_METHOD(d,s) FTOUCH_CLEAR(d,s)
#define WIPE_METHOD(d) FTOUCH_WIPE(d)
#endif

// Convenience macros for allocate+zero, with & without halo
#define CLEAR(dom) do {                         \
    dom = malloc ( Nx*Ny*Nz*sizeof(real_t) );   \
    CLEAR_METHOD(dom,Nx*Ny*Nz);                 \
} while ( false )

#define CLEAR_HALO(dom) do {                        \
    dom = malloc ( HNx*HNy*HNz*sizeof(real_t) );    \
    CLEAR_METHOD(dom,HNz*HNy*HNx);                  \
} while ( false )

/* Numerical types */
typedef int64_t int_t;
typedef float real_t;

/* Magic values
 * TODO: Would be nice to remove number values from source types, but
 *       retaining scheme for backwards compatibility, these things
 *       come from the command line, and would require a bit of scanning
 *       to remain compat with both number and symbolic name
 */
typedef enum { STRESS=1, F_MONOP=2, F_DIP=3 } source_t;
typedef enum { XDIR, YDIR, ZDIR } direction_t;


/* Model parameter struct */
typedef struct {
    real_t
        * restrict input,
        * restrict vp,
        * restrict vs,
        * restrict rho,
        * restrict lambda,
        * restrict mu,
        * restrict l,
        * restrict m;
} model_t;

void model_init ( model_t * );
void model_set_uniform ( model_t * );
void model_destroy ( model_t * );

/* Mesh struct */
typedef struct {
    real_t dx, dy, dz, dt;
    real_t
        * restrict sxx ,
        * restrict syy ,
        * restrict szz ,
        * restrict sxy ,
        * restrict syz ,
        * restrict sxz ,
        * restrict vx  ,
        * restrict vy  ,
        * restrict vz  ,
        * restrict del1,
        * restrict del2,
        * restrict del3;
} mesh_t;

void mesh_init ( mesh_t * );
void mesh_destroy ( mesh_t * );

/* Receiver struct */
typedef struct {
    int_t n, *x, *y, *z;
    real_t *p, *vx, *vy, *vz;
} receiver_t;

void receiver_init (
    receiver_t *recv, int_t n, bool press, bool vel_x, bool vel_y, bool vel_z
);
void receiver_setup ( receiver_t *recv );
void receiver_save ( int_t ts, receiver_t *recv );
void receiver_write ( receiver_t *recv, double elapsed_time );
void receiver_destroy ( receiver_t *recv );

void receiver_save_MPI ( int_t ts, receiver_t *recv );
void receiver_write_MPI ( receiver_t *recv, double elapsed_time );