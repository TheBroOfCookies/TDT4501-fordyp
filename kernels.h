// Differentiators
#define XF(dst,src) do { \
                for ( int_t l=0; l<HALF; l++ )                              \
                    dst(k,j,i) += W[l] * (src(k,j,i+l+1) - src(k,j,i-l) );  \
                dst(k,j,i) *= (1.0 / kDx);                                  \
} while ( false )

#define XB(dst,src) do { \
                for ( int_t l=0; l<HALF; l++ )                              \
                    dst(k,j,i) += W[l] * (src(k,j,i+l) - src(k,j,i-l-1) );  \
                dst(k,j,i) *= (1.0 / kDx);                                  \
} while ( false )

#define YF(dst,src) do { \
                for ( int_t l=0; l<HALF; l++ )                              \
                    dst(k,j,i) += W[l] * (src(k,j+l+1,i) - src(k,j-l,i) );  \
                dst(k,j,i) *= (1.0 / kDy);                                  \
} while ( false )

#define YB(dst,src) do { \
                for ( int_t l=0; l<HALF; l++ )                              \
                    dst(k,j,i) += W[l] * (src(k,j+l,i) - src(k,j-l-1,i) );  \
                dst(k,j,i) *= (1.0 / kDy);                                  \
} while ( false )

#define ZF(dst,src) do { \
                for ( int_t l=0; l<HALF; l++ )                              \
                    dst(k,j,i) += W[l] * (src(k+l+1,j,i) - src(k-l,j,i) );  \
                dst(k,j,i) *= (1.0 / kDz);                                  \
} while ( false )

#define ZB(dst,src) do { \
                for ( int_t l=0; l<HALF; l++ )                              \
                    dst(k,j,i) += W[l] * (src(k+l,j,i) - src(k-l-1,j,i) );  \
                dst(k,j,i) *= (1.0 / kDz);                                  \
} while ( false )

#define FDD(dst,src,krn) do {                                   \
    WIPE_METHOD(dst);                                           \
    for ( int_t k=HALF; k<Nz-HALF; k++ )                        \
        for ( int_t j=HALF; j<Ny-HALF; j++ )                    \
            for ( int_t i=HALF; i<Nx-HALF; i++ )                \
                krn(dst,src);                                   \
} while ( false )
