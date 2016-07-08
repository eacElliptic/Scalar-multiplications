#ifndef STRUCTS_GLOBD

#define STRUCTS_GLOBD 

#define BIT_SIZE 256
#define SCAL_SIZE 129  // (ceil part of log_2 (efp)/2) + 1 

typedef struct AffPoint{
	mpz_t X;
	mpz_t Y;
} AffPoint;

typedef struct JacPoint{
	mpz_t X;
	mpz_t Y;
	mpz_t Z;
} JacPoint;

typedef struct SGLVScalar{
	char tk1[SCAL_SIZE];
	char tk2[SCAL_SIZE];
	int j;
	char isEven;
} SGLVScalar;

typedef struct SGLVData{
	AffPoint PP[2];
	mpz_t phi, efp, refp, aa, bb, Na;
} SGLVData;


static mpz_t p, ca, cb, beta, A, B, C, D;

static gmp_randstate_t rand_gen;


#endif



