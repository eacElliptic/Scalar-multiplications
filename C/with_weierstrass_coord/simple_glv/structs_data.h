#ifndef STRUCTS_GLOBD

#define STRUCTS_GLOBD 

#define BIT_SIZE 256

typedef struct AffPoint{
	mpz_t X;
	mpz_t Y;
} AffPoint;

typedef struct JacPoint{
	mpz_t X;
	mpz_t Y;
	mpz_t Z;
} JacPoint;

typedef struct GLVScalar{
	char tk1[BIT_SIZE];
	char tk2[BIT_SIZE];
	int j;
} GLVScalar;

typedef struct GLVData{
	AffPoint PP[4];
	mpz_t phi, efp, refp, aa, bb, Na;
} GLVData;


static mpz_t p, ca, cb, beta, A, B, C, D;


#endif



