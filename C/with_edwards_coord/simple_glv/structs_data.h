#ifndef STRUCTS_GLOBD

#define STRUCTS_GLOBD 

#define BIT_SIZE 253   // size of base point order.

typedef struct ExtAffPoint{
	mpz_t x;
	mpz_t y;
	mpz_t t;
} ExtAffPoint;

typedef struct ExtProjPoint{
	mpz_t X;
	mpz_t Y;
	mpz_t Z;
	mpz_t T;
} ExtProjPoint;

typedef struct GLVScalar{
	char tk1[BIT_SIZE];
	char tk2[BIT_SIZE];
	int j;
} GLVScalar;

typedef struct GLVData{
	ExtAffPoint PP[4];
	mpz_t phi, efp, refp, aa, bb, Na;
} GLVData;


static mpz_t p, curve_a, curve_d, beta, A, B;


#endif



