#include "eac_timing.h"
#include <stdint.h> 
#define BILLION 1000000000L





void init_datas(Point *P){
	
	mpz_t X, Y, tmp_exp;
	mpz_inits (X, Y, tmp_exp, NULL);
	
	//~ p = 2^358-36855
	mpz_setbit (p, 358);
	mpz_sub_ui (p, p, 36855);
	
	//~ Curve is y^2 = x^3 + 17
	mpz_set_ui (ca, 0);
	mpz_set_ui (cb, 17);
	
	//~ 2 is a generator of Fp and beta is an element of order 3
	mpz_set_ui (beta, 2);
	mpz_sub_ui (tmp_exp, p, 1); 
	mpz_divexact_ui (tmp_exp, tmp_exp, 3);
	mpz_powm (beta, beta, tmp_exp, p);

	mpz_set_str (X, "568554169514108215082878628903438431409037009413741676480372717420637861491120753156713065337030704866093074", 10);
	mpz_set_str (Y, "125804816918392729528454862544298396901812597725830867877741173065843595271902914550421712710898294754082775", 10);
	//~ P represents point P itself and the point Phi_P see point.h for it definition 
	mpz_set (P->PQ[0], X);
	mpz_set (P->PQ[2], Y);
	mpz_set_ui (P->Z, 1);
	
	if (is_on_curve(P, 0, ca, cb)) 
		printf("P OK\n");
	else
		printf("P not OK\n");
	
	
	mpz_clears (X, Y, tmp_exp, NULL);
}




void go(Point *P, int nbiter){
	
	int i,j;
	int eac_len = 256; // We are going to use eac of length 256 
	uchar eac[256]; 
	unsigned long seed = time(NULL); 
	mpz_t ZZ, tmp;
	gmp_randstate_t r; 
	
	struct timespec start1;
    struct timespec end1;
    uint64_t diff1, diff2, diff3;
	
	mpz_inits (ZZ, tmp, NULL);
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);

	printf("EAC 358 bits benchmark, running...\n\n");

    
    diff1 = 0;
    
	for(i=0; i<nbiter; i++){
		//~ Generation of an eac
		mpz_urandomb(tmp, r, eac_len);    
		memset(eac, 0, eac_len);
		j=0;
		while(mpz_cmp_ui (tmp, 0) && j<eac_len){ // Check if tmp is not equal to zero
			eac[j] = mpz_tstbit (tmp, j);
			mpz_clrbit (tmp, j);
			j++;
		}
		
		mpz_set_ui (P->Z, 1);
		//~ We compute PHI_P 
        mpz_mul (P->PQ[1], P->PQ[0], beta);
		mpz_mod (P->PQ[1], P->PQ[1], p);
        mpz_set (P->PQ[3], P->PQ[2]);
        
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start1);
        
        //~ Coordinates of P are updated so that at the end of the method P contains ((k-t)P, kP) for some integer t 
		point_from_eac(P, eac, eac_len);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end1); 
		diff1 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec;

		mpz_invert (ZZ, P->Z, p);
		mpz_mul (tmp, ZZ, ZZ); // tmp = ZZ^2
		//~ P is transformed in affine coordinates for next iteration
		mpz_mul (P->PQ[0], P->PQ[0], tmp);
		mpz_mod (P->PQ[0], P->PQ[0], p);
		
		mpz_mul (P->PQ[2], P->PQ[2], tmp);
		mpz_mul (P->PQ[2], P->PQ[2], ZZ);
		mpz_mod (P->PQ[2], P->PQ[2], p);
	}
	
	if (is_on_curve(P, 1, ca, cb)) {
		printf("MULT: kP OK\n");
		print_point(P, 1);
		printf("Elapsed time = %llu  nanoseconds\n", (long long unsigned int) diff1);
	}
	else{
		printf("MULT: kP not OK\n");
		print_point(P, 1);
		printf("Elapsed time = %llu  nanoseconds\n", (long long unsigned int) diff1);
	}

	mpz_clears(ZZ, tmp, NULL);
	gmp_randclear(r);

}


void print_eac(uchar *eac, int eac_len){
	int i;
	for(i=0; i<eac_len; i++){
		printf("%d", eac[i]);
	}
	printf("\n");
}








