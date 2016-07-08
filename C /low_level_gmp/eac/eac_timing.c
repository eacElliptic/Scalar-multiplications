#include "eac_timing.h"
#include <stdint.h> 
#define BILLION 1000000000L



void init_datas(Point *p0p1){
	int j;
	mpz_t X, Y, tmp_exp;
	mp_limb_t *tx_limbs, *ty_limbs;
	
	mpz_inits (X, Y, tmp_exp, NULL);
	
	//~ p = 2^358-36855
	mpz_setbit (p, 358);
	mpz_sub_ui (p, p, 36855);
	p_limbs = mpz_limbs_modify (p, NB_LIMBS);
	
	//~ add_overflow = 2^384-p
	mpz_setbit (add_overflow, 384); // 384 = 64*NB_LIMBS
	mpz_sub (add_overflow, add_overflow, p);
	add_overflow_limbs = mpz_limbs_modify (add_overflow, NB_LIMBS);
	
	//~ Curve is y^2 = x^3 + 17
	ca = 0;
	cb = 17;
	
	//~ 2 is a generator of Fp and beta is an element of order 3
	mpz_set_ui (beta, 2);
	mpz_sub_ui (tmp_exp, p, 1); 
	mpz_divexact_ui (tmp_exp, tmp_exp, 3);
	mpz_powm (beta, beta, tmp_exp, p);
	beta_limbs = mpz_limbs_modify (beta, NB_LIMBS);

	for(j=0; j<NB_LIMBS; j++)
		Z_limbs[j] = 0;
	init_point(p0p1);
	init_point((p0p1+1));
	
	mpz_set_str (X, "568554169514108215082878628903438431409037009413741676480372717420637861491120753156713065337030704866093074", 10);
	mpz_set_str (Y, "125804816918392729528454862544298396901812597725830867877741173065843595271902914550421712710898294754082775", 10);
	tx_limbs = mpz_limbs_modify (X, NB_LIMBS);
	ty_limbs = mpz_limbs_modify (Y, NB_LIMBS);
	
	//~ the point P
	mpn_copyi (p0p1->X_limbs, tx_limbs, NB_LIMBS);
	mpn_copyi (p0p1->Y_limbs, ty_limbs, NB_LIMBS);
	Z_limbs[0] = 1;
	
	if (is_on_curve(p0p1, &ca, 1, &cb, 1)) {
		printf("P OK\n");
		print_point(p0p1);
	}
	else
		printf("P not OK\n");
	
	mpz_clears (X, Y, tmp_exp, NULL);
}




void go(Point *p0p1, int nbiter){
	
	int i,j;
	int eac_len = 256; // We are going to use eac of length 256 
	uchar eac[256]; 
	
	
	unsigned long seed = time(NULL); 
	gmp_randstate_t r; 
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	struct timespec start1;
    struct timespec end1;
    uint64_t diff1;
	
	
	mpz_t ZZ, tmp_inv, tmp_rand;
	mp_limb_t *zz_limbs;
	mp_limb_t tmp_limbs[(NB_LIMBS*2)];
	mp_limb_t tmp_prod_limbs[(NB_LIMBS*2)];
	mpz_inits (ZZ, tmp_rand, NULL);
	

	printf("\nEAC 358 bits benchmark, running...\n\n");

    
    diff1 = 0;
    
	for(i=0; i<nbiter; i++){
		//~ Generation of an eac
		mpz_urandomb(tmp_rand, r, eac_len);    
		memset(eac, 0, eac_len);
		j=0;
		while(mpz_cmp_ui (tmp_rand, 0) && j<eac_len){
			eac[j] = mpz_tstbit (tmp_rand, j);
			mpz_clrbit (tmp_rand, j);
			j++;
		}

        //~ We compute PHI_P
        apply_endo((p0p1+1), p0p1, beta_limbs); // Point PHI_P 
		
        for(j=1; j<NB_LIMBS; j++)
			Z_limbs[j] = 0;
        Z_limbs[0] = 1;
        
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start1);
        
        //~ Coordinates of P are updated so that at the end of the method P contains ((k-t)P, kP) for some integer t 
		point_from_eac(p0p1, eac, eac_len);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end1); // CLOCK_MONOTONIC CLOCK_MONOTONIC_RAW CLOCK_PROCESS_CPUTIME_ID CLOCK_REALTIME
		diff1 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec;


		mpz_roinit_n (tmp_inv, Z_limbs, NB_LIMBS);
		mpz_invert (ZZ, tmp_inv, p);
		zz_limbs = mpz_limbs_modify (ZZ, NB_LIMBS);
		mpn_sqr (tmp_limbs, zz_limbs, NB_LIMBS);  // tmp = ZZ^2
		mpn_tdiv_qr (useless_q, tmp_limbs, 0, tmp_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
       
		//~ P is transformed in affine coordinates for next iteration
		mpn_mul_n (tmp_prod_limbs, p0p1->X_limbs, tmp_limbs, NB_LIMBS);
		mpn_tdiv_qr (useless_q, p0p1->X_limbs, 0, tmp_prod_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
       
		mpn_mul_n (tmp_prod_limbs, p0p1->Y_limbs, tmp_limbs, NB_LIMBS);
		mpn_tdiv_qr (useless_q, p0p1->Y_limbs, 0, tmp_prod_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
		mpn_mul_n (tmp_prod_limbs, p0p1->Y_limbs, zz_limbs, NB_LIMBS);
		mpn_tdiv_qr (useless_q, p0p1->Y_limbs, 0, tmp_prod_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	}
	
	
	if (is_on_curve((p0p1+1), &ca, 1, &cb, 1)) {
		printf("MULT: kP OK\n");
		print_point((p0p1+1));
		printf("Elapsed time = %llu  nanoseconds\n", (long long unsigned int) diff1);
	}
	else{
		printf("MULT: kP not OK\n");
		print_point((p0p1+1));
		printf("Elapsed time = %llu  nanoseconds\n", (long long unsigned int) diff1);
	}

	mpz_clears (ZZ, tmp_rand, NULL);
	gmp_randclear(r);
}


void print_eac(uchar *eac, int eac_len){
	int i;
	for(i=0; i<eac_len; i++){
		printf("%d", eac[i]);
	}
	printf("\n");
}








