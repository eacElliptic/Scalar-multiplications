#include "sglv_timing.h"
#include <stdint.h> 
#define BILLION 1000000000L



void init_glvObjects(SGLVData *op){
	mpz_inits (op->phi, op->efp, op->refp, op->aa, op->bb, op->Na, NULL);
	
	init_affPoint(op->PP);
	init_affPoint((op->PP)+1);
}

void free_glvObjects(SGLVData *op){
	mpz_clears (op->phi, op->efp, op->refp, op->aa, op->bb, op->Na, NULL);

	free_affPoint(op->PP);
	free_affPoint((op->PP)+1);
}

void init_data(SGLVData *op){
	//~ Note: 'op' is supposed to be already initialised (see main.c)
	mpz_t  X, Y, tmp_exp;
	mpz_inits (X, Y, tmp_exp, NULL);
	
	//~ p = 2^256-1539
	mpz_setbit (p, 256);
	mpz_sub_ui (p, p, 1539);
	//~ Curve is y^2 = x^3 + 5 
	//~ mpz_set_ui (ca, 0);  // Already set to ZERO through the initialisation
	mpz_set_ui (cb, 5);
	mpz_set_str (X, "66043340678279369258981193985450715448372246692524118399919799979198175326756", 10);
	mpz_set_str (Y, "52931614173969837860927627620185828149756459753631369661231760168280520196499", 10);
	//~ Point P
	mpz_set (op->PP[0].X , X);
	mpz_set (op->PP[0].Y , Y);
	//~ 5 is a generator of Fp and beta is an element of order 3
	mpz_set_ui (beta, 5);
	mpz_sub_ui (tmp_exp, p, 1); 
	mpz_divexact_ui (tmp_exp, tmp_exp, 3);
	mpz_powm (beta, beta, tmp_exp, p);
	
	//~ Order of P
	mpz_set_str (op->efp, "115792089237316195423570985008687907852920869663551026044856920202296625078603", 10);
	//~ Integer part of square root of efp
	mpz_set_str (op->refp, "340282366920938463463374607431768211455", 10);
	//~ Root of  phi^2+phi+1 = 0 mod #E(Fp)
	mpz_set_str (op->phi, "86115113571596370384202877098459685148505633832452548870570058599719872579953", 10);
	
	init_glvRData(op);
	 
	if (is_on_curve_aff(op->PP, ca, cb)) {
		printf("P OK\n");
		print_affPoint(op->PP);
	}
	else
		printf("P not OK\n");
	
	mpz_clears (X, Y, tmp_exp, NULL);
}

//~ Initialises remainning data of 'op'
void init_glvRData(SGLVData *op){
	
	mpz_t  u, v, x1, x2, y1, y2, q, r, x, y, tmp;
	mpz_inits (u, v, x1, x2, y1, y2, q, r, x, y, tmp, NULL);
	
	mpz_set (u, op->efp);
	mpz_set (v, op->phi);
	
	mpz_set_ui (x1, 1);
	mpz_set_ui (y2, 1);
	
	while(mpz_cmp(u, op->refp) > 0){
		mpz_fdiv_qr (q, r, v, u);
		
		mpz_mul (tmp, q, x1);
		mpz_sub (x, x2, tmp);
		mpz_mul (tmp, q, y1);
		mpz_sub (y, y2, tmp);
		
		mpz_set (v, u);
		mpz_set (u, r);
		
		mpz_set (x2, x1);
		mpz_set (x1, x);
		mpz_set (y2, y1);
		mpz_set (y1, y);
	}
	
	mpz_add (op->aa, u, y);		// aa = u + y
	mpz_neg (op->bb, y);	 	// bb = -y
	
	mpz_mul (x1, op->bb, op->bb);		// x1 = (bb)^2
	mpz_mul (x2, u, op->bb);			// x2 = u * bb
	mpz_mul (op->Na, u, u);				// Na = u^2
	mpz_add (op->Na, op->Na, x1);		// Na = u^2 + (bb)^2 
	mpz_sub (op->Na, op->Na, x2);		// Na = u^2 + (bb)^2 - (u * bb)
	
	mpz_clears (u, v, x1, x2, y1, y2, q, r, x, y, tmp, NULL);
}



void go(SGLVData *data, int nbiter){
	int i;
	mpz_t  k_rd, ZZ, tmp;
	mpz_inits (k_rd, ZZ, tmp, NULL);
	
	gmp_randstate_t rand_gen;
	unsigned long seed = time(NULL); 
	gmp_randinit_default(rand_gen);
	gmp_randseed_ui(rand_gen, seed);
	
	
	SGLVScalar k;
	JacPoint jp;
	
	struct timespec start1,start2;
    struct timespec end1,end2;
    uint64_t diff1,diff2, diff3;
	
	
	printf("\nSGLV 256 bits benchmark, running...\n\n");
	
	
	diff1 = 0;
	diff2 = 0;
	
	init_jacPoint(&jp);
	
	
	
	for (i = 0; i < nbiter; i++) {
		
		mpz_urandomb (k_rd, rand_gen, BIT_SIZE);
		mpz_mod (k_rd, k_rd,  data->efp);
		
		
		apply_endo(((data->PP)+1), (data->PP), beta); /* Point PHI_P */
		add_aff_aff(((data->PP)+1), (data->PP));      /* Point P+PHI_P */
		
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
		
		init_glvScalar(&k);
		build_glvScalar(&k, k_rd, data->aa, data->bb, data->Na);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
		diff1 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec; 
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
		
		point_from_SGLVScalar(&jp, &k, data->PP);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
		diff2 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec; 
		
		mpz_invert (ZZ, jp.Z, p);
		mpz_mul (tmp, ZZ, ZZ); // tmp = ZZ^2
		//~ P is transformed in affine coordinates for next iteration
		mpz_mul (data->PP[0].X, jp.X, tmp);
		mpz_mod (data->PP[0].X, data->PP[0].X, p);
		
		mpz_mul (data->PP[0].Y, jp.Y, tmp);
		mpz_mul (data->PP[0].Y, data->PP[0].Y, ZZ);
		mpz_mod (data->PP[0].Y, data->PP[0].Y, p);
	}
	diff3 = diff1 + diff2;
	
	
	if (is_on_curve_jac(&jp, ca, cb)){ 
		printf("\nMULT: kP OK\n");
		print_jacPoint(&jp);
		printf("Elapsed time = %llu %llu %llu nanoseconds\n", (long long unsigned int) diff1, (long long unsigned int) diff2, (long long unsigned int) diff3);
	}
	else
		printf("MULT: kP non OK\n");
	
	
	mpz_clears (k_rd, ZZ, tmp, NULL);
	free_jacPoint(&jp);
	gmp_randclear(rand_gen);
}




