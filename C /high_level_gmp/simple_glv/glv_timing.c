#include "glv_timing.h"
#include <stdint.h> 
#define BILLION 1000000000L



void init_glvObjects(GLVData *op){
	mpz_inits (op->phi, op->efp, op->refp, op->aa, op->bb, op->Na, NULL);
	
	init_affPoint(op->PP);
	init_affPoint((op->PP)+1);
	init_affPoint((op->PP)+2);
	init_affPoint((op->PP)+3);

}

void free_glvObjects(GLVData *op){
	mpz_clears (op->phi, op->efp, op->refp, op->aa, op->bb, op->Na, NULL);

	free_affPoint(op->PP);
	free_affPoint((op->PP)+1);
	free_affPoint((op->PP)+2);
	free_affPoint((op->PP)+3);
}

void init_data(GLVData *op){
	//~ Note: 'op' is supposed to be already initialised (see main.c)
	mpz_t  X, Y, tmp_exp;
	mpz_inits (X, Y, tmp_exp, NULL);
	
	//~ p = 2^256-1539
	mpz_setbit (p, 256);
	mpz_sub_ui (p, p, 1539);
	//~ Curve is y^2 = x^3 + 5 
	//~ mpz_set_ui (ca, 0);
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
		printf("P not OK\n");;
	
	mpz_clears (X, Y, tmp_exp, NULL);
}

void init_glvRData(GLVData *op){
	//~ Initialises remainning data of 'op'
	
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


void go(GLVData *data, int nbiter){
	int i;
	mpz_t  ZZ, tmp;
	mpz_inits (ZZ, tmp, NULL);
	
	GLVScalar k;
	JacPoint jp;
	
	struct timespec start1;
    struct timespec end1;
    uint64_t diff1,diff2,diff3;
	
	init_jacPoint(&jp);
	
	printf("\nGLV 256 bits benchmark, running...\n");
	
	diff1 = 0;
	diff2 = 0;
	
	for (i = 0; i < nbiter; i++) {
		apply_endo(((data->PP)+2), (data->PP), beta); /* Point PHI_P */
		
		mpz_set (data->PP[3].X, data->PP[2].X);
		mpz_set (data->PP[3].Y, data->PP[2].Y);
		add_aff_aff(((data->PP)+3), (data->PP));     /* Point P+PHI_P */
		
		mpz_set (data->PP[1].X, data->PP[0].X);
		mpz_neg (tmp, data->PP[0].Y);
		mpz_set (data->PP[1].Y, tmp);
		add_aff_aff(((data->PP)+1), ((data->PP)+2));  /* Point -P+PHI_P */
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
		
		init_glvScalar(&k);
		build_glvScalar(&k, BIT_SIZE, data->aa, data->bb, data->Na, data->efp);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);   
		diff1 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec; 
		
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
		
		point_from_GLVScalar(&jp, &k, data->PP);
		
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
	
	diff3=diff2+diff1;
	
	
	if (is_on_curve_jac(&jp, ca, cb)){ 
		printf("\nMULT: kP OK\n");
		print_jacPoint(&jp);
		printf("Elapsed time = %llu %llu %llu nanoseconds\n", (long long unsigned int) diff1, (long long unsigned int) diff2, (long long unsigned int) diff3);
	}
	else{
		printf("MULT: kP not OK\n");
		print_jacPoint(&jp);
		printf("Elapsed time = %llu %llu %llu nanoseconds\n", (long long unsigned int) diff1, (long long unsigned int) diff2, (long long unsigned int) diff3);
	}
	
	
	mpz_clears (ZZ, tmp, NULL);
	free_jacPoint(&jp);
}




