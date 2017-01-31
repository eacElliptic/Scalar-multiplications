#include "glv_timing.h"
#include <stdint.h> 
#define BILLION 1000000000L



void init_glvObjects(GLVData *op){
	mpz_inits (op->phi, op->efp, op->refp, op->aa, op->bb, op->Na, NULL);
	
	init_affPoint(op->PP);
	init_affPoint((op->PP)+1);
}

void free_glvObjects(GLVData *op){
	mpz_clears (op->phi, op->efp, op->refp, op->aa, op->bb, op->Na, NULL);

	free_affPoint(op->PP);
	free_affPoint((op->PP)+1);
}

void init_data(GLVData *op){
	//~ Note: 'op' is supposed to be already initialised (see main.c)
	mpz_t  X, Y, tmp_exp;
	mpz_inits (X, Y, tmp_exp, NULL);
	
	//~ p = 2^256-43443
	mpz_setbit (p, 256);
	mpz_sub_ui (p, p, 43443);
	//~ Curve is : -(x^2) + y^2 = 1 + (x^2)*(y^2)
	mpz_set_si (curve_a, -1);
	mpz_set_ui (curve_d, 1);
	mpz_set_str (X, "22174792725782664025844270156401847138026886831119649840138701280259288954423", 10);
	mpz_set_str (Y, "87737099222887081277295330229931648697203685991117339269318803967832062780143", 10);
	//~ Point P
	mpz_set (op->PP[0].x , X);
	mpz_set (op->PP[0].y , Y);
	mpz_mul (op->PP[0].t, op->PP[0].x, op->PP[0].y);
	mpz_mod (op->PP[0].t, op->PP[0].t, p);
	//~ beta is an element of order 4
	mpz_set_str (beta, "9058966984510276008007633023999327412385577717625298889713350828587238576126", 10);
	
	//~ Order of P
	mpz_set_str (op->efp, "14474011154664524427946373126085988481582509411733023890661985475889632981401", 10);
	//~ Integer part of square root of efp
	mpz_set_str (op->refp, "120307984584002255772516886238812528463", 10);
	//~ Root of  phi^2 + 1 = 0 mod #E(Fp)
	mpz_set_str (op->phi, "10749260569431236026102217475317958236773166001345671471970756507155106163337", 10);
	
	init_glvRData(op);
	 
	if (is_on_curve_aff(op->PP)) {
		printf("P OK\n");
		print_affPoint(op->PP);
	}
	else
		printf("P non OK\n");;
	
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
	
	mpz_set (op->aa, u);		// aa = u 
	mpz_neg (op->bb, y);	 	// bb = -y
	
	mpz_mul (x1, op->bb, op->bb);		// x1 = (bb)^2
	mpz_mul (op->Na, u, u);				// Na = u^2
	mpz_add (op->Na, op->Na, x1);		// Na = u^2 + (bb)^2 
	
	mpz_clears (u, v, x1, x2, y1, y2, q, r, x, y, tmp, NULL);
}


void go(GLVData *data, int nbiter){
	int i;
	mpz_t  k_rd, ZZ;
	mpz_inits (k_rd, ZZ, NULL);
	
	gmp_randstate_t rand_gen;
	unsigned long seed = time(NULL); 
	gmp_randinit_default(rand_gen);
	gmp_randseed_ui(rand_gen, seed);

	
	GLVScalar k;
	ExtProjPoint jp;
	
	struct timespec start1;
    struct timespec end1;
    uint64_t diff1,diff2,diff3;
    
    
	printf("\nSGLV 256/253 bits benchmark, running...\n\n");
	
	diff1 = 0;
	diff2 = 0;
	
	init_extProjPoint(&jp);
	
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
		
		point_from_GLVScalar(&jp, &k, data->PP);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
		diff2 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec;
		
		
		mpz_invert (ZZ, jp.Z, p);
		//~ //P is transformed in affine coordinates for next iteration
		mpz_mul (data->PP[0].x, jp.X, ZZ);
		mpz_mod (data->PP[0].x, data->PP[0].x, p);
		
		mpz_mul (data->PP[0].y, jp.Y, ZZ);
		mpz_mod (data->PP[0].y, data->PP[0].y, p);
		
		mpz_mul (data->PP[0].t, data->PP[0].x, data->PP[0].y);
		mpz_mod (data->PP[0].t, data->PP[0].t, p);
	}
	
	diff3=diff2+diff1;
	
	
	if (is_on_curve_extProj(&jp)){ 
		printf("\nMULT: kP OK\n");
		print_extProjPoint(&jp);
		printf("Elapsed time = %llu %llu %llu nanoseconds\n", (long long unsigned int) diff1, (long long unsigned int) diff2, (long long unsigned int) diff3);
	}
	else
		printf("MULT: kP non OK\n");
	
	
	mpz_clears (k_rd, ZZ, NULL);
	free_extProjPoint(&jp);
	gmp_randclear(rand_gen);
}




