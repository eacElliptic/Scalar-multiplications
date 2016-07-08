#include "glv_timing.h"
#include <stdint.h> 
#define BILLION 1000000000L




void init_data(GLVData *op){
	mpz_t  tmp1, tmp2, tmp_exp;
	mp_limb_t *t1_limbs, *t2_limbs;
	mpz_inits (tmp1, tmp2, tmp_exp, NULL);
	
	mpz_t  tmp_efp, tmp_refp, tmp_phi;
	mp_limb_t *tmp_efpl, *tmp_refpl, *tmp_phil;
	mpz_inits (tmp_efp, tmp_refp, tmp_phi, NULL);
	
	//~ p = 2^256-1539
	mpz_setbit (p, 256);
	mpz_sub_ui (p, p, 1539);
	p_limbs = mpz_limbs_modify (p, NB_LIMBS);
	
	//~ add_overflow = 2^256-p
	mpz_setbit (add_overflow, 256); // 256 = 64*NB_LIMBS
	mpz_sub (add_overflow, add_overflow, p);
	add_overflow_limbs = mpz_limbs_modify (add_overflow, NB_LIMBS);
	
	//~ Curve is y^2 = x^3 + 5 
	ca = 0;
	cb = 5;
	
	//~ 5 is a generator of Fp and beta is an element of order 3
	mpz_set_ui (beta, 5);
	mpz_sub_ui (tmp_exp, p, 1); 
	mpz_divexact_ui (tmp_exp, tmp_exp, 3);
	mpz_powm (beta, beta, tmp_exp, p);
	beta_limbs = mpz_limbs_modify (beta, NB_LIMBS);
	
	init_aff_point(op->PP);
	init_aff_point((op->PP +1));
	
	
	mpz_set_str (tmp1, "66043340678279369258981193985450715448372246692524118399919799979198175326756", 10);
	mpz_set_str (tmp2, "52931614173969837860927627620185828149756459753631369661231760168280520196499", 10);
	t1_limbs = mpz_limbs_modify (tmp1, NB_LIMBS);
	t2_limbs = mpz_limbs_modify (tmp2, NB_LIMBS);
	
	//~ the point P
	mpn_copyi (op->PP[0].X_limbs, t1_limbs, NB_LIMBS);
	mpn_copyi (op->PP[0].Y_limbs, t2_limbs, NB_LIMBS);
	
	
	mpn_zero (op->efp_limbs, NB_LIMBS);
	mpn_zero (op->refp_limbs, NB_LIMBS);
	mpn_zero (op->phi_limbs, NB_LIMBS);
	mpn_zero (op->aa_limbs, NB_LIMBS);
	mpn_zero (op->neg_bb_limbs, NB_LIMBS);
	mpn_zero (op->Na_limbs, NB_LIMBS);
	
	//~ Order of P
	mpz_set_str (tmp_efp, "115792089237316195423570985008687907852920869663551026044856920202296625078603", 10);
	tmp_efpl = mpz_limbs_modify (tmp_efp, mpz_size(tmp_efp));
	mpn_copyi (op->efp_limbs, tmp_efpl, mpz_size(tmp_efp));
	//~ Integer part of square root of efp
	mpz_set_str (tmp_refp, "340282366920938463463374607431768211455", 10);
	tmp_refpl = mpz_limbs_modify (tmp_refp, mpz_size(tmp_refp));
	mpn_copyi (op->refp_limbs, tmp_refpl, mpz_size(tmp_refp));
	//~ Root of  phi^2+phi+1 = 0 mod #E(Fp)
	mpz_set_str (tmp_phi, "86115113571596370384202877098459685148505633832452548870570058599719872579953", 10);
	tmp_phil = mpz_limbs_modify (tmp_phi, mpz_size(tmp_phi));
	mpn_copyi (op->phi_limbs, tmp_phil, mpz_size(tmp_phi));
	
	init_glvRData(op);
	 
	 
	if (is_on_curve_aff(op->PP, &ca, 1, &cb, 1)) {
		printf("P OK\n");
		print_aff_point(op->PP);
	}
	else
		printf("P not OK\n");
	
	mpz_clears (tmp1, tmp2, tmp_exp, NULL);
	mpz_clears (tmp_efp, tmp_refp, tmp_phi, NULL);
}

//~ Initialises remainning data of 'op'
void init_glvRData(GLVData *op){
	
	mpz_t  u, v, x1, x2, y1, y2, q, r, x, y;
	mpz_t tmp, laa, lbb, lNa, lefp, lrefp, lphi;
	mpz_inits (u, v, x1, x2, y1, y2, q, r, x, y, tmp, laa, lbb, lNa, NULL);
	
	mp_limb_t *laa_limbs, *neg_lbb_limbs, *lNa_limbs;
	
	mpz_roinit_n (lefp, op->efp_limbs, NB_LIMBS);
	mpz_roinit_n (lrefp, op->refp_limbs, NB_LIMBS);
	mpz_roinit_n (lphi, op->phi_limbs, NB_LIMBS);
	
	mpz_set (u, lefp);
	mpz_set (v, lphi);
	
	mpz_set_ui (x1, 1);
	mpz_set_ui (y2, 1);
	
	while(mpz_cmp(u, lrefp) > 0){
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
	
	mpz_add (laa, u, y);		// aa = u + y
	mpz_neg (lbb, y);	 	// bb = -y
	
	mpz_mul (x1, lbb, lbb);		// x1 = (bb)^2
	mpz_mul (x2, u, lbb);			// x2 = u * bb
	mpz_mul (lNa, u, u);			// Na = u^2
	mpz_add (lNa, lNa, x1);		// Na = u^2 + (bb)^2 
	mpz_sub (lNa, lNa, x2);		// Na = u^2 + (bb)^2 - (u * bb)
	
	laa_limbs = mpz_limbs_modify (laa, mpz_size(laa));
	mpn_copyi (op->aa_limbs, laa_limbs, mpz_size(laa));
	neg_lbb_limbs = mpz_limbs_modify (lbb, mpz_size(lbb));
	mpn_copyi (op->neg_bb_limbs, neg_lbb_limbs, mpz_size(lbb));
	lNa_limbs = mpz_limbs_modify (lNa, mpz_size(lNa));
	mpn_copyi (op->Na_limbs, lNa_limbs, mpz_size(lNa));
	
	mpz_clears (u, v, x1, x2, y1, y2, q, r, x, y, tmp, laa, lbb, lNa, NULL);
}



void go(GLVData *data, int nbiter){
	int i;
	mpz_t  ZZ, tmp_inv, k, l_efp;
	mp_limb_t *zz_limbs, *k_limbs;
	mp_limb_t tmp_limbs[(NB_LIMBS*2)];
	mp_limb_t tmp_prod_limbs[(NB_LIMBS*2)];
	mpz_inits (ZZ, k, NULL);
	
	gmp_randstate_t rand_gen;
	unsigned long seed = time(NULL); 
	gmp_randinit_default(rand_gen);
	gmp_randseed_ui(rand_gen, seed);
	
	GLVScalar scl;
	JacPoint jp;
	
	struct timespec start1;
    struct timespec end1;
    uint64_t diff1,diff2, diff3;
	
	
	printf("\nGLV 256 bits benchmark, running...\n\n");
	
	
	diff1 = 0;
	diff2 = 0;
	
	init_jac_point(&jp);
	
	mpz_roinit_n (l_efp, data->efp_limbs, NB_LIMBS);
	
	for (i = 0; i < nbiter; i++) {
		
		mpz_urandomb (k, rand_gen, BIT_SIZE);
		mpz_mod (k, k, l_efp);
		k_limbs = mpz_limbs_modify (k, NB_LIMBS);

		
		apply_endo(((data->PP)+2), (data->PP), beta_limbs); // Point PHI_P 
		
		mpn_copyi (data->PP[3].X_limbs, data->PP[2].X_limbs, NB_LIMBS);
		mpn_copyi (data->PP[3].Y_limbs, data->PP[2].Y_limbs, NB_LIMBS);
		add_aff_aff(((data->PP)+3), (data->PP));      // Point P+PHI_P
		
		mpn_copyi (data->PP[1].X_limbs, data->PP[0].X_limbs, NB_LIMBS);
		mpn_sub_n (data->PP[1].Y_limbs, p_limbs, data->PP[0].Y_limbs, NB_LIMBS); 
		add_aff_aff(((data->PP)+1), ((data->PP)+2));      // Point -P+PHI_P
		
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
		
		init_glvScalar(&scl);
		build_glvScalar(&scl, k_limbs, data->aa_limbs, data->neg_bb_limbs, data->Na_limbs);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
		diff1 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec; 
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start1);
		
		point_from_GLVScalar(&jp, &scl, data->PP);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
		diff2 += BILLION * (end1.tv_sec - start1.tv_sec) + end1.tv_nsec - start1.tv_nsec; 
		
		mpz_roinit_n (tmp_inv, jp.Z_limbs, NB_LIMBS);
		mpz_invert (ZZ, tmp_inv, p);
		zz_limbs = mpz_limbs_modify (ZZ, NB_LIMBS);
		mpn_sqr (tmp_limbs, zz_limbs, NB_LIMBS);  // tmp = ZZ^2
		mpn_tdiv_qr (useless_q, tmp_limbs, 0, tmp_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
       
		//~ P is transformed in affine coordinates for next iteration
		mpn_mul_n (tmp_prod_limbs, jp.X_limbs, tmp_limbs, NB_LIMBS);
		mpn_tdiv_qr (useless_q, data->PP[0].X_limbs, 0, tmp_prod_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
       
		mpn_mul_n (tmp_prod_limbs, jp.Y_limbs, tmp_limbs, NB_LIMBS);
		mpn_tdiv_qr (useless_q, data->PP[0].Y_limbs, 0, tmp_prod_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
		mpn_mul_n (tmp_prod_limbs, data->PP[0].Y_limbs, zz_limbs, NB_LIMBS);
		mpn_tdiv_qr (useless_q, data->PP[0].Y_limbs, 0, tmp_prod_limbs, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	}
	diff3 = diff1 + diff2;
	
	
	if (is_on_curve_jac(&jp, &ca, 1, &cb, 1)){ 
		printf("\nMULT: kP OK\n");
		print_jac_point(&jp);
		printf("Elapsed time = %llu %llu %llu nanoseconds\n", (long long unsigned int) diff1, (long long unsigned int) diff2, (long long unsigned int) diff3);
	}
	else
		printf("MULT: kP not OK\n");
	
	
	mpz_clears (ZZ, k, NULL);
	gmp_randclear(rand_gen);
}




