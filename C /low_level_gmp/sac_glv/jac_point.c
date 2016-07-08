#include "jac_point.h"


void print_jac_point(JacPoint *pc){
	gmp_printf ("X : %Nu\n", pc->X_limbs, NB_LIMBS);
	gmp_printf ("Y : %Nu\n", pc->Y_limbs, NB_LIMBS);
	gmp_printf ("Z : %Nu\n", pc->Z_limbs, NB_LIMBS);
}

void init_jac_point(JacPoint *pc){
	int j;
	for(j=0; j<NB_LIMBS; j++){
		pc->X_limbs[j] = 0;
		pc->Y_limbs[j] = 0;
		pc->Z_limbs[j] = 0;
	}
	pc->Z_limbs[0] = 1;
}

int is_on_curve_jac(JacPoint *pc, const mp_limb_t *a, int size_a, const mp_limb_t *b, int size_b){
	int rep;
	mp_limb_t c[(NB_LIMBS*2)];
	mp_limb_t d[(NB_LIMBS*2)];
	mp_limb_t e[(NB_LIMBS*2)];
	mp_limb_t f[(NB_LIMBS*2)];
	
	//~ c= Z^2 (p)
	mpn_sqr (c, pc->Z_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, c, 0, c, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ d= Z^4 (p)
	mpn_sqr (d, c, NB_LIMBS);
	mpn_tdiv_qr (useless_q, d, 0, d, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ c = bZ^6 (p)
	mpn_mul_n (e, c, d, NB_LIMBS);
	mpn_tdiv_qr (useless_q, e, 0, e, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul (c, e, NB_LIMBS, b, size_b);
	mpn_tdiv_qr (useless_q, c, 0, c, (NB_LIMBS+size_b), p_limbs, NB_LIMBS);
	//~ d=aXZ^4
	mpn_mul (e, d, NB_LIMBS, a, size_a);
	mpn_tdiv_qr (useless_q, e, 0, e, (NB_LIMBS+size_a), p_limbs, NB_LIMBS);
	mpn_mul_n (d, e, pc->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, d, 0, d, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ e= X^3 (p)
	mpn_sqr (f, pc->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, f, 0, f, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (e, f, pc->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, e, 0, e, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	
	//~ e = X^3 +aXZ^4 +bZ^6 (p)
	add_limbs_n(e, e, c);
	add_limbs_n(e, e, d);
	mpn_tdiv_qr (useless_q, e, 0, e, NB_LIMBS, p_limbs, NB_LIMBS);
	
	//~ c= Y^2 (p)
	mpn_sqr (c, pc->Y_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, c, 0, c, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	
	rep = mpn_cmp (c, e, NB_LIMBS);
	
	return (rep==0);
}


void dblu(JacPoint *jp){
    //~ "dbl-2009-l" hyperelliptic.org formulas 
    
    //~ Z3 <- 2*Y1*Z1
    mpn_mul_n (op_tmp_mult, jp->Z_limbs, jp->Y_limbs, NB_LIMBS);
    op_tmp_mult[(NB_LIMBS*2)] = mpn_lshift (op_tmp_mult, op_tmp_mult, (NB_LIMBS*2), 1);
    mpn_tdiv_qr (useless_q, jp->Z_limbs, 0, op_tmp_mult, ((NB_LIMBS*2)+1), p_limbs, NB_LIMBS);
	//~ A <- X1^2
	mpn_sqr (op_tmp_mult, jp->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, A, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ B <- Y1^2
	mpn_sqr (op_tmp_mult, jp->Y_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, B, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ C <- B^2
	mpn_sqr (op_tmp_mult, B, NB_LIMBS);
	mpn_tdiv_qr (useless_q, C, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ D <- X+B
	add_limbs_n(B, B, jp->X_limbs);
	//~ D <- 2*(B^2-A-C)
	mpn_sqr (op_tmp_mult, B, NB_LIMBS);
	mpn_tdiv_qr (useless_q, B, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	sub_limbs_n(B, B, A);
	sub_limbs_n(B, B, C);
    lshift_limbs_n(B, B);
	mpn_tdiv_qr (useless_q, B, 0, B, NB_LIMBS, p_limbs, NB_LIMBS);
	//~ E <- 3A
	lshift_limbs_n(jp->Y_limbs, A);
	add_limbs_n(jp->Y_limbs, jp->Y_limbs, A);
	//~ F <- E^2
	mpn_sqr (op_tmp_mult, jp->Y_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, jp->X_limbs, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ X3 <- F -2*D
	lshift_limbs_n(A, B);
	sub_limbs_n(jp->X_limbs, jp->X_limbs, A);
	mpn_tdiv_qr (useless_q, jp->X_limbs, 0, jp->X_limbs, NB_LIMBS, p_limbs, NB_LIMBS);
	
	//~ Y3 <- E*(D-X3)-8*C
	sub_limbs_n(A, B, jp->X_limbs);
	mpn_mul_n (op_tmp_mult, A, jp->Y_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, jp->Y_limbs, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	
	lshift_limbs_n(A, C);
	lshift_limbs_n(A, A);
	lshift_limbs_n(A, A);  // A <- 8*C
	sub_limbs_n(jp->Y_limbs, jp->Y_limbs, A);
	mpn_tdiv_qr (useless_q, jp->Y_limbs, 0, jp->Y_limbs, NB_LIMBS, p_limbs, NB_LIMBS);
}


//~ Jacobian <-- Jacobian + Affine
//~ Note: Z coordinate of Q must be 1 (affine coordinate) 
//~ Note: jp will be updated so that it will contain the result
void add_jac_aff(JacPoint *jp, mp_limb_t *QX, mp_limb_t *QY){
	//~ "madd-2004-hmv" hyperelliptic.org formulas

	mpn_sqr (op_tmp_mult, jp->Z_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, A, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (op_tmp_mult, A, jp->Z_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, B, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (op_tmp_mult, A, QX, NB_LIMBS);
	mpn_tdiv_qr (useless_q, A, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (op_tmp_mult, B, QY, NB_LIMBS);
	mpn_tdiv_qr (useless_q, B, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	sub_limbs_n(A, A, jp->X_limbs);
	sub_limbs_n(B, B, jp->Y_limbs);

	mpn_mul_n (op_tmp_mult, jp->Z_limbs, A, NB_LIMBS);
	mpn_tdiv_qr (useless_q, jp->Z_limbs, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	
	mpn_sqr (op_tmp_mult, A, NB_LIMBS);
	mpn_tdiv_qr (useless_q, C, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (op_tmp_mult, C, A, NB_LIMBS);
	mpn_tdiv_qr (useless_q, D, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (op_tmp_mult, C, jp->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, C, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	lshift_limbs_n(A, C);
	mpn_sqr (op_tmp_mult, B, NB_LIMBS);
	mpn_tdiv_qr (useless_q, jp->X_limbs, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	sub_limbs_n(jp->X_limbs, jp->X_limbs, A);
	sub_limbs_n(jp->X_limbs, jp->X_limbs, D);
	mpn_tdiv_qr (useless_q, jp->X_limbs, 0, jp->X_limbs, NB_LIMBS, p_limbs, NB_LIMBS);
	
	sub_limbs_n(C, C, jp->X_limbs);
	mpn_mul_n (op_tmp_mult, C, B, NB_LIMBS);
	mpn_tdiv_qr (useless_q, C, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (op_tmp_mult, D, jp->Y_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, D, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	sub_limbs_n(jp->Y_limbs, C, D);
	mpn_tdiv_qr (useless_q, jp->Y_limbs, 0, jp->Y_limbs, NB_LIMBS, p_limbs, NB_LIMBS);
}



//~ Note: PP contains P and P+PHI_P
void point_from_SGLVScalar(JacPoint *rop, SGLVScalar *k, AffPoint *pp){
	
	int j,u;
	mp_limb_t neg_pp_Y[2][NB_LIMBS];
	
	j = k->j;
	
	if(k->tk2[j] < 0)
		u = - (k->tk2[j]);
	else
		u = k->tk2[j];
	
	mpn_copyi (rop->X_limbs, pp[u].X_limbs, NB_LIMBS);
	mpn_zero (rop->Z_limbs, NB_LIMBS);
	rop->Z_limbs[0] = 1;
	
	mpn_sub_n (neg_pp_Y[0], p_limbs, pp[0].Y_limbs, NB_LIMBS); // rop <- -pp[0].Y (p)
	mpn_sub_n (neg_pp_Y[1], p_limbs, pp[1].Y_limbs, NB_LIMBS); // rop <- -pp[1].Y (p)
	
	if(k->tk1[j] < 0){
		mpn_copyi (rop->Y_limbs, neg_pp_Y[u], NB_LIMBS);
	}
	else{
		mpn_copyi (rop->Y_limbs, pp[u].Y_limbs, NB_LIMBS);
	}
	
	j = j - 1;
	while (j >= 0){
		dblu(rop);
		
		if(k->tk2[j] < 0)
			u = - (k->tk2[j]);
		else
			u = k->tk2[j];
			
		if (k->tk1[j] < 0)
			add_jac_aff(rop, pp[u].X_limbs, neg_pp_Y[u]);
		else
			add_jac_aff(rop, pp[u].X_limbs, pp[u].Y_limbs);
		
		j = j - 1;
	}
	if(k->isEven == 0)
		add_jac_aff(rop, pp[0].X_limbs, pp[0].Y_limbs);
}








