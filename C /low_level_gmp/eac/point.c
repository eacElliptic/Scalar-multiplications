#include "point.h"


void print_point(Point *pc){
	gmp_printf ("X : %Nu\n", pc->X_limbs, NB_LIMBS);
	gmp_printf ("Y : %Nu\n", pc->Y_limbs, NB_LIMBS);
	gmp_printf ("Z : %Nu\n", Z_limbs, NB_LIMBS);
}

void init_point(Point *pc){
	int j;
	for(j=0; j<NB_LIMBS; j++){
		pc->X_limbs[j] = 0;
		pc->Y_limbs[j] = 0;
	}
}

//~ Note: 'e' is supposed to have NB_LIMBS limbs.
void apply_endo(Point *rop, Point *op, mp_limb_t *e){
	mpn_copyi (rop->Y_limbs, op->Y_limbs, NB_LIMBS);
	
	mpn_mul_n (zaddu_tmp_mult, op->X_limbs, e, NB_LIMBS);
	mpn_tdiv_qr (useless_q, rop->X_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
}


void point_from_eac(Point *p0p1, uchar *eac, uint eac_len){
	uint j;
	for(j=0; j<eac_len; j++){
		zaddu((p0p1+(1-eac[j])), (p0p1+eac[j]));
	}
	zaddu(p0p1, (p0p1+1));
}

//~ Always computes (p1, p1+p2) with p1 updated.
void zaddu(Point *p1, Point *p2){
	//~ A <- X2-X1 
	sub_limbs_n(a_limbs, p2->X_limbs, p1->X_limbs);
	//~ Z <- Z(X2-X1) 
	mpn_mul_n (zaddu_tmp_mult, Z_limbs, a_limbs , NB_LIMBS);
	mpn_tdiv_qr (useless_q, Z_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ A <- (X2-X1)^2 
	mpn_sqr (zaddu_tmp_mult, a_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, a_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ X1 <- X1*(X2-X1)^2 
	mpn_mul_n (zaddu_tmp_mult, p1->X_limbs, a_limbs , NB_LIMBS);
	mpn_tdiv_qr (useless_q, p1->X_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ A <- X2*(X2-X1)^2 
	mpn_mul_n (zaddu_tmp_mult, p2->X_limbs, a_limbs , NB_LIMBS);
	mpn_tdiv_qr (useless_q, a_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ B <- (Y2-Y1)^2 
	sub_limbs_n(p2->Y_limbs, p2->Y_limbs, p1->Y_limbs); // Y2 <- (Y2-Y1)
	mpn_sqr (zaddu_tmp_mult, p2->Y_limbs, NB_LIMBS);  // (Y2-Y1)^2 
	mpn_tdiv_qr (useless_q, b_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ X2 <- B - X1 - A 
	sub_limbs_n(p2->X_limbs, b_limbs, p1->X_limbs); 
	sub_limbs_n(p2->X_limbs, p2->X_limbs, a_limbs); 
	mpn_tdiv_qr (useless_q, p2->X_limbs, 0, p2->X_limbs, NB_LIMBS, p_limbs, NB_LIMBS);
	//~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
	sub_limbs_n(a_limbs, a_limbs, p1->X_limbs);  // A <- (A - X1)
	mpn_mul_n (zaddu_tmp_mult, p1->Y_limbs, a_limbs , NB_LIMBS);
	mpn_tdiv_qr (useless_q, p1->Y_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	//~ Y2 <- Y2*(X1-X2)-Y1
	sub_limbs_n(b_limbs, p1->X_limbs, p2->X_limbs);  // B <- (X1 - X2)
	mpn_mul_n (zaddu_tmp_mult, p2->Y_limbs, b_limbs , NB_LIMBS);
	sub_limbs_2n(zaddu_tmp_mult, zaddu_tmp_mult, p1->Y_limbs);
	mpn_tdiv_qr (useless_q, p2->Y_limbs, 0, zaddu_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
}





int is_on_curve(Point *pc, const mp_limb_t *a, int size_a, const mp_limb_t *b, int size_b){
	int rep;
	mp_limb_t c[(NB_LIMBS*2)];
	mp_limb_t d[(NB_LIMBS*2)];
	mp_limb_t e[(NB_LIMBS*2)];
	mp_limb_t f[(NB_LIMBS*2)];
	
	//~ c= Z^2 (p)
	mpn_sqr (c, Z_limbs, NB_LIMBS);
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




//~ Compute (op1 + op2) and put the result in rop
//~ Note: rop, op1 and op2 are supposed to have NB_LIMBS limbs.
void add_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2){
	add_flag = mpn_add_n (rop, op1, op2, NB_LIMBS);
	while(add_flag){
		add_flag = mpn_add_n (rop, rop, add_overflow_limbs, NB_LIMBS);
	}
}

//~ Compute (op1 - op2) and put the result in rop
//~ Note: rop, op1 and op2 are supposed to have NB_LIMBS limbs.
void sub_limbs_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2){
	sub_flag = mpn_sub_n (rop, op1, op2, NB_LIMBS);
	while(sub_flag){
		mpn_com (rop, rop, NB_LIMBS);
		mpn_add_1 (rop, rop, NB_LIMBS, 1);
		sub_flag = mpn_sub_n (rop, p_limbs, rop, NB_LIMBS);
	}
}

//~ Compute (op1 - op2) and put the result in rop
//~ Note: rop and op1 are supposed to have (NB_LIMBS*2) limbs while op2 is supposed to have NB_LIMBS limbs.
void sub_limbs_2n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2){
	sub_flag = mpn_sub (rop, op1, (NB_LIMBS*2), op2, NB_LIMBS);
	while(sub_flag){
		mpn_com (rop, rop, (NB_LIMBS*2));
		mpn_add_1 (rop, rop, (NB_LIMBS*2), 1);
		sub_flag = mpn_sub_n (rop, p_limbs, rop, NB_LIMBS);
	}
}












