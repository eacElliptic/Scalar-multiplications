#include "aff_point.h"


void print_aff_point(AffPoint *pc){
	gmp_printf ("X : %Nu\n", pc->X_limbs, NB_LIMBS);
	gmp_printf ("Y : %Nu\n", pc->Y_limbs, NB_LIMBS);
}

void init_aff_point(AffPoint *pc){
	int j;
	for(j=0; j<NB_LIMBS; j++){
		pc->X_limbs[j] = 0;
		pc->Y_limbs[j] = 0;
	}
}

//~ Note: size_a <= NB_LIMBS and size_b <= NB_LIMBS
int is_on_curve_aff(AffPoint *pc, const mp_limb_t *a, int size_a, const mp_limb_t *b, int size_b){
	int rep;
	mp_limb_t c[(NB_LIMBS*2)];
	mp_limb_t d[(NB_LIMBS*2)];
	mp_limb_t e[(NB_LIMBS*2)];
	mp_limb_t f[(NB_LIMBS*2)];
	
	//~ d=aX
	mpn_mul (d, pc->X_limbs, NB_LIMBS, a, size_a);
	mpn_tdiv_qr (useless_q, d, 0, d, (NB_LIMBS+size_a), p_limbs, NB_LIMBS);
	//~ e= X^3 (p)
	mpn_sqr (f, pc->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, f, 0, f, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	mpn_mul_n (e, f, pc->X_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, e, 0, e, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	
	//~ e = X^3 + aX + b (p)
	add_limbs_n(e, e, d);
	mpn_zero (c, (NB_LIMBS*2));
	mpn_copyi (c, b, size_b);  // c <- (b, size_b)
	add_limbs_n(e, e, c);
	mpn_tdiv_qr (useless_q, e, 0, e, NB_LIMBS, p_limbs, NB_LIMBS);
	
	//~ c= Y^2 (p)
	mpn_sqr (c, pc->Y_limbs, NB_LIMBS);
	mpn_tdiv_qr (useless_q, c, 0, c, (NB_LIMBS*2), p_limbs, NB_LIMBS);
	
	rep = mpn_cmp (c, e, NB_LIMBS);
	
	return (rep==0);
}

//~ Note: e is supposed to have NB_LIMBS limbs.
void apply_endo(AffPoint *rop, AffPoint *op, mp_limb_t *e){
	mpn_copyi (rop->Y_limbs, op->Y_limbs, NB_LIMBS);
	
	mpn_mul_n (op_tmp_mult, op->X_limbs, e, NB_LIMBS);
	mpn_tdiv_qr (useless_q, rop->X_limbs, 0, op_tmp_mult, (NB_LIMBS*2), p_limbs, NB_LIMBS);
}

//~ For GLV initialization
//~ Affine <-- Affine + Affine
//~ 'op1' will be updated so that it will contain the result
void add_aff_aff(AffPoint *op1, AffPoint *op2){
	
	mpz_t l1X, l1Y, l2X, l2Y;
	mpz_roinit_n (l1X, op1->X_limbs, NB_LIMBS);
	mpz_roinit_n (l1Y, op1->Y_limbs, NB_LIMBS);
	mpz_roinit_n (l2X, op2->X_limbs, NB_LIMBS);
	mpz_roinit_n (l2Y, op2->Y_limbs, NB_LIMBS);
	
	mpz_t lA, lB, lC, lop1X, lop1Y;
	mp_limb_t *lop1X_limbs, *lop1Y_limbs;
	mpz_inits (lA, lB, lC, lop1X, lop1Y, NULL);
	
	mpz_set(lop1X, l1X);
	mpz_set(lop1Y, l1Y);
	
	
	mpz_sub (lA, l2X, l1X);
	mpz_invert (lA, lA, p);
	mpz_sub (lB, l2Y, l1Y);
	mpz_mul (lB, lB, lB);
	mpz_mod (lB, lB, p); 
	mpz_mul (lB, lB, lA);
	mpz_mod (lB, lB, p); 
	mpz_mul (lB, lB, lA);
	mpz_mod (lB, lB, p); 
	mpz_set (lC, l1X);
	mpz_sub (lop1X, lB, l1X);
	mpz_sub (lop1X, lop1X, l2X);
	mpz_mod (lop1X, lop1X, p); 
	
	mpz_sub (lB, l2Y, l1Y);
	mpz_mul (lB, lB, lA);
	mpz_mod (lB, lB, p); 
	mpz_sub (lA, lC, lop1X);
	mpz_mul (lA, lB, lA);
	mpz_mod (lA, lA, p); 
	mpz_sub (lop1Y, lA, l1Y);
	mpz_mod (lop1Y, lop1Y, p); 
	
	
	lop1X_limbs = mpz_limbs_modify (lop1X, NB_LIMBS);
	lop1Y_limbs = mpz_limbs_modify (lop1Y, NB_LIMBS);
	mpn_copyi (op1->X_limbs, lop1X_limbs, NB_LIMBS);
	mpn_copyi (op1->Y_limbs, lop1Y_limbs, NB_LIMBS);
	
	
	mpz_clears(lA, lB, lC, lop1X, lop1Y, NULL);
}













