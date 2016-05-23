#include "aff_point.h"


void init_affPoint(AffPoint *ap){
	mpz_init (ap->X);
	mpz_init (ap->Y);
}

void free_affPoint(AffPoint *ap){
	mpz_clear (ap->X);
	mpz_clear (ap->Y);
}

void print_affPoint(AffPoint *ap){
	gmp_printf("X : %Zd\nY : %Zd\n", ap->X, ap->Y);
}

void add_aff_aff(AffPoint *op1, AffPoint *op2){
	//~ For GLV initialization
	//~ Note: op1 will be updated so that it will contain the result
	//~ Affine <-- Affine + Affine
	
	mpz_t lA, lB, lC;
	mpz_inits (lA, lB, lC, NULL);
	
	mpz_sub (lA, op2->X, op1->X);
	mpz_invert (lA, lA, p);
	mpz_sub (lB, op2->Y, op1->Y);
	mpz_mul (lB, lB, lB);
	mpz_mod (lB, lB, p); 
	mpz_mul (lB, lB, lA);
	mpz_mod (lB, lB, p); 
	mpz_mul (lB, lB, lA);
	mpz_mod (lB, lB, p); 
	mpz_set (lC, op1->X);
	mpz_sub (op1->X, lB, op1->X);
	mpz_sub (op1->X, op1->X, op2->X);
	mpz_mod (op1->X, op1->X, p); 
	
	mpz_sub (lB, op2->Y, op1->Y);
	mpz_mul (lB, lB, lA);
	mpz_mod (lB, lB, p); 
	mpz_sub (lA, lC, op1->X);
	mpz_mul (lA, lB, lA);
	mpz_mod (lA, lA, p); 
	mpz_sub (op1->Y, lA, op1->Y);
	mpz_mod (op1->Y, op1->Y, p); 
	
	mpz_clears(lA, lB, lC, NULL);
}

void apply_endo(AffPoint *rop, AffPoint *op, mpz_t e){
	mpz_set (rop->X, op->X);
	mpz_set (rop->Y, op->Y);
	
	mpz_mul (rop->X, rop->X, e);
	mpz_mod (rop->X, rop->X, p); 
}


int is_on_curve_aff(AffPoint *ap, const mpz_t a, const mpz_t b){
	int rep;
	mpz_t c, e;
	mpz_inits (c, e, NULL);
	
	//~ e= X^3 (p)
	mpz_mul (e, ap->X, ap->X);
	mpz_mul (e, e, ap->X);
	mpz_mod (e, e, p); 
	//~ c= aX+b (p)
	mpz_mul (c, ap->X, a);
	mpz_add (c, c, b);
	mpz_mod (c, c, p); 
	//~ e= X^3+ aX+ b (p)
	mpz_add (e, e, c);
	mpz_mod (e, e, p); 
	//~ c= Y^2 (p)
	mpz_mul (c, ap->Y, ap->Y);
	mpz_mod (c, c, p); 
	
	
	rep = mpz_cmp (c, e);
	
	mpz_clears(c, e, NULL);
	
	return (rep==0);
}








