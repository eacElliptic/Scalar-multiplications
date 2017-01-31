#include "point.h"


void init_point(Point *pc){
	mpz_init (pc->X);
	mpz_init (pc->Y);
}

void free_point(Point *pc){	
	mpz_clear (pc->X);
	mpz_clear (pc->Y);
}

void print_point(Point *pc){
	gmp_printf("X : %Zd\n", pc->X);
	gmp_printf("Y : %Zd\n", pc->Y);
	gmp_printf("Z : %Zd\n", Z);
}

void apply_endo(Point *rop, Point *op, mpz_t e){
	mpz_mul (rop->X, op->X, e);
	mpz_mod (rop->X, rop->X, p); 
	
	mpz_set (rop->Y, op->Y);
}



void point_from_eac(Point *p0p1, uchar *eac, uint eac_len){
	uint j;
	for(j=0; j<eac_len; j++){
		safe_perm(p0p1, (p0p1+1), (mp_limb_t) (- eac[j]));
		
		zaddu(p0p1, (p0p1+1));
	}
	zaddu(p0p1, (p0p1+1));
}


/**
 * Permuts p1 and p2 according the value of 'mask' 
 * Assumes elements (p1->X, p1->Y, p2->X and p1->Y) have NB_LIMBS limbs.
 */
void safe_perm(Point *p1, Point *p2, mp_limb_t mask){
	int j;
	mp_limb_t l_xor;
	
	/* for X coordinates */
	//~ IMPORTANT : including '_mp_alloc' fields bring errors and incoherences, as we don't make any allocation.
	l_xor = p1->X->_mp_size ^ p2->X->_mp_size;
	p2->X->_mp_size = (mask & l_xor) ^ p1->X->_mp_size;
	p1->X->_mp_size = l_xor ^ p2->X->_mp_size;
	for(j=0; j<NB_LIMBS; j++) {
		l_xor = p1->X->_mp_d[j] ^ p2->X->_mp_d[j];
		p2->X->_mp_d[j] = (mask & l_xor) ^ p1->X->_mp_d[j];
		p1->X->_mp_d[j] = l_xor ^ p2->X->_mp_d[j];
	}
	
	/* for Y coordinates */
	//~ IMPORTANT : idem as for X coordinates about '_mp_alloc'.
	l_xor = p1->Y->_mp_size ^ p2->Y->_mp_size;
	p2->Y->_mp_size = (mask & l_xor) ^ p1->Y->_mp_size;
	p1->Y->_mp_size = l_xor ^ p2->Y->_mp_size;
	for(j=0; j<NB_LIMBS; j++) {
		l_xor = p1->Y->_mp_d[j] ^ p2->Y->_mp_d[j];
		p2->Y->_mp_d[j] = (mask & l_xor) ^ p1->Y->_mp_d[j];
		p1->Y->_mp_d[j] = l_xor ^ p2->Y->_mp_d[j];
	}
}



//~ Always computes (p1, p1+p2) with p1 updated.
void zaddu(Point *p1, Point *p2){
	//~ A <- X2-X1 
	mpz_sub (A, p2->X, p1->X);
	//~ Z <- Z(X2-X1) 
	mpz_mul (Z, Z, A);
	mpz_mod (Z, Z, p); 
	//~ A <- (X2-X1)^2 
	mpz_mul (A, A, A);
	mpz_mod (A, A, p); 
	//~ X1 <- X1*(X2-X1)^2 
	mpz_mul (p1->X, p1->X, A);
	mpz_mod (p1->X, p1->X, p); 
	//~ A <- X2*(X2-X1)^2 
	mpz_mul (A, p2->X, A);
	mpz_mod (A, A, p);
	//~ Y2 <- (Y2-Y1)^2 
	mpz_sub (p2->Y, p2->Y, p1->Y);  	// Y2 <- (Y2-Y1)
	//~ X2 <- (Y2-Y1)^2  - X1 - A 
	mpz_mul (p2->X, p2->Y, p2->Y);		// X2 <- (Y2-Y1)^2
	mpz_mod (p2->X, p2->X, p);
	mpz_sub (p2->X, p2->X, p1->X);
	mpz_sub (p2->X, p2->X, A);
	mpz_mod (p2->X, p2->X, p); 
	//~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
	mpz_sub (A, A, p1->X);  		// A <- (A-X1)
	mpz_mul (p1->Y, p1->Y, A);
	mpz_mod (p1->Y, p1->Y, p);
	//~ Y2 <- Y2*(X1-X2)-Y1
	mpz_sub (A, p1->X, p2->X);  	// A <- (X1-X2)
	mpz_mul (p2->Y, p2->Y, A);
	mpz_sub (p2->Y, p2->Y, p1->Y);
	mpz_mod (p2->Y, p2->Y, p);
}



int is_on_curve_aff(Point *pc, const mpz_t a, const mpz_t b){
	int rep;
	mpz_t c, e;
	mpz_inits (c, e, NULL);
	
	//~ e= X^3 (p)
	mpz_mul (e, pc->X, pc->X);
	mpz_mul (e, e, pc->X);
	mpz_mod (e, e, p); 
	//~ c= aX+b (p)
	mpz_mul (c, pc->X, a);
	mpz_add (c, c, b);
	mpz_mod (c, c, p); 
	//~ e= X^3+ aX+ b (p)
	mpz_add (e, e, c);
	mpz_mod (e, e, p); 
	//~ c= Y^2 (p)
	mpz_mul (c, pc->Y, pc->Y);
	mpz_mod (c, c, p); 
	
	
	rep = mpz_cmp (c, e);
	
	mpz_clears(c, e, NULL);
	
	return (rep==0);
}


int is_on_curve(Point *pc, const mpz_t a, const mpz_t b){
	int rep;
	mpz_t c, d, e;
	mpz_inits (c, d, e, NULL);

	//~ c= Z^2 (p)
	mpz_mul (c, Z, Z);
	mpz_mod (c, c, p);
	//~ d= Z^4 (p)
	mpz_mul (d, c, c);
	mpz_mod (d, d, p);
	//~ c = bZ^6 (p)
	mpz_mul (c, c, d);
	mpz_mul (c, c, b);
	mpz_mod (c, c, p);
	//~ d=aXZ^4
	mpz_mul (d, d, a);
	mpz_mul (d, d, pc->X);
	mpz_mod (d, d, p);
	//~ e= X^3 (p)
	mpz_mul (e, pc->X, pc->X);
	mpz_mul (e, e, pc->X);
	mpz_mod (e, e, p);
	//~ e = X^3 +aXZ^4 +bZ^6 (p)
	mpz_add (e, e, c);
	mpz_add (e, e, d);
	mpz_mod (e, e, p);
	//~ c= Y^2 (p)
	mpz_mul (c, pc->Y, pc->Y);
	mpz_mod (c, c, p);
	
	rep = mpz_cmp (c, e);
	
	mpz_clears(c, d, e, NULL);
	
	return (rep==0);
}


