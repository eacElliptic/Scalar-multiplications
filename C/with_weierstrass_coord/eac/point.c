#include "point.h"


void init_point(Point *pc){
	mpz_init_set_ui (pc->Z, 1);
	mpz_init (pc->PQ[0]);
	mpz_init (pc->PQ[1]);
	mpz_init (pc->PQ[2]);
	mpz_init (pc->PQ[3]);
}

void free_point(Point *pc){	
	mpz_clear (pc->Z);
	mpz_clear (pc->PQ[0]);
	mpz_clear (pc->PQ[1]);
	mpz_clear (pc->PQ[2]);
	mpz_clear (pc->PQ[3]);
}

void print_point(Point *pc, char j){
	gmp_printf("X : %Zd\nY : %Zd\nZ : %Zd\n", pc->PQ[j], pc->PQ[j+2], pc->Z);
}

void zaddu(Point *pc, uchar b){
	//~ C <- X1-X2 
	mpz_sub (C, pc->PQ[1-b], pc->PQ[b]); 
	//~ Z3 <- Z(X1-X2) 
	mpz_mul (pc->Z, pc->Z, C);
	mpz_mod (pc->Z, pc->Z, p); 
	//~ C <- (X1-X2)^2 
	mpz_mul (C, C, C);
	mpz_mod (C, C, p); 
	//~ W2 <- X2C
	mpz_mul (W, pc->PQ[b], C);
	mpz_mod (W, W, p); 
	//~ W1 <- X1C
	mpz_mul (pc->PQ[0], pc->PQ[1-b], C);
	mpz_mod (pc->PQ[0], pc->PQ[0], p); 
	//~ Y1-Y2
	mpz_sub (C, pc->PQ[3-b], pc->PQ[2+b]); 
	//~ D
	mpz_mul (pc->PQ[1], C, C);
	mpz_mod (pc->PQ[1], pc->PQ[1], p); 
	mpz_sub (pc->PQ[1], pc->PQ[1], W); 
	//~ X3 <- D -W1 -W2 
	mpz_sub (pc->PQ[1], pc->PQ[1], pc->PQ[0]); 
	mpz_mod (pc->PQ[1], pc->PQ[1], p);
	
	mpz_sub (W, pc->PQ[0], W); 
	//~ A1 <- Y1(W1-W2) 
	mpz_mul (pc->PQ[2], pc->PQ[3-b], W);
	mpz_mod (pc->PQ[2], pc->PQ[2], p);
	
    mpz_sub (W, pc->PQ[0], pc->PQ[1]); 
    mpz_mul (C, C, W);
	mpz_mod (C, C, p);
	//~ Y3 <- (Y1-Y2)(W1-X3) -A1
    mpz_sub (pc->PQ[3], C, pc->PQ[2]); 
	mpz_mod (pc->PQ[3], pc->PQ[3], p);
}

void point_from_eac(Point *pc, uchar *eac, uint eac_len){
	uint j;
	for(j=0; j<eac_len; j++){
		zaddu(pc, eac[j]);
	}
	zaddu(pc, 1);
}


int is_on_curve(Point *pc, char j, const mpz_t a, const mpz_t b){
	int rep;
	mpz_t c, d, e;
	mpz_inits (c, d, e, NULL);

	//~ c= Z^2 (p)
	mpz_mul (c, pc->Z, pc->Z);
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
	mpz_mul (d, d, pc->PQ[j]);
	mpz_mod (d, d, p);
	//~ e= X^3 (p)
	mpz_mul (e, pc->PQ[j], pc->PQ[j]);
	mpz_mul (e, e, pc->PQ[j]);
	mpz_mod (e, e, p);
	//~ e = X^3 +aXZ^4 +bZ^6 (p)
	mpz_add (e, e, c);
	mpz_add (e, e, d);
	mpz_mod (e, e, p);
	//~ c= Y^2 (p)
	mpz_mul (c, pc->PQ[j+2], pc->PQ[j+2]);
	mpz_mod (c, c, p);
	
	rep = mpz_cmp (c, e);
	
	mpz_clears(c, d, e, NULL);
	
	return (rep==0);
}


