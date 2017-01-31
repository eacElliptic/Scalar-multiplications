#include "aff_point.h"


void init_affPoint(ExtAffPoint *ap){
	mpz_init (ap->x);
	mpz_init (ap->y);
	mpz_init (ap->t);
}

void free_affPoint(ExtAffPoint *ap){
	mpz_clear (ap->x);
	mpz_clear (ap->y);
	mpz_clear (ap->t);
}

void print_affPoint(ExtAffPoint *ap){
	gmp_printf("X : %Zd\nY : %Zd\nT : %Zd\n", ap->x, ap->y, ap->t);
}


/* 
 * For GLV initialization
 * Note: op1 will be updated so that it will contain the result
 * ExtAffine <-- ExtAffine + ExtAffine
 */
void add_aff_aff(ExtAffPoint *op1, ExtAffPoint *op2){
	mpz_t laa, lbb, lcc;
	mpz_inits (laa, lbb, lcc, NULL);
	
	mpz_mul (laa, op1->x, op2->y);
	mpz_mul (lbb, op1->y, op2->x);
	
	mpz_set (lcc, op1->x);			// lcc <-- x1
	mpz_add (op1->x, laa, lbb);	
	mpz_mod (op1->x, op1->x, p);	// x1  <-- (x1*y2 + y1*x2)
	mpz_mul (lbb, laa, lbb);		
	mpz_mul (lbb, lbb, curve_d);		
	mpz_mod (lbb, lbb, p);			// lbb <-- (d*x1*x2*y1*y2)

	mpz_add_ui (laa, lbb, 1);	
	mpz_invert (laa, laa, p);	
	mpz_mul (op1->x, op1->x, laa);		
	mpz_mod (op1->x, op1->x, p);	// x1 <-- (x1*y2+y1*x2)/(1+d*x1*x2*y1*y2)
	
	
	mpz_mul (laa, lcc, op2->x);		// laa <-- x1*x2
	mpz_mul (laa, laa, curve_a);
	mpz_mul (lcc, op1->y, op2->y);
	mpz_sub (op1->y, lcc, laa);
	mpz_mod (op1->y, op1->y, p);	// y1 <-- (y1*y2 - a*x1*x2)
	mpz_ui_sub (laa, 1, lbb);	
	mpz_invert (laa, laa, p);	
	mpz_mul (op1->y, op1->y, laa);		
	mpz_mod (op1->y, op1->y, p);	// y1 <-- (y1*y2-a*x1*x2)/(1-d*x1*x2*y1*y2)
	
	
	mpz_mul (op1->t, op1->x, op1->y);		
	mpz_mod (op1->t, op1->t, p);	// t1 <-- x1 * y1
	
	
	mpz_clears(laa, lbb, lcc, NULL);
}

/*
 * Computes : phi(x ,y) = (beta*x ,1/y)
 * We assume that : a = -1 and d = 1, so eq. : -x^2 + y^2 = 1 + (x^2)*(y^2)
 */
void apply_endo(ExtAffPoint *rop, ExtAffPoint *op, mpz_t beta){
	mpz_mul (rop->x, beta, op->x);
	mpz_mod (rop->x, rop->x, p);
	
	mpz_invert (rop->y, op->y, p);
	
	mpz_mul (rop->t, rop->x, rop->y);		
	mpz_mod (rop->t, rop->t, p);	// t1 <-- x1 * y1
}

//~ eq : a*x^2 + y^2 = 1 + d*(x^2)*(y^2)
int is_on_curve_aff(ExtAffPoint *ap){
	int rep;
	mpz_t leq, req, lx, ly;
	mpz_inits (leq, req, lx, ly, NULL);
	
	mpz_mul (lx, ap->x, ap->x);
	mpz_mod (lx, lx, p); 
	mpz_mul (ly, ap->y, ap->y);
	mpz_mod (ly, ly, p);
	
	mpz_mul (leq, lx, curve_a);
	mpz_add (leq, leq, ly);
	mpz_mod (leq, leq, p);		// leq <-- (a*x^2 + y^2)
	
	mpz_mul (req, lx, ly);
	mpz_mul (req, req, curve_d);
	mpz_add_ui (req, req, 1);
	mpz_mod (req, req, p);		// req <-- (1 + d*(x^2)*(y^2))
	
	
	rep = mpz_cmp (leq, req);
	
	mpz_clears (leq, req, lx, ly, NULL);
	
	return (rep==0);
}








