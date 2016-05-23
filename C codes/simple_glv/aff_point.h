#ifndef AFF_POINT

#define AFF_POINT 


static void init_affPoint(AffPoint *ap);
static void free_affPoint(AffPoint *ap);

static void print_affPoint(AffPoint *ap);

static void add_aff_aff(AffPoint *op1, AffPoint *op2);

static void apply_endo(AffPoint *rop, AffPoint *op, mpz_t e);


static int is_on_curve_aff(AffPoint *ap, const mpz_t a, const mpz_t b);



#endif


