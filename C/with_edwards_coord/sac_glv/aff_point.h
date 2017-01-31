#ifndef AFF_POINT

#define AFF_POINT 


static void init_affPoint(ExtAffPoint *ap);
static void free_affPoint(ExtAffPoint *ap);

static void print_affPoint(ExtAffPoint *ap);

static void add_aff_aff(ExtAffPoint *op1, ExtAffPoint *op2);

static void apply_endo(ExtAffPoint *rop, ExtAffPoint *op, mpz_t e);


static int is_on_curve_aff(ExtAffPoint *ap);



#endif


