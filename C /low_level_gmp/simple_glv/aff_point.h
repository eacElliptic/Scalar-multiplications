#ifndef AFF_POINT

#define AFF_POINT 


static void init_aff_point(AffPoint *pc);

static void print_aff_point(AffPoint *pc);


static int is_on_curve_aff(AffPoint *pc, const mp_limb_t *a, int size_a, const mp_limb_t *b, int size_b);



static void apply_endo(AffPoint *rop, AffPoint *op, mp_limb_t *e);

static void add_aff_aff(AffPoint *op1, AffPoint *op2);



#endif


