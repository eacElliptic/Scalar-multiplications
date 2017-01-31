#ifndef JAC_POINT

#define JAC_POINT 


static void init_jacPoint(JacPoint *jp);
static void free_jacPoint(JacPoint *jp);

static void print_jacPoint(JacPoint *jp);

static void dblu(JacPoint *jp);

static void add_jac_aff(JacPoint *jp, mpz_t QX, mpz_t QY);

static void point_from_GLVScalar(JacPoint *rop, GLVScalar *k, AffPoint *pp);

static void double_and_add(JacPoint *rop, AffPoint *op, mpz_t k);

static int is_on_curve_jac(JacPoint *jp, const mpz_t a, const mpz_t b);


#endif


