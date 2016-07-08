
#ifndef EAC_TIMING

#define EAC_TIMING


static mpz_t beta;
static mp_limb_t *beta_limbs;

static mp_limb_t ca, cb;


static void init_datas(Point *p0p1);

static void go(Point *p0p1, int nbiter);

static void print_eac(uchar *eac, int eac_len);

#endif



