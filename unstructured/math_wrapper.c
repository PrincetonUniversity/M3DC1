/* Self-contained special functions for math.f90 — no external library
 * and no POSIX/libm special functions required (replaces the former
 * GSL/Boost wrappers).  Only ISO C99 <math.h> primitives are used
 * (exp, log, sin, cos, sqrt, fabs, ceil).
 *
 *   bessel_J, bessel_Y : Cephes j0/y0/j1/y1 kernels + stable recurrences
 *                        (upward for Y_n and J_n with n < x; Miller
 *                        downward with sum normalization for n >= x)
 *   erf                : Cephes erf/erfc rational approximations
 *   elliptic K, E      : arithmetic-geometric mean (AGM) iteration,
 *                        modulus-k convention (same as GSL Kcomp/Ecomp)
 *   bessel_I           : Miller downward recurrence, normalized with
 *                        exp(x) = I_0(x) + 2*sum_{k>=1} I_k(x)
 *
 * The J0/Y0/J1/Y1 and erf/erfc kernels are from the Cephes Math Library
 * (Release 2.8, Copyright 1984-2000 Stephen L. Moshier), as vendored by
 * SciPy; peak errors ~5e-16 (bessel), ~6e-14 (erfc far tail).
 * https://moshier.net/#Cephes
 */

#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif
#ifndef M_2_PI
#define M_2_PI 0.63661977236758134308
#endif

#define SQ2OPI 7.9788456080286535587989E-1 /* sqrt(2/pi)   */
#define THPIO4 2.35619449019234492885      /* 3*pi/4       */
#define MAXLOG 7.09782712893383996732E2    /* log(DBL_MAX) */

/*==================================================================
 * Polynomial evaluation (Cephes polevl.c)
 *==================================================================*/
static double polevl(double x, const double coef[], int N)
{
  double ans;
  int i;
  const double *p = coef;

  ans = *p++;
  i = N;
  do
    ans = ans * x + *p++;
  while (--i);
  return ans;
}

/* Same, with an implied leading coefficient of 1.0 */
static double p1evl(double x, const double coef[], int N)
{
  double ans;
  int i;
  const double *p = coef;

  ans = x + *p++;
  i = N - 1;
  do
    ans = ans * x + *p++;
  while (--i);
  return ans;
}

/*==================================================================
 * Bessel J0, Y0 (Cephes j0.c)
 *==================================================================*/
static const double J0_PP[7] = {
  7.96936729297347051624E-4, 8.28352392107440799803E-2,
  1.23953371646414299388E0,  5.44725003058768775090E0,
  8.74716500199817011941E0,  5.30324038235394892183E0,
  9.99999999999999997821E-1,
};
static const double J0_PQ[7] = {
  9.24408810558863637013E-4, 8.56288474354474431428E-2,
  1.25352743901058953537E0,  5.47097740330417105182E0,
  8.76190883237069594232E0,  5.30605288235394617618E0,
  1.00000000000000000218E0,
};
static const double J0_QP[8] = {
  -1.13663838898469149931E-2, -1.28252718670509318512E0,
  -1.95539544257735972385E1,  -9.32060152123768231369E1,
  -1.77681167980488050595E2,  -1.47077505154951170175E2,
  -5.14105326766599330220E1,  -6.05014350600728481186E0,
};
static const double J0_QQ[7] = {
  /* 1.0 */
  6.43178256118178023184E1, 8.56430025976980587198E2,
  3.88240183605401609683E3, 7.24046774195652478189E3,
  5.93072701187316984827E3, 2.06209331660327847417E3,
  2.42005740240291393179E2,
};
static const double J0_YP[8] = {
  1.55924367855235737965E4,  -1.46639295903971606143E7,
  5.43526477051876500413E9,  -9.82136065717911466409E11,
  8.75906394395366999549E13, -3.46628303384729719441E15,
  4.42733268572569800351E16, -1.84950800436986690637E16,
};
static const double J0_YQ[7] = {
  /* 1.0 */
  1.04128353664259848412E3,  6.26107330137134956842E5,
  2.68919633393814121987E8,  8.64002487103935000337E10,
  2.02979612750105546709E13, 3.17157752842975028269E15,
  2.50596256172653059228E17,
};
static const double J0_DR1 = 5.78318596294678452118E0;
static const double J0_DR2 = 3.04712623436620863991E1;
static const double J0_RP[4] = {
  -4.79443220978201773821E9,  1.95617491946556577543E12,
  -2.49248344360967716204E14, 9.70862251047306323952E15,
};
static const double J0_RQ[8] = {
  /* 1.0 */
  4.99563147152651017219E2,  1.73785401676374683123E5,
  4.84409658339962045305E7,  1.11855537045356834862E10,
  2.11277520115489217587E12, 3.10518229857422583814E14,
  3.18121955943204943306E16, 1.71086294081043136091E18,
};

static double cephes_j0(double x)
{
  double w, z, p, q, xn;

  if (x < 0)
    x = -x;

  if (x <= 5.0) {
    z = x * x;
    if (x < 1.0e-5)
      return 1.0 - z / 4.0;
    p = (z - J0_DR1) * (z - J0_DR2);
    p = p * polevl(z, J0_RP, 3) / p1evl(z, J0_RQ, 8);
    return p;
  }

  w = 5.0 / x;
  q = 25.0 / (x * x);
  p = polevl(q, J0_PP, 6) / polevl(q, J0_PQ, 6);
  q = polevl(q, J0_QP, 7) / p1evl(q, J0_QQ, 7);
  xn = x - M_PI_4;
  p = p * cos(xn) - w * q * sin(xn);
  return p * SQ2OPI / sqrt(x);
}

static double cephes_y0(double x)
{
  double w, z, p, q, xn;

  if (x <= 5.0) {
    if (x == 0.0)
      return -HUGE_VAL;
    if (x < 0.0)
      return NAN;
    z = x * x;
    w = polevl(z, J0_YP, 7) / p1evl(z, J0_YQ, 7);
    w += M_2_PI * log(x) * cephes_j0(x);
    return w;
  }

  w = 5.0 / x;
  z = 25.0 / (x * x);
  p = polevl(z, J0_PP, 6) / polevl(z, J0_PQ, 6);
  q = polevl(z, J0_QP, 7) / p1evl(z, J0_QQ, 7);
  xn = x - M_PI_4;
  p = p * sin(xn) + w * q * cos(xn);
  return p * SQ2OPI / sqrt(x);
}

/*==================================================================
 * Bessel J1, Y1 (Cephes j1.c)
 *==================================================================*/
static const double J1_RP[4] = {
  -8.99971225705559398224E8, 4.52228297998194034323E11,
  -7.27494245221818276015E13, 3.68295732863852883286E15,
};
static const double J1_RQ[8] = {
  /* 1.0 */
  6.20836478118054335476E2,  2.56987256757748830383E5,
  8.35146791431949253037E7,  2.21511595479792499675E10,
  4.74914122079991414898E12, 7.84369607876235854894E14,
  8.95222336184627338078E16, 5.32278620332680085395E18,
};
static const double J1_PP[7] = {
  7.62125616208173112003E-4, 7.31397056940917570436E-2,
  1.12719608129684925192E0,  5.11207951146807644818E0,
  8.42404590141772420927E0,  5.21451598682361504063E0,
  1.00000000000000000254E0,
};
static const double J1_PQ[7] = {
  5.71323128072548699714E-4, 6.88455908754495404082E-2,
  1.10514232634061696926E0,  5.07386386128601488557E0,
  8.39985554327604159757E0,  5.20982848682361821619E0,
  9.99999999999999997461E-1,
};
static const double J1_QP[8] = {
  5.10862594750176621635E-2, 4.98213872951233449420E0,
  7.58238284132545283818E1,  3.66779609360150777800E2,
  7.10856304998926107277E2,  5.97489612400613639965E2,
  2.11688757100572135698E2,  2.52070205858023719784E1,
};
static const double J1_QQ[7] = {
  /* 1.0 */
  7.42373277035675149943E1, 1.05644886038262816351E3,
  4.98641058337653607651E3, 9.56231892404756170795E3,
  7.99704160447350683650E3, 2.82619278517639096600E3,
  3.36093607810698293419E2,
};
static const double J1_YP[6] = {
  1.26320474790178026440E9,  -6.47355876379160291031E11,
  1.14509511541823727583E14, -8.12770255501325109621E15,
  2.02439475713594898196E17, -7.78877196265950026825E17,
};
static const double J1_YQ[8] = {
  /* 1.0 */
  5.94301592346128195359E2,  2.35564092943068577943E5,
  7.34811944459721705660E7,  1.87601316108706159478E10,
  3.88231277496238566008E12, 6.20557727146953693363E14,
  6.87141087355300489866E16, 3.97270608116560655612E18,
};
static const double J1_Z1 = 1.46819706421238932572E1;
static const double J1_Z2 = 4.92184563216946036703E1;

static double cephes_j1(double x)
{
  double w, z, p, q, xn;

  w = x;
  if (x < 0)
    return -cephes_j1(-x);

  if (w <= 5.0) {
    z = x * x;
    w = polevl(z, J1_RP, 3) / p1evl(z, J1_RQ, 8);
    w = w * x * (z - J1_Z1) * (z - J1_Z2);
    return w;
  }

  w = 5.0 / x;
  z = w * w;
  p = polevl(z, J1_PP, 6) / polevl(z, J1_PQ, 6);
  q = polevl(z, J1_QP, 7) / p1evl(z, J1_QQ, 7);
  xn = x - THPIO4;
  p = p * cos(xn) - w * q * sin(xn);
  return p * SQ2OPI / sqrt(x);
}

static double cephes_y1(double x)
{
  double w, z, p, q, xn;

  if (x <= 5.0) {
    if (x == 0.0)
      return -HUGE_VAL;
    if (x <= 0.0)
      return NAN;
    z = x * x;
    w = x * (polevl(z, J1_YP, 5) / p1evl(z, J1_YQ, 8));
    w += M_2_PI * (cephes_j1(x) * log(x) - 1.0 / x);
    return w;
  }

  w = 5.0 / x;
  z = w * w;
  p = polevl(z, J1_PP, 6) / polevl(z, J1_PQ, 6);
  q = polevl(z, J1_QP, 7) / p1evl(z, J1_QQ, 7);
  xn = x - THPIO4;
  p = p * sin(xn) + w * q * cos(xn);
  return p * SQ2OPI / sqrt(x);
}

/*==================================================================
 * Bessel J_n, integer order.
 * n < |x| : upward recurrence from J0, J1 (stable there).
 * n >= |x|: Miller downward recurrence normalized with
 *           1 = J_0 + 2*sum_{k>=1} J_{2k}.
 *==================================================================*/
static double bessel_Jn(int n, double x)
{
  int k, m, sign = 1;
  double ax, r;

  if (n < 0) {
    n = -n;                             /* J_{-n} = (-1)^n J_n */
    if (n & 1) sign = -sign;
  }
  if (x < 0.0 && (n & 1)) sign = -sign; /* J_n(-x) = (-1)^n J_n(x) */
  ax = fabs(x);

  if (n == 0) return sign * cephes_j0(ax);
  if (n == 1) return sign * cephes_j1(ax);
  if (ax == 0.0) return 0.0;

  if ((double)n < ax) {
    /* upward recurrence */
    double bjm = cephes_j0(ax), bj = cephes_j1(ax), bjp;
    for (k = 1; k < n; k++) {
      bjp = (2.0 * k / ax) * bj - bjm;
      bjm = bj;
      bj  = bjp;
    }
    r = bj;
  }
  else {
    /* Miller downward recurrence */
    double bjp = 0.0, bj = 1e-300, bjm, sum = 0.0, val = 0.0;
    m = n + (int)ceil(ax) + 40;
    for (k = m; k >= 1; k--) {
      bjm = (2.0 * k / ax) * bj - bjp;  /* J_{k-1} */
      bjp = bj;
      bj  = bjm;
      if (fabs(bj) > 1e250) {           /* rescale to avoid overflow */
        bj  *= 1e-250;
        bjp *= 1e-250;
        sum *= 1e-250;
        val *= 1e-250;
      }
      if (((k - 1) & 1) == 0 && k - 1 > 0) sum += bj;  /* even orders >= 2 */
      if (k - 1 == n) val = bj;
    }
    sum = 2.0 * sum + bj;               /* + J_0 */
    r = val / sum;
  }
  return sign * r;
}

/*==================================================================
 * Bessel Y_n, integer order: upward recurrence from Y0, Y1
 * (always stable, since Y_n grows with n).  Domain: x > 0.
 *==================================================================*/
static double bessel_Yn(int n, double x)
{
  int k, sign = 1;
  double bym, by, byp;

  if (n < 0) {
    n = -n;                             /* Y_{-n} = (-1)^n Y_n */
    if (n & 1) sign = -sign;
  }
  if (x == 0.0) return -HUGE_VAL;
  if (x < 0.0) return NAN;

  if (n == 0) return sign * cephes_y0(x);
  if (n == 1) return sign * cephes_y1(x);

  bym = cephes_y0(x);
  by  = cephes_y1(x);
  for (k = 1; k < n; k++) {
    byp = (2.0 * k / x) * by - bym;
    bym = by;
    by  = byp;
  }
  return sign * by;
}

/*==================================================================
 * erf, erfc (Cephes ndtr.c)
 *==================================================================*/
static const double E_P[9] = {
  2.46196981473530512524E-10, 5.64189564831068821977E-1,
  7.46321056442269912687E0,   4.86371970985681366614E1,
  1.96520832956077098242E2,   5.26445194995477358631E2,
  9.34528527171957607540E2,   1.02755188689515710272E3,
  5.57535335369399327526E2,
};
static const double E_Q[8] = {
  /* 1.0 */
  1.32281951154744992508E1, 8.67072140885989742329E1,
  3.54937778887819891062E2, 9.75708501743205489753E2,
  1.82390916687909736289E3, 2.24633760818710981792E3,
  1.65666309194161350182E3, 5.57535340817727675546E2,
};
static const double E_R[6] = {
  5.64189583547755073984E-1, 1.27536670759978104416E0,
  5.01905042251180477414E0,  6.16021097993053585195E0,
  7.40974269950448939160E0,  2.97886665372100240670E0,
};
static const double E_S[6] = {
  /* 1.0 */
  2.26052863220117276590E0, 9.39603524938001434673E0,
  1.20489539808096656605E1, 1.70814450747565897222E1,
  9.60896809063285878198E0, 3.36907645100081516050E0,
};
static const double E_T[5] = {
  9.60497373987051638749E0, 9.00260197203842689217E1,
  2.23200534594684319226E3, 7.00332514112805075473E3,
  5.55923013010394962768E4,
};
static const double E_U[5] = {
  /* 1.0 */
  3.35617141647503099647E1, 5.21357949780152679795E2,
  4.59432382970980127987E3, 2.26290000613890934246E4,
  4.92673942608635921086E4,
};

static double cephes_erf(double x);

static double cephes_erfc(double a)
{
  double p, q, x, y, z;

  x = (a < 0.0) ? -a : a;

  if (x < 1.0)
    return 1.0 - cephes_erf(a);

  z = -a * a;
  if (z < -MAXLOG)                      /* underflow */
    return (a < 0) ? 2.0 : 0.0;
  z = exp(z);

  if (x < 8.0) {
    p = polevl(x, E_P, 8);
    q = p1evl(x, E_Q, 8);
  }
  else {
    p = polevl(x, E_R, 5);
    q = p1evl(x, E_S, 6);
  }
  y = (z * p) / q;

  if (a < 0)
    y = 2.0 - y;
  if (y == 0.0)
    return (a < 0) ? 2.0 : 0.0;
  return y;
}

static double cephes_erf(double x)
{
  double z;

  if (x < 0.0)
    return -cephes_erf(-x);
  if (x > 1.0)
    return 1.0 - cephes_erfc(x);

  z = x * x;
  return x * polevl(z, E_T, 4) / p1evl(z, E_U, 5);
}

/*==================================================================
 * Modified Bessel function of the first kind, integer order n.
 * Miller's algorithm: downward recurrence I_{k-1} = I_{k+1} + (2k/x) I_k
 * from a trial value high above max(n,|x|), normalized by
 *   exp(|x|) = I_0 + 2 * sum_{k>=1} I_k.
 *==================================================================*/
static double bessel_In(int n, double x)
{
  int k, m;
  double ax, bi, bip, bim, sum, val, r;

  if (n < 0) n = -n;                    /* I_{-n} = I_n */
  ax = fabs(x);
  if (ax == 0.0) return (n == 0) ? 1.0 : 0.0;

  m   = n + (int)ceil(ax) + 40;         /* start order: safely above n, |x| */
  bip = 0.0;                            /* trial I_{m+1} */
  bi  = 1e-300;                         /* trial I_m (arbitrary small)      */
  sum = 0.0;
  val = 0.0;

  for (k = m; k >= 1; k--) {
    bim = bip + (2.0 * k / ax) * bi;    /* I_{k-1} */
    bip = bi;
    bi  = bim;
    if (fabs(bi) > 1e250) {             /* rescale to avoid overflow */
      bi  *= 1e-250;
      bip *= 1e-250;
      sum *= 1e-250;
      val *= 1e-250;
    }
    sum += 2.0 * bip;                   /* accumulates 2*I_k, k = m..1 */
    if (k - 1 == n) val = bi;
  }
  sum += bi;                            /* + I_0 */

  r = (val / sum) * exp(ax);
  if (x < 0.0 && (n & 1)) r = -r;       /* I_n(-x) = (-1)^n I_n(x) */
  return r;
}

/*==================================================================
 * Complete elliptic integrals K(k), E(k), modulus-k convention:
 *   K(k) = int_0^{pi/2} dt / sqrt(1 - k^2 sin^2 t)
 * AGM iteration (Abramowitz & Stegun 17.6); quadratic convergence.
 *==================================================================*/
static void ellip_KE(double k, double *K, double *E)
{
  double k2 = k * k;
  double a, b, c, t, fac, sum;

  if (k2 >= 1.0) {                      /* K diverges at |k| = 1 */
    *K = HUGE_VAL;
    *E = (k2 == 1.0) ? 1.0 : NAN;
    return;
  }

  a   = 1.0;
  b   = sqrt(1.0 - k2);
  sum = 0.5 * k2;                       /* 2^{-1} c_0^2, c_0 = k */
  fac = 1.0;

  while (fabs(a - b) > DBL_EPSILON * a) {
    c = 0.5 * (a - b);
    t = 0.5 * (a + b);
    b = sqrt(a * b);
    a = t;
    fac *= 2.0;
    sum += 0.5 * fac * c * c;           /* 2^{n-1} c_n^2 */
  }

  *K = M_PI / (2.0 * a);
  *E = *K * (1.0 - sum);
}

/*------------------------------------------------------------------
 * Fortran-callable wrappers (trailing underscore, args by reference)
 *------------------------------------------------------------------*/
void m3d_bessel_i_(const int *n, const double *x, double *out)
{
  *out = bessel_In(*n, *x);
}

void m3d_bessel_j_(const int *n, const double *x, double *out)
{
  *out = bessel_Jn(*n, *x);
}

void m3d_bessel_y_(const int *n, const double *x, double *out)
{
  *out = bessel_Yn(*n, *x);
}

void m3d_elliptic_e_(const double *k, double *out)
{
  double K, E;
  ellip_KE(*k, &K, &E);
  *out = E;
}

void m3d_elliptic_k_(const double *k, double *out)
{
  double K, E;
  ellip_KE(*k, &K, &E);
  *out = K;
}

void m3d_erf_(const double *x, double *out)
{
  *out = cephes_erf(*x);
}
