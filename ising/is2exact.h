/* stand-alone code for the free energy of the Ising model */
#include "lnum.h"

/* log(exp(a) + exp(b)) */
__inline static double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = b - a) < ln0) ? a : a + log(1 + exp(c));
}

/* log(exp(a) - exp(b)), only works for a > b */
__inline static double lndif(double a, double b)
{
  double c;
  return ((c = b - a) < ln0) ? a : a + log(1 - exp(c));
}

/* log(exp(a)+b) */
__inline static double lnaddn(double a, double b)
{
  return (a > -ln0) ? a : a + log(1 + b*exp(-a));
}

/* the exact logarithm of the partition function, Z
 * the average energy and heat capacity of the Ising model */
__inline static double is2exact(int n, int m, double beta, double *eav, double *cv)
{
  double mh, nm, ex, f, th, sech, bet2, ln2, x;
  double lnz, lnz1, lnz2, lnz3, lnz4, dz, ddz;
  double z21, z31, z41, za;
  double dr1, dr2, dr3, dr4, ddr1, ddr2, ddr3, ddr4;
  double g, g0, dg, ddg, dg0;
  double xn2b, sh2b, coth2b;
  double lnch2b, lncc2b, lncl, lnsl, cd, cdsqr, lnddcl;
  int k, sgn4 = 1;

  mh = .5*m;
  nm = m * n;
  ln2 = log(2.0);
  bet2 = 2.*beta;
  xn2b = exp(-bet2);

  lnz1 = lnz2 = lnz3 = lnz4 = 0;
  dr1 = dr2 = dr3 = dr4 = 0;
  ddr1 = ddr2 = ddr3 = ddr4 = 0;
  lnch2b = lnadd(bet2, -bet2) - ln2;
  coth2b = 2./(1. - xn2b*xn2b) - 1.;
  lncc2b = lnch2b + log(coth2b); /* ln[ cosh(2b) * coth(2b) ] */
  g0 = bet2 + log(2./(1. + xn2b) - 1.);
  sgn4 = (g0 >= 0) ? 1 : -1;

  sh2b = 0.5*(1./xn2b - xn2b);
  dg0 = 2. + 2./sh2b;
  x = sh2b*sh2b;
  cd = 2. - 2./x; /* cl' = cd * cosh(2b) */
  cdsqr = cd*cd;
  lnddcl = lnaddn(lncc2b, 2.0/(x * sh2b)) + 2.*ln2; /* log(cl'') */

  for ( k = 0; k < n; k++ ) { /* for odd number */
    /* cosh(gamma_k) = cosh^2(2K)/sinh(2K) - cos(pi*k/n) */
    lncl = lnaddn(lncc2b, -cos((2*k + 1)*M_PI/n));
    lnsl = lncl + 0.5*log(1. - exp(-2.*lncl));
    g = lnadd(lncl, lnsl); /* gamma_{2k+1} */
    f = mh * g;
    lnz1 += lnadd(f, -f);
    lnz2 += lndif(f, -f);

    dg = exp(lnch2b - lnsl)*cd; /* g' = cl'/sl; */
    ex = exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    dr1 += mh * dg * th;
    dr2 += mh * dg / th;

    /* g''=cl''/sl - cl' ^2 *cl/sl^3; */
    ddg = exp(lnddcl - lnsl);
    ddg -= exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    sech = 2.0*dg/(ex + 1.0/ex); /* g' * sech(0.5*m*g) */
    ddr1 += mh*(ddg*th + mh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex); /* g' * csch(0.5*m*g) */
    ddr2 += mh*(ddg/th - mh*(sech*sech));

    if ( k == 0 ) {
      g = g0;
    } else {
      lncl = lnaddn(lncc2b, -cos(2.0*M_PI*k/n));
      lnsl = lncl + 0.5*log(1 - exp(-2*lncl));
      g = lnadd(lncl, lnsl);
    }
    f = mh * g;
    lnz3 += lnadd(f, -f); /* ln[2 cosh(f)] */
    lnz4 += (f < 0) ? lndif(-f, f) : lndif(f, -f); /* avoid neg. g0 */

    ex = exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    dg = (k == 0) ? dg0 : exp(lnch2b - lnsl)*cd;
    dr3 += mh * dg * th;
    dr4 += mh * dg / th;

    if ( k == 0 ) {
      ddg = -4*coth2b*coth2b*exp(-lnch2b);
    } else {
      ddg = exp(lnddcl - lnsl);
      ddg -= exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    }
    sech = 2.0*dg/(ex + 1.0/ex);
    ddr3 += mh*(ddg*th + mh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex);
    ddr4 += mh*(ddg/th - mh*(sech*sech));
  }

  z21 = exp(lnz2 - lnz1);
  z31 = exp(lnz3 - lnz1);
  z41 = sgn4*exp(lnz4 - lnz1);
  za = 1.0 + z21 + z31 + z41;
  lnz = lnz1 + log(za);
  lnz += .5*nm*log(2.*sh2b) - ln2;
  dz = (dr1 + z21*dr2 + z31*dr3 + z41*dr4)/za;
  if ( eav != NULL ) {
    *eav = - nm*coth2b - dz;
  }
  ddr1 += dr1 * dr1;
  ddr2 += dr2 * dr2;
  ddr3 += dr3 * dr3;
  ddr4 += dr4 * dr4;
  ddz = (ddr1 + z21*ddr2 + z31*ddr3 + z41*ddr4)/za;
  if ( cv != NULL ) {
    *cv = beta * beta * (-2.*nm/(sh2b*sh2b) + ddz - dz*dz);
  }
  return lnz;
}



/* return the exact logarithmic partition function
 * also compute the average energy and heat capacity */
__inline static double is2_exact(int n, int m, double beta, double *eav, double *cv)
{
  lnum_t a, b, c, d, ch, sh, ccs, dg, ddg, cg, sg, z[4];
  double mh, xp, den, th, gam, dgam, ddgam, tg, ln2, ish2;
  double dlnz[4], ddlnz[4], r[4], rt = 0, dz = 0, ddz = 0;
  int k, i;

  mh = m * 0.5;
  ln2 = log(2);
  /* sh = sinh(2*K) = exp(2*beta)(1-exp(-4*beta)/2); */
  xp = exp(-2*beta);
  lnum_set(&a, (1 - xp*xp)/2);
  lnum_imul(lnum_setln(&sh, 2*beta), &a);
  lnum_copy(&c, &sh);
  c.ln = -c.ln; /* c = 1/sinh(2*K); */
  lnum_add(&ccs, &sh, &c); /* ccs = sinh(2*K) + 1/sinh(2*K) */
  ish2 = exp(2*c.ln); /* 1/sinh^2(2*K) */
  /* dg = ccs' = 2*cosh(2*K)*(1 - 1/sinh^2(2*K)) */
  lnum_set(&d, (1 + xp*xp)*(1 - ish2));
  lnum_imul(lnum_setln(&dg, 2*beta), &d);
  /* ddg = ccs'' */
  lnum_setln(&ddg, sh.ln + ln2*2 + log(1 + ish2 + 2*ish2*ish2));
  for ( i = 0; i < 4; i++ ) {
    lnum_set(&z[i], 1);
    dlnz[i] = ddlnz[i] = 0;
  }
  for ( k = 0; k < n * 2; k++ ) {
    if ( k == 0 ) {
      gam = 2*beta + log((1 - xp)/(1 + xp));
      den = 1 - xp*xp;
      dgam = 2 + 4*xp/den;
      ddgam = -8*xp*(1+xp*xp)/(den*den);
    } else {
      lnum_add(&ch, &ccs, lnum_set(&d, -cos(M_PI*k/n))); /* ch = cosh(gam) */
      /* gam = arccosh(ch) = ln(ch + sh) + ln(ch) + ln(1+th) */
      th = sqrt(1 - exp(-2*ch.ln)); /* th = tanh(gam) */
      gam = ch.ln + log(1 + th);
      /* sinh(gam) = cosh(gam) * tanh(gam); */
      lnum_mul(&sh, &ch, lnum_set(&a, th));
      dgam = lnum_get( lnum_div(&a, &dg, &sh) ); /* dg/sinh(gam) */
      lnum_imul(lnum_set(&c, dgam*dgam), &ch);
      lnum_sub(&d, &ddg, &c); /* d = ddg - cosh(gam)*dgam^2 */
      ddgam = lnum_get( lnum_div(&b, &d, &sh) ); /* dg/sinh(gam) */
    }
    lnum_setln(&a,  m * gam * 0.5); /* a = exp(m*gam/2) */
    lnum_setln(&b, -m * gam * 0.5); /* b = exp(-m*gam/2) */
    lnum_add(&cg, &a, &b); /* cg = 2 cosh(m*gam/2) */
    lnum_sub(&sg, &a, &b); /* sg = 2 sinh(m*gam/2) */
    tg = lnum_get( lnum_div(&c, &sg, &cg) );
    i = (k % 2) * 2;
    lnum_imul(&z[i],   &cg);
    lnum_imul(&z[i+1], &sg);
    dlnz[i]   += dgam * tg * mh;
    dlnz[i+1] += dgam / tg * mh;
    ddlnz[i]   += ( dgam*dgam*4*exp(-2*cg.ln)*mh + ddgam*tg) * mh;
    ddlnz[i+1] += (-dgam*dgam*4*exp(-2*sg.ln)*mh + ddgam/tg) * mh;
  }

  lnum_add(&d, lnum_add(&b, lnum_add(&a, &z[0], &z[1]), &z[2]), &z[3]);
  lnum_setln(&a, 2*beta);
  lnum_sub(&sh, &a, lnum_setln(&b, -2*beta));
  for ( i = 0; i < 4; i++ ) {
    r[i] = lnum_get( lnum_div(&a, &z[i], &z[0]) );
    rt += r[i];
    dz += r[i] * dlnz[i];
    ddz += r[i] * (ddlnz[i] + dlnz[i] * dlnz[i]);
  }
  dz /= rt;
  *eav = n*m*(1 - 2/(1 - xp*xp)) - dz;
  *cv = beta * beta * (ddz/rt - dz * dz - 2*n*m*ish2);
  return d.ln + sh.ln*n*m*0.5 - ln2;
}

/* exact density of states */
__inline static double *is2dos(int n, int m)
{
  lpoly_t *a, *b, *c, *s, *bm, *u, *v, *w, *xp, *xm;
  lpoly_t **p, **q, *Z1, *Z2, *Z3, *Z4;
  lnum_t tmp;
  double *lndos, lnmfac, lnimax;
  int i, k;

  a = lpoly_open();
  b = lpoly_open();
  c = lpoly_open();
  s = lpoly_open();
  u = lpoly_open();
  v = lpoly_open();
  w = lpoly_open();
  bm = lpoly_open();
  xp = lpoly_open();
  xm = lpoly_open();
  Z1 = lpoly_open();
  Z2 = lpoly_open();
  Z3 = lpoly_open();
  Z4 = lpoly_open();
  p = calloc(m + 1, sizeof(*p));
  q = calloc(m + 1, sizeof(*q));
  for ( i = 0; i <= m; i++ ) {
    p[i] = lpoly_open();
    q[i] = lpoly_open();
  }
  lpoly_set(b, 4, 0.0, 2.0, 0.0, -2.0); /* beta = 2x(1-x^2) */
  lpoly_pow(u, b, m, v); /* beta^m */
  lnum_setln(&tmp, (m - 1)*log(0.5)); /* tmp = 0.5^(m-1) */
  lpoly_mullnum(bm, u, &tmp); /* bm = 0.5^(m-1) beta^m */

  for ( lnmfac = 0, i = 2; i <= m; i++ )
    lnmfac += log(i * 0.5); /* lnmfac = ln(m!/2^(m-1)) */

  lpoly_set(u, 2, 1.0, 1.0); /* u = 1 + x */
  lpoly_pow(xp, u, m, v); /* xp = (1 + x)^m */
  lpoly_set(u, 2, 1.0, -1.0); /* u = 1 - x */
  lpoly_pow(xm, u, m, v); /* xm = (1 - x)^m */
  lpoly_resize(u, m + 1);
  lnum_set(u->a + m, 1); /* u = x^m */
  lpoly_mul(v, u, xp); /* v = x^m (1 + x)^m */
  lpoly_mul(w, u, xm); /* w = x^m (1 - x)^m */
  lpoly_add(Z3, xm, v); /* Z3 = c0 = (1 - x)^m + x^m (1 + x)^m */
  lpoly_sub(Z4, xm, v); /* Z4 = s0 = (1 - x)^m - x^m (1 + x)^m */
  lpoly_add(Z1, xp, w); /* Z1 = cn = (1 + x)^m + x^m (1 - x)^m */
  lpoly_sub(Z2, xp, w); /* Z2 = sn = (1 + x)^m - x^m (1 - x)^m */

  if ( n % 2 == 0 ) { /* n is even */
    lpoly_imul(Z3, Z1, u);
    lpoly_imul(Z4, Z2, u);
    lpoly_set(Z1, 1, 1.0);
    lpoly_set(Z2, 1, 1.0);
  }

  for ( k = 1; k < n; k++ ) {
    lnum_set(&tmp, -cos(M_PI*k/n));
    lpoly_mullnum(u, b, &tmp); /* u = beta*cos(pi*k/n) */
    lpoly_set(v, 5, 1.0, 0.0, 2.0, 0.0, 1.0); /* v = (1+x^2)^2 */
    lpoly_add(a, v, u); /* ak = (1+x^2)^2 - beta*cos(pi*k/n) */

    /* p[i] = (ak^2 - beta^2)^i/(2i)! */
    lpoly_pow(u, a, 2, w);
    lpoly_pow(v, b, 2, w);
    lpoly_sub(w, u, v); /* w = ak^2 - beta^2 */
    lpoly_set(p[0], 1, 1.0);
    for ( i = 1; i <= m/2; i++ ) {
      lpoly_mul(p[i], p[i-1], w);
      lpoly_imulnum(p[i], 1.0/(2.*i*(2*i-1)));
    }
    /* q[i] = ak^i m!/i!/2^(m-1) */
    lpoly_set(q[0], 1, 1.0);
    lnum_setln(&q[0]->a[0], lnmfac);
    for ( i = 1; i <= m; i++ ) {
      lpoly_mul(q[i], q[i-1], a);
      lpoly_imulnum(q[i], 1./i);
    }

    /* ck^2 = bm + Sum_{i=0 to m/2} p[i]*q[m-2i] */
    lpoly_copy(c, bm);
    for ( i = 0; i <= m/2; i++ ) {
      lpoly_mul(u, p[i], q[m-2*i]);
      lpoly_iadd(c, u, v); /* c += u */
    }
    /* sk^2 = ck^2 - 2*(beta^m/2^(m-1)) */
    lpoly_sadd(s, c, bm, -2.0);
    if ( k % 2 ) {
      lpoly_imul(Z1, c, u);
      lpoly_imul(Z2, s, u);
    } else {
      lpoly_imul(Z3, c, u);
      lpoly_imul(Z4, s, u);
    }
  }
  lpoly_iadd(Z1, Z2, u);
  lpoly_iadd(Z1, Z3, u);
  lpoly_iadd(Z1, Z4, u);
  lpoly_imulnum(Z1, 0.5);

  /* export to a 1D array */
  lnimax = log((double) ULONG_MAX);
  lndos = calloc(n*m + 1, sizeof(*lndos));
  for ( i = 0; i < Z1->n; i+= 2 ) {
    double y = Z1->a[i].ln;
    if ( y < 0 ) y = ln0;
    else if ( y < lnimax ) /* round to nearest integer */
      y = log((double)((unsigned long) (exp(y)+.5)));
    lndos[i/2] = y;
  }

  lpoly_close(a);
  lpoly_close(b);
  lpoly_close(c);
  lpoly_close(s);
  lpoly_close(u);
  lpoly_close(v);
  lpoly_close(w);
  lpoly_close(bm);
  lpoly_close(xp);
  lpoly_close(xm);
  lpoly_close(Z1);
  lpoly_close(Z2);
  lpoly_close(Z3);
  lpoly_close(Z4);
  for ( i = 0; i <= m; i++ ) {
    lpoly_close(p[i]);
    lpoly_close(q[i]);
  }
  free(p);
  free(q);
  return lndos;
}

__inline static int is2dos_save(double *lndos, int n, int m)
{
  char fn[32];
  FILE *fp;
  int i, digs;

  sprintf(fn, "is2lndos%dx%d.dat", n, m);
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  digs = (int)(-log10(lndos[n*m/2]*DBL_EPSILON));
  if ( digs < 0 ) digs = 0;
  for ( i = 0; i <= n*m; i++ ) {
    fprintf(fp, "%d %.*f\n", -2*n*m+4*i, digs, lndos[i]);
  }
  fclose(fp);
  return 0;
}


