#ifndef LJ3D_H__
#define LJ3D_H__



/* Three-dimensional Lennard-Jones fluid */



#ifndef D
#define D 3
#endif



/* initialize a fcc lattice */
static void lj_initfcc(lj_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  double a, noise;

  n1 = (int) (pow(2*n, 1.0/D) + .999999); /* # of particles per side */
  a = lj->l / n1;
  noise = a * 1e-5;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++)
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 != 0) continue;
        /* add some noise to prevent two atoms happened to
         * be separated by precisely some special cutoff distance,
         * which might be half of the box */
        lj->x[id][0] = (i + .5) * a + noise * (2*rand01() - 1);
        lj->x[id][1] = (j + .5) * a + noise * (2*rand01() - 1);
        lj->x[id][2] = (k + .5) * a + noise * (2*rand01() - 1);
        id++;
      }
}



/* get the tail correction */
static double lj_gettail(lj_t *lj, double rho, int n, double *ptail)
{
  double irc, irc3, irc6, utail;

  irc = 1/lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  utail = 8*M_PI*rho*n/9*(irc6 - 3)*irc3;
  if (ptail != NULL)
    *ptail = 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3;
  return utail;
}



#endif /* LJ_H__ */
