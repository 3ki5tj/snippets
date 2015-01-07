/* K-means clustering */



"use strict";



// small amount to be added to the variance
// to ensure the stability of cholesky decomposition
var KMEANS_EPSILON = 1e-8;



function Kmeans(dim, K, dat)
{
  this.dim = dim;
  this.n = dat.length;
  this.K = K;
  this.x = dat;
  this.av = newarr2d(K, dim);
  this.var = newarr2d(K, dim * dim);
  this.lndetL = newarr(K);
  this.lnp = newarr2d(this.n, K);
  this.psum = newarr(k);
  for ( var k = 0; k < K; k++ ) {
    var i = Math.floor( this.n * k / K );
    if ( i >= this.n ) i = this.n - 1;
    for ( var d = 0; d < dim; d++ ) {
      this.av[k][d] = dat[i][d] + 1e-3 * (2.*Math.random() - 1);
    }
    for ( i = 0; i < dim; i++ )
      this.var[k][i*dim+i] = KMEANS_EPSILON;
    this.psum[k] = 1./K;
  }
}



/* normalize lna such that Sum { i = 1 to n } exp(lna) = 1 */
function lnnorm(lna, n)
{
  var max = -1e300, sum = 0, c;
  var im, i;

  for ( i = 0; i < n; i++ )
    if ( lna[i] > max )
      max = lna[im = i];
  for ( sum = 0, i = 0; i < n; i++ )
    sum += Math.exp(lna[i] - max);
  // we want max + log(sum) + c = log(1.0) = 0
  c = -max - Math.log(sum);
  for ( i = 0; i < n; i++ )
    lna[i] += c;
  return c;
}



/* given the means and variance of each cluster,
 * compute the probability of each frame belonging to the cluster */
Kmeans.prototype.estep = function()
{
  var k, i, d, K = this.K, n = this.n, dim = this.dim;
  var L, dx, dy, s, c;

  L = newarr(dim * dim);
  dx = newarr(dim);
  dy = newarr(dim);
  for ( k = 0; k < K; k++ ) {
    // perform the cholesky decomposition of the covariance matrix
    //    var = L.L^T
    for ( d = 0; d < dim * dim; d++ ) // copy the matrix
      L[d] = this.var[k][d];
    choldecomp(L, dim); // var = L . L^T
    this.lndetL[k] = chollogdetL(L, dim, KMEANS_EPSILON);
    c = -this.lndetL[k] + Math.log(this.psum[k]);

    for ( i = 0; i < n; i++ ) {
      for ( d = 0; d < dim; d++ )
        dy[d] = dx[d] = this.x[i][d] - this.av[k][d];
      // we need s = dx^T var^(-1) dx
      // now var = L.L^T, so that var^(-1) = L^(-1)^T . L^(-1)
      // and s = |L^(-1) dx|^2 = |dy|^2
      // here L^(-1) dx = dy, or L dy = dx
      cholsolveL(L, dy, dim);
      for ( s = 0, d = 0; d < dim; d++ )
        s += dy[d] * dy[d];
      this.lnp[i][k] = -0.5 * s + c;
    }
  }

  // normalize p[i][k] for each frame i, such that
  // Sum {k = 0 to K - 1} p[i][k] = 1
  for ( i = 0; i < n; i++ )
    lnnorm(this.lnp[i], K);
}



/* update the averages and covariance according to the current p */
Kmeans.prototype.mstep = function()
{
  var k, i, d, d2, K = this.K, n = this.n, dim = this.dim;
  var s, p, dx;

  dx = newarr(dim);
  for ( k = 0; k < K; k++ ) {
    // compute the mean of cluster k
    for ( d = 0; d < dim; d++ )
      this.av[k][d] = 0;
    s = 0;
    for ( i = 0; i < n; i++ ) {
      p = Math.exp(this.lnp[i][k]);
      for ( d = 0; d < dim; d++ )
        this.av[k][d] += p * this.x[i][d];
      s += p;
    }
    for ( d = 0; d < dim; d++ )
      this.av[k][d] /= s;
    this.psum[k] = s;

    /* compute the variance of cluster k */
    for ( d = 0; d < dim; d++ )
      for ( d2 = 0; d2 < dim; d2++ )
        this.var[k][d*dim + d2] = 0;

    for ( i = 0; i < n; i++ ) {
      p = Math.exp(this.lnp[i][k]) / s;
      for ( d = 0; d < dim; d++ )
        dx[d] = this.x[i][d] - this.av[k][d];

      for ( d = 0; d < dim; d++ )
        for ( d2 = 0; d2 < dim; d2++ )
          this.var[k][d*dim + d2] += p * dx[d] * dx[d2];
    }
  }
}



Kmeans.prototype.print = function()
{
  var k, d, d2, K = this.K, dim = this.dim;

  for ( k = 0; k < K; k++ ) {
    console.log("cluster " + k + ", pop " + this.psum[k]);
    console.log("ave:", this.av[k]);

    var s = "var:\n";
    for ( d = 0; d < dim; d++ ) {
      for ( d2 = 0; d2 < dim; d2++ )
        s += this.var[k][d*dim + d2] + "\t";
      s += "\n";
    }
    s += "\n";
    console.log(s);
  }
}



