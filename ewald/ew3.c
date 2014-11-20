/* Calculate the electrostatic energy of NaCl using Ewald sum */
#include <stdio.h>
#include <math.h>
#define N 2
#define M 5
int main(void){
  int i,j,k,zeros;
  double x,multi,k0,cwe,ewc=2.0; // ewc=sqrt(alpha)
  double phi_real,phi_self,phi_recip;
  
  // self energy,  charge and its own cloud
  phi_self=-2.0*ewc/sqrt(M_PI); 
  
  // real space sum
  phi_real=0;  
  for(i=0; i<=N; i++)
  for(j=0; j<=N; j++)
  for(k=0; k<=N; k++){
    zeros=(i==0)+(j==0)+(k==0);
    if(zeros==3)continue;
    else if(zeros==2) multi=2;  // 6 axes, but I'll cover 3, +x,+y,+z
    else if(zeros==1) multi=4;  // 12 faces, but I'll cover 3, +xy,+yz,+zx
    else multi=8; // 8 octants
    if((i+j+k)%2 != 0) multi=-multi; // sign of the charge
    
    x=sqrt(i*i+j*j+k*k);
    phi_real += multi*erfc(ewc*x)/x;
  }

  // we choose a 2x2x2 unit cell, so there are 4 Na+ and 4 Cl- inside
  k0=(2.0*M_PI)/2.0;  
  cwe=0.25*k0*k0/(ewc*ewc);
  // Here we calculate sum phi_i q_i over the 8 ions, 
  // the eight terms should be the same, 
  // we'll include the factor 8 at the end
  phi_recip=0;
  // only sum over odd wave vectors due to cancellation between Na+ and Cl-
  for(i=1; i<=M; i+=2)
  for(j=1; j<=M; j+=2)
  for(k=1; k<=M; k+=2){
    x=i*i+j*j+k*k;
    phi_recip += exp(-x*cwe)/x;
  }
  // the normalization factor
  // 1. |rho_k|^2 is always 8.0 because 4 Na+ & 4 Cl- are supposed to 
  //    reinforce each other, so |rho_k|^2 = 64,
  // 2. we are supposed to divide it by k^2 = k0^2 * index^2
  //    = pi^2 index^2 instead of by index^2
  // 3. here we assume 4 pi eps0 = 1.0
  //    so 1.0/eps0 = 4 pi
  // 4. each cell has 4 Na+ & 4 Cl- ions, or 8 equivalent ions,
  //    multiply it by 1/8.
  // 5. eight octants in the reciprocal k space, 
  //    +/-kx, +/-ky, +/-kz, so multiply it by 8
  // 6. divided by the volume 8.0
  // so the normalization is 64*4*pi/(pi^2)*(1/8)*8/8.0=32.0/pi;
  phi_recip*=32.0/M_PI;

  printf("error in realspace=%g, in k-space=%g\n", 
      erfc(ewc*(N+1))/(N+1), exp(-.25*(M+2)*(M+2)*k0*k0/ewc/ewc)/((M+2)*(M+2)*k0*k0));
  printf("self   =%+17.14lf\n"
         "real   =%+17.14lf\n"
         "recip  =%+17.14lf\n"
         "phi    =%+17.14lf\n", 
         phi_self,phi_real,phi_recip,phi_recip+phi_real+phi_self);
  return 0;
}
