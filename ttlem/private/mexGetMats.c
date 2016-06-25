// mexGetMats.c
// [Ar Ac Av Af B] = mexGetMats(z,bdy,dt,dx,dy,Nx,Ny,a,b,U)

#include "mex.h"
#include "matrix.h"


void GetfF(double z[], double bdy[], double *dx, double *dy, const int K, const int J, double *aconst, double *bconst, double *U, double f[], double F[])
{
   
int i, j, KJ, c, dc[8];
double a, b, ab;
double invdx, invdx2, invdy, invdy2;
double zijp1, zim1jp1, zim1j, zim1jm1, zijm1, zip1jm1, zip1j, zip1jp1, zij;
double zx, zy, zxx, zyy, zxy, Lap, invden, invden2, invden3;
double zxzx, zyzy, zxzxx, zyzyy, zxzxy, zyzxy, q1, q1x, q1y, q2, q3, q3x, q3y, q4, q4x, q4y, q5, q5x, q5y;

// create arrays of RHS and partial derivatives.
// Columns run CCW from E:
// 0 (i,j+1) 
// 1 (i-1,j+1)
// 2 (i-1,j)
// 3 (i-1,j-1)
// 4 (i,j-1)
// 5 (i+1,j-1)
// 6 (i+1,j)
// 7 (i+1,j+1)
// 8 (i,j)

KJ = K*J;
a = *aconst;
b = *bconst;
ab = a*b;
invdx = 1/(*dx);
invdx2 = invdx*invdx;
invdy = 1/(*dy);
invdy2 = invdy*invdy;

dc[0]=K; dc[1]=K-1; dc[2]=-1; dc[3]=-K-1; dc[4]=-K; dc[5]=-K+1; dc[6]=1; dc[7]=K+1;
                          
for (i=0; i<K; i++) {
   for (j=0; j<J; j++) {
       if (i>0 && i<K-1 && j>0 && j<J-1) { // this version doesn't handle periodic or non-pinned boundaries
           c=K*j+i; // vector index of element (i,j)
           if (bdy[c]!=1) { // leave f and F zero if it's a pinned element
                
                // elevations in 9-point neighborhood
                // 0 [i,j+1]
                zijp1=z[c+dc[0]];
                // 1 [i-1,j+1]
                zim1jp1=z[c+dc[1]];
                // 2 [i-1,j]
                zim1j=z[c+dc[2]];
                // 3 [i-1,j-1]
                zim1jm1=z[c+dc[3]];
                // 4 [i,j-1]
                zijm1=z[c+dc[4]];
                // 5 [i+1,j-1]
                zip1jm1=z[c+dc[5]];
                // 6 [i+1,j]
                zip1j=z[c+dc[6]];
                // 7 [i+1,j+1]
                zip1jp1=z[c+dc[7]];
                // 8 [i,j]
                zij=z[c];
                
                
                // derivatives
                zx=0.5*(zijp1-zijm1)*invdx;
                zy=0.5*(zip1j-zim1j)*invdy;
                zxx=(zijp1-2*zij+zijm1)*invdx*invdx;
                zyy=(zip1j-2*zij+zim1j)*invdy*invdy;
                zxy=0.25*(zip1jp1-zip1jm1-zim1jp1+zim1jm1)*invdx*invdy;
                Lap=zxx+zyy;
                invden = 1/(1-b*(zx*zx+zy*zy));
                invden2 = invden*invden;
                invden3 = invden2*invden;
                
                // calculate f
                f[c] = *U + a*Lap*invden + 2*ab*invden2*(zx*zx*zxx + zy*zy*zyy + 2*zx*zy*zxy);

                
                // calculate partials
                
                // calculate repeated quantities up front
                zxzx = zx*zx;
                zyzy = zy*zy;
                zxzxx = zx*zxx;
                zyzyy = zy*zyy;
                zxzxy = zx*zxy;
                zyzxy = zy*zxy;

                q1 = 4*ab*b*invden3*(zxzx*zxx+zyzy*zyy+2*zx*zyzxy);
                q1x = q1*zx*invdx;
                q1y = q1*zy*invdy;
                
                q2 = ab*invden2*zx*invdx*zy*invdy;
                
                q3 = ab*Lap*invden2;
                q3x = q3*zx*invdx;
                q3y = q3*zy*invdy;

                q4 = a*invden;
                q4x = q4*invdx2;
                q4y = q4*invdy2;
                
                q5 = 2*ab*invden2;
                q5x = q5*invdx;
                q5y = q5*invdy;
                
                
                // 0 (i,j+1) 
                F[c] = q4x + q3x + q1x + q5x*(zxzxx + zxzx*invdx + zyzxy);

                // 1 (i-1,j+1)
                F[KJ+c] = -q2;
                
                // 2 (i-1,j)
                F[2*KJ+c] = q4y - q3y - q1y + q5y*(-zyzyy + zyzy*invdy - zxzxy);
                
                // 3 (i-1,j-1)
                F[3*KJ+c] = q2;
                
                // 4 (i,j-1)
                F[4*KJ+c] = q4x - q3x - q1x + q5x*(-zxzxx + zxzx*invdx - zyzxy);
                
                // 5 (i+1,j-1)
                F[5*KJ+c] = -q2;
                
                // 6 (i+1,j)
                F[6*KJ+c] = q4y + q3y + q1y + q5y*(zyzyy + zyzy*invdy + zxzxy);
                
                // 7 (i+1,j+1)
                F[7*KJ+c] = q2;
                
                // 8 (i,j)
                F[8*KJ+c] = -2*(q4x + q4y) - 2*q5*(zxzx*invdx2 + zyzy*invdy2);
                
           }
       }
   }
}
    
} // end GetfF


int BuildMats(double z[], double *deltat, const int K, const int J, double f[], double F[], double Ar[], double Ac[], double Av[], double B[])
{

int c, a=0, k, n, dc[8], KJ;
double dt;

dt = (*deltat);
dc[0]=K; dc[1]=K-1; dc[2]=-1; dc[3]=-K-1; dc[4]=-K; dc[5]=-K+1; dc[6]=1; dc[7]=K+1;
KJ = K*J;


// B is RHS vector (present time). The full expression is 
// B = (1 - dt*Fi)*zi^n + dt*f - dt*sum(Fk*zk^n)

// A is matrix that operates on future time. This will give us the quantity
// (1 - dt*Fi)*zi^n+1 - dt*sum(Fk*zk^n+1)


for (c=0; c<KJ; c++) { // for each element


    // diagonal elements of A
    Ar[a] = c+1;
    Ac[a] = c+1;
    Av[a] = 1 - dt*F[KJ*8+c];
    a++; // move to the next element in the A vectors

    // first 2 terms of B
    B[c] = (1 - dt*F[KJ*8+c])*z[c] + dt*f[c];

    
    for (k=0; k<8; k++) { // loop through neighbors
        
        n = c + dc[k]; // vector index (and matrix column) of neighbor element

        if (n>-1 && n<KJ) { // don't go off boundaries
            
            
            // off-diagonal elements of A
            Ar[a] = c+1;
            Ac[a] = n+1;
            Av[a] = -dt*F[KJ*k+c];
            a++; // move to the next element in the A vectors

            // add neighbor terms of B
            B[c] = B[c] - dt*F[KJ*k+c]*z[n];
            
        }
    }
}

return a-1;
    
} // end BuildMats


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *fptr, *Fptr;
  double *z,*bdy,*dt,*dx,*dy,*Nx,*Ny,*a,*b,*U,*f,*F,*Ar,*Ac,*Av,*Af,*B;
  int K, J, nz;

  // Get pointers to inputs
  // syntax is [Ar Ac Av Af B] = mexGetMats(z,bdy,dt,dx,dy,Nx,Ny,a,b,U);
  z = (double *)mxGetPr(prhs[0]); // elevations
  bdy = (double *)mxGetPr(prhs[1]); // 1 if no uplift, 2 if periodic, zero if free to uplift and erode
  dt = (double *)mxGetPr(prhs[2]);   
  dx = (double *)mxGetPr(prhs[3]); 
  dy = (double *)mxGetPr(prhs[4]); 
  Nx = (double *)mxGetPr(prhs[5]); 
  Ny = (double *)mxGetPr(prhs[6]); 
  a = (double *)mxGetPr(prhs[7]); // K/(rho_r/rho_s)
  b = (double *)mxGetPr(prhs[8]); // 1/Sc^2
  U = (double *)mxGetPr(prhs[9]);

  // convert dimensions to integers
  K = (int)(*Ny+0.5);
  J = (int)(*Nx+0.5);
  
  // Create arrays for return arguments
  Ar = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K*J*9, 1, mxREAL)); // 1-BASED (Matlab) row indices of the LHS matrix operator. 
  Ac = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K*J*9, 1, mxREAL)); // 1-BASED (Matlab) column indices of the LHS matrix operator. 
  Av = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(K*J*9, 1, mxREAL)); // values of the LHS matrix operator. 
  Af = (double *)mxGetPr(plhs[3]= mxCreateDoubleMatrix(1, 1, mxREAL)); // index of last nonzero element in the A vectors 
  B = (double *)mxGetPr(plhs[4]= mxCreateDoubleMatrix(K*J, 1, mxREAL)); // The RHS vector
  
  // Create internally used arrays
  f = (double *)mxGetPr(fptr= mxCreateDoubleMatrix(K*J, 1, mxREAL)); // RHS of time derivative at present time
  F = (double *)mxGetPr(Fptr= mxCreateDoubleMatrix(K*J, 9, mxREAL)); // partial derivatives of RHS of time derivative at present time

  GetfF(z,bdy,dx,dy,K,J,a,b,U,f,F); // get values of finite difference approximation and its partial derivatives
  
  nz=BuildMats(z,dt,K,J,f,F,Ar,Ac,Av,B); // calculate values of the nonzero elements in the matrix equation
  
  *Af = nz+1; // now *Af is the 1-BASED (Matlab) index of the last nonzero element that was assigned to the A vectors
    
  // free memory
  mxDestroyArray(fptr);
  mxDestroyArray(Fptr);
  
} // End mexFunction
