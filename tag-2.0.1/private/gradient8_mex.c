// TopoToolbox
// file: gradient8_mex.c
// steepest gradient
// input:   DEM [m], cellsize [m]
// output:  G [m/m]

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize nrow, ncol; // use mwIndex
	mwIndex ix, ixn, row, col, r, c, ixnn;
	int coln[4] = {1, 1, 1, 0};
	int rown[4] = {-1, 0, 1, 1};
	double z, zn, g;
	double *cs;
	double *dem;
	double *G;
	double dist[4] = {sqrt(2), 1, sqrt(2), 1 };

    
    if( nrhs > 2 || nlhs > 1 ) {
        mexErrMsgTxt("Need 2 inputs and no more than 1 output");
    }
    
    dem    = mxGetPr(prhs[0]);
	nrow   = mxGetM(prhs[0]);
	ncol   = mxGetN(prhs[0]);
	cs     = mxGetPr(prhs[1]);
	
	/* Create an m-by-n mxArray and either copy initial values into it or create a new one with ones*/
	plhs[0] = mxCreateDoubleMatrix((mwSize) nrow,(mwSize) ncol, mxREAL);
	/* pointer into drainage basins*/
	G       = mxGetPr(plhs[0]);
	
	// set cell distances
	for( ixn = 0; ixn < 4; ixn++) {
		dist[ixn] = cs[0] * dist[ixn];
	}
    
    for( row= 0; row<nrow; row++ ) {
	
		for( col= 0; col<ncol; col++) {
			
			ix = (col*nrow) + row;
			z = dem[ix];
			
			// nan?
			if (mxIsNaN(z)) {
				G[ix] = dem[ix];
				continue;
			}
		
			for ( ixn = 0; ixn < 4; ixn++ ) {
				// row index of neighbor
				r = row + rown[ixn];
				// column index of neighbor
				c = col + coln[ixn];
				
				// check if indices may be outside the grid 
				if (r < 0) { 
					continue;
				}
				if (r >= nrow) { 
					continue;
				}
				if (c < 0) { 
					continue;
				}
				if (c >= ncol) { 
					continue;
				}
				ixnn = (c*nrow) + r;
				// calculate slope
				g = (z - dem[ixnn])/ dist[ixn];
				// find maximum slope
				if ( g > G[ix] ) {
				   G[ix] = g;
				   } 
				else if (g < 0) {
					g = -g;
					if ( g > G[ixnn] ) {
						G[ixnn] = g;
						}
					}
			}
		}
		
    }
}