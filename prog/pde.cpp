/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include "crash.h"
#include "parameter.h"
#include "ca.h"
#include "pde.h"
#include "conrec.h"
#include "random.h"

/* STATIC DATA MEMBER INITIALISATION */
const int PDE::nx[9] = {0, 1, 1, 1, 0,-1,-1,-1, 0 };
const int PDE::ny[9] = {0, 1, 0,-1,-1,-1, 0, 1, 1 };

extern Parameter par;

/** PRIVATE **/

PDE::PDE(const int l, const int sx, const int sy) {
  
//  sigma=0;	//sigma not scalar
  thetime=0;
  sizex=sx;
  sizey=sy;
  layers=l;
  
  sigma=AllocateSigma(l,sx,sy);
  alt_sigma=AllocateSigma(l,sx,sy);
  sigmaExt=AllocateSigma(l,sx+200,sy+200);
  alt_sigmaExt=AllocateSigma(l,sx+200,sy+200);

}


PDE::PDE(void) {

//  sigma=0;
//  alt_sigma=0;
  sizex=0; sizey=0; layers=0;
  thetime=0;
  
}

// destructor (virtual)
PDE::~PDE(void) {
  if (sigma) {
    for (int i=0;i<layers;i++){
      for (int j=0;j<sizex;j++){
        delete [] sigma[i][j];
      };
      delete [] sigma[i];
    };
  }
  
  if (alt_sigma) {
    for (int i=0;i<layers;i++){
      for (int j=0;j<sizex;j++){
        delete [] alt_sigma[i][j];
      };
      delete [] alt_sigma[i];
    };
  }
  
}

double ***PDE::AllocateSigma(const int layers, const int sx, const int sy) {
  
  double ***mem;
  int lsizex=sx, lsizey=sy;
    
  mem=new double** [layers];
   
  for (int i=0;i<layers;i++){
    mem[i]=new double* [lsizex];
  
    for (int j=0;j<lsizex;j++) {
      mem[i][j]=new double [lsizey];
    };
  };
  
  
  /* Initialize PDE plane */
  
  // empty dish
  if (par.pde_initialization == 0) {
    for (int i=0;i<layers;i++){
      for (int j=0;j<lsizex;j++){
        for (int k=0;k<lsizey;k++) 
          mem[i][j][k]=0.0; 
      };
    };
  }
  

   return mem;
}

int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}


void PDE::Plot(Graphics *g,const int l) {
  // l=layer: default layer is 0
  for (int x=0;x<sizex;x++){
    for (int y=0;y<sizey;y++) {
      // Make the pixel four times as large
      // to fit with the CPM plane
      g->Point(MapColour(sigma[l][x][y]),2*x,2*y);
      g->Point(MapColour(sigma[l][x][y]),2*x+1,2*y);
      g->Point(MapColour(sigma[l][x][y]),2*x,2*y+1);
      g->Point(MapColour(sigma[l][x][y]),2*x+1,2*y+1);
    } 
  };
}

// Plot the value of the PDE only in the medium of the CPM
void PDE::Plot(Graphics *g, CellularPotts *cpm, const int l) {
  
  // suspend=true suspends calling of DrawScene
  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++) 
      if (cpm->Sigma(x,y)==0) {
	// Make the pixel four times as large
	// to fit with the CPM plane
	g->Point(MapColour(sigma[l][x][y]),2*x,2*y);
	g->Point(MapColour(sigma[l][x][y]),2*x+1,2*y);
	g->Point(MapColour(sigma[l][x][y]),2*x,2*y+1);
	g->Point(MapColour(sigma[l][x][y]),2*x+1,2*y+1);
      }
 
}

void PDE::PlotGradientMagnitude(Graphics *g, const int l, const int k){
// Plot the magnitude of the grandient in layer l
  double mag, magx, magy;
  int c;

  for (int y=2;y<sizey-2;y++) {
    for (int x=2;x<sizex-2;x++) {
      magx = (sigma[l][x+1][y]-sigma[l][x-1][y])/2.;
      magy = (sigma[l][x][y+1]-sigma[l][x][y-1])/2.;
      mag = sqrt( magx*magx + magy*magy );
      // Make the pixel four times as large
      // to fit with the CPM plane
      c = (((int)((mag/((mag)+k))*100))%100)+155;
      if (c>156){
        g->Point(c,2*x,2*y);
        g->Point(c,2*x+1,2*y);
        g->Point(c,2*x,2*y+1);
        g->Point(c,2*x+1,2*y+1);
      }
    } 
  }
}

void PDE::SaveGradientMagnitude(const int l, const char * fname){
// Save the magnitude of the grandient in layer l
  double mag, magx, magy;
  FILE *s;
  s=fopen(fname,"wt");
  
  for (int y=2;y<sizey-2;y++) {
    for (int x=2;x<sizex-2;x++) {
      magx = (sigma[l][x+1][y]-sigma[l][x-1][y])/2.;
      magy = (sigma[l][x][y+1]-sigma[l][x][y-1])/2.;
      mag = sqrt( magx*magx + magy*magy );
      fprintf (s,"%i\t%i\t%lf\n",x, y, mag);
    } 
  }
  fprintf(s,"\n");
  fclose (s);
}

void PDE::SaveConcentration(const int l, const char * fname){
// Save concentrations in layer l
  FILE *s;
  s=fopen(fname,"wt");
  
  for (int y=1;y<sizey-1;y++) {
    for (int x=1;x<sizex-1;x++) {
      fprintf (s,"%i\t%i\t%lf\n",x, y, sigma[l][x][y]);
    } 
  }
  fprintf(s,"\n");
  fclose (s);
}

void PDE::ContourPlot(Graphics *g, int l, int colour) {
  
  // calls "conrec" routine by Paul Bourke, as downloaded from
  // http://astronomy.swin.edu.au/~pbourke/projection/conrec

  // number of contouring levels
  int nc = 10;

  // A one dimensional array z(0:nc-1) that saves as a list of the contour levels in increasing order.   
  double *z=(double *)malloc(nc*sizeof(double));
  double min=Min(l), max=Max(l);
  double step=(max-min)/nc;
  {for (int i=0;i<nc;i++)
    z[i]=(i+1)*step;}
  
  double *x=(double *)malloc(sizex*sizeof(double));
  {for (int i=0;i<sizex;i++)
    x[i]=i;}
  
  double *y=(double *)malloc(sizey*sizeof(double));
  {for (int i=0;i<sizey;i++)
    y[i]=i;}
  
  conrec(sigma[l],0,sizex-1,0,sizey-1,x,y,nc,z,g,colour);
  
  free(x);
  free(y);
  free(z);
  
 
}



// public
void PDE::Secrete(CellularPotts *cpm, int ipmcs) {

  //const double dt=par.dt;
  double dt=1.0/(double)(ipmcs);
  int s,tau;

  for (int x=0;x<sizex;x++){
    for (int y=0;y<sizey;y++) {
	
	s=cpm->Sigma(x,y);
	if (s>0){
	  tau=(cpm->getCell(s)).getTau();
	  
	  if (tau==1) {					// neural crest
            sigma[0][x][y]+=par.secr_rate[0]*dt;	// secrete CoA
	    sigma[1][x][y]-=par.Sdf1UptakeNC*dt;	// take up Sdf1
	  }
	
	  if (tau==2) {					// placode 
            sigma[1][x][y]+=par.secr_rate[1]*dt;	// secrete Sd1
            sigma[2][x][y]+=par.secr_rate[2]*dt;	// secrete inhibitory molecule
	  }
	} else {
          // Decay only outside cellular regions:
          for (int l=0;l<layers;l++) {
            sigma[l][x][y]-=par.decay_rate[l]*dt*sigma[l][x][y];
          }
        }
        // Decay everywhere
	//for (int l=0;l<layers;l++) {
        //  sigma[l][x][y]-=par.decay_rate[l]*dt*sigma[l][x][y];
        //}
	
	if (sigma[0][x][y]<0) sigma[0][x][y]=0;	// boundary condition
	if (sigma[1][x][y]<0) sigma[1][x][y]=0;
	if (sigma[2][x][y]<0) sigma[2][x][y]=0;
    }
  }
  if (par.pdeExtendedField){
    // decay outside the cells in the extended area:
    for (int x=1;x<(sizex+200-1);x++){
      for (int y=1;y<101;y++) {
        for (int l=0;l<layers;l++) {
          sigmaExt[l][x][y]-=par.decay_rate[l]*dt*sigmaExt[l][x][y];
	}
      }
      for (int y=sizey;y<sizey+100-1;y++) {
        for (int l=0;l<layers;l++) {
          sigmaExt[l][x][y]-=par.decay_rate[l]*dt*sigmaExt[l][x][y];
	}
      }
    }
    for (int y=101;y<(sizey+101);y++) {
      for (int x=1;x<101;x++){
        for (int l=0;l<layers;l++) {
          sigmaExt[l][x][y]-=par.decay_rate[l]*dt*sigmaExt[l][x][y];
	}
      }
      for (int x=sizex;x<sizex+100-1;x++) {
        for (int l=0;l<layers;l++) {
          sigmaExt[l][x][y]-=par.decay_rate[l]*dt*sigmaExt[l][x][y];
	}
      }
    }
  }

}
                                            
void PDE::ClearMatrix(CellularPotts *cpm, int layer) {

  for (int x=0;x<sizex;x++) {
    for (int y=0;y<sizey;y++) {
	if (cpm->Sigma(x,y)) {
	    sigma[layer][x][y]=0;
	}
    }
  }
}

void PDE::ClearMatrix(int x, int y) {
  for (int i=0;i<layers;i++) sigma[i][x][y]=0;
}

void PDE::Diffuse(int ipmcs) {
 
  // Just diffuse everywhere (cells are transparent), using finite difference
  // (We're ignoring the problem of how to cope with moving cell
  // boundaries right now)
  
  double dt=1.0/(double)(ipmcs);
  int sw=0;

  if (par.pdeExtendedField){
    // solve diffusion on a larger grid to avoid boundary effects
    for (int l=0; l<layers; l++){
      for (int x=0;x<sizex; x++){
        for (int y=0;y<sizey; y++){
          sigmaExt[l][x+100][y+100]=sigma[l][x][y];
	}
      }
    }
    
    AbsorbingBoundariesExt();
    for (int l=0;l<layers;l++) {
      for (int x=1;x<(sizex+200-1);x++){
  	for (int y=1;y<(sizey+200-1);y++) {
  	  
  	  double sum=0.;
  	  sum+=sigmaExt[l][x+1][y];
  	  sum+=sigmaExt[l][x-1][y];
  	  sum+=sigmaExt[l][x][y+1];
  	  sum+=sigmaExt[l][x][y-1];
        
  	  sum-=4*sigmaExt[l][x][y];
  	  alt_sigmaExt[l][x][y]=sigmaExt[l][x][y]+sum*dt*par.diff_coeff[l];
  	  if ( fabs(alt_sigmaExt[l][x][y]) > 1.7e307 ){
  	    alt_sigmaExt[l][x][y]=alt_sigmaExt[l][x][y]/2.;
  	    sw=1;
  	  };
        }
      }
      if (sw==1) fprintf(stderr,"Diffusion overflow\n");
    }
    
    // Swap alt_sigmaExt and sigmaExt around:
    double ***tmp;
    tmp=sigmaExt;
    sigmaExt=alt_sigmaExt;
    alt_sigmaExt=tmp;
    
    // Copy the correct region of sigmaExt into the real sigma matrix:
    for (int l=0; l<layers; l++){
      for (int x=0;x<sizex; x++){
        for (int y=0;y<sizey; y++){
          sigma[l][x][y]=sigmaExt[l][x+100][y+100];
        }
      }
    }

  } else {
    //NoFluxBoundaries();
    if (par.periodic_boundaries) {
      PeriodicBoundaries();
    } else {
      AbsorbingBoundaries();
      //NoFluxBoundaries();
    }
    
    for (int l=0;l<layers;l++) {
      for (int x=1;x<sizex-1;x++) {
	for (int y=1;y<sizey-1;y++) {
	  
	  double sum=0.;
	  sum+=sigma[l][x+1][y];
	  sum+=sigma[l][x-1][y];
	  sum+=sigma[l][x][y+1];
	  sum+=sigma[l][x][y-1];
      
	  sum-=4*sigma[l][x][y];
	  alt_sigma[l][x][y]=sigma[l][x][y]+sum*dt*par.diff_coeff[l];
	  if ( fabs(alt_sigma[l][x][y]) > 1.7e307 ){
	    alt_sigma[l][x][y]=alt_sigma[l][x][y]/2.;
	    sw=1;
	  };
	}
      }
    }
    if (sw==1) fprintf(stderr,"Diffusion overflow\n");
    
    for (int l=0; l<layers; l++){
      for (int x=0;x<sizex; x++){
        for (int y=0;y<sizey; y++){
          sigma[l][x][y]=alt_sigma[l][x][y];
        }
      }
    }
  }
  
  
  thetime+=dt;
}

//public
void PDE::IteratePDE(CellularPotts *cpm) {
    int a, b=1;

    for (int l=0; l<layers; l++){
      a=1+int(par.diff_coeff[l] / 0.1);
      if (b<a) b=a;
    }
    for (int r=0; r<b; r++){
      Secrete(cpm, b);
      //ClearMatrix(cpm,0);
      Diffuse(b);
      //ClearMatrix(cpm,0);
    }

}

double PDE::GetChemAmount(const int layer) {
  
  // Sum the total amount of chemical in the lattice
  // in layer l
  // (This is useful to check particle conservation)
  double sum=0.;
  if (layer==-1) { // default argument: sum all chemical species
    for (int l=0;l<layers;l++) {
      for (int x=1;x<sizex-1;x++)
	for (int y=1;y<sizey-1;y++) {
	  
	  sum+=sigma[l][x][y];
	}
      
    }
  } else {
    for (int x=1;x<sizex-1;x++)
      for (int y=1;y<sizey-1;y++) {
	
	sum+=sigma[layer][x][y];
      }
  } 
  return sum;
  
}

// private
void PDE::NoFluxBoundaries(void) {

  // all gradients at the edges become zero, 
  // so nothing flows out
  // Note that four corners points are not defined (0.)
  // but they aren't used in the calculations
  
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      sigma[l][x][0]=sigma[l][x][1];
      sigma[l][x][sizey-1]=sigma[l][x][sizey-2];
    }
  
    for (int y=0;y<sizey;y++) {
      sigma[l][0][y]=sigma[l][1][y];
      sigma[l][sizex-1][y]=sigma[l][sizex-2][y];
    }
  }
}


// private
void PDE::AbsorbingBoundaries(void) {

  // all boundaries are sinks, 
  
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      sigma[l][x][0]=0.;
      sigma[l][x][sizey-1]=0.;
    }
  
    for (int y=0;y<sizey;y++) {
      sigma[l][0][y]=0.;
      sigma[l][sizex-1][y]=0.;
    }
  }
}

// private
void PDE::AbsorbingBoundariesExt(void) {

  // all boundaries are sinks, 
  
  for (int l=0;l<layers;l++) {
    for (int x=0;x<(sizex+200);x++) {
      sigmaExt[l][x][0]=0.;
      sigmaExt[l][x][sizey+200-1]=0.;
    }
  
    for (int y=0;y<(sizey+200);y++) {
      sigmaExt[l][0][y]=0.;
      sigmaExt[l][sizex+200-1][y]=0.;
    }
  }
}

// private
void PDE::GradientBoundaries(void) {

  // all boundaries are sinks, 
  
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      sigma[l][x][0]=2. * sigma[l][x][1] - sigma[l][x][2];
      if (sigma[l][x][0] < 0) sigma[l][x][0]=0;
      sigma[l][x][sizey-1]=2. * sigma[l][x][sizey-2] - sigma[l][x][sizey-3];
      if (sigma[l][x][sizey-1] < 0) sigma[l][x][sizey-1]=0;
    }
  
    for (int y=0;y<sizey;y++) {
      sigma[l][0][y]=2. * sigma[l][1][y] - sigma[l][2][y];
      if (sigma[l][0][y] < 0) sigma[l][0][y]=0;
      sigma[l][sizex-1][y]=2. * sigma[l][sizex-2][y] - sigma[l][sizey-3][y];
      if (sigma[l][sizex-1][y] < 0) sigma[l][sizex-1][y]=0;
    }
  }
}

// private
void PDE::PeriodicBoundaries(void) {

  // periodic...
  
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      sigma[l][x][0]=sigma[l][x][sizey-2];
      sigma[l][x][sizey-1]=sigma[l][x][1];
    }
    for (int y=0;y<sizey;y++) {
      sigma[l][0][y]=sigma[l][sizex-2][y];
      sigma[l][sizex-1][y]=sigma[l][1][y];
    }
  }
}

void PDE::GradC(int layer, int first_grad_layer) {
  
  // calculate the first and second order gradients and put
  // them in the next chemical fields
  if (par.n_chem<5) {
    throw("PDE::GradC: Not enough chemical fields");
  }
  
  // GradX
  for (int y=0;y<sizey;y++) {
    for (int x=1;x<sizex-1;x++) {
      sigma[first_grad_layer][x][y]=(sigma[layer][x+1][y]-sigma[layer][x-1][y])/2.;
    } 
  }
  
  // GradY
  for (int x=0;x<sizex;x++) {
    for (int y=1;y<sizey-1;y++) {
      sigma[first_grad_layer+1][x][y]=(sigma[layer][x][y+1]-sigma[layer][x][y-1])/2.;
    } 
  }

  // GradXX
  for (int y=0;y<sizey;y++) {
    for (int x=1;x<sizex-1;x++) {
      sigma[first_grad_layer+2][x][y]=sigma[layer][x+1][y]-sigma[layer][x-1][y]-2*sigma[layer][x][y];
    } 
  }

  // GradYY
  for (int x=0;x<sizex;x++) {
    for (int y=1;y<sizey-1;y++) {
      sigma[first_grad_layer+3][x][y]=sigma[layer][x][y-1]-sigma[layer][x][y+1]-2*sigma[layer][x][y];
    } 
  }
}

void PDE::PlotVectorField(Graphics &g, int stride, int linelength, int first_grad_layer) {
  
  // Plot vector field assuming it's in layer 1 and 2
  for (int x=1;x<sizex-1;x+=stride) {
    for (int y=1;y<sizey-1;y+=stride) {
      
      // calculate line
      int x1,y1,x2,y2;
      
      x1=(int)(x-linelength*sigma[first_grad_layer][x][y]);
      y1=(int)(y-linelength*sigma[first_grad_layer+1][x][y]);
      x2=(int)(x+linelength*sigma[first_grad_layer][x][y]);
      y2=(int)(y+linelength*sigma[first_grad_layer+1][x][y]);
      
      if (x1<0) x1=0;
      if (x1>sizex-1) x1=sizex-1;
      if (y1<0) y1=0;
      if (y1>sizey-1) y1=sizey-1;
      
      if (x2<0) x2=0;
      if (x2>sizex-1) x2=sizex-1;
      if (y2<0) y2=0;
      if (y2>sizey-1) y2=sizey-1;

      // And draw it :-)
      // perhaps I can add arrowheads later to make it even nicer :-)
      
      g.Line(2*x1,2*y1,2*x2,2*y2,1);
                  
    }
  }
  
}
