/* 

Copyright 1995-2006 Roeland Merks, Nick Savill

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


/* CA.cpp: implementation of Glazier & Graner's Cellular Potts Model */

// This code derives from a Cellular Potts implementation written around 1995
// by Nick Savill

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include "sticky.h"
#include "random.h"
#include "ca.h"
#include "parameter.h"
#include "dish.h"
#include "sqr.h"
#include "crash.h"

#define PI (4.0*atan2(1.0,1.0))

/* STATIC DATA MEMBER INITIALISATION */
double copyprob[BOLTZMANN]; 
double disc_fact=100;
int **countnum;
int **countnumA;
int NCType=1;
int PlacodeType=2;

const int CellularPotts::nx[25] = {0, 0, 1, 0,-1, 1, 1,-1,-1, 0, 2, 0, -2, 1, 2, 2, 1,-1,-2,-2,-1, 0, 2, 0,-2 };
const int CellularPotts::ny[25] = {0,-1, 0, 1, 0,-1, 1, 1,-1,-2, 0, 2,  0,-2,-1, 1, 2, 2, 1,-1,-2,-2, 0, 2, 0 };

const int CellularPotts::nbh_level[5] = { 0, 4, 8, 20, 24 };
int CellularPotts::shuffleindex[9]={0,1,2,3,4,5,6,7,8};
double CellularPotts::Tscale;

extern Parameter par;

/* global variables */


/** PRIVATE **/

using namespace std;
void CellularPotts::BaseInitialisation(vector<Cell> *cells) {
  CopyProb(disc_fact);
  cell=cells;
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).";
  
}

CellularPotts::CellularPotts(vector<Cell> *cells,
			     const int sx, const int sy) {


  
  sigma=0;
  frozen=false;
  thetime=0;

  
  BaseInitialisation(cells);
  sizex=sx;
  sizey=sy;

  AllocateSigma(sx,sy);
									
  
  // fill borders with special border state
  for (int x=0;x<sizex;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }
  
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";

}

CellularPotts::CellularPotts(void) {



  sigma=0;
  sizex=0; sizey=0;
  frozen=false;
  thetime=0;

  CopyProb(disc_fact);

  // fill borders with special border state
  for (int x=0;x<sizex;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
}

// destructor (virtual)
CellularPotts::~CellularPotts(void) {
  if (sigma) {
    free(sigma[0]);
    free(sigma);
    sigma=0;
  }
}


void CellularPotts::AllocateSigma(int sx, int sy) {

  sizex=sx; sizey=sy;
  
  sigma=(int **)malloc(sizex*sizeof(int *));
  if (sigma==NULL)
    MemoryWarning();
  
  sigma[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (sigma[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<sizex;i++) 
    sigma[i]=sigma[i-1]+sizey;}
  
  /* Clear CA plane */
   {for (int i=0;i<sizex*sizey;i++) 
     sigma[0][i]=0; }

}

void CellularPotts::IndexShuffle() {

  int i;
  int temp;
  int index1,index2;
  
  for (i=0;i<9;i++) {
    
    index1=RandomNumber(8);
    index2=RandomNumber(8);

    temp=shuffleindex[index1];
    shuffleindex[index1]=shuffleindex[index2];
    shuffleindex[index2]=temp;

  }
}


double sat(double x) {
  
  return x/(par.saturation*x+1.);

}

double CellularPotts::DeltaH(int x,int y, int xp, int yp, PDE *PDEfield, int active)       
{
//
// Get energy-difference for copying sigma[xp][yp] into sigma[x][y].
// If (x,y) is the same as (xp,yp), then it means that a spin is not 
// copied from one site to another but a spin=0 (cell-free area) is 
// inserted at that spot. 
// Switch 'active' is controlling which components to include in the
// energy calculation. By default, active=1, which is a normal MCS.
// If active==0, this is a relaxation step, so 'active' cell 
// behaviour is excluded. 
//

  double dH = 0, dDH = 0, H;
  int i, sxy, sxyp;
  int neighsite;
  double dx,dy;


  //--- Compute energydifference *IF* the copying were to occur ---//
  sxy = sigma[x][y];
  sxyp = sigma[xp][yp];





  //--- Terms of the energy-function ---//

  // COMPONENTS FOR BOTH ACTIVE AND PASSIVE MODES //
  //--- Cell adhesion due to surface tension ---//
  for (i=1;i<=n_nb;i++) {
    int xp2,yp2;
    xp2=x+nx[i]; yp2=y+ny[i];
    if (par.periodic_boundaries) {
      
      // since we are asynchronic, we cannot just copy the borders once 
      // every MCS
      
      if (xp2<=0) xp2=sizex-2+xp2;
      if (yp2<=0) yp2=sizey-2+yp2;
      if (xp2>=sizex-1) xp2=xp2-sizex+2;
      if (yp2>=sizey-1) yp2=yp2-sizey+2;
    
      neighsite=sigma[xp2][yp2];
    } else {
      if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1) neighsite=-1;
      else neighsite=sigma[xp2][yp2];
    }
    
    if (x==xp && y==yp) {				//inserting medium
      if (neighsite==-1) {
        H = par.border_energy;			// border
        dH-= H; 
      } else {
        H = (*cell)[MEDIUM].EnergyDifference((*cell)[neighsite])-(*cell)[sxy].EnergyDifference((*cell)[neighsite]);
        dH+= H;
      }
    } else {						//spin copy
      if (neighsite==-1){ 				// border 
        H = (sxyp==MEDIUM?0:par.border_energy)-(sxy==MEDIUM?0:par.border_energy);
	dH += H;
      } else {  
	H = (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) -  (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
	dH += H;
      }
    }
  }
  
  
  //--- Area constraint ---//
  if (par.lambda != 0){
    if ( sxyp == MEDIUM ) {			// retraction from MEDIUM
      H = par.lambda * (1. - 2. * (double) ( (*cell)[sxy].Area() - (*cell)[sxy].TargetArea()) );
      dH += H;
    } else if ( sxy == MEDIUM ) {			// extension into MEDIUM
      H = par.lambda * (1. + 2. * (double) ( (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) );
      dH += H;
    } else if (sxy==sxyp) { 			// inserting medium (same as retraction from MEDIUM)
      H = par.lambda * (1. - 2. * (double) ( (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) );
      dH += H;
    } else {					// two cell border
      H = par.lambda * (2. + 2. * (double)( ((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) 
  	    				    - ((*cell)[sxy].Area()  - (*cell)[sxy].TargetArea())) );
      dH += H;
    };
  }

  
  //--- Length constraint ---//
  // sp is expanding cell, s is retracting cell
 
  if (par.lambda2 != 0){ 
    if ( sxyp == MEDIUM ) {	// retraction
      double a=DSQR(((*cell)[sxy].Length()-(*cell)[sxy].TargetLength()));
      double b=DSQR(((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) - (*cell)[sxy].TargetLength()));
      H = par.lambda2*(a-b);
    } else if ( sxy == MEDIUM ) {	// extension
      double a=DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength());
      double b=DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y) - (*cell)[sxyp].TargetLength());
      H=par.lambda2*(a-b);
    } else {			// cell-cell border
      double a=DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength());
      double b=DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength());
      double c=DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength());
      double d=DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) - (*cell)[sxy].TargetLength());
      H=par.lambda2*(a-b+c-d);
    }
    dH -= H;
  }


  //--- Cytoskeletal bonds ---//
  if (par.cytosk_bonds){
    if ( (sxy > MEDIUM && (*cell)[sxy].AliveP())
      || (sxyp> MEDIUM && (*cell)[sxyp].AliveP()) ){
	// Check only if either of the two sites belongs to a living cell

      double sumfR_x, sumfR_y;
      double sumfE_x, sumfE_y;
      double cxR=0., cyR=0.; 
      double cxE=0., cyE=0.; 
      double cx2=0., cy2=0.; 
      double dcx,dcy,dlen,ddlen;
      int tauR=0, tauE=0, tau2=0;
      bool searchR = false, searchE = false;
    
      if (sxy > MEDIUM && (*cell)[sxy].AliveP()){		
							// retracting cell (pulling)
        searchR = true;
	tauR = (*cell)[sxy].tau;
        sumfR_x=0; sumfR_y=0;
        (*cell)[sxy].CenterOfMass(&cxR, &cyR);
      }
      if (sxyp > MEDIUM && (*cell)[sxyp].AliveP() && !(x==xp && y==yp)){		
							// extending cell (pushing)
        searchE = true;				// only consider if not inserting medium: (x,y)!=(xp,yp)
	tauE = (*cell)[sxyp].tau;
        sumfE_x=0; sumfE_y=0;
        (*cell)[sxyp].CenterOfMass(&cxE, &cyE);
      }

      for (i=1; i< (int)(cell->size()); i++){		// go through cells to find bound neighbours 
							// and get sum force from connections
	if (searchR){
	  if ((*cell)[sxy].cytosk_bonds[i] >= 0 && (*cell)[i].AliveP()){
            tau2 = (*cell)[i].tau;
            (*cell)[i].CenterOfMass(&cx2, &cy2);
	    dlen = sqrt((cx2-cxR)*(cx2-cxR) + (cy2-cyR)*(cy2-cyR));
 	    dcx = cx2 - cxR; if (dlen != 0) {dcx /= dlen;}
 	    dcy = cy2 - cyR; if (dlen != 0) {dcy /= dlen;}
 
 	    ddlen = (dlen - (*cell)[sxy].cytosk_bonds[i]);
 
 	    sumfR_x += dcx * ddlen * par.SC[tauR][tau2];
 	    sumfR_y += dcy * ddlen * par.SC[tauR][tau2];
          }
        }
	if (searchE){
	  if ((*cell)[sxyp].cytosk_bonds[i] >= 0 && (*cell)[i].AliveP()){
            tau2 = (*cell)[i].tau;
            (*cell)[i].CenterOfMass(&cx2, &cy2);
	    dlen = sqrt((cx2-cxE)*(cx2-cxE) + (cy2-cyE)*(cy2-cyE));
 	    dcx = cx2 - cxE; if (dlen != 0) {dcx /= dlen;}
 	    dcy = cy2 - cyE; if (dlen != 0) {dcy /= dlen;}
 
 	    ddlen = (dlen - (*cell)[sxyp].cytosk_bonds[i]);
 
 	    sumfE_x += dcx * ddlen * par.SC[tauE][tau2];
 	    sumfE_y += dcy * ddlen * par.SC[tauE][tau2];
	  }      
	}
      }
							// Calculate energy bias for the cells;
      if (searchR){					// for the retracting cell
        // Movement direction:
	(*cell)[sxy].DecrementArea();
        (*cell)[sxy].RemoveSiteFromMoments(x,y);
        (*cell)[sxy].CenterOfMass(&cx2,&cy2);
        (*cell)[sxy].IncrementArea();
        (*cell)[sxy].AddSiteToMoments(x,y);
        dx=cx2-cxR;
        dy=cy2-cyR;
        double norm = sqrt((dx*dx)+(dy*dy));
        if (norm != 0) { dx=dx/norm; dy=dy/norm;}
        else {dx=0; dy=0;}
      
	// Get energy change from (direction of movement)*(sum force):
        dH -= (dx * sumfR_x + dy * sumfR_y);
      }
      if (searchE){					// for the expanding cell
	// Movement direction:
        (*cell)[sxyp].IncrementArea();
        (*cell)[sxyp].AddSiteToMoments(x,y);
        (*cell)[sxyp].CenterOfMass(&cx2,&cy2);
        (*cell)[sxyp].DecrementArea();
        (*cell)[sxyp].RemoveSiteFromMoments(x,y);
        dx=cx2-cxE;
        dy=cy2-cyE;
        double norm = sqrt((dx*dx)+(dy*dy));
        if (norm != 0) { dx=dx/norm; dy=dy/norm;}
        else {dx=0; dy=0;}
      
	// Get energy change from (direction of movement)*(sum force):
        dH -= (dx * sumfE_x + dy * sumfE_y);
      }
    }
  }



  //--- ACTIVE COMPONENTS ---//
  if (active==1){
    //--- Chemotaxis ---//
    if (PDEfield){
      dDH=0;
      if ( (par.vecadherinknockout || (sxyp==MEDIUM || sxy==MEDIUM))
           && (!(par.extension_only_chemotaxis && sxyp == MEDIUM))
  	){
  		// Conditions for chemotaxis:
  		// 1. chemotax only at free boundaries, except if vecadherinknockout is activated
		// 2. only allow extension, except if extension_only_chemotaxis is not true
  		//
  	double c=0, nc=0;
  	for (int l=0;l<par.n_chem;l++){			// chemical 0: CoA, chemical 1: Sdf1
    	  nc=0.0; c=0.0;
  	  if (sxy>MEDIUM)  {c+=(*cell)[sxy].getChem(l);  nc++;}	// just in case extension only is not wanted
  	  if (sxyp>MEDIUM) {c+=(*cell)[sxyp].getChem(l); nc++;}
  	  if (nc>0) c /= nc;
  	  dDH+=(double)(c*(sat(PDEfield->Sigma(l,x,y))-sat(PDEfield->Sigma(l,xp,yp))));
  	}
  	if (dDH > 0.0) {
  	  // Make a note that the cell is chemotaxing:
	  if (sxy>MEDIUM)  (*cell)[sxy].chemotaxing  = true;		// just in case extension only is not wanted
  	  if (sxyp>MEDIUM) (*cell)[sxyp].chemotaxing = true;
  	}
      }
      dH-=dDH;
    }


    //--- Polarisation vector ---//
    if ((sxyp && ((*cell)[sxyp].getPolStr() != 0)) || (sxy && ((*cell)[sxy].getPolStr() != 0))){ 
      
      double px=0.0,py=0.0,norm=0.0,p=0.0;
      double cx1=0.0,cy1=0.0,cx2=0.0,cy2=0.0, sss=0;
      int tc;
  
  
      if (					// expanding cell
      	(tc = sxyp) && 				// valid cell
  	(*cell)[sxyp].AliveP() && 		// alive cell
  	(sxyp!=sxy)  				// not inserting medium (in which case there is no extending cell)
  	){				
  					
        (*cell)[tc].CenterOfMass(&cx1,&cy1);	// get current and proposed CoM of cell
        (*cell)[tc].IncrementArea();
        (*cell)[tc].AddSiteToMoments(x,y);
        (*cell)[tc].CenterOfMass(&cx2,&cy2);
        (*cell)[tc].DecrementArea();
        (*cell)[tc].RemoveSiteFromMoments(x,y);
        dx=cx2-cx1;
        dy=cy2-cy1;
        sss = (dx*dx)+(dy*dy);
	if (sss > 0){
		norm = sqrt(sss);
        	dx=dx/norm; dy=dy/norm;		// displacement vector
	}
        
        p=(*cell)[tc].getPolStr();		// polarity
        (*cell)[tc].getPolarisation (&px, &py);
	sss = (px * px) + (py * py);
	if (sss>0){
          norm = sqrt(sss);
  	  px = px/norm; py = py/norm;	
        } else {norm=0; px = 0; py = 0;};	// polarity vector
        H =  p*(dx*px + dy*py);			// polarity bias: dr * P
        dH -= H;
      }
      
      if ((tc = sxy)) {				//retracting cell
        if ((*cell)[tc].Area() > 1){		// get current and proposed CoM of cell
  	  (*cell)[tc].CenterOfMass(&cx1,&cy1);
          (*cell)[tc].DecrementArea();
          (*cell)[tc].RemoveSiteFromMoments(x,y);
          (*cell)[tc].CenterOfMass(&cx2,&cy2);
          (*cell)[tc].IncrementArea();
          (*cell)[tc].AddSiteToMoments(x,y);
          dx=cx2-cx1;
          dy=cy2-cy1;
        } else {
          dx=x-xp; 
          dy=y-yp;
        };
        sss = (dx*dx)+(dy*dy);
	if (sss > 0){
        	norm = sqrt(sss);
		dx/=norm; dy/=norm;		// displacement vector
        }

        p=(*cell)[tc].getPolStr();		// polarity
        (*cell)[tc].getPolarisation (&px, &py);
        sss = px*px + py*py;
        if (sss>0) { 
  	  norm = sqrt(sss); 
	  px = px/norm; py = py/norm;		
        } else {norm=0; px = 0; py = 0;};		// polarity vector
        H =  p*(dx*px + dy*py) ;		// polarity bias: dr * P
        dH -= H;
      }
    }
  }
  //--- END OF ACTIVE COMPONENTS ---//


  return dH;
}




void CellularPotts::ConvertSpin(int x,int y,int xp,int yp)
{

  if (x==xp && y==yp) {		//inserting medium
    (*cell)[(sigma[x][y])].DecrementArea();
    (*cell)[(sigma[x][y])].RemoveSiteFromMoments(x,y);
    Displace (sigma[x][y],x,y);
    
    if (!(*cell)[(sigma[x][y])].Area()) {
      (*cell)[(sigma[x][y])].Apoptose();
      cerr << "Cell " << sigma[x][y] << " apoptosed\n";
    }
    sigma[x][y]=0;
  
  }else{			//regular conversion
    if (sigma[x][y]) { 		// if (x,y) is not MEDIUM
      (*cell)[sigma[x][y]].DecrementArea();
      (*cell)[sigma[x][y]].RemoveSiteFromMoments(x,y);
      Displace (sigma[x][y],x,y);

      if (!(*cell)[sigma[x][y]].Area()) {
        (*cell)[sigma[x][y]].Apoptose();
        cerr << "Cell " << sigma[x][y] << " apoptosed\n";
      }
    }
  
    if (sigma[xp][yp]) {		// if (xp,yp) is not MEDIUM
      (*cell)[sigma[xp][yp]].IncrementArea();
      (*cell)[sigma[xp][yp]].AddSiteToMoments(x,y);
      Displace (sigma[xp][yp],x,y);
    }
    sigma[x][y] = sigma[xp][yp];
  }
};

/** PUBLIC **/

int CellularPotts::CopyvProb(double DH,  double stiff) {

  double dd; 
  int s;

//  s=(int)stiff;
  if (DH<=-stiff) return 2;
  else {
    s=(int)((DH+stiff) * disc_fact / (double)(par.T));
    // if DH becomes extremely large, probability becomes zero
    if ((s) > (BOLTZMANN-1)) return 0;
    else {
      dd=copyprob[s]; 
      if (RANDOM()<dd) return 1; 
      else return 0;
    }
  }
} 

void CellularPotts::CopyProb(double f) {
  int i;
  for ( i = 0; i < BOLTZMANN; i++ )
    copyprob[i] = exp( -( (double)(i)/f ) );
}

void CellularPotts::FreezeAmoebae(void) 
{
  if (frozen) 
    frozen=FALSE;
  else
    frozen=TRUE;
}

#include <fstream>
void CellularPotts::AmoebaeMove(PDE *PDEfield)
{
  int loop,p;
  //int updated=0;
  thetime++;
  
  
  if (frozen) 
    return;
  if (par.print_sigmaB == true) PrintSigma ();
  if (par.print_CMB == true) PrintCM ();
  if (par.print_AreasB == true) PrintAreas ();

  if (par.time==0)	//initial arrangements
    SpecInitConditions();
   


  loop=(sizex-2)*(sizey-2);
  //Stochastic timestep:
  for (int i=0;i<loop;i++) {  		//beginning of MCS

    
    // take a random site
    int xy = (int)(RANDOM()*(sizex-2)*(sizey-2));
    int x = xy%(sizex-2)+1;
    int y = xy/(sizex-2)+1; 
   
    // take a random neighbour
    int xyp=(int)((n_nb)*RANDOM()+1);
    int xp = nx[xyp]+x;
    int yp = ny[xyp]+y;
     
    //--- See if we want to insert medium ---//
    if (RANDOM() < par.insertMediumProb) {	// Insert medium with probability
      xyp=0; xp=x; yp=y;
    }	
    
    if (xyp==0) {
      // Inserting medium: 
      // Calling the functions below with the same coordinates for target and source means inserting medium
      double H_diss=0;
      if (CheckConnectivityConstraints(x,y,x,y)){
        double D_H=DeltaH(x,y,x,y,PDEfield,1);
        if ((p=CopyvProb(D_H,H_diss))>0) {
          ConvertSpin (x,y,x,y);
        }
      }
    } else {
      // Normal sigma-copy attempt
      int sxy=sigma[x][y];
      int sxyp;
      if (par.periodic_boundaries) {
        // since we are asynchronic, we cannot just copy the borders once 
        // every MCS
        if (xp<=0) xp=sizex-2+xp;
        if (yp<=0) yp=sizey-2+yp;
        if (xp>=sizex-1) xp=xp-sizex+2;
        if (yp>=sizey-1) yp=yp-sizey+2;
        sxyp=sigma[xp][yp];
      } else {
        if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1) sxyp=-1;
        else sxyp=sigma[xp][yp];
      }
      // test for border state (relevant only if we do not use periodic boundaries)
      if ((sxyp!=-1) && (sxy!=sxyp)){
        // Don't copying border state, try to copy if sites do not belong to the same cell
        // connectivity dissipation:
        double H_diss=0;
	if (CheckConnectivityConstraints(xp,yp,x,y)) {
	  double D_H=DeltaH(x,y,xp,yp,PDEfield,1);
	  if ((p=CopyvProb(D_H,H_diss))>0) {
            ConvertSpin (x,y,xp,yp);
	  }
        }
      }
    }
  }						//end of MCS
  

  
  setColours();

  if (!par.internalrelaxation){			//if relax. period is over update cell properties & output
    InsertNCCellsAtTop();				// check if new cells will be entered in the simulation
    if (par.takeoutCellsBelow >= 0) takeoutCellsBelow(par.takeoutCellsBelow);
    UpdateNeighbours();					// update links and contact information
    PolUpdate ();					// update polarity (using new contact information)
    if (par.PassiveRelaxationSteps) PassiveRelaxation(PDEfield);
							// insert relaxation steps to allow passive equilibration
    if (par.print_asParticles>0 && (par.time % par.print_asParticles) == 0) PrintAsParticles();
    if (par.print_PIF>0 && (par.time % par.print_PIF) == 0) PrintPIF();
  }
  SyncCellTime();

}


void CellularPotts::takeoutCellsBelow(double limit){
  double x,y;
  for (int i=1;i<(int)(cell->size());i++){
    (*cell)[i].CenterOfMass(&x, &y);
    if (y > (sizey-limit)){
      (*cell)[i].Apoptose();
      int ta = (*cell)[i].target_area;
      
      if (ta > 0) (*cell)[i].SetTargetArea(ta-1);
    }
  }
}


void CellularPotts::PassiveRelaxation(PDE *PDEfield){
// Run passive MCSteps until the system is relaxed
// Criterion for being relaxed is that the average (sum) of energy change (dH)
// in X steps (t-X --> t) is less than the dH in the X steps before (t-2X --> t-X).
// 

  double dH=0.0, prev_dH=0.0;
  int n=0, X=par.PassiveRelaxationWindow;
  bool cycle = true, just_started = true;
  char * fname;

  if (!(par.time % par.storage_stride) && par.PassiveRelaxationMonitorDH){
    fname=(char*)calloc(100,sizeof(char));
    sprintf(fname,"%s/cells_%05i-A.pif",par.datadir,par.time);
    PrintPIF(fname);
    free (fname);
  }
  

  while (cycle){
    dH += PassiveRelaxationStep(PDEfield);
    n+=1;
    if (!(n%X)){
      // Stop relaxation steps if the energy change from the X steps is less than 
      // the energy change in the X steps before this. This will stop the relaxation
      // process when the energy cannot be relaxed further with these steps. 

      if (!just_started && dH <= prev_dH){
        cycle = false;
      }

      if (par.PassiveRelaxationMonitorDH){
        // To monitor the relaxation process
        std::cout << "Relaxation-dH_" << X << "(cycle= " << n << " )= " << dH << "\n";
      }
      prev_dH = dH;
      dH = 0.0;
      just_started = false;
    }
  }
  
  if (!(par.time % par.storage_stride) && par.PassiveRelaxationMonitorDH){
    fname=(char*)calloc(100,sizeof(char));
    sprintf(fname,"%s/cells_%05i-B.pif",par.datadir,par.time);
    PrintPIF(fname);
    free (fname);
  }
}


void CellularPotts::getCellCellAdhesionTensions(double *tension, int *n_links){
// Return the value of tension in the cytoskeletal links
// 
  double cx1, cy1, cx2, cy2, dlen, sumf=0.0;
  int n = 0;

  for (int i=1; i<(int)(cell->size()); i++) {
    if ((*cell)[i].AliveP()){
      int tau1=(*cell)[i].tau;
      (*cell)[i].CenterOfMass(&cx1, &cy1);
      
      for (int j=1;j<i;j++) {
        if ((*cell)[i].cytosk_bonds[j] >= 0.0 && (*cell)[j].AliveP()){
	  int tau2=(*cell)[j].tau;
          (*cell)[j].CenterOfMass(&cx2, &cy2);

   	  dlen = sqrt(((cx2-cx1)*(cx2-cx1)) + ((cy2-cy1)*(cy2-cy1)));
   	  sumf += fabs( (dlen - (*cell)[i].cytosk_bonds[j]) * par.SC[tau1][tau2] );
	  n += 1;
        }
      }
    }
  }
  *tension = sumf;
  *n_links = n;
}


double CellularPotts::PassiveRelaxationStep(PDE *PDEfield){
// One MCS for relaxing cell-cell adhesion tension affecting only CIL-active cell-cell boundaries
// 
  double dH=0;

  int loop=(sizex-2)*(sizey-2);
  for (int i=0;i<loop;i++) {  		//beginning of MCS

    
    // take a random site
    int xy = (int)(RANDOM()*(sizex-2)*(sizey-2));
    int x = xy%(sizex-2)+1, y = xy/(sizex-2)+1; 
   
    // take a random neighbour
    int xyp=(int)((n_nb)*RANDOM()+1);
    int xp = nx[xyp]+x, yp = ny[xyp]+y;
     
    // (No medium insertion here)
    // determine the cells on the two selected sites, considering boundary conditions
    int sxy=sigma[x][y];
    int sxyp;
    if (par.periodic_boundaries) {
      // since we are asynchronic, we cannot just copy the borders once 
      // every MCS
      if (xp<=0) xp=sizex-2+xp;
      if (yp<=0) yp=sizey-2+yp;
      if (xp>=sizex-1) xp=xp-sizex+2;
      if (yp>=sizey-1) yp=yp-sizey+2;
      sxyp=sigma[xp][yp];
    } else {
      if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1) sxyp=-1;
      else sxyp=sigma[xp][yp];
    }
    
    // Attempt spin-copy only if not copying boundary state; the boundary is between 
    // two (different) CIL-active cells (not bound to underlying substrate). 
    if ((sxyp>0) && (sxy>0) && (sxy!=sxyp) && (*cell)[sxy].CILstate && (*cell)[sxyp].CILstate){
      // connectivity dissipation:
      double H_diss=0;
      if (CheckConnectivityConstraints(xp,yp,x,y)) {
	double D_H=DeltaH(x,y,xp,yp,PDEfield,0);
	if ((CopyvProb(D_H,H_diss))>0) {
          ConvertSpin (x,y,xp,yp);
	  dH += D_H;
	}
      }
    }
  }						//end of MCS
  
  for (int i=1;i<(int)(cell->size());i++){	// clear cell displacements so these passive movements 
    (*cell)[i].setDisplacement (0,0);		// are not registered in the dynamics of the polarisation
  }						// vector
  return dH; 
}



void CellularPotts::InsertNCCellsAtTop(){
// Check if cell-free space is available at the top of the simulation area to insert "virtual" NC cells.
//

  double areaFactor = 0.5;
  int cs = (int)(sqrt((double)(par.size_init_cells))), csl = (int)(areaFactor * cs)+1;
  int freearea = 0, startx = 1+par.sidePadding;

  // Check which cells are available (Apoptosed):
  std::vector <int> deadcells;

  for (int i=(int)(cell->size())-1; i>0;i--){
    if ((*cell)[i].Area() == 0){
      (*cell)[i].Apoptose();
      deadcells.push_back(i);
    }
  }


  // Search the top of the simulation for free space:
  for (int x=1+par.sidePadding;x<sizex-par.sidePadding;x++){
    for (int y=1;y<csl; y++){


      if (sigma[x][y] != MEDIUM) {
        // Found a cell in the area. Start search for empty space in the next column
        x++; y=1; freearea = 0; startx = x; if (x>=sizex-par.sidePadding) y=csl+1;
      
      } else {
        // This is a free area
        freearea+=1; 

        if (freearea >= (int)(areaFactor * par.size_init_cells)){
	 
          // Found sufficient free space, add a new cell if still available.
	  if (deadcells.size()>0){
	    int newcellid = deadcells.back();
	    deadcells.pop_back();

            (*cell)[newcellid].alive = true;
            // Assign index (sigma) and area of cell: 
	    for (int xc = startx; xc<=x; xc++){
	      for (int yc = 1; yc<csl; yc++){
	        if (freearea >0){
	          sigma[xc][yc] = newcellid; 
                  (*cell)[newcellid].IncrementArea();
                  (*cell)[newcellid].AddSiteToMoments(xc,yc);
		  freearea -= 1;
	        }
	      }
	    }

	    // Initialise other parameters of the new cell
	    (*cell)[newcellid].setTau(NCType);
	    (*cell)[newcellid].setPolarisation (0.0,0.0);
            (*cell)[newcellid].setDisplacement (0.0,0.0);
            (*cell)[newcellid].initCellParameters();
            (*cell)[newcellid].initContacts();
          }
	  // Clear counters and prepare for next possible cell
          x++; y=1; freearea = 0; startx = x; if (x>=sizex-par.sidePadding) y=csl+1;
	}
      }
    }
  }
  std::vector<int>().swap(deadcells);

}

void CellularPotts::UpdateNeighbours(){
    int tau1, tau2;
    double cm1_x, cm1_y, cm2_x, cm2_y;
    double dlen = 0;
    double eq_dist = 2.0 * sqrt(((double)par.target_area)/PI);
    std::vector<int> foundNeighbor[par.n_init_cells+1];

    // clear contact information
    for (int i=1;i<(int)(cell->size());i++) (*cell)[i].inContact=0;
    
    
    // creating new bonds
    for (int x = 1;x<sizex-1;x++) {
	for (int y = 1;y<sizey-1;y++) {
	    for (int n = 1;n<=n_nb;n++) {
		int xn = x+nx[n];
		int yn = y+ny[n];
		
		if (xn > 0 && xn < sizex-1 && yn > 0 && yn < sizey-1 && (xn > x || (xn == x && yn > y))) {
		    int sxy = sigma[x][y]; 
		    int sxyn = sigma[xn][yn]; 
		    if (sxy != MEDIUM && sxyn != MEDIUM && sxy != sxyn){
		        // sxy and sxyn are touching cells
			(*cell)[sxy].inContact=1;
		    	if ((*cell)[sxy].cytosk_bonds[sxyn] < 0 && RANDOM() < par.CC[(*cell)[sxy].tau][(*cell)[sxyn].tau]) {
  			  // sxy and sxyn were not connected, link them together
			  // NOTE: Use this below to set the equilibrium distance for the links 
			  //       equal to the initial distance of the cells when they come in 
			  //       contact. 
			  //(*cell)[sxy].CenterOfMass(&cm1_x, &cm1_y);
  			  //(*cell)[sxyn].CenterOfMass(&cm2_x, &cm2_y);
  			  //dlen = sqrt(((cm2_x-cm1_x)*(cm2_x-cm1_x)) + ((cm2_y-cm1_y)*(cm2_y-cm1_y)));
			  // NOTE: Use this below to set the equilibrium distance for the links
			  //       equal to a given constant distance eq_dist. 
			  dlen = eq_dist;

  			  (*cell)[sxy].cytosk_bonds[sxyn] = dlen;
  			  (*cell)[sxyn].cytosk_bonds[sxy] = dlen;
  			  //cout << "new bond: " << sxy << " + " << sxyn << endl;
			}
			if ((*cell)[sxy].cytosk_bonds[sxyn]>=0){ 
			  foundNeighbor[sxy].push_back(sxyn);
			}
		    }
		}
	    }
	}
    }
    
    //uncouple bonds

    for (int i=1; i<(int)(cell->size()); i++) {
	tau1 = (*cell)[i].tau;
	for (int j=1;j<i;j++) {
	    if (((*cell)[i].cytosk_bonds[j]>=0.) && (*cell)[i].AliveP() && (*cell)[j].AliveP()) {
		tau2 = (*cell)[j].tau;
		// Check if the two cells are still in contact
		bool inContact=false;
		for (std::vector<int>::iterator it=foundNeighbor[i].begin(); it!=foundNeighbor[i].end(); it++){

		  if (*it == j) {inContact=true; (*cell)[i].inContact=1;}
		}
		if (par.breakBondIfNoContact && !inContact){
		  // i and j are not in contact anymore
		  (*cell)[i].cytosk_bonds[j] = -1;
		  (*cell)[j].cytosk_bonds[i] = -1;
		  (*cell)[i].contactDurationWith[j] = 0;
		  (*cell)[j].contactDurationWith[i] = 0;
		} else {
		  // i and j are in contact
		  // Check distance of pair: if their centres of mass are too far, break link with probability
		  (*cell)[i].CenterOfMass(&cm1_x, &cm1_y);
		  (*cell)[j].CenterOfMass(&cm2_x, &cm2_y);
		  dlen = sqrt(((cm2_x-cm1_x)*(cm2_x-cm1_x)) + ((cm2_y-cm1_y)*(cm2_y-cm1_y)));
		  if ((dlen-(*cell)[i].cytosk_bonds[j]) > (RANDOM() * par.UCC[tau1][tau2])) {
		    (*cell)[i].cytosk_bonds[j] = -1;
		    (*cell)[j].cytosk_bonds[i] = -1;
		    //cout << "uncoupled bond: " << i+1 << " + " << j+1 << endl;
		    (*cell)[i].contactDurationWith[j] = 0;
		    (*cell)[j].contactDurationWith[i] = 0;
		  } else {
		    // i and j can remain in contact
		    (*cell)[i].contactDurationWith[j] += 1;
		    (*cell)[j].contactDurationWith[i] += 1;
		  }
		}
	    }
	}
    }
}



void CellularPotts::PlotCytoskBonds(Graphics *g, int colour) {
  double x1, y1, x2, y2;

  if (g && par.cytosk_bonds) {
    for (int i=1;i<=par.n_init_cells;i++) {
      for (int j=1;j<=par.n_init_cells;j++) {
        if (i > j && ((*cell)[i].cytosk_bonds[j]>=0) && (*cell)[i].AliveP() && (*cell)[j].AliveP()) {
          (*cell)[i].CenterOfMass(&x1,&y1);
          (*cell)[j].CenterOfMass(&x2,&y2);
          g->Line(2*x1,2*y1,2*x2,2*y2,colour);
        }
      }
    }
  }
}

void CellularPotts::PlotCILcells(Graphics *g, int colour) {
  double x1, y1;
  int size=2;

  if (g) {
    for (int i=1;i<=par.n_init_cells;i++) {
      if ((*cell)[i].CILstate) {
        (*cell)[i].CenterOfMass(&x1,&y1);
	for (int x=-size;x<size+1;x++){
	for (int y=-size;y<size+1;y++){
        g->Point(colour, 2*x1+x,2*y1+y);
        }}
      }
    }
  }
}

void CellularPotts::PlotPolarityVectors(Graphics *g, double size, int colour) {
  double x, y, dx, dy, dl;

  if (g) {
    for (int i=1;i<=par.n_init_cells;i++) {
        if ((*cell)[i].AliveP()) {
          (*cell)[i].CenterOfMass(&x,&y);
          (*cell)[i].getPolarisation (&dx,&dy);
	  dl=sqrt( dx*dx + dy*dy );
	  if (dl > 0) {dx *= size/dl; dy *= size/dl;}
	  dx = (x+dx); dy = (y+dy);
	  if (x<0) x=0; if (x>sizex) x=sizex; 
	  if (y<0) y=0; if (y>sizey) y=sizey; 
	  if (dx<0) dx=0; if (dx>sizex) dx=sizex; 
	  if (dy<0) dy=0; if (dy>sizey) dy=sizey; 
          g->Line(2*x, 2*y, 2*dx, 2*dy, colour);
          g->Line(2*x, 2*y+1, 2*dx, 2*dy, colour);
          g->Line(2*x+1, 2*y, 2*dx, 2*dy, colour);
          g->Line(2*x+1, 2*y+1, 2*dx, 2*dy, colour);
        }
    }
  }
}

double CellularPotts::getEnergy () {
  int x,y,xx,yy,i,s,sn;
  double E=0,k;

  for (x=1;x<sizex-1;x++){
    for (y=1;y<sizey-1;y++){
      for (i=1;i<=n_nb;i++){
        xx=x+nx[i];
	yy=y+ny[i];

        s=sigma[x][y];
	if (xx<=0 || xx>=sizex-1 || yy<=0 || y>=sizey-1) sn=-1;
	else sn=sigma[xx][yy];

        if (s==0) {}			//medium has no energy
	else if (s>0){
	/* contact energies */
	  if (sn==-1) {			//cell-border interface
		E+=par.border_energy;}
	  else if (sn==0 && s!=sn){	//cell-medium interface
		E+= (*cell)[s].EnergyDifference((*cell)[sn]);}
	  else if (sn!=0 && s!=sn){	//cell-cell interface
		  			//don't count twice --> 0.5 multiplier
		E+= 0.5* (*cell)[s].EnergyDifference((*cell)[sn]);}
        }
      }
    }
  }
 
  // area constraint
  s=par.target_area;
  for (i=1; i<par.n_init_cells; i++){
    k=(double)((*cell)[i].Area()-s);
    E+=(double) (par.lambda * k * k);
  }

  return E;
}


void CellularPotts::SyncCellTime (){
//signal the end of MCS to every cell and do the neccessary accounting:
//zero the displacement counter for cells at the end of each MCS
  
  for (int i=1;i<(int)(cell->size());i++){
    (*cell)[i].setDisplacement (0,0);
  }
};


void CellularPotts::PolUpdate (){
// Update polarization vectors for all cell once after every MCS
// CIL is handled here
//


  int i;
  double px,py,dx,dy,decay,rx,ry;
  double newpx[(int)(cell->size())], newpy[(int)(cell->size())];
  double CCMx=0, CCMy=0;

  //--- Polarisation vector update ---//
  for (i=1;i<(int)(cell->size());i++){
    newpx[i] = 0.0;
    newpy[i] = 0.0;
    if ((*cell)[i].AliveP()){
      (*cell)[i].CenterOfMass(&rx, &ry);
      CCMx+=rx;
      CCMy+=ry;
      (*cell)[i].getPolarisation (&px,&py);
      (*cell)[i].getDisplacement (&dx,&dy);
      
      newpx[i]=px;
      newpy[i]=py;
      
      
      (*cell)[i].CILstate = 0;
      if (  (*cell)[i].inContact && \
            (par.CILprob[(*cell)[i].tau] > RANDOM())  ){
  	// Cell is in contact and is subject to CIL
	
	(*cell)[i].CILstate = 1;
      
        //--- Decay ---//
        decay = (*cell)[i].getPolDecayContact();
        newpx[i] += (-1.0 * decay * px);
        newpy[i] += (-1.0 * decay * py);
        //---
       
         
      	//--- Increase with cell displacement ---//
        if (par.PolDynamicUpdateContact){
	  newpx[i] += (dx);
	  newpy[i] += (dy);
        }
        //---

	
        //--- Repolarisation away from contacting neighbours ---//
        if (par.repolarisationPM != 0){
	  double rnx, snx=0, rny, sny=0, norm;
	  int d=0, nn=0;

	  if (par.repolarisationDelay > 0) d=par.repolarisationDelay;

	  //--- Contact inhibition of locomotion
	  // ... 1. Get direction of CIL-neighbours:
	  for (int j=1;j<(int)(cell->size());j++){
  	    if ( (*cell)[i].contactDurationWith[j] > d ){
	      // i and j are touching cells, for the previous 'd' MCSteps
		
	      (*cell)[j].CenterOfMass(&rnx,&rny);
	      norm=sqrt( (rx-rnx)*(rx-rnx) + (ry-rny)*(ry-rny) );
	      if (norm > 0){
	        snx+=((rx-rnx)/norm);
	        sny+=((ry-rny)/norm);
	        nn += 1;
	      }
	    }
	  }
	
	  // ... 2. Apply repolarisation:
  	  norm=sqrt( (snx*snx) + (sny*sny) );
  	  if (norm > 0 && nn > 0 && norm > 0.2*nn){
	    snx=(snx/norm);
	    sny=(sny/norm);
	  } else {
	    // No CIL effect because the neighbours are symmetrically arranged around the cell
	    snx=0;
	    sny=0;
	  }
	  newpx[i] += par.repolarisationPM * snx;
	  newpy[i] += par.repolarisationPM * sny;
        } 
        //---
    
      } else {
      	// Cell does not undergo CIL
  
  	//--- Decay ---//
  	decay = (*cell)[i].getPolDecayFree();
  	newpx[i] += (-1.0 * decay * px);
  	newpy[i] += (-1.0 * decay * py);
  	//---
  	
  	
  	//--- Increase with cell displacement ---//
  	if (par.PolDynamicUpdateFree){
  	  newpx[i] += (dx);
  	  newpy[i] += (dy);
  	}
  	//---
	
      }
    }
  }
  
  CCMx /= (int)(cell->size());				// Displacement of the whole cell mass
  CCMy /= (int)(cell->size());				//
    
  
  for (i=1;i<(int)(cell->size());i++){
    (*cell)[i].setPolarisation (newpx[i],newpy[i]);	// Update polarisation vector from the immediately passed MCS
    (*cell)[i].chemotaxing=false;			// Reset chemotaxis indicator
  }
};


void CellularPotts::Displace (int k, int x, int y){
//record the elementary step (xp,yp)-->(x,y) for cell k
//increment the displacement bookkeeping of cell

  double dx,dy;
  double cdx,cdy;
  double cx1,cx2,cy1,cy2;


  if (sigma[x][y] == k) {       	//the retracting cell
    (*cell)[k].CenterOfMass(&cx1,&cy1);
    (*cell)[k].DecrementArea();
    (*cell)[k].RemoveSiteFromMoments(x,y);
    (*cell)[k].CenterOfMass(&cx2,&cy2);
    (*cell)[k].IncrementArea();
    (*cell)[k].AddSiteToMoments(x,y);
    dx=cx2-cx1;
    dy=cy2-cy1;
  } else {                      	//the expanding cell
    (*cell)[k].CenterOfMass(&cx1,&cy1);
    (*cell)[k].IncrementArea();
    (*cell)[k].AddSiteToMoments(x,y);
    (*cell)[k].CenterOfMass(&cx2,&cy2);
    (*cell)[k].DecrementArea();
    (*cell)[k].RemoveSiteFromMoments(x,y);
    dx=cx2-cx1;
    dy=cy2-cy1;
  };
  if (dx > (sizex*0.5)){ 
    do 
    dx-=sizex;
    while (dx > (sizex*0.5));
  }
  if (dx < -1*(sizex*0.5)){ 
    do 
    dx+=sizex;
    while (dx < -1*(sizex*0.5));
  }
  if (dy > (sizey*0.5)){ 
    do 
    dy-=sizey;
    while (dy > (sizey*0.5));
  }
  if (dy < -1*(sizey*0.5)){ 
    do 
    dy+=sizey;
    while (dy < -1*(sizey*0.5));
  }

  (*cell)[k].getDisplacement (&cdx,&cdy);
  (*cell)[k].setDisplacement (cdx+dx,cdy+dy);
};


void CellularPotts::setColours (){
  // Set cell colours
  int col;

  if (par.colour_cells == 0){
    // Update the cell colours only if they are not supposed to be fixed
    for (int c=1;c<(int)(cell->size());c++){
      col = (*cell)[c].tau+1;
  
      if (par.plotCILcells > 0 && (*cell)[c].CILstate){
        col = par.plotCILcells;
      } else if (par.plotChemotaxCells > 0 && (*cell)[c].chemotaxing) {
        col = par.plotChemotaxCells;
      } 
      
      (*cell)[c].SetColour(col);
    };
  };
};

void CellularPotts::SpecInitConditions (){
  
  vector<Cell>::iterator c=cell->begin(); ++c;
  for (;c!=cell->end();c++) {
    c->setPolarisation (0.0,0.0);
    c->setDisplacement (0.0,0.0);
    c->initCellParameters();
    c->initContacts();
  }
  
  setColours ();
  // set cells to different colors
  if (par.colour_cells > 0){
    int i,j,c;
    j=(int)(par.n_init_cells * par.colour_cells);
    for (i=1;i<j;i++){
      c = 3 + (i % 10);
      (*cell)[i].SetColour (c);
    };
  };
  
};


void CellularPotts::PrintPIF(){
  // Output in pixel image format, a format introduced by CompuCell3D
  int i,t,x,y,z=0;
  FILE *s;
  char *fname, *cellType;
  fname=(char*)calloc(100,sizeof(char));
  cellType=(char*)calloc(100,sizeof(char));
  sprintf(fname,"%s/cells_%05i.pif",par.datadir,par.time);
  s=fopen(fname,"wt");
  
  for (x=0;x<sizex;x++){  
    for (y=0;y<sizey;y++){
      i=sigma[x][y];
      if (i>0){
        t=(*cell)[i].tau;
	if (t==NCType) sprintf(cellType,"NC");
	if (t==PlacodeType) sprintf(cellType,"Placode");
        fprintf (s,"%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\n",i, cellType, x, x, y, y, z, z);
      }
    }
  }
  fclose (s);
  free(fname);
  free(cellType);
}

void CellularPotts::PrintPIF(const char* fname){
  // Output in pixel image format, a format introduced by CompuCell3D
  int i,t,x,y,z=0;
  FILE *s;
  char *cellType;
  cellType=(char*)calloc(100,sizeof(char));
  s=fopen(fname,"wt");
  
  for (x=0;x<sizex;x++){  
    for (y=0;y<sizey;y++){
      i=sigma[x][y];
      if (i>0){
        t=(*cell)[i].tau;
	if (t==NCType) sprintf(cellType,"NC");
	if (t==PlacodeType) sprintf(cellType,"Placode");
        fprintf (s,"%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\n",i, cellType, x, x, y, y, z, z);
      }
    }
  }
  fclose (s);
  free(cellType);
}


void CellularPotts::PrintAsParticles(){
  double x,y,vx,vy,px,py;
  FILE *s;
  char *fname;
  fname=(char*)calloc(5000,sizeof(char));
  sprintf(fname,"%s/particles.dat",par.datadir);
  if (par.time == 0) s=fopen(fname,"wt");
  else s=fopen(fname,"at");
  for (int i=1;i<(int)(cell->size());i++){
    if ( (*cell)[i].AliveP() ){
      (*cell)[i].CenterOfMass(&x, &y);
      (*cell)[i].getDisplacement(&vx, &vy);
      (*cell)[i].getPolarisation(&px, &py);
      // Print format:
      // time ID type colour x y vx vy px py
      fprintf(s,"%i\t%i\t%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", par.time, i, (*cell)[i].tau, (*cell)[i].CILstate, x, y, vx, vy, px, py );
    }
  }
  fprintf(s,"\n\n");
  fclose (s);
  free(fname);
}

void CellularPotts::PrintCM (){

  int x,fc;
  double cx,cy;
  FILE *s;
  char *fname;
  fname=(char*)calloc(50,sizeof(char));
  sprintf(fname,"%s/cm_%05i.out",par.datadir,par.time);
  s=fopen(fname,"wt");
  if (par.print_CMarea == true)
    fprintf(s,"#%i %i\n#t\tID\tx\t\ty\t\tstate\tarea\n",sizex,sizey);
  else 
    fprintf(s,"#%i %i\n#t\tID\tx\t\ty\tstate\n",sizex,sizey);
  
  for (x=0;x<(par.n_init_cells+1);x++){  
    (*cell)[x].CenterOfMass(&cx,&cy);
    fc=1;
    if (cx>0 && cy>0){ //if cell still exists, print coordinates
      if (par.print_CMarea == true){
        fprintf (s,"%i\t%i\t%lf\t%lf\t%i\t%i\n",par.time,x,cx,cy,fc,
					(*cell)[x].Area());
      }else{
        fprintf (s,"%i\t%i\t%lf\t%lf\t%i\n",par.time,x,cx,cy,fc);
      }
    }
  };
  fprintf(s,"\n");
  fclose (s);
  free(fname);
}


void CellularPotts::PrintSigma (){

  int x,y;
  FILE *s;
  char *fname;
	        
  fname=(char*)calloc(50,sizeof(char));
  sprintf(fname,"%s/sigma_%05i.out",par.datadir,par.time);  
  s=fopen(fname,"wt");
  fprintf(s,"#%i %i\n",sizex,sizey);
  for (y=0;y<sizey;y++){
    for (x=0;x<sizex;x++){
      fprintf (s,"%i ", sigma[x][y]);
    };
    fprintf(s,"\n");
  };
  fprintf(s,"\n");
  fclose (s);
  free(fname);
};

void CellularPotts::PrintAreas(){
  FILE *fp;
  char *fn;
  int i,a=0,aa=0,n=0,dummy;
  double r,sd;

  fn=(char*)calloc(50,sizeof(char));
  sprintf(fn,"%s/areas.out",par.datadir);  
  fp=fopen(fn,"at");
  for (i=1; i<par.n_init_cells+1; i++){
    dummy = (*cell)[i].Area();
    a  += dummy;
    aa += dummy*dummy;
    n++;
  };

  if (n>0){ 
    r = (double)(a) /( (double)(par.target_area) * (double)(n) );
    sd = sqrt((double)(aa)/(double)(n) - DSQR((double)(a)/(double)(n)));
    fprintf (fp,"%i %i %i %lf %lf\n",par.time,a,n,r,sd);
  };


  fclose (fp);
  free (fn);
};



void CellularPotts::PlotSigma(Graphics *g, int mag) {
/** A simple method to plot all sigma's in window
    without the black lines */
  
  for (int x=0;x<sizex;x++) 
    for (int y=0;y<sizey;y++) {
      for (int xm=0;xm<mag;xm++)
	for (int ym=0;ym<mag;ym++)
      g->Point( sigma[x][y], mag*x+xm, mag*y+ym);
  }
  
}

void CellularPotts::PlotNC(Graphics *g, int mag, int colour) {
/** A simple method to plot NCType cells only
    without the black lines */
  int s;

  for (int x=0;x<sizex;x++){
    for (int y=0;y<sizey;y++){
      s = sigma[x][y];
      if (s>MEDIUM){
        if ((*cell)[s].tau == NCType){
          for (int xm=0;xm<mag;xm++){
	    for (int ym=0;ym<mag;ym++){
	      g->Point( colour, mag*x+xm, mag*y+ym);
	    }
	  }
	}
      }
    }
  }
}

void CellularPotts::PlotPL(Graphics *g, int mag, int colour) {
/** A simple method to plot PlacodeType cells only
    without the black lines */
  int s;

  for (int x=0;x<sizex;x++){
    for (int y=0;y<sizey;y++){
      s = sigma[x][y];
      if (s>MEDIUM){
        if ((*cell)[s].tau == PlacodeType){
          for (int xm=0;xm<mag;xm++){
	    for (int ym=0;ym<mag;ym++){
	      g->Point( colour, mag*x+xm, mag*y+ym);
	    }
	  }
	}
      }
    }
  }
}

int **CellularPotts::SearchNandPlot(Graphics *g, bool get_neighbours)
{
  int i, j,q;
  int **neighbours=0;
  
  
  /* Allocate neighbour matrix */
  if (get_neighbours) {
    neighbours=(int **)malloc(((int)(cell->size())+1)*sizeof(int *));
    if (neighbours==NULL) 
      MemoryWarning();
    
    neighbours[0]=(int *)malloc(((int)(cell->size())+1)*((int)(cell->size())+1)*sizeof(int));
    if (neighbours[0]==NULL)
      MemoryWarning();
   
    for (i=1;i<(int)(cell->size())+1;i++)
      neighbours[i]=neighbours[i-1]+((int)(cell->size())+1);
    
    /* Clear this matrix */
    for (i=0;i<((int)(cell->size())+1)*((int)(cell->size())+1);i++)
      neighbours[0][i]=EMPTY;  
  }

  for ( i = 0; i < sizex-1; i++ )
    for ( j = 0; j < sizey-1; j++ ) {
      

      int colour;
      if (sigma[i][j]<=0) {
	colour=0;
      } else {
	colour = (*cell)[sigma[i][j]].Colour();
	//colour = sigma[i][j];
      }
      
      if (g && sigma[i][j]>0)  /* if draw */ 
        g->Point( colour, 2*i, 2*j);
      
      if ( sigma[i][j] != sigma[i+1][j] )  /* if cellborder */ /* etc. etc. */
	{
	  if (g) 
	    g->Point( BLACK, 2*i+1, 2*j );
	  if (get_neighbours) {
	    if (sigma[i][j]>0) {
	      for (q=0;q<(int)(cell->size());q++)
		if (neighbours[sigma[i][j]][q]==EMPTY) { 
		  neighbours[sigma[i][j]][q]=sigma[i+1][j];  
		  break;
		}
		else
		  if (neighbours[sigma[i][j]][q]==sigma[i+1][j]) 
		    break;
	    }
	    if (sigma[i+1][j]>0) {
	      for (q=0;q<(int)(cell->size());q++)
		if (neighbours[sigma[i+1][j]][q]==EMPTY) { 
		  neighbours[sigma[i+1][j]][q]=sigma[i][j]; 
		  break;
		}
		else
		  if (neighbours[sigma[i+1][j]][q]==sigma[i][j]) 
		    break;
	    }
	  }
	} 
      else
        if (g && sigma[i][j]>0) 
          g->Point( colour, 2*i+1, 2*j );
      
      
      if ( sigma[i][j] != sigma[i][j+1] ) {
	
        if (g) 
	  g->Point( BLACK, 2*i, 2*j+1 );
	
	if (get_neighbours) {
	  if (sigma[i][j]>0) {
	    for (q=0;q<(int)(cell->size());q++)
	      if (neighbours[sigma[i][j]][q]==EMPTY) { 
		neighbours[sigma[i][j]][q]=sigma[i][j+1];  
		break; 
	      }
	      else
		if (neighbours[sigma[i][j]][q]==sigma[i][j+1]) 
		  break;
	  }
	  
	  if (sigma[i][j+1]>0) {
	    
	    for (q=0;q<(int)(cell->size());q++)
	      if (neighbours[sigma[i][j+1]][q]==EMPTY) { 
		neighbours[sigma[i][j+1]][q]=sigma[i][j]; 
		break;
	      }
	      else
		if (neighbours[sigma[i][j+1]][q]==sigma[i][j]) 
		  break;
	  }
	}
      } 
      else
        if (g && sigma[i][j]>0) 
          g->Point( colour, 2*i, 2*j+1 );
      
      /* Cells that touch eachother's corners are NO neighbours */ 
      
      if (sigma[i][j]!=sigma[i+1][j+1] 
	  || sigma[i+1][j]!=sigma[i][j+1] ) { 
        if (g) 
          g->Point( BLACK, 2*i+1, 2*j+1 ); 
      }
      else
        if (g && sigma[i][j]>0) 
          g->Point( colour, 2*i+1, 2*j+1 );
    }
  
  if (get_neighbours)
    return neighbours;
  else 
    return 0;

}


void CellularPotts::ConstructInitCells (Dish &beast) {
  
  // construct enough cells for the zygote.  "cells", contains the
  // number of colours (excluding background).
  { for (int i=0; i<par.n_init_cells; i++) {
      cell->push_back(Cell(beast));
    }}
  
  // Set the area and target area of the cell
  // makes use of the pointer to the Cell pointer of Dish
  // which is a member of CellularPotts 
  MeasureCellSizes();
  
}

void CellularPotts::MeasureCellSizes(void) {
  
  // Clean areas of all cells, including medium
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->SetTargetArea(0);
    c->area = 0;
  }
  
  // calculate the area of the cells
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	(*cell)[sigma[x][y]].IncrementTargetArea();
	(*cell)[sigma[x][y]].IncrementArea();
	(*cell)[sigma[x][y]].AddSiteToMoments(x,y);

      }
    }
  }
  
  // set the actual area to the target area
  {
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->SetAreaToTarget();

  }
  }
}

void CellularPotts::MeasureCellSize(Cell &c) {
  
  c.CleanMoments();
  
  // calculate the area of the cell
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y] == c.sigma) {
	(*cell)[sigma[x][y]].IncrementTargetArea();
	(*cell)[sigma[x][y]].IncrementArea();
	(*cell)[sigma[x][y]].AddSiteToMoments(x,y);

      }
    }
  }
  
//   // set the actual area to the target area
//   {
//   for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
//     c->SetAreaToTarget();

//   }

}

Dir *CellularPotts::FindCellDirections(void) const
{ 
  
  double *sumx=0,*sumy=0;
  double *sumxx=0,*sumxy=0,*sumyy=0;
  double *n=0;  

  double xmean=0,ymean=0,sxx=0,sxy=0,syy=0;
  double D,lb1=0,lb2=0;

  Dir *celldir;

  /* Allocation of sufficient memory space */
  if( (sumx= (double *)malloc((int)(cell->size())*sizeof(double)))==NULL)
    MemoryWarning(); 
  else
    if( (sumy= (double *)malloc((int)(cell->size())*sizeof(double)))==NULL) 
      MemoryWarning();
    else
      if ((sumxx=(double *)malloc((int)(cell->size())*sizeof(double)))==NULL) 
	MemoryWarning();
      else
	if((sumxy=(double *)malloc((int)(cell->size())*sizeof(double)))==NULL) 
	  MemoryWarning();
	else
	  if((sumyy=(double *)malloc((int)(cell->size())*sizeof(double)))==NULL)
	    MemoryWarning();
	  else
	    if((n=(double *)malloc((int)(cell->size())*sizeof(double)))==NULL) 
	      MemoryWarning();
  
  
  if ( !(celldir=new Dir[(int)(cell->size())]) )
    MemoryWarning();

  	
  /* Initialization of the variables */
   
  for (int i=0;i<(int)(cell->size());i++) {
    
    sumx[i]=0.;
    sumy[i]=0.;
    sumxx[i]=0.;
    sumxy[i]=0.;
    sumyy[i]=0.;
    n[i]=0L;

  }


  /* Find sumx, sumy, sumxx and sumxy for all cells */

  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++) 
      if (sigma[x][y]>0) {
	sumx[0]+=(double)x;
	sumy[0]+=(double)y;
	sumxx[0]+=(double)x*x;
	sumxy[0]+=(double)x*y;
	sumyy[0]+=(double)y*y;
	
	n[0]++;
	
	sumx[sigma[x][y]]+=(double)x;
	sumy[sigma[x][y]]+=(double)y;
	
	sumxx[sigma[x][y]]+=(double)x*x;
	sumxy[sigma[x][y]]+=(double)x*y;
	sumyy[sigma[x][y]]+=(double)y*y;
	
	n[sigma[x][y]]++;
	
      }
  
  /* Compute the principal axes for all cells */
  
  {
    for (int i=0;i<(int)(cell->size());i++) {
    
      if (n[i]>10) {
      
	xmean=((double)sumx[i])/((double)n[i]);
	ymean=((double)sumy[i])/((double)n[i]);

	sxx=(double)(sumxx[i])-((double)(sumx[i]*sumx[i]))/(double)n[i];
	sxx=sxx/(double)(n[i]-1);

	sxy=(double)(sumxy[i])-((double)(sumx[i]*sumy[i]))/(double)n[i];
	sxy=sxy/(double)(n[i]-1);

	syy=(double)(sumyy[i])-((double)(sumy[i]*sumy[i]))/(double)n[i];
	syy=syy/(double)(n[i]-1);

	D=sqrt( (sxx+syy)*(sxx+syy)-4.*(sxx*syy-sxy*sxy) );
	lb1=(sxx+syy+D)/2.;lb2=(sxx+syy-D)/2.;
	celldir[i].lb1=lb1; celldir[i].lb2=lb2; 
      }
      if (sxy==0.0)
	celldir[i].bb1=1.; 
      else
	celldir[i].bb1=sxy/(lb1-syy);
    
      if (fabs(celldir[i].bb1)<.00001) {
	if (celldir[i].bb1>0.) 
	  celldir[i].bb1=.00001;
	else 
	  celldir[i].bb1=-.00001;
      }
 
      celldir[i].aa1=ymean-xmean*celldir[i].bb1;
      celldir[i].bb2= (-1.)/celldir[i].bb1;
    
      celldir[i].aa2=ymean-celldir[i].bb2*xmean;     
    }
		  
  }

  /* bevrijd gealloceerd geheugen */
  free(sumx);
  free(sumy);
  free(sumxx);
  free(sumxy);
  free(sumyy);
  free(n);

  return celldir;
 
}

void CellularPotts::ShowDirections(Graphics &g, const Dir *celldir) const
{
  int i;
  
  if ((int)(cell->size())>1) 
    for (i=1;i<(int)(cell->size());i++)
      g.Line(0,(int)(2*celldir[i].aa1),sizex*2,(int)((celldir[i].aa1+celldir[i].bb1*sizey)*2),2);
  
}

void CellularPotts::DivideCells(vector<bool> which_cells)
{
  
  // for the cell directions
  Dir *celldir=0;
  
  /* Allocate space for divisionflags */
  int *divflags=(int *)malloc(((int)(cell->size())*2+5)*sizeof(int));
  
  /* Clear divisionflags */
  for (int i=0;i<(int)((int)(cell->size())*2+5);i++) 
    divflags[i]=0;
  
  
  if ( !((int)(which_cells.size())==0 || (int)(which_cells.size())>=(int)(cell->size())) ) {
    throw "In CellularPotts::DivideCells, Too few elements in vector<int> which_cells.";
  }
  
  /* division */
  {for (int i=0;i<sizex;i++)
      for (int j=0;j<sizey;j++) 
	if (sigma[i][j]>0) // i.e. not medium and not border state (-1)
	  { 
      
      
	    // Pointer to mother. Warning: Renew pointer after a new
	    // cell is added (push_back). Then, the array *cell is relocated and
	    // the pointer will be lost...
      
	    Cell *mother=&((*cell)[sigma[i][j]]);
	    Cell *daughter;
      
	    /* Divide if NOT medium and if DIV bit set or divide_always is set */
	    // if which_cells is given, divide only if the cell
	    // is marked in which_cells.
	    if  ( !which_cells.size() || which_cells[mother->sigma] )    {
	
	      if (!(divflags[ mother->Sigma() ]) ) {
	  
		// add daughter cell, copying states of mother
		daughter=new Cell(*mother);
		cell->push_back(*daughter);
	  
		// renew pointer to mother
		mother=&((*cell)[sigma[i][j]]);

		divflags[ mother->Sigma() ]=daughter->Sigma();
		delete daughter;
	  
		// array may be relocated after "push_back"
	  
		// renew daughter pointers
		daughter=&(cell->back());
	  
		/* administration on the onset of mitosis */
	  
		/* Ancestry is taken care of in copy constructor of Cell 
		   see cell.hh: Cell(const Cell &src, bool newcellP=false) : Cytoplasm(src) {} */
	  
		/* inherit  polarity of mother */
		// All that needs to be copied is copied in the copy constructor
		// of Cell and in the default copy constr. of its base class Cytoplasm
		// note: also the celltype is inherited
	  


	      } else {
		daughter=&((*cell)[ divflags[mother->Sigma()] ]);
	      }
	    
	    
	      /* Now the actual division takes place */
	    
	      /* If celldirections where not yet computed: do it now */
	      if (!celldir) 
		celldir=FindCellDirections();
	
	      /* if site is below the minor axis of the cell: sigma of new cell */
	      if (j>((int)(celldir[mother->sigma].aa2+
			   celldir[mother->sigma].bb2*(double)i))) { 

		mother->DecrementArea();
		mother->DecrementTargetArea();
		mother->RemoveSiteFromMoments(i,j);
		sigma[i][j]=daughter->Sigma();
		daughter->AddSiteToMoments(i,j);
		daughter->IncrementArea();
		daughter->IncrementTargetArea();
	  
	      } 


	    }
      
	  }
  }  
  if (celldir) 
    delete[] (celldir);
  
  if (divflags)
    free(divflags);
}        


/**! Fill the plane with initial cells 
 \return actual amount of cells (some are not draw due to overlap) */
int CellularPotts::ThrowInCells(int n,int cellsize) {
  
  //  int gapx=(sizex-nx*cellsize)/(nx+1);
  //int gapy=(sizey-ny*cellsize)/(ny+1);
  
  int cellnum=1;

  for (int i=0;i<n;i++) {
    
    // draw a circle at x0, y0
    int x0=RandomNumber(sizex);
    int y0=RandomNumber(sizey);
   
    bool overlap=false;
    
    // check overlap
    for (int x=0;x<cellsize;x++)
      for (int y=0;y<cellsize;y++)
	if ( ( 
	      ( (x-cellsize/2)*(x-cellsize/2)+(y-cellsize/2)*(y-cellsize/2) )<
	      ( (cellsize/2)*(cellsize/2))) &&
	     ( x0+x<sizex && y0+y<sizey ) )
	  if (sigma[x0+x][y0+y]) {
	    overlap=true;
	    break;
	  }
    
    if (!overlap) {
      for (int x=0;x<cellsize;x++)
	for (int y=0;y<cellsize;y++)
	  if ( ( 
		( (x-cellsize/2)*(x-cellsize/2)+(y-cellsize/2)*(y-cellsize/2) )<
		( (cellsize/2)*(cellsize/2))) &&
	       ( x0+x<sizex && y0+y<sizey ) )
	    sigma[x0+x][y0+y]=cellnum;
      
      cellnum++;
    }
  }
  cerr << "[ cellnum = " << cellnum << "]";

  // repair borders
  // fill borders with special border state
  for (int x=0;x<sizex-1;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey-1;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }

  {for (int x=1;x<sizex-2;x++) {
      sigma[x][1]=0;
      sigma[x][sizey-2]=0;
    }}
  {for (int y=1;y<sizey-2;y++) {
      sigma[1][y]=0;
      sigma[sizex-2][y]=0;
    }}
  return cellnum;
} 

  
int CellularPotts::GrowInCells(int n_cells, int cell_size, double subfield) {

  
  int sx = (int)((sizex-2)/subfield);
  int sy = (int)((sizey-2)/subfield);
  
  int offset_x = (sizex-2-sx)/2;
  int offset_y = (sizey-2-sy)/2;

  if (n_cells==1) {
    return GrowInCells(1, cell_size, sizex/2, sizey/2, 0, 0);
  } else {
    return GrowInCells(n_cells, cell_size, sx, sy, offset_x, offset_y);
  }
}

void CellularPotts::EdenGrowth(int cell_size, int cellnum, int n_cells, int** new_sigma){

  int cs = (int)(log((double)(cell_size)) / log(2.));
  int clock=1;
  for (int i=0;i<cs+1;i++) {
    for (int x=1;x<sizex-1;x++){
      for (int y=1;y<sizey-1;y++) {

        if (sigma[x][y]==0) {
          // take a random neighbour
          int xyp= clock; //(int)(8*RANDOM()+1);
          int xp = nx[xyp]+x;
          int yp = ny[xyp]+y;
          int kp;
          // NB removing this border test yields interesting effects :-)
          // You get a ragged border, which you may like!
          if ((kp=sigma[xp][yp])!=-1){
	    if (kp>(cellnum-n_cells))
	      new_sigma[x][y]=kp;
   	    else
	      new_sigma[x][y]=0;
          }else
	    new_sigma[x][y]=0;
  
        } else {
          new_sigma[x][y]=sigma[x][y];
        };
      };
    };
    if (clock > n_nb - 1) clock = 1;
    else clock++;
    // copy sigma to new_sigma, but do not touch the border!
    for (int x=1;x<sizex-1;x++) {
      for (int y=1;y<sizey-1;y++) {
        sigma[x][y]=new_sigma[x][y];
      };
    };
  };
}


int CellularPotts::GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y) {
  
  // make initial cells using Eden Growth
  
  int **new_sigma=(int **)malloc(sizex*sizeof(int *));
  if (new_sigma==NULL)
    MemoryWarning();
  
  new_sigma[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (new_sigma[0]==NULL)  
    MemoryWarning();
  
  for (int i=1;i<sizex;i++) 
    new_sigma[i]=new_sigma[i-1]+sizey;
  
  /* Clear CA plane */
  { for (int i=0;i<sizex*sizey;i++) 
     new_sigma[0][i]=0; 
  }

  int cellnum=(int)(cell->size())-1;
  if (par.init_cell_distribution == 0){
    // scatter initial points, or place a cell in the middle 
    // if only one cell is desired
    if (n_cells>1) {    
      for (int i=0;i<n_cells;i++) {
      
        sigma[RandomNumber(sx)+offset_x][RandomNumber(sy)+offset_y]=++cellnum;
      
      }
    } else {
      sigma[sx][sy]=++cellnum;
    };
  }else if (par.init_cell_distribution == 1){
    //set cells in middle in a cluster
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int il = (int)(sqrt(n_cells) );
    double spacing = 1.5*sqrt((double)(par.target_area)/PI);
    double step;

    for (j=0;j<=il;j++){
      for (i=0;i<il;i++) {
        if (k<n_cells){
	  step = spacing * (double) ( i-((double)(il-1)/2) );
	  xx=(int)((double)(sizex/2-1) + step);
	  step = spacing * (double) ( j-((double)(il-1)/2) );
          yy=(int)((double)(sizey/2-1) + step);
	  sigma[xx][yy]=++cellnum;
	  //for (int xxx=xx; xxx<xx+step; xxx++)
	  //  for (int yyy=yy; yyy<yy+step; yyy++)
	  //    sigma[xxx][yyy]=sigma[xx][yy];
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 2){
    //set cells in middle in a cluster with 
    //size_init_cells sized cells
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int xmin, xmax, ymin, ymax;
    int il = (int)(sqrt(n_cells) );
    double spacing = 1.5*sqrt((double)(par.size_init_cells)/PI);
    double stepx, stepy;

    for (i=0;i<=il;i++) {
      for (j=0;j<il;j++){
        if (k<n_cells){
	  stepx = spacing * (double) ( i-((double)(il-1)/2) );
	  stepy = spacing * (double) ( j-((double)(il-1)/2) );
	  xmin=(int)((double)(sizex/2-1) + stepx - (spacing/2));
	  xmax=(int)((double)(sizex/2-1) + stepx + (spacing/2));
	  ymin=(int)((double)(sizey/2-1) + stepy - (spacing/2));
	  ymax=(int)((double)(sizey/2-1) + stepy + (spacing/2));
	  cellnum++;
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      sigma[xx][yy]=cellnum;
	    }
	  }
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 3){
    //set cells in center in a disk with 
    //size_init_cells sized cells
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int xmin, xmax, ymin, ymax;
    int il = (int)(2*sqrt(n_cells/3) +1);
    double spacing = 1.5*sqrt((double)(par.size_init_cells)/PI);
    double stepx, stepy, rr;

    for (i=0;i<=il;i++) {
      for (j=0;j<il;j++){
        if (k<n_cells){
	  stepx = spacing * (double) ( i-((double)(il-1)/2) );
	  stepy = spacing * (double) ( j-((double)(il-1)/2) );
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt(n_cells/3)+1)){
            xmin=(int)((double)(sizex/2-1) + stepx - (spacing/2));
	    xmax=(int)((double)(sizex/2-1) + stepx + (spacing/2));
	    ymin=(int)((double)(sizey/2-1) + stepy - (spacing/2));
	    ymax=(int)((double)(sizey/2-1) + stepy + (spacing/2));
	    cellnum++;
	    for (xx=xmin; xx<xmax; xx++){
	      for (yy=ymin;yy<ymax;yy++){
		sigma[xx][yy]=cellnum;
	      }
	    }
	    k++;
	  }
	}
      }
    }
  } else if (par.init_cell_distribution == 4){
    //fill cells from top, downwards and left
    
    double spacing = sqrt((double)(par.size_init_cells));
    double stepx, stepy;
    //int cellnum=(int)(cell->size())-1;
    int cellnum = 1;
    int xx,yy,i,j, ync;
    int xmin, xmax, ymin, ymax;
    int imx = (int)(sizex / (spacing + 0.1));
    int imy = (int)(sizey / (spacing + 0.1));
    int cs = (int)(spacing);
    bool addedCell=false;
  
    ync = -1; ymax = -1;
    if (cs == 0) cs = 1;
    for (j=1;j<imy;j++){
      for (i=1;i<imx;i++) {
        if (cellnum <= par.n_NC){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + cs);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + cs);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {
	  cellnum++; addedCell=false;}
	} else if (ync == -1){ ync = ymax; }
      }
    }
    if (ync == -1) ync = ymax;
    
    // Placodes
    if (par.init_PL_distance >= 0){
      spacing = par.init_PL_distance;
    } else {
      spacing = sqrt((double)(par.size_init_cells));
    }
    imx = (int)(sizex / (spacing + 0.1));
    imy = (int)((sizey-ync) / (spacing + 0.1));

    for (j=0;j<imy;j++){
      for (i=1;i<imx;i++) {
        if (cellnum <= n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + cs);
	  ymin=(int)(double)(stepy + ync);
	  ymax=(int)(double)(stepy + ync + cs);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {
	  cellnum++; addedCell=false;}
	}
      }
    }
    if (cellnum <= n_cells) {
      std::cerr << "Warning! Not enough space to initialise all cells! Starting with " << (cellnum-par.n_NC) << " instead of " << (n_cells-par.n_NC) << ".\n";
      cellnum = n_cells;
    }
  } else if (par.init_cell_distribution == 5){
    //fill cells from top, downwards and left
    // leave padding at sides
    
    double spacing = sqrt((double)(par.size_init_cells));
    double stepx, stepy;
    //int cellnum=(int)(cell->size())-1;
    int cellnum = 1;
    int xx,yy,i,j,k=0;
    int xmin, xmax, ymin, ymax;
    int imx = (int)(sizex / (spacing + 0.1));
    int imy = (int)(sizey / (spacing + 0.1));
    bool addedCell=false;
  
    for (j=1;j<imy;j++){
      for (i=1;i<imx;i++) {
        if (k<par.n_NC){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + spacing);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + spacing);
	  if (xmin>par.sidePadding && xmax<(sizex-par.sidePadding)){
	    for (xx=xmin; xx<xmax; xx++){
	      for (yy=ymin;yy<ymax;yy++){
	        if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	      }
	    }
	    if (addedCell) {cellnum++; addedCell=false;}
	    k++;
	  }
	}else if (k<n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + spacing);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + spacing);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {cellnum++; addedCell=false;}
	  k++;
	}
      }
    }
  }else if (par.init_cell_distribution == 6){
    //set cells in middle in two clusters, one for NC and one for Placode
    int cellnum=0;
    int xx,yy,i,j,k=0;
    int iln = (int)(sqrt(par.n_NC));
    int ilp = (int)(sqrt(par.n_init_cells - par.n_NC));
    int clusterDistance=par.init_cluster_distance;
    double spacing = 1.5*sqrt((double)(par.target_area)/PI);
    double step;
    int cxn=(int)(sizex/2-1 - (iln/2)*spacing - clusterDistance);
    int cyn=(int)(sizey/2-1);
    int cxp=(int)(sizex/2-1 + (ilp/2)*spacing + clusterDistance);
    int cyp=(int)(sizey/2-1);

    for (j=0;j<=iln;j++){
      for (i=0;i<iln;i++) {
        if (k<par.n_NC){
	  step = spacing * (double) ( i-((double)(iln-1)/2) );
	  xx=(int)((double)(cxn) + step);
	  step = spacing * (double) ( j-((double)(iln-1)/2) );
          yy=(int)((double)(cyn) + step);
	  sigma[xx][yy]=++cellnum;
	  k++;
	}
      }
    }
    for (j=0;j<=ilp;j++){
      for (i=0;i<ilp;i++) {
	if (k<par.n_init_cells){
	  step = spacing * (double) ( i-((double)(ilp-1)/2) );
	  xx=(int)((double)(cxp) + step);
	  step = spacing * (double) ( j-((double)(ilp-1)/2) );
          yy=(int)((double)(cyp) + step);
	  sigma[xx][yy]=++cellnum;
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 7){
    //set cells in center two disks with 
    //size_init_cells sized cells
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int ila = (int)(2*sqrt(par.n_init_cells/(2.*PI)) +1);
    int ilb = (int)(2*sqrt(par.n_init_cells/(2.*PI)) +1);
    double spacing = 1.5*sqrt((double)(par.size_init_cells)/PI);
    double stepx, stepy, rr;
    int clusterDistance=par.init_cluster_distance;
    int cxa=(int)(sizex/2-1 - (ila/2)*spacing - clusterDistance/2.);
    int cya=(int)(sizey/2-1);
    int cxb=(int)(sizex/2-1 + (ilb/2)*spacing + clusterDistance/2.);
    int cyb=(int)(sizey/2-1);

    for (j=0;j<=ila;j++){
      for (i=0;i<ila;i++) {
        if (k<par.n_init_cells/2){
	  stepx = spacing * (double) ( i-((double)(ila-1)/2) );
	  xx=(int)((double)(cxa) + stepx);
	  stepy = spacing * (double) ( j-((double)(ila-1)/2) );
          yy=(int)((double)(cya) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt(par.n_init_cells/(2.*PI))+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
    for (j=0;j<=ilb;j++){
      for (i=0;i<ilb;i++) {
	if (k<par.n_init_cells){
	  stepx = spacing * (double) ( i-((double)(ilb-1)/2) );
	  xx=(int)((double)(cxb) + stepx);
	  stepy = spacing * (double) ( j-((double)(ilb-1)/2) );
          yy=(int)((double)(cyb) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt(par.n_init_cells/(2.*PI))+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
  } else if (par.init_cell_distribution == 8){
    //set cells in center two disks with 
    //size_init_cells sized cells
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int iln = (int)(2*sqrt(par.n_NC/PI) +1);
    int ilp = (int)(2*sqrt((par.n_init_cells - par.n_NC)/PI) +1);
    double spacing = 1.5*sqrt((double)(par.size_init_cells)/PI);
    double stepx, stepy, rr;
    int clusterDistance=par.init_cluster_distance;
    int cxn=(int)(sizex/2-1 - (iln/2)*spacing - clusterDistance/2.);
    int cyn=(int)(sizey/2-1);
    int cxp=(int)(sizex/2-1 + (ilp/2)*spacing + clusterDistance/2.);
    int cyp=(int)(sizey/2-1);

    for (j=0;j<=iln;j++){
      for (i=0;i<iln;i++) {
        if (k<par.n_NC){
	  stepx = spacing * (double) ( i-((double)(iln-1)/2) );
	  xx=(int)((double)(cxn) + stepx);
	  stepy = spacing * (double) ( j-((double)(iln-1)/2) );
          yy=(int)((double)(cyn) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt(par.n_NC/PI)+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
    for (j=0;j<=ilp;j++){
      for (i=0;i<ilp;i++) {
	if (k<par.n_init_cells){
	  stepx = spacing * (double) ( i-((double)(ilp-1)/2) );
	  xx=(int)((double)(cxp) + stepx);
	  stepy = spacing * (double) ( j-((double)(ilp-1)/2) );
          yy=(int)((double)(cyp) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt((n_cells-par.n_NC)/PI)+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
  } else if (par.init_cell_distribution == 9){
    // set Placode cells in center in a disk and two NC cell clusters on both sides,
    // with size_init_cells sized cells
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int iln = (int)(2*sqrt(par.n_NC/(2*PI)) +1);
    int ilp = (int)(2*sqrt((par.n_init_cells - par.n_NC)/PI) +1);
    double spacing = 1.5*sqrt((double)(par.size_init_cells)/PI);
    double stepx, stepy, rr;
    int clusterDistance=par.init_cluster_distance;
    int cxp=(int)(sizex/2-1);	// Placodes go in the centre
    int cyp=(int)(sizey/2-1);
    int cxn1=(int)(sizex/2-1 - (ilp/2+iln/2)*spacing - clusterDistance);
    int cyn1=(int)(sizey/2-1);
    int cxn2=(int)(sizex/2-1 + (ilp/2+iln/2)*spacing + clusterDistance);
    int cyn2=(int)(sizey/2-1);

    // NC cluster 1
    for (j=0;j<=iln;j++){
      for (i=0;i<iln;i++) {
        if (k<(par.n_NC/2)){
	  stepx = spacing * (double) ( i-((double)(iln-1)/2) );
	  xx=(int)((double)(cxn1) + stepx);
	  stepy = spacing * (double) ( j-((double)(iln-1)/2) );
          yy=(int)((double)(cyn1) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt(par.n_NC/(2*PI))+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
    // NC cluster 2
    for (j=0;j<=iln;j++){
      for (i=0;i<iln;i++) {
        if (k<par.n_NC){
	  stepx = spacing * (double) ( i-((double)(iln-1)/2) );
	  xx=(int)((double)(cxn2) + stepx);
	  stepy = spacing * (double) ( j-((double)(iln-1)/2) );
          yy=(int)((double)(cyn2) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt(par.n_NC/(2*PI))+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
    // Placode cluster
    for (j=0;j<=ilp;j++){
      for (i=0;i<ilp;i++) {
	if (k<par.n_init_cells){
	  stepx = spacing * (double) ( i-((double)(ilp-1)/2) );
	  xx=(int)((double)(cxp) + stepx);
	  stepy = spacing * (double) ( j-((double)(ilp-1)/2) );
          yy=(int)((double)(cyp) + stepy);
	  rr=sqrt( stepx*stepx + stepy*stepy );
	  if (rr < (spacing*sqrt((n_cells-par.n_NC)/PI)+1)){
	    cellnum++;
	    sigma[xx][yy]=cellnum;
	    k++;
	  }
	}
      }
    }
  } else if (par.init_cell_distribution == 10){
    //set cells on the left side of the system 
    
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int xmin, xmax, ymin, ymax;
    int il = (int)(sizey/ (sqrt(par.size_init_cells)));
    double spacing = sqrt((double)(par.size_init_cells));
    double stepx, stepy;
  
    for (i=0;i<=il;i++) {
      for (j=0;j<il;j++){
        if (k<n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + spacing);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + spacing);
	  cellnum++;
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) sigma[xx][yy]=cellnum;
	    }
	  }
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 11){
    //fill cells from middle up, downwards and left
    //useful for piston testing configuration
    
    int cellnum=0;
    //int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int xmin, xmax, ymin, ymax;
    int il = (int)(sizey/ (sqrt(par.size_init_cells)));
    double spacing = sqrt((double)(par.size_init_cells));
    double stepx, stepy;
  
    for (i=il/2;i>0;i--) {
      for (j=0;j<il;j++){
        if (k<n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + spacing);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + spacing);
	  cellnum++;
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) sigma[xx][yy]=cellnum;
	    }
	  }
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 12){
    //set cells on the left side of the system 
    //same as above, but with slightly increased spacing between cells
    //to avoid initial pressure
    
    int cellnum=(int)(cell->size())-1;
    int xx,yy,i,j,k=0;
    int xmin, xmax, ymin, ymax;
    int il = (int)(sizey/ (sqrt(par.size_init_cells)));
    double spacing = sqrt((double)(par.size_init_cells)) + 0.5;
    double stepx, stepy;
  
    for (i=1;i<=il;i++) {
      for (j=1;j<il;j++){
        if (k<n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + spacing);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + spacing);
	  cellnum++;
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {
	        sigma[xx][yy]=cellnum;
              }
	    }
	  }
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 13){
    //fill cells from top, NC cells in the middle top as a cluster, Placodes as a spaced block below them
    
    int iln = (int)(sqrt(par.n_NC)), k=0;
    double spacing = 1.5*sqrt((double)(par.target_area)/PI);
    double step, stepx, stepy;
    int cxn=(int)(sizex/2-1);
    int cyn=(int)((iln/2.0)*spacing+1.0);
    int xx,yy,i,j;
    int xmin, xmax, ymin, ymax;
    int imx = (int)(sizex / (spacing + 0.1));
    int imy = (int)(sizey / (spacing + 0.1));
    bool addedCell=false;
    // The NC cluster

    for (j=0;j<=iln;j++){
      for (i=0;i<=iln;i++) {
        if (k<par.n_NC){
	  cellnum++;
	  step = spacing * (double) ( i-((double)(iln-1)/2) );
	  xx=(int)((double)(cxn) + step);
	  step = spacing * (double) ( j-((double)(iln-1)/2) );
          yy=(int)((double)(cyn) + step);
	  sigma[xx][yy]=cellnum;
	  k++;
	}
      }
    }
    
    // Placodes
    if (par.init_PL_distance >= 0){
      spacing = par.init_PL_distance;
    } else {
      spacing = sqrt((double)(par.size_init_cells));
    }
    cellnum++;
    imx = (int)(sizex / (spacing + 0.1));
    imy = (int)(sizey / (spacing + 0.1));
  
    for (j=iln;j<imy;j++){
      for (i=1;i<imx;i++) {
        if (k<n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + spacing);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + spacing);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {cellnum++; addedCell=false;}
	  k++;
	}
      }
    }
  } else if (par.init_cell_distribution == 14){
    // Fill cells from top, downwards and left as in "4" but only one row of NC is actually present
    // the rest is just initialised, but without area at this point. They will be introduced later, but 
    // here their cell numbers are still reserved. 
    
    double spacing = sqrt((double)(par.size_init_cells));
    double stepx, stepy, ncymax=0;
    int cellnum = 1;
    int xx,yy,i,j;
    int xmin, xmax, ymin, ymax;
    int imx = (int)(sizex / (spacing + 0.1));
    int imy = (int)(sizey / (spacing + 0.1));
    int cs = (int)(spacing);
    bool addedCell=false;
  
    if (cs == 0) cs = 1;
    for (j=1;j<2;j++){
      for (i=1;i<imx;i++) {
        if (cellnum <= par.n_NC){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + cs);
	  ymin=(int)(double)(stepy);
	  ymax=(int)(double)(stepy + cs);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {cellnum++; addedCell=false; ncymax=ymax;}
	  else {(*cell)[cellnum].Apoptose(); cellnum++;};
        }
      }
    }
    
    // Placodes
    cellnum=par.n_NC+1;
    if (par.init_PL_distance >= 0){
      spacing = par.init_PL_distance;
    } else {
      spacing = sqrt((double)(par.size_init_cells));
    }
    imx = (int)(sizex / (spacing + 0.1));
    imy = (int)((sizey-ncymax) / (spacing + 0.1));

    for (j=0;j<imy;j++){
      for (i=1;i<imx;i++) {
        if (cellnum <= n_cells){
	  stepx = spacing * (double) ( i );
	  stepy = spacing * (double) ( j );
	  xmin=(int)(double)(stepx);
	  xmax=(int)(double)(stepx + cs);
	  ymin=(int)(double)(stepy + ncymax);
	  ymax=(int)(double)(stepy + ncymax + cs);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {
	  cellnum++; addedCell=false;}
	}
      }
    }
    if (cellnum <= n_cells) {
      std::cerr << "Warning! Not enough space to initialise all cells! Starting with " << (cellnum-par.n_NC) << " instead of " << (n_cells-par.n_NC) << ".\n";
      cellnum = n_cells;
    }
  } else if (par.init_cell_distribution == 15){
    // Fill cells from top, downwards and left as in "4" but only one row of NC is actually present
    // the rest is just initialised, but without area at this point. They will be introduced later, but 
    // here their cell numbers are still reserved. Difference from 14 is that the remainder spaces 
    // around and inside the population are distributed more equally. 
    
    double spacing = sqrt((double)(par.size_init_cells)), ncymax=0, gap_xspacing, gap_yspacing;
    int nx, ny, nx_gaps, ny_gaps;
    nx = (int)((sizex-2) / (spacing));
    nx_gaps = fmod(sizex, spacing);
    ny_gaps = fmod(sizey, spacing);
    
    if (nx_gaps > 0){
      gap_xspacing = (double)(sizex-2) / (double)(nx_gaps);
    } else {
      gap_xspacing = (double)(2*sizex);
    }
    if (ny_gaps > 0){
      gap_yspacing = (double)(sizey-2) / (double)(ny_gaps);
    } else {
      gap_yspacing = (double)(2*sizey);
    }

    int cellnum = 1;
    int xx,yy,i,j;
    int xmin, xmax, ymin, ymax, xgap, ygap;
    bool addedCell=false;
  
    for (j=0;j<1;j++){
      for (i=0;i<nx;i++) {
        if (cellnum <= par.n_NC){
	  xgap=(int)((double)(i)/gap_xspacing + 0.5);
	  xmin=(int)(2 + i * spacing + xgap);
	  xmax=(int)(xmin + spacing);
	  ygap=(int)((double)(j)/gap_yspacing + 0.5);
	  ymin=(int)(2 + j * spacing + ygap);
	  ymax=(int)(ymin + spacing);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {cellnum++; addedCell=false; ncymax=ymax;}
	  else {(*cell)[cellnum].Apoptose(); cellnum++;};
        }
      }
    }
    
    // Placodes
    cellnum=par.n_NC+1;
    if (par.init_PL_distance >= 0){
      spacing = par.init_PL_distance;
    } else {
      spacing = sqrt((double)(par.size_init_cells));
    }
    nx = (int)((sizex-2) / (spacing));
    ny = (int)((sizey-2) / (spacing));
    nx_gaps = fmod(sizex, spacing);
    ny_gaps = fmod(sizey, spacing);
    
    if (nx_gaps > 0){
      gap_xspacing = (double)(sizex-2) / (double)(nx_gaps);
    } else {
      gap_xspacing = (double)(2*sizex);
    }
    if (ny_gaps > 0){
      gap_yspacing = (double)(sizey-2) / (double)(ny_gaps);
    } else {
      gap_yspacing = (double)(2*sizey);
    }

    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++) {
        if (cellnum <= n_cells){
	  xgap=(int)((double)(i)/gap_xspacing + 0.5);
	  xmin=(int)(2 + i * spacing + xgap);
	  xmax=(int)(xmin + spacing);
	  ygap=(int)((double)(j)/gap_yspacing + 0.5);
	  ymin=(int)(ncymax + j * spacing + ygap);
	  ymax=(int)(ymin + spacing);
	  for (xx=xmin; xx<xmax; xx++){
	    for (yy=ymin;yy<ymax;yy++){
	      if (xx>0 && xx<sizex && yy>0 &&yy<sizey) {sigma[xx][yy]=cellnum; addedCell=true;}
	    }
	  }
	  if (addedCell) {cellnum++; addedCell=false;}
        }
      }
    }
    if (cellnum <= n_cells) {
      std::cerr << "Warning! Not enough space to initialise all cells! Starting with " << (cellnum-par.n_NC) << " instead of " << (n_cells-par.n_NC) << ".\n";
      cellnum = n_cells;
    }
  }

  if (par.EdenGrowth == false) {
    // Do Eden growth for a number of time steps
    EdenGrowth(cell_size, cellnum, n_cells, new_sigma);
  };

  free(new_sigma[0]);
  free(new_sigma);
  
  return cellnum;
}
  

// Predicate returns true when connectivity is locally preserved
// if the value of the central site would be changed
bool CellularPotts::ConnectivityPreservedP(int x, int y) {
  
  // Use local nx and ny in a cyclic order (starts at upper left corner)
  // first site is repeated, for easier looping
  const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1 };
  const int cyc_ny[10] = {0, -1,-1,-1, 0, 1, 1,  1,  0, -1 };
  
  int sxy=sigma[x][y]; // the central site
  if (sxy==0) return true;

  int n_borders=0; // to count the amount of sites in state sxy bordering a site !=sxy

  static int stack[8]; // stack to count number of different surrounding cells
  int stackp=-1;
  bool one_of_neighbors_medium=false;
  
  for (int i=1;i<=8;i++) {
    
    int s_nb=sigma[x+cyc_nx[i]][y+cyc_ny[i]];
    int s_next_nb=sigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]];
    
    if ((s_nb==sxy || s_next_nb==sxy) && (s_nb!=s_next_nb)) {
      
      // check whether s_nb is adjacent to non-identical site,
      // count it
      n_borders++;
    }
    int j;
    bool on_stack_p=false;
    
    // we need the next heuristic to prevent stalling at
    // cell-cell borders
    // do not enforce constraint at two cell interface(no medium)
    if (s_nb) {
      for (j=stackp;j>=0;j--) {
	if (s_nb==stack[j]) {
	  on_stack_p=true;
	  break;
	}
      }
      if (!on_stack_p) {
	if (stackp>6) {
	  cerr << "Stack overflow, stackp=" << stackp << "\n";
	}
	stack[++stackp]=s_nb;
      }
    } else {
      one_of_neighbors_medium=true;
    }
  }
  
  // number of different neighbours is stackp+1;
  if (n_borders>2 && ( (stackp+1)>2 || one_of_neighbors_medium) ) {
    return false;
  }
  else 
    return true;

}

bool CellularPotts::CheckConnectivityConstraints(int xp, int yp, int x, int y){
// Examines that the expanding cell remains connected and simply connected,
// and the rectracting cell remains connected
// when copy sigma[xp][yp] into sigma[x][y]
// Returns true, if the above are fulfilled (connectivity preserved), returns false if not.

    //Retracting cell remains connected:
    if (sigma[x][y] != 0) {
	int sxy = sigma[x][y]; // retracting cell's id
	int empty_connected = 0;
	const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1};
	const int cyc_ny[10] = {0, -1,-1,-1, 0, 1, 1,  1,  0, -1};

	if (sigma[x+cyc_nx[0]][y+cyc_ny[0]] != sxy && sigma[x+cyc_nx[7]][y+cyc_ny[7]] != sxy) {
	    empty_connected -= 1;
	}

	for (int i=0;i<8;i++) {
	    if (sigma[x+cyc_nx[i]][y+cyc_ny[i]] != sxy) {
		if (empty_connected == 1) {return false;}
		else {
		    if((cyc_nx[i] != 0 && cyc_ny[i] != 0) &&
			sigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]] == sxy ) {
			empty_connected+=1;
		    }
		    else if ((cyc_nx[i]==0 || cyc_ny[i]==0) &&
			sigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]] == sxy &&
			sigma[x+cyc_nx[i+2]][y+cyc_ny[i+2]] == sxy) {
			    empty_connected+=1;
		    }
		}
	    }
	}
    }
    // Do not check expanding cell if inserting medium

    //Expanding cell remains connected:
    int sxyp=sigma[xp][yp];  //expanding cell's id
    int c_nx[4]={-1,0,1,0};
    int c_ny[4]={0,-1,0,1};
    if (sigma[xp][yp]!=0 && !(xp==x && yp==y)) {
	bool empty_neigh = false;
	for (int i=0;i<4;i++){
	    if (sigma[x+c_nx[i]][y+c_ny[i]] == sxyp) {
		empty_neigh = true;
	    }
	}
	if (empty_neigh == false) return false;
    }

    //Expanding cell remains simply connected:
    if (sigma[xp][yp]!=0 && !(xp==x && yp==y)) {

        int empty_connected=0;
        int sx = xp-x;
        int sy = yp-y;
        int m_nx[10];
        int m_ny[10];
        if ( sx == 0 || sy == 0 ){
	    m_nx[0]=sx;m_nx[1]=sx-sy;m_nx[2]=-sy;m_nx[3]=-sy-sx;m_nx[4]=-sx;m_nx[5]=sy-sx;m_nx[6]=sy;m_nx[7]=sx+sy;m_nx[8]=sx;m_nx[9]=sx-sy;
	    m_ny[0]=sy;m_ny[1]=sy-sx;m_ny[2]=-sx;m_ny[3]=-sy-sx;m_ny[4]=-sy;m_ny[5]=sx-sy;m_ny[6]=sx;m_ny[7]=sx+sy;m_ny[8]=sy;m_ny[9]=sy-sx;
	} else {
	    m_nx[0]=sx;m_nx[1]=0;m_nx[2]=-sx;m_nx[3]=-sx;m_nx[4]=-sx;m_nx[5]=0;m_nx[6]=sx;m_nx[7]=sx;m_nx[8]=sx;m_nx[9]=0;
	    m_ny[0]=sy;m_ny[1]=sy;m_ny[2]=sy;m_ny[3]=0;m_ny[4]=-sy;m_ny[5]=-sy;m_ny[6]=-sy;m_ny[7]=0;m_ny[8]=sy;m_ny[9]=sy;
	    if (sigma[x+m_nx[1]][y+m_ny[1]] != sxyp && sigma[x+m_nx[7]][y+m_ny[7]] != sxyp) {
		empty_connected -= 1;
	    }
	}

	for (int i=0;i<8;i++){
	    if ( sigma[x+m_nx[i]][y+m_ny[i]] != sxyp ){
		if (empty_connected == 1){
		    return false;
		}
		else {
		    if((m_nx[i] != 0 && m_ny[i] != 0) &&
			sigma[x+m_nx[i+1]][y+m_ny[i+1]] == sxyp ) {
			empty_connected+=1;
		    }
		    else if ((m_nx[i]==0 || m_ny[i]==0) &&
			sigma[x+m_nx[i+1]][y+m_ny[i+1]] == sxyp &&
			sigma[x+m_nx[i+2]][y+m_ny[i+2]] == sxyp) {
			    empty_connected+=1;
		    }
		}
	    }
	}
    }

    return true;
}

double CellularPotts::CellDensity(void) const {
  
  // return the density of cells
  int sum=0;
  for (int i=0;i<sizex*sizey;i++) {
    if (sigma[0][i]) {
      sum++;
    }
  }
  return (double)sum/(double)(sizex*sizey);

}

double CellularPotts::MeanCellArea(void) const {
  
  int sum_area=0, n=0;
  double sum_length=0.;
  vector<Cell>::iterator c=cell->begin(); ++c;
  
  for (; 
	c!=cell->end();
	c++) {
    
    sum_area+=c->Area();
    sum_length+=c->Length();
    n++;    
  }
  
  cerr << "Mean cell length is " << sum_length/((double)n) << endl;
  return (double)sum_area/(double)n;
}

void CellularPotts::ResetTargetLengths(void)  {
   vector<Cell>::iterator c=cell->begin(); ++c;

   for (;
        c!=cell->end();
        c++) {

     c->SetTargetLength(par.target_length);

} 

}

void CellularPotts::SetRandomTypes(void) {
  
  // each cell gets a random type 1..maxtau
  
  vector<Cell>::iterator c=cell->begin(); ++c;
  
  for (;
       c!=cell->end();
       c++) {
    
    int celltype = RandomNumber(Cell::maxtau);
    c->setTau(celltype);
    
  };
  
}

void CellularPotts::SetNCandPlacodeTypes(int n_NC) {
  
  // assign the first n_NC cell type=1, and the rest is type=2
  
  vector<Cell>::iterator c=cell->begin(); ++c;
  int n=1;
  
  for (;c!=cell->end();c++) {
    if (n<=n_NC) c->setTau(NCType);
    else c->setTau(PlacodeType);
    n++;
  };
}

void CellularPotts::GrowAndDivideCells(int growth_rate) {

  vector<Cell>::iterator c=cell->begin(); ++c;
  vector<bool> which_cells((int)(cell->size()));

  for (;
       c!=cell->end();
       c++) {

    // only tumor cells grow and divide
    if (c->getTau()==2) {
     
      c->SetTargetArea(c->TargetArea()+growth_rate);
    
      if (c->Area()>par.target_area) {
	which_cells[c->Sigma()]=true;
      } else {
	which_cells[c->Sigma()]=false;
      }

      if (c->chem[1]<0.9) { //arbitrary oxygen threshold for the moment
	c->setTau(3);
      }
    } else {
      which_cells[c->Sigma()]=false;
    }

  }

  DivideCells(which_cells);

}


