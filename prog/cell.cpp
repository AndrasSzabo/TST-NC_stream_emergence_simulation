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
#include <list>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include "cell.h"
#include "sticky.h"
#include "parameter.h"
#include "dish.h"

#define HASHCOLNUM 255

extern Parameter par;

double **Cell::J=0;
int Cell::amount=0;
int Cell::capacity=0;
int Cell::maxsigma=0;
int Cell::maxtau=0;

//Cell::Cell(const Dish &who) : Cytoplasm(who);
// Note: g++ wants to have the body of this constructor in cell.hh
// body is defined in "ConstructorBody" below
class Dish;

using namespace std;


Cell::~Cell(void) {

  amount--;
  if (amount==0) {
    // clear J if last cell has been destructed
    for (int i=0;i<par.nJ;i++){ 
      if (J[i]) {free(J[i]);}
    }
    if (J) free(J);
    capacity=0;
    maxsigma=0;
    J=0;
  }
  delete[] chem;
  delete[] contactDurationWith;
  delete[] cytosk_bonds;
}

Cell::Cell(Cell &mother, int settau) {
  
  // settau is optional
  ConstructorBody(settau);
  
  // Administrate ancestry
  mother.daughter=sigma;
  this->mother=mother.sigma;
  times_divided=++mother.times_divided;
  owner=mother.owner;
  
  date_of_birth=owner->Time();

  colour_of_birth=mother.colour;
  colour=mother.colour;
  
  alive=mother.alive;
  
  tau=mother.tau;
  target_length = mother.target_length;

  n_copies=0;

  grad[0]=0.;
  grad[1]=0.;
  polx=0;
  poly=0;

}

void Cell::ConstructorBody(int settau) {
  
  // Note: Constructor of Cytoplasm will be called first
  alive=true;
  colour=1; // undifferentiated
  
  colour_of_birth=1;
  date_of_birth=0;
  times_divided=0;
  mother=0;
  daughter=0;
    
  // add new elements to each of the dimensions of "J"
  
  // amount gives the total number of Cell instantiations (including copies)
  amount++;
  
  // maxsigma keeps track of the last cell identity number given out to a cell
  sigma=maxsigma++;
  
  if (!J) {
    ReadStaticJTable();
  }
  
  tau=settau;
  area=0;
  target_area=0;
  length=0;
  target_length=par.target_length;
  sum_x=0;
  sum_y=0;
  sum_xx=0;
  sum_yy=0;
  sum_xy=0;

  //  growth_threshold=par.dthres;
  growth_threshold=0;
  n_copies=0;

  chem = new double[par.n_chem];
  contactDurationWith = new int[par.n_init_cells+1];
  cytosk_bonds = new double[par.n_init_cells+1];

}

void Cell::CenterOfMass(double *cx, double *cy) {
    double a;
    
    a=(double)area;
    if (a>0){
      *cx = (double)sum_x/a;
      *cy = (double)sum_y/a;
    }
    else {
    //probably new cell, has no area yet 
      *cx = -1.0;
      *cy = -1.0;
    };
}

void Cell::Theta(double* t){

  if (area>0){
    double iyy=(double)sum_xx-(double)sum_x*sum_x/(double)area;
    double ixx=(double)sum_yy-(double)sum_y*sum_y/(double)area;
    double ixy=-(double)sum_xy+(double)sum_x*sum_y/(double)area;
          
    double rhs1=(ixx+iyy)/2., rhs2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )/2.;
              
    double lambda_b=rhs1+rhs2;
    double lambda_a=rhs1-rhs2;
    if (lambda_a>lambda_b){
      if (lambda_b!=0.) *t=sqrt(lambda_a/lambda_b);
    }
    else {
      if (lambda_a!=0.) *t=sqrt(lambda_b/lambda_a);  
    };                    
  }
  else {
    *t=-1.0;
  };
};

void Cell::getLongAxisDir (double* vx,double* vy){
  if (area>0){
    double iyy=(double)sum_xx-(double)sum_x*sum_x/(double)area;
    double ixx=(double)sum_yy-(double)sum_y*sum_y/(double)area;
    double ixy=-(double)sum_xy+(double)sum_x*sum_y/(double)area;
    double rhs1=(ixx+iyy)/2., rhs2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )/2.;
    double lambda_b=rhs1+rhs2;
    double lambda_a=rhs1-rhs2;

    if (lambda_a>lambda_b){
      *vx=ixy;   //the smaller eigenvalue corresponds with the longer axis!
      *vy=lambda_b-ixx;
    }else{
      *vx=ixy;
      *vy=lambda_a-ixx;
    };
  }else{
    *vx=0.0;
    *vy=0.0;
  };
};

void Cell::setDisplacement (double dx, double dy){
//set the displacement of the cell

  displacement[0] = dx;
  displacement[1] = dy;
};

void Cell::getDisplacement (double *dx, double *dy){
//returns the displacement of the cell

  *dx=displacement[0];
  *dy=displacement[1];
};


void Cell::setPolarisation (double px, double py){
//set the polarisation vector to the give vales 
  
  polx = px;
  poly = py;

};

void Cell::getPolarisation (double *px, double *py){
//read the polarisation vector of the cell

  *px=polx;
  *py=poly;

};

void Cell::setPolDecay (double p){
//set the polarisation strength of cell
  
  PolDecay = p;
};

void Cell::setPolDecayFree (double p){
//set the polarisation strength of cell
  
  PolDecayFree = p;
};

void Cell::setPolDecayContact (double p){
//set the polarisation strength of cell
  
  PolDecayContact = p;
};
double Cell::getPolDecay (void){
//read the polarisation strength of the cell

  return (PolDecay);
};

double Cell::getPolDecayFree (void){
//read the polarisation strength of the cell

  return (PolDecayFree);
};

double Cell::getPolDecayContact (void){
//read the polarisation strength of the cell

  return (PolDecayContact);
};

void Cell::setPolStr (double p){
//set the polarisation strength of cell
  
  PolStr = p;
};

double Cell::getPolStr (void){
//read the polarisation strength of the cell

  return (PolStr);
};

void Cell::setEcmContact (int e){
  ecmcontact = e;
};

int Cell::getEcmContact (void){
  return (ecmcontact);
};


void Cell::ReadStaticJTable() {

  // Allocate
  int n; // number of taus
  n=par.nJ;
  maxtau=n-1;
  if (J) { 
    for (int i=0; i<n; i++) {if (J[i]) {free(J[i]);};}; 
    free(J); 
  }
  J=(double **)malloc(n*sizeof(double *));
  J[0]=(double *)malloc(n*n*sizeof(double));
  for (int i=1;i<n;i++) {
    J[i]=J[i-1]+n;
  }
  
  capacity = n;
  {for (int i=0;i<n;i++) {
    for (int j=0;j<=i;j++) {
      J[i][j] = par.J[i][j];
      //symmetric:
      J[j][i] = J[i][j];
    }
  
  }}
}

double Cell::EnergyDifference(const Cell &cell2) const
{ 
  if (sigma==cell2.sigma) return 0;
  
  if (par.internalrelaxation == 0) return J[tau][cell2.tau];
  else return par.relaxEM[tau][cell2.tau];
  
}

void Cell::ClearJ(void) {

  for (int i=0;i<capacity*capacity;i++) {
    J[0][i]=EMPTY;
  }
}



void Cell::initCellParameters(){
  // Initialise cellular parameters
  int NCType=1, PlacodeType=2;


  if (tau==NCType){
                                                      // NC cell
    for (int l=0;l<par.n_chem;l++){                   // chemotaxis
      if (par.saturation>0)                           // renormalise with saturation if there is any
        setChem(l, par.chemotaxisNC[l]*par.saturation);   // par.chemotaxis will be the max. value of chemotaxis energy
      else 
        setChem(l, par.chemotaxisNC[l]);
    }

    setPolStr(par.Polarisation);           // polarisation parameters
    setPolDecay(par.Pol_decay);
    setPolDecayContact(par.Pol_decay);
    setPolDecayFree(par.Pol_decayFree);
  } else if (tau==PlacodeType){
                                                        // Placode cell
    for (int l=0;l<par.n_chem;l++){                   // chemotaxis
      if (par.saturation>0)                           // renormalise with saturation if there is any
        setChem(l, par.chemotaxisPlacode[l]*par.saturation);   // par.chemotaxis will be the max. value of chemotaxis energy
      else 
        setChem(l, par.chemotaxisPlacode[l]);
    }
      
    setPolStr(par.PolarisationPlacode);    // polarisation parameters
    setPolDecay(par.Pol_decayPlacode);
    setPolDecayContact(par.Pol_decayPlacode);
    setPolDecayFree(par.Pol_decayFreePlacode);
  }
}  

void Cell::initContacts(){
  // Initiate contact durations
  for (int j=1;j<= par.n_init_cells; j++){
    contactDurationWith[j]=0;
    cytosk_bonds[j] = -1;
  }
}
