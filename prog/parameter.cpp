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


#include "parameter.h"
#include <cstdio>
#include <math.h>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include "output.h"
#include "parse.h"

Parameter::Parameter() {

  T = 50;
  target_area = 100;
  target_length = 60;
  lambda = 50;
  lambda2 = 5.0;
  colour_cells = 0.0;
  colour_active_cells = 0.0;
  meas_dens = false;
  show_time = false;
  
  Polarisation = 0;
  Pol_decay = 1.0;
  Pol_decayFree = 1.0;
  PolarisationPlacode = 0;
  Pol_decayPlacode = 1.0;
  Pol_decayFreePlacode = 1.0;
  PolDynamicUpdateFree = true;
  plotPolarityVectors = 0;
  polVecPlotLen = 0;
  
  EdenGrowth = false;
  nJ = 2;
  J = (double**) NULL;
  relaxEM = (double**) NULL;
  
  SC = (double**) NULL;
  CC = (double**) NULL;
  UCC = (double**) NULL;
  cytosk_bonds = false;
  plotCytoskBonds = 0;
  breakBondIfNoContact = false;
  
  CILprob = (double*) NULL;
  repolarisationPM = 0.0;
  repolarisationDelay = 0;
  plotCILcells = 0;
  plotChemotaxCells = 0;

  PassiveRelaxationSteps = false;
  PassiveRelaxationWindow = 5;
  PassiveRelaxationMonitorDH = false;
  
  insertMediumProb = 0;
  conn_diss = 2000;
  vecadherinknockout = true;
  extension_only_chemotaxis = true;
  chemotaxis = 1000;
  border_energy = 100;
  neighbours = 2;
  periodic_boundaries = false;
  n_chem = 1;
  saturation = 0;
  Sdf1UptakeNC = 0;
  pdeExtendedField = false;
  pde_initialization = 0;
  n_init_cells = 100;
  n_NC = 0;
  size_init_cells = 10;
  sizex = 200;
  sizey = 200;
  divisions = 0;
  mcs = 10000;
  rseed = -1;
  subfield = 1.0;
  relaxation = 0;
  internalrelaxation = 0;
  init_cell_distribution = 0;
  sidePadding = 0;
  init_cluster_distance = 10;
  init_PL_distance = -1.0;
  takeoutCellsBelow = -1;
  storage_stride = 10;
  print_Areas = 0;
  print_AreasB = false;
  print_CM = 0;
  print_CMB = false; 
  print_CMsample = 50;
  print_CMfreq = 1000;
  print_CMarea = false;
  print_sigma = 0;
  print_asParticles = 0;
  print_PIF = 0;
  time = 0;
  print_sigmaB = false;
  graphics = true;
  printcc = false;
  store = false;
  store_start = 0;
  storeNCimage = false;
  storePLimage = false;
  //datadir = strdup("data_film");
}

Parameter::~Parameter() {
  // destruct parameter object
  // free string parameter
  CleanUp();
}

void Parameter::CleanUp(void) {

  if (diff_coeff) 		free(diff_coeff);
  if (chemotaxisNC) 		free(chemotaxisNC);
  if (chemotaxisPlacode) 	free(chemotaxisPlacode);
  if (decay_rate) 		free(decay_rate);
  if (secr_rate) 		free(secr_rate);
  if (plotcc)			free(plotcc);
  if (plotccCol)		free(plotccCol);
  if (plotccgrad)		free(plotccgrad);
  if (storeccgrad)		free(storeccgrad);
  if (storecc)			free(storecc);
  if (datadir) 			free(datadir);
  if (CILprob)			free(CILprob);
  if (J){
     for (int i=0;i<nJ;i++){
       if (J[i]) free (J[i]);};
     free (J);
  };
  if (relaxEM){
     for (int i=0;i<nJ;i++){
       if (relaxEM[i]) free (relaxEM[i]);};
     free (relaxEM);
  };
  
  if (CC){
     for (int i=0;i<nJ;i++){
       if (CC[i]) free (CC[i]);};
     free (CC);
  };
  if (UCC){
     for (int i=0;i<nJ;i++){
       if (UCC[i]) free (UCC[i]);};
     free (UCC);
  };
  if (SC){
     for (int i=0;i<nJ;i++){
       if (SC[i]) free (SC[i]);};
     free (SC);
  };
  
}

void Parameter::Read(const char *filename) {
  
  static bool ReadP=false;

  if (ReadP) {

    //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
    CleanUp();
    
  } else
    ReadP=true;

  FILE *fp=OpenReadFile(filename);

  T = fgetpar(fp, "T", 50, true);
  target_area = igetpar(fp, "target_area", 100, true);
  target_length = igetpar(fp, "target_length", 60, true);
  lambda = fgetpar(fp, "lambda", 50, true);
  lambda2 = fgetpar(fp, "lambda2", 5.0, true);
  colour_cells = fgetpar(fp, "colour_cells", 0.0, true);
  colour_active_cells = igetpar(fp, "colour_active_cells", 0, true);
  show_time = bgetpar(fp, "show_time", false, true);
  
  Polarisation = fgetpar(fp, "Polarisation", 0.0, true);
  Pol_decay = fgetpar(fp, "Pol_decay", 1.0, true);
  Pol_decayFree = fgetpar(fp, "Pol_decayFree", 1.0, true);
  PolarisationPlacode = fgetpar(fp, "PolarisationPlacode", 0.0, true);
  Pol_decayPlacode = fgetpar(fp, "Pol_decayPlacode", 1.0, true);
  Pol_decayFreePlacode = fgetpar(fp, "Pol_decayFreePlacode", 1.0, true);
  plotPolarityVectors = igetpar(fp, "plotPolarityVectors", 0, true);
  polVecPlotLen = fgetpar(fp, "polVecPlotLen", 0, true);
  PolDynamicUpdateFree = bgetpar(fp, "PolDynamicUpdateFree", true, true);
  PolDynamicUpdateContact = bgetpar(fp, "PolDynamicUpdateContact", true, true);
  EdenGrowth = bgetpar(fp, "EdenGrowth", false, true);
  nJ = igetpar(fp, "nJ", 2, true);
  {
    int i,j;
    char s[10];
    
    J = (double**)calloc(nJ,sizeof(double*));
    for (i=0;i<nJ;i++){
      J[i]=(double*)calloc(nJ,sizeof(double)); };  
    for (i=0;i<nJ;i++){
      for (j=0;j<=i;j++){
        sprintf(s,"J[%i][%i]",i,j);
        J[i][j]=(double) fgetpar(fp,s, 20, true);      
  };};};
  {
    int i,j;
    char s[100];

    relaxEM = (double**)calloc(nJ,sizeof(double*));
    for (i=0;i<nJ;i++){
      relaxEM[i]=(double*)calloc(nJ,sizeof(double)); };  
    for (i=0;i<nJ;i++){
      for (j=0;j<=i;j++){
        sprintf(s,"relaxEM[%i][%i]",i,j);
        relaxEM[i][j]=(double) fgetpar(fp,s, 20, true);      
        relaxEM[j][i]=relaxEM[i][j];      
  };};};
								 
  
  int i,j;
  char s[200];
              
  CC  = (double**)calloc(nJ,sizeof(double*));
  SC  = (double**)calloc(nJ,sizeof(double*));
  UCC = (double**)calloc(nJ,sizeof(double*));
  CILprob = (double*)calloc(nJ,sizeof(double));
  plotCILcells = igetpar(fp, "plotCILcells", 0, true);
  plotChemotaxCells = igetpar(fp, "plotChemotaxCells", 0, true);
  for (i=0;i<nJ;i++){
    CC[i] =(double*)calloc(nJ,sizeof(double)); 
    SC[i] =(double*)calloc(nJ,sizeof(double)); 
    UCC[i]=(double*)calloc(nJ,sizeof(double)); 
    if (i==0){
      CILprob[i]=0;
    } else {
      sprintf(s,"CILprob[%i]",i);
      CILprob[i]=fgetpar(fp,s, 0, true);
    }
  }
  for (i=1;i<nJ;i++){
    for (j=1;j<=i;j++){
      if (i==0 || j==0){
        CC[i][j]=0;
        SC[i][j]=0;
        UCC[i][j]=0;
        CC[j][i]=0;
        SC[j][i]=0;
        UCC[j][i]=0;
      } else {
        sprintf(s,"CC[%i][%i]",i,j);
        CC[i][j]=fgetpar(fp,s, 0.0, true);
        sprintf(s,"SC[%i][%i]",i,j);
        SC[i][j]=fgetpar(fp,s, 0.0, true);
        sprintf(s,"UCC[%i][%i]",i,j);
        UCC[i][j]=fgetpar(fp,s, 0.0, true);
        CC[j][i]=CC[i][j];
        SC[j][i]=SC[i][j];
        UCC[j][i]=UCC[i][j];
      }
    }
  }

  repolarisationPM = fgetpar(fp, "repolarisationPM", 0, true);
  repolarisationDelay = igetpar(fp, "repolarisationDelay", 0, true);
  PassiveRelaxationSteps = bgetpar(fp, "PassiveRelaxationSteps", false, true);
  PassiveRelaxationWindow = igetpar(fp, "PassiveRelaxationWindow", 5, true);
  PassiveRelaxationMonitorDH = bgetpar(fp, "PassiveRelaxationMonitorDH", false, true);
  cytosk_bonds = bgetpar(fp, "cytosk_bonds", false, true);
  plotCytoskBonds = fgetpar(fp, "plotCytoskBonds", 0, true);
  breakBondIfNoContact = bgetpar(fp, "breakBondIfNoContact", false, true);
  insertMediumProb = fgetpar(fp, "insertMediumProb", 0, true);
  conn_diss = fgetpar(fp, "conn_diss", 2000, true);
  border_energy = fgetpar(fp, "border_energy", 100, true);
  neighbours = igetpar(fp, "neighbours", 2, true);
  periodic_boundaries = bgetpar(fp, "periodic_boundaries", false, true);
  n_chem = igetpar(fp, "n_chem", 1, true);
  diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, "0.0", true);
  decay_rate = dgetparlist(fp, "decay_rate", n_chem, "0.0", true);
  secr_rate = dgetparlist(fp, "secr_rate", n_chem, "0.0", true);
  saturation = fgetpar(fp, "saturation", 0, true);
  Sdf1UptakeNC = fgetpar(fp, "Sdf1UptakeNC", 0, true);
  vecadherinknockout = bgetpar(fp, "vecadherinknockout", true, true);
  extension_only_chemotaxis = bgetpar(fp, "extension_only_chemotaxis", true, true);
  chemotaxis = fgetpar(fp, "chemotaxis", 1000, true);
  chemotaxisNC = dgetparlist(fp, "chemotaxisNC", n_chem, "0.0", true);
  chemotaxisPlacode = dgetparlist(fp, "chemotaxisPlacode", n_chem, "0.0", true);
  pdeExtendedField = bgetpar(fp, "pdeExtendedField", false, true);
  pde_initialization = igetpar(fp, "pde_initialization", 0, true);
  n_init_cells = igetpar(fp, "n_init_cells", 100, true);
  n_NC = igetpar(fp, "n_NC", 0, true);
  size_init_cells = igetpar(fp, "size_init_cells", 10, true);
  sizex = igetpar(fp, "sizex", 200, true);
  sizey = igetpar(fp, "sizey", 200, true);
  divisions = igetpar(fp, "divisions", 0, true);
  mcs = igetpar(fp, "mcs", 10000, true);
  rseed = igetpar(fp, "rseed", -1, true);
  subfield = fgetpar(fp, "subfield", 1.0, true);
  relaxation = igetpar(fp, "relaxation", 0, true);
  init_cell_distribution = igetpar(fp, "init_cell_distribution", 0, true);
  sidePadding = igetpar(fp, "sidePadding", 0, true);
  init_cluster_distance = igetpar(fp, "init_cluster_distance", 10, true);
  init_PL_distance = fgetpar(fp, "init_PL_distance", 0, true);
  takeoutCellsBelow = fgetpar(fp, "takeoutCellsBelow", -1, true);
  storage_stride = igetpar(fp, "storage_stride", 10, true);
  print_sigma = igetpar(fp, "print_sigma", 0, true);
  print_CMstart = igetpar(fp, "print_CMstart", 0, true);
  print_CM = igetpar(fp, "print_CM", 0, true);
  print_CMsample = igetpar(fp, "print_CMsample", 0, true);
  print_CMfreq = igetpar(fp, "print_CMfreq", 0, true);
  print_CMarea = bgetpar(fp, "print_CMarea", false, true);
  print_Areas = igetpar(fp, "print_areas", 0, true);
  print_asParticles = igetpar(fp, "print_asParticles", 0, true);
  print_PIF = igetpar(fp, "print_PIF", 0, true);
  graphics = bgetpar(fp, "graphics", true, true);
  plotcc = dgetparlist(fp, "plotcc", n_chem, "0.0", true);
  plotccCol = dgetparlist(fp, "plotccCol", n_chem, "0.0", true);
  plotccgrad = dgetparlist(fp, "plotccgrad", n_chem, "0.0", true);
  storeccgrad = dgetparlist(fp, "storeccgrad", n_chem, "0.0", true);
  storecc = dgetparlist(fp, "storecc", n_chem, "0.0", true);
  printcc = bgetpar(fp, "printcc", false, true);
  store = bgetpar(fp, "store", false, true);
  store_start = igetpar(fp, "store_start", 0, true);
  storeNCimage = bgetpar(fp, "storeNCimage", false, true);
  storePLimage = bgetpar(fp, "storePLimage", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);

  
  if (target_length == -1) 
    target_length = (int)(2.0 * sqrt((double)(target_area)/3.14));
}


const char *sbool(const bool &p) {

  const char *true_str="true";
  const char *false_str="false";
  if (p)
    return true_str;
  else
    return false_str;
}

void Parameter::Write(ostream &os) const {

  os << " T = " << T << endl;
  os << " target_area = " << target_area << endl;
  os << " target_length = " << target_length << endl;
  os << " lambda = " << lambda << endl;
  os << " lambda2 = " << lambda2 << endl;
  os << " colour_cells = " << colour_cells << endl;
  os << " colour_active_cells = " << colour_active_cells << endl;
  os << " meas_dens = " << meas_dens << endl;
  os << " show_time = " << sbool(show_time) << endl;
  os << " Polarisation = " << Polarisation << endl;
  os << " Pol_decay = " << Pol_decay << endl;
  os << " EdenGrowth = " << sbool(EdenGrowth) << endl;
          
  os << " insertMediumProb = " << insertMediumProb << endl;
  os << " conn_diss = " << conn_diss << endl;
  os << " vecadherinknockout = " << sbool(vecadherinknockout) << endl;
  os << " extension_only_chemotaxis = " << sbool(extension_only_chemotaxis) << endl;
  os << " chemotaxis = " << chemotaxis << endl;
  os << " border_energy = " << border_energy << endl;
  os << " neighbours = " << neighbours << endl;
  os << " periodic_boundaries = " << sbool(periodic_boundaries) << endl;
  os << " n_chem = " << n_chem << endl;
  os << " diff_coeff = "<< diff_coeff[0] << endl;
  os << " decay_rate = "<< decay_rate[0] << endl;
  os << " secr_rate = "<< secr_rate[0] << endl;
  os << " saturation = " << saturation << endl;
  os << " Sdf1UptakeNC = " << Sdf1UptakeNC << endl;
  os << " pde_initialization = " << pde_initialization << endl;
  os << " n_init_cells = " << n_init_cells << endl;
  os << " n_NC = " << n_NC << endl;
  os << " size_init_cells = " << size_init_cells << endl;
  os << " sizex = " << sizex << endl;
  os << " sizey = " << sizey << endl;
  os << " divisions = " << divisions << endl;
  os << " mcs = " << mcs << endl;
  os << " rseed = " << rseed << endl;
  os << " subfield = " << subfield << endl;
  os << " relaxation = " << relaxation << endl;
  os << " init_cell_distribution = " << init_cell_distribution << endl;
  os << " sidePadding = " << sidePadding << endl;
  os << " init_cluster_distance = " << init_cluster_distance << endl;
  os << " init_PL_distance = " << init_PL_distance << endl;
  os << " takeoutCellsBelow = " << takeoutCellsBelow << endl;
  os << " storage_stride = " << storage_stride << endl;
  os << " graphics = " << sbool(graphics) << endl;
  os << " store = " << sbool(store) << endl;
  os << " storeNCimage = " << sbool(storeNCimage) << endl;
  os << " storePLimage = " << sbool(storePLimage) << endl;

  if (datadir) 
    os << " datadir = " << datadir << endl;
}


ostream &operator<<(ostream &os, Parameter &p) {
  p.Write(os);
  return os;
}

Parameter par;
