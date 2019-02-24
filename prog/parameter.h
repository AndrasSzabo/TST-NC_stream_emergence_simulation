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
#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <iostream>
using namespace std;
class Parameter {
  
 public: 
  Parameter();
  ~Parameter();
  void CleanUp(void);
  void Read(const char *filename);
  void Write(ostream &os) const;
  double T;
  int target_area;
  int target_length;
  double lambda;
  double lambda2;
  double colour_cells;
  int colour_active_cells;
  bool show_time;

  double Polarisation;
  double Pol_decay;
  double Pol_decayFree;
  bool PolDynamicUpdateFree, PolDynamicUpdateContact;
  double PolarisationPlacode;
  double Pol_decayPlacode;
  double Pol_decayFreePlacode;
  
  double * CILprob;
  double repolarisationPM;
  int repolarisationDelay;
  int plotCILcells;
  int plotChemotaxCells;

  bool PassiveRelaxationSteps;
  int PassiveRelaxationWindow;
  bool PassiveRelaxationMonitorDH;
  
  bool EdenGrowth;
  int nJ;
  double ** J;
  double ** relaxEM;
  double ** CC;
  double ** UCC;
  double ** SC;
  bool cytosk_bonds;
  int plotCytoskBonds;
  bool breakBondIfNoContact;
  
  int plotPolarityVectors;
  double polVecPlotLen;
  double insertMediumProb;
  double conn_diss;
  bool vecadherinknockout;
  bool extension_only_chemotaxis;
  double chemotaxis;
  double * chemotaxisNC;
  double * chemotaxisPlacode;
  double border_energy;
  int neighbours;
  bool periodic_boundaries;
  int n_chem;
  double * diff_coeff;
  double * decay_rate;
  double * secr_rate;
  double saturation;
  double Sdf1UptakeNC;
  bool pdeExtendedField;
  int pde_initialization;
  int n_init_cells;
  int n_NC;
  int size_init_cells;
  int sizex;
  int sizey;
  int divisions;
  int mcs;
  int rseed;
  double subfield;
  bool meas_dens;
  int relaxation;
  int internalrelaxation;
  int init_cell_distribution;
  int sidePadding;
  int init_cluster_distance;
  double init_PL_distance;
  double takeoutCellsBelow;
  int storage_stride;
  int print_sigma;
  int print_CM;
  int print_CMstart;
  int print_CMsample;
  int print_CMfreq;
  int print_asParticles;
  int print_PIF;
  bool print_CMarea;
  int print_Areas;
  bool print_AreasB;
  int time;
  bool print_sigmaB;
  bool print_CMB;
  bool graphics;
  double * plotcc;
  double * plotccCol;
  double * plotccgrad;
  double * storeccgrad;
  double * storecc;
  bool printcc;
  bool store;
  bool storeNCimage;
  bool storePLimage;
  int store_start;
  char * datadir;
 private:
};

ostream &operator<<(ostream &os, Parameter &p);
const char *sbool(const bool &p);


#endif
