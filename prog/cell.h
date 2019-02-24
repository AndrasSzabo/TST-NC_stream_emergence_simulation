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
#ifndef _CELL_HH_
#define _CELL_HH_

#include "parameter.h"
//#define EMPTY -1
#include <math.h>

extern Parameter par;
class Dish;

class Cell
{ 

  friend class Dish;
  friend class CellularPotts;
  friend class Info;

public:
  /*! \brief Constructor to insert a cell into Dish "who"
    
  Used to add a new Cell to the dish: new Cell(dish,
  celtype).
  */
  
  void CenterOfMass(double *cx, double *cy);
  void Theta(double *theta); 
  void getLongAxisDir (double *vx, double *vy);
  void setDisplacement (double dx, double dy);
  void getDisplacement (double *dx, double *dy);
  void getPolarisation (double *px, double *py);
  void setPolarisation (double px, double py);
  double getPolStr (void);
  void setPolStr (double p);
  double getPolDecay (void);
  double getPolDecayContact (void);
  double getPolDecayFree (void);
  void setPolDecay (double p);
  void setPolDecayContact (double p);
  void setPolDecayFree (double p);
  int getEcmContact (void);
  void setEcmContact (int e);
  void initCellParameters (void);
  void initContacts (void);

  Cell(const Dish &who, int settau=1) 
  {    
    owner=&who;
    ConstructorBody(settau);
  }

  Cell(void) {
    chem = new double[par.n_chem];
    contactDurationWith = new int[par.n_init_cells+1];
    cytosk_bonds = new double[par.n_init_cells+1];
  };
  
  ~Cell(void);
  
  //! Default copy constructor.
  Cell(const Cell &src) {
        
    // make an exact copy (for internal use)
    sigma=src.sigma;
    amount++;
    area=src.area;
    target_area=src.target_area;
    length=src.length;
    target_length=src.target_length;
    growth_threshold=src.growth_threshold;
    mother=src.mother;
    daughter=src.daughter;
    times_divided=src.times_divided;
    date_of_birth=src.date_of_birth;
    colour_of_birth=src.colour_of_birth;
    tau=src.tau;
    alive=src.alive;
    polx=src.polx;
    poly=src.poly;
    PolStr=src.PolStr;
    PolDecay=src.PolDecay;
    PolDecayFree=src.PolDecayFree;
    PolDecayContact=src.PolDecayContact;
    displacement[0]=src.displacement[0];
    displacement[1]=src.displacement[1];
    n_copies=src.n_copies;
    sum_x=src.sum_x;
    sum_y=src.sum_y;
    sum_xx=src.sum_xx;
    sum_yy=src.sum_yy;
    sum_xy=src.sum_xy;
    owner=src.owner;
    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];    
    ecmcontact = src.ecmcontact;
    contactDurationWith = new int[par.n_init_cells+1];
    cytosk_bonds = new double[par.n_init_cells+1];
    for (int ch=0;ch<=par.n_init_cells;ch++){
      contactDurationWith[ch]=src.contactDurationWith[ch];    
      cytosk_bonds[ch] = src.cytosk_bonds[ch];
    }
  }
  
    /*! \brief Special copy constructor to add a new cell to the dish.
   
       Call it as: new Cell(parent, true); mother will be modified for
       ancestry administration!  

       \param settau 
       Cell type of daughter cell.
    */
  Cell(Cell &mother, int settau=1);
 
    /*! \brief Assignment operator.

    Called if one cell is assigned to another. Remember to change both
    assignment operator and copy constructor when adding new attributes
    to Cell.
    */
  inline Cell &operator=(const Cell &src) {
    colour=src.colour;
    alive=src.alive;
    sigma=src.sigma;
    area=src.area;
    tau=src.tau;
    target_area=src.target_area;
    polx=src.polx;
    poly=src.poly;
    PolStr=src.PolStr;
    PolDecay=src.PolDecay;
    PolDecayFree=src.PolDecayFree;
    PolDecayContact=src.PolDecayContact;
    displacement[0]=src.displacement[0];
    displacement[1]=src.displacement[1];
    n_copies=src.n_copies;
    
    sum_x=src.sum_x;
    sum_y=src.sum_y;
    sum_xx=src.sum_xx;
    sum_yy=src.sum_yy;
    sum_xy=src.sum_xy;
    
    length=src.length;
    target_length=src.target_length;
    amount++;
fprintf(stderr, "# Cell assignment\n");
    owner=src.owner;
    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];
    ecmcontact = src.ecmcontact;
    contactDurationWith = new int[par.n_init_cells+1];
    cytosk_bonds = new double[par.n_init_cells+1];
    for (int ch=0;ch<par.n_init_cells+1;ch++){
      contactDurationWith[ch]=src.contactDurationWith[ch];    
      cytosk_bonds[ch] = src.cytosk_bonds[ch];
    }
    
    return *this;
  }

    /*! \brief Returns false if Cell has apoptosed (vanished). */
  inline bool AliveP(void) const {
    return alive;
  }
  
    //! Returns the cell colour.
  inline int Colour(void) const {   
    return colour;
  };

    //! Set cell type of this Cell.
  inline void setTau(int settau) {
    tau=settau;
  }
  
    //! Get cell type of this Cell.
  inline int getTau(void) {
    return tau;
  }
    //! Set color of this cell to new_colour, irrespective of type.
  inline int SetColour(const int new_colour) {
    return colour=new_colour;
  }

    /* \brief Returns the energy between this cell and cell2. 

    Called from CellularPotts::DeltaH.
    */
  double EnergyDifference(const Cell &cell2) const;

    //! Return Cell's actual area.
  inline int Area() const {
    return area;
  }
  
    //! Return Cell's target area.
  inline int TargetArea() const {
    return target_area;
  }
  
    /*! \brief Return Cell's target length.
    
    Length constraint is documented in Merks et al. 2006, Dev. Biol. 
    */
  inline double TargetLength() const {
    return target_length;
  }

    //! Set the Cell's target length
  inline double SetTargetLength(double l) {
    return target_length=l;
  }
  

    //! Debugging function used to print the cell's current inertia tensor (as used for calculations of the length )
  inline void PrintInertia(void) {
    
    double ixx=(double)sum_xx-(double)sum_x*sum_x/(double)area;
    double iyy=(double)sum_yy-(double)sum_y*sum_y/(double)area;
    double ixy=(double)sum_xy-(double)sum_x*sum_y/(double)area;

    cerr << "ixx = " << ixx << "\n";
    cerr << "iyy = " << iyy << "\n";
    cerr << "ixy = " << ixy << "\n";
    
  }

    // return the current length
  inline double Length(void) {
    return length;
  }

    /*! \brief Clears the table of J's.  

    This is only important for a
    feature called "DynamicJ's", where J-values depend on internal states
    of the cells (such as a genetic network; see e.g. Hogeweg et
    al. 2000). The current version of TST does not include such functionality.
    */
  static void ClearJ(void);
  double polarvec[9]; 
  void RenormPolarVec(void);
  
     /*! \brief Returns the maximum cell identity number in the Dish.
       This would normally be the number of cells in the Dish, although
       the number includes apoptosed cells.
     */
  static inline int MaxSigma() {
    return maxsigma;
  }

     //! Returns the cell's cell identity number.
  inline int Sigma() const {
    return sigma;
  }
  
     //! Sets the target area of the cell.
  inline int SetTargetArea(const int new_area) {
    return target_area=new_area;
  }
  
     //! Sends the current cell into apoptosis
  inline void Apoptose() {
    alive=false;
  }
     //! Increment the cell's target area by one unit.
  inline int IncrementTargetArea() {
    return ++target_area;
  }
  
     //! Decrement the cell's target area by one unit.
  inline int DecrementTargetArea() {
    return --target_area;
  }

    //! Cell lineage tracking: get the cell's parent
  inline int Mother(void) const { return mother; }
  
    //! Cell lineage tracking: get the cell's daughter
  inline int Daughter(void) const { return daughter; }
  
    //! Returns a counter keeping track of the number of divisions
  inline int TimesDivided(void) const { return times_divided; }
  
    //! Returns Monte Carlo Step (MCS) when this cell originated.
  inline int DateOfBirth(void) const { return date_of_birth; }
  
    //! Returns the cell type at the time of birth. 
  inline int ColourOfBirth(void) const { return colour_of_birth; }

    //! Returns the bond energy J between this cell and cell c2.
  inline double GetJ(const Cell &c2) const {
    return J[sigma][c2.sigma];
  }

    //! Sets bond energy J between cell type t1 and t2 to val
  inline static double SetJ(int t1,int t2, int val) {
    return J[t2][t1]=J[t1][t2]=(double)val;
  }

  
  //! Return (or set) chemotaxis coefficient 'chem' for the cell
  inline double getChem(int l){
    return chem[l];
  }

  inline double setChem(int l, double v){
    chem[l]=v;
    return chem[l];
  }

    
     /*! After introducing a new Cell (e.g. with GrowInCell)
         call this function to set the moments and areas right.
      */
  void MeasureCellSize(Cell &c);
     //! End of public declarations.

private:
     /*! \brief Read a table of static Js.
         First line: number of types (including medium)
         Next lines: diagonal matrix, starting with 1 element (0 0)
         ending with n elements */
  static void ReadStaticJTable();

     // used internally by dish in "CellGrowthAndDivision"
  inline int GrowthThreshold(void) const { return growth_threshold; }
  
     // used internally by class CellularPotts
  inline void CleanMoments(void) {
    sum_x = sum_y = sum_xx = sum_xy = sum_yy = area = target_area  = 0;
  }
     // used internally by class CellularPotts
  inline double AddSiteToMoments(int x,int y, double new_l=-1. ) {
    
     // Add a site to the raw moments, then update and return the
     // length of the cell
    
     // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
     // Eventually this function may be used to carry
     // out all necessary adminstration at once
    sum_x+=x;
    sum_y+=y;
    sum_xx+=x*x;
    sum_yy+=y*y;
    sum_xy+=x*y;
    
       // update length (see appendix. A, Zajac.jtb03), if length is not given
       // NB. 24 NOV 2004. Found mistake in Zajac's paper. See remarks in
       // method "Length(..)". 
    if (new_l<0.) {
      length=Length(sum_x,sum_y,sum_xx,sum_yy,sum_xy,area);
    } else {
      length=new_l;
    }
    return length;
  }

    // used internally by class CellularPotts
  inline double RemoveSiteFromMoments(int x,int y, double new_l=-1.) {
    
      // Remove a site from the raw moments, then update and return the
      // length of the cell
    
      // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
      // Eventually this function may be used to carry
      // out all necessary adminstration at once
    sum_x-=x;
    sum_y-=y;
    sum_xx-=x*x;
    sum_yy-=y*y;
    sum_xy-=x*y;
    
      // update length (see app. A, Zajac.jtb03), if length is not given
    if (new_l<0.) {
      length=Length(sum_x,sum_y,sum_xx,sum_yy,sum_xy,area);
    } else {
      length=new_l;
    }
    return length;
  }

    //! \brief Calculates the length based on the given inertia tensor
    //components (used internally)
  inline double Length(int s_x,int s_y,int s_xx,
				int s_yy,int s_xy,int n) {
    
    if (n==0) return (0);
    else{
        // inertia tensor (constructed from the raw momenta, see notebook)
      double iyy=(double)s_xx-(double)s_x*s_x/(double)n;
      double ixx=(double)s_yy-(double)s_y*s_y/(double)n;
      double ixy=-(double)s_xy+(double)s_x*s_y/(double)n;
        
      double rhs1=(ixx+iyy)/2., rhs2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )/2.;

      double lambda_b=rhs1+rhs2;
        //double lambda_a=rhs1-rhs2;
    
        // according to Zajac et al. 2003:
        //return 2*sqrt(lambda_b);
        // Grumble, this is not right!!!
        // Must divide by mass!!!!!!

        // see: http://scienceworld.wolfram.com/physics/MomentofInertiaEllipse.html
        //    cerr << "n = " << n << "\n";
      return 4*sqrt(lambda_b/n);

        // 2*sqrt(lambda_b/n) give semimajor axis. We want the length.
    };
  }
  
    // return the new length that the cell would have
    // if site (x,y) were added.
    // used internally by CellularPotts
  inline double GetNewLengthIfXYWereAdded(int x, int y) {
    
    return Length(sum_x+x,sum_y+y,
			   sum_xx+x*x,sum_yy+y*y,sum_xy+x*y,area+1);    
  }

    // return the new length that the cell would have
    // if site (x,y) were removed
    // used internally by CellularPotts
  inline double GetNewLengthIfXYWereRemoved(int x, int y) {

        return Length(sum_x-x,sum_y-y,
			       sum_xx-x*x,sum_yy-y*y,sum_xy-x*y,area-1);
  }
    //! End of private declarations.
private:
    //! Increments the cell's actual area by 1 unit.
  inline int IncrementArea() {
    return ++area;
  }

    //! Decrements the cell's actual area by 1 unit.
  inline int DecrementArea() {
    return --area;
  }

  
    /*! \brief Sets target area to actual area, to remove "pressure".
        This is useful when reading an initial condition from an image.
    */
  inline int SetAreaToTarget(void) {
    return area=target_area;
  }

    //! Called whenever a cell is constructed, from constructor
  void ConstructorBody(int settau=1);
  
    // returns the maximum cell type index
  static int MaxTau(void) {
    return maxtau;
  }


    //! End of private declarations.
protected:
  int colour;
  bool alive;
  int sigma;   // cell identity, 0 if medium
  int tau;     // Cell type, when dynamicJ's are not used
  

  double length;         // length of the cell;
  double target_length;
    // Two dimensional (square) array of ints, containing the J's.
    // Dynamically increased when cells are added to the system
    // unless a static Jtable is used (currently this is the default situation)
  static double **J;
  static int maxtau;

    // Amount: the number of Cell instantations, INCLUDING copies
    // For internal use only.
    // Reading amount is NOT the way to get the number of cells!!
  static int amount;
  static int capacity;
  static int maxsigma;  // the last cell identity number given out
  

    // indices of mother and daughter
    // (Note: no pointers, cells may be relocated)
  int mother;
  int daughter;
  int times_divided;
  int date_of_birth;
  int colour_of_birth;
  
  int area;
  double polx, poly, PolStr, PolDecay, PolDecayFree, PolDecayContact;
  double displacement[2];
  int inContact;
  int * contactDurationWith;
  double * cytosk_bonds;
  int CILstate;
  int ecmcontact;
  int ActivelyMoving;
  bool chemotaxing;
   
  int target_area;
  int growth_threshold;

  int n_copies;   // number of expansions of this cell
      // gradient of a chemical (to be extended to the total number chemicals)
  double grad[2];  
  double *chem;
    // Raw moments of the cells
    // Are used to calculate minor and major axes
    // and center of mass
    // are locally adjusted, so axes are easily
    // and quickly calculated!  
    // N.B: N is area!  
  int sum_x;
  int sum_y;
  int sum_xx;
  int sum_yy;
  int sum_xy;
  
  const Dish *owner; // pointer to owner of cell
};

#endif
