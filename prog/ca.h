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

// mainpage.h contains no C++ code, it is for the main page of the
// documentation
#include "mainpage.h"

/*! Implementation of the Glazier & Graner cellular Potts model **/
#ifndef _CA_HH_
#define _CA_HH_
#include <vector>
#include <stdio.h>
#include "graph.h"
#include "pde.h"
//#include "dish.h"
#include "cell.h"

class Dish;

class Dir {
  
  /* To store a celldirection matrix */
  friend class CellularPotts;
public:

  Dir() {
    aa1=0.; aa2=0.;
    bb1=0.; bb2=0.;
    lb1=0.; lb2=0.;
  }

  double  aa1,aa2;
  double  bb1,bb2;
  double  lb1,lb2;
};

class CellularPotts {

  friend class Info;
  friend class Morphometry;
    
public:
  //! \brief Constructs a CA field. This should be done in "Dish".
  CellularPotts(std::vector<Cell> *cells, const int sizex=200, 
		const int sizey=200 );
  // empty constructor
  // (necessary for derivation)
  CellularPotts(void);

  // Keyword virtual means, that derived classed (cppvmCellularPotts) can override
  // this function and carry out the memory allocation in their preferred way
  // Every time AllocateSigma is called in the base class methods
  // the function belonging the actual type will be called
  virtual void AllocateSigma(int sx, int sy);
  
  // destructor must also be virtual
  virtual ~CellularPotts();

  /*! \brief Plots the dish to the screen or to a movie and searches the
   neighbours. 

   These distinct tasks have been lumped together in the
   same method because both for drawing the black lines between the
   cells and for searching the neighbours the cell borders have to be
   determined. */
  int **SearchNandPlot(Graphics *g=0, bool get_neighbours=true);

  /*! \brief Counts the neighbouring sites between every cell and, in a 
	very similar manner to SearchNandPlot function */
  int **CountN(Graphics *g=0, bool get_neighbours=true);
  
  //! Plot the dish to Graphics window g
  inline void Plot(Graphics *g) {
    SearchNandPlot(g, false);
  }
  
  //! Searches the cells' neighbors without plotting
  inline int **SearchNeighbours(void) {
    return SearchNandPlot(0, true);
  }

  //! Return the total area occupied by the cells
  inline int Mass(void) {
    int mass=0;
    for (int i=0;i<sizex*sizey;i++) {
      if (sigma[0][i]>0) mass++;
    }
    return mass;
  }

  /*! Plot the cells according to their cell identity, not their type.
    
  The black lines are omitted.
  */
  void PlotSigma(Graphics *g, int mag=2);
  
  /*! Plot only NCType cells without borders */
  void PlotNC(Graphics *g, int mag=2, int colour=2);
  /*! Plot only PLType cells without borders */
  void PlotPL(Graphics *g, int mag=2, int colour=2);

  
  //! Divide all cells.
    void DivideCells(void) {
	  std::vector<bool> tmp;
    DivideCells(tmp);
  }
  
    /*! Divide all cells marked "true" in which_cells.
      
    \param which_cells is a vector<bool> with the same number of
    elements as the number of cells. It is a mask indicating which
    cells should be divided; each cell marked true will be divided.
      
     If which_cells is empty, this method divides all cells.
    */
    void DivideCells(std::vector<bool> which_cells);
    
    //! Implements the core CPM algorithm
    void AmoebaeMove(PDE *PDEfield=0);
    void PlotCytoskBonds(Graphics* g, int colour);
    void PlotCILcells(Graphics* g, int colour);
    void PlotPolarityVectors(Graphics* g, double size, int colour);
  
    
    void ConstructInitCells(Dish &beast);
    
    //! Returns the number of completed Monte Carlo steps.
    inline int Time() const {
      return thetime;
    }
  

  //! \brief Return the horizontal size of the CA plane.
  inline int SizeX() const {
    return sizex;
  }
  
  //! \brief Return the vertical size of the CA plane.
  inline int SizeY() const {
    return sizey;
  }
  
  /*! \brief Return the value of lattice site (x,y).

  i.e. This will return the index of the cell which occupies site (x,y). */
  inline int Sigma(const int x, const int y) const {
    return sigma[x][y];
  }
  
  // Was used to make it possible to enlarge the Graphics window in
  // X11 and replace the contents interactively. Not currently supported.
  void Replace(Graphics *g);

  /*! In this method the principal axes of the cells are computed using
   the method described in "Biometry", box 15.5 
   \return a pointer to a "new[]"ed array containing the directions.
   The memory has to be freed afterwards using the delete[] operator
  */
  Dir *FindCellDirections(void) const;

  /*! \brief Initialize the CA plane with n circular cells fitting in
    a cellsize^2 square.*/
  int ThrowInCells(int n, int cellsize);

  /*! \brief Initialize the CA plane with n cells using an Eden growth algorithm.

  \param n: Number of cells.
  \param cellsize: Number of Eden growth iterations.
  \param subfield: Defines a centered frame of size (size/subfield)^2 in which all cell will be positioned. 
  \return Index of last cell inserted.
  */
  int GrowInCells(int n_cells, int cellsize, double subfield=1.);
  void EdenGrowth(int cell_size, int cellnum, int n_cells, int ** new_sigma);
  int GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y);
  
  //! \brief Adds a new Cell and returns a reference to it.
  inline Cell &AddCell(Dish &beast) {
    cell->push_back(Cell(beast));
    return cell->back();
  }
  /*! \brief Display the division planes returned by FindCellDirections.
    
  \param g: Graphics window
  \param celldir: cell axes as returned by FindCellDirections.
  */
  void ShowDirections(Graphics &g, const Dir *celldir) const;
  
  //! \brief Returns the mean area of the cells. 
  double MeanCellArea(void) const;
  
  /*! \brief Returns the cell density.

  Cell density is defined as the area occupied by cells divided by the size of the field.
  */
  double CellDensity(void) const; 
  
  //! \brief Set target lengths of all cells to the value given in parameter file.
  void ResetTargetLengths(void);
 
  int spins_converted;
  
  /*! \brief Give each cell a random cell type.
  */
  void SetRandomTypes(void);
  void SetNCandPlacodeTypes(int nNC);

  /*! Cells grow until twice their original target_length, then
    divide, with rate "growth_rate"
  */
  void GrowAndDivideCells(int growth_rate);
  
  inline Cell &getCell(int c) {
    return (*cell)[c];
  }
  
  double MatrixDecay(double cc);

private:
  
  void SpecInitConditions (void);
  void setColours (void);
  void setPolarisationParameters (void);
  void PolUpdate (void);
  void SyncCellTime (void);
  void Displace (int k, int x, int y);
  void PrintCM(void);
  void PrintAsParticles(void);
  void PrintAnis(void);
  void PrintECM(PDE * PDEfield);
  void PrintTAU(void);
  void PrintAreas(void);
  void PrintSigma(void);
  void PrintPIF(void);
  void PrintPIF(const char * fname);
  void IndexShuffle(void);
  double DeltaH(int x,int y, int xp, int yp, PDE *PDEfield=0, int i=1);
  void ConvertSpin(int x,int y,int xp,int yp);
  int CopyvProb(double DH,  double stiff);
  void FreezeAmoebae(void);
  void MeasureCellSizes(void);
  void MeasureCellSize(Cell &c);
  void CopyProb(double T);
  bool ConnectivityPreservedP(int x, int y);
  bool CheckConnectivityConstraints(int xp, int yp, int x, int y);
  double getEnergy ();
  void PassiveRelaxation(PDE *PDEfield=0);
  double PassiveRelaxationStep(PDE *PDEfield=0);
  void getCellCellAdhesionTensions(double *tension, int *n);
  void UpdateNeighbours(void);
  void takeoutCellsBelow(double limit);
  void InsertNCCellsAtTop();


  // little debugging function to print the site and its neighbourhood
  inline void PrintSite(int x,int y) {
	  std::cerr << "--------\n";
	  std::cerr << "[" << sigma[x-1][y-1] << " " << sigma[x][y-1] << " " << sigma[x+1][y-1] << "]\n";
	  std::cerr << "[" << sigma[x-1][y] << " " << sigma[x][y] << " " << sigma[x+1][y] << "]\n";
	  std::cerr << "[" << sigma[x-1][y+1] << " " << sigma[x][y+1] << " " << sigma[x+1][y+1] << "]\n";
  }

protected:
	void BaseInitialisation(std::vector<Cell> *cell);
  
protected:
//  double **cytosk_bonds;
  int **sigma;
  int sizex;
  int sizey;

  double *Epass;
  int *EpassN;
  double Edetail[100], EdetailSum[100];
  

private:
  bool frozen;
  static const int nx[25], ny[25];
  static const int nbh_level[5];
  static int shuffleindex[9];
  static double Tscale;
  std::vector<Cell> *cell;
  int thetime;
  int n_nb;
};


#endif
