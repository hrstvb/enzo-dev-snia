/// @file      AMRsolve_grid.h
/// @author    James Bordner
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Declaration of the AMRsolve_Grid class

#ifndef AMRSOLVE_GRID_H
#define AMRSOLVE_GRID_H

class AMRsolve_Grid
{

  /// @class    AMRsolve_Grid
  /// @brief    AMRsolve_Grid class for representing a grid patch in an AMR hierarchy
  ///
  /// A AMRsolve_Grid is an orthogonal grid of zones in \f$ R^d, 1 \le d \le 3 \f$.  Zones
  /// are congruent, and each zone contains a single double value at its
  /// center
  ///
  /// @verbatim
  ///         o---o---o---o---o---X
  ///         |   |   |   |   |   |
  ///         | * | * | * | * | * |
  ///         |   |   |   |   |   |
  ///         o---o---o---o---o---o
  ///         |   |   |   |   |   |
  ///         | * | * | * | * | * |
  ///         |   |   |   |   |   |
  ///         o---o---o---o---o---o
  ///         |   |   |   |   |   |
  ///         | * | * | * | * | * |
  ///         |   |   |   |   |   |
  ///         X---o---o---o---o---o
  ///  
  ///   o vertices
  ///   X extreme vertices xl_ and xu_
  ///   * unknowns
  ///   @endverbatim 

  friend class AMRsolve_Hierarchy;
  friend class ItGridNeighbors;
  friend class ItGridChildren;

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  int id_;         /// Unique id > 0
  int id_parent_;  /// Parent grid

  // data defined at AMRsolve_Grid creation (in AMRsolve_Grid())

  int     ip_;     /// Owning processor
  Scalar  xl_[3];  /// Position of lowest vertex
  Scalar  xu_[3];  /// Position of highest vertex
  int     il_[3];  /// Global index of lowest vertex
  int     n_[3];   /// Number of zones
  AMRsolve_Faces* faces_;  /// AMRsolve_Faces class associated with AMRsolve_Grid

  // data computed at hierarchy creation (in input())

  int     level_;           /// Containing Level number (0 = root)

  int     nu_[3];           /// Allocated size of solution u_
  Scalar* u_;               /// Unknowns (single cell-centered variable) 
                            ///   Stored as 3D fortran-style array
  int     offset_u_;        /// Offset of u_, such that u_=new[]+offset, or 
                            ///   u_-offset=new[]
  bool    is_u_allocated_;  /// Whether u_ was allocated by hypre-solve, or
                            ///   attached to an external array

  int     nf_[3];           /// Allocated size of right-hand side f_
  Scalar* f_;               /// Right-hand side (single cell-centered variable) 
                            ///   Stored as 3D fortran-style array
  int     offset_f_;        /// Offset of f_, such that f_=new[]+offset, or 
                            ///   f_-offset=new[]
  bool    is_f_allocated_;  /// Whether f_ was allocated by hypre-solve, or 
                            ///   attached to an external array

  int* counters_;           /// Counters for nonstencil entries.  Required by
                            ///   hypre to maintain state between nonstencil 
                            ///   and matrix initialization.
  
  // optional pointers to Enzo data arrays (if used; set using set_* routines)
  Scalar* E_;               /// ptr to Radiation energy density array
  Scalar* eta_;             /// ptr to emissivity array
  Scalar* HI_;              /// ptr to HI density array
  Scalar* HeI_;             /// ptr to HeI density array
  Scalar* HeII_;            /// ptr to HeII density array

  // optional integers describing grid relationship with Enzo data arrays
  int Ghosts_[3][2];        /// number of ghost zones per grid face

  //--------------------------------------------------------------------
  // STATIC MEMBER DATA
  //--------------------------------------------------------------------

  static AMRsolve_Mpi    mpi_;
  static AMRsolve_Domain domain_; // Only used by geomview stuff

  //--------------------------------------------------------------------
  // PROTECTED MEMBER DATA
  //--------------------------------------------------------------------

 protected:

  /// Vector of neighboring Grids, with extra sentinal 0-pointer at the end
  std::vector<AMRsolve_Grid *> neighbors0_; 
  /// Vector of child Grids
  std::vector<AMRsolve_Grid *> children0_;   


 public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  /// Create routine given the string "<id> <id_parent> <ip> <xl>[3] <xu>[3], <n>[3]"
  AMRsolve_Grid(std::string parms) throw();

  /// Create routine given the args <id> <id_parent> <ip> <xl>[3] <xu>[3], <n>[3]
  AMRsolve_Grid(int id, int ip_parent, int ip, Scalar* xl, Scalar* xu, 
		int* il, int* n) throw();

  /// Create a AMRsolve_Grid given a AMRsolve_Grid file as written by write()
  AMRsolve_Grid(std::string field, FILE*) throw();

  ~AMRsolve_Grid() throw();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  // IO

  /// Write the Grid definition to standard out in human-readable format
  void print() throw();

  /// Write the Grid to a file in compact format
  void write(std::string field, std::string name) throw()
  { FILE* fp = fopen(name.c_str(),"w"); this->write(field, fp); fclose(fp); };

  /// Write the Grid to a file in compact format
  void write(std::string field, FILE* fp, bool brief=0) throw();

  /// Read the Grid from a file written using write()
  void read(std::string field, std::string name) throw()
  { FILE* fp = fopen(name.c_str(),"r"); this->read(field, fp); fclose(fp); };

  /// Read the Grid from a file written using write()
  void read(std::string field, FILE* fp=0, bool brief=0) throw();

  /// Write the Grid outline to a geomview file 
  void geomview_grid(FILE* fpr, bool full=true) throw();

  /// Write the Grid's face data to a geomview file
  void geomview_face(FILE* fpr, bool full=true) throw();

  /// Write the Grid's face data of specified types to a geomview file
  void geomview_face_type(FILE* fpr, AMRsolve_Faces::Label* types, 
			  int num_types, bool full=true) throw();

  /// Return a pointer to the solution array associated with the Grid
  Scalar* get_u(int* nu0, int* nu1, int* nu2) throw();

  /// Set the solution array associated with the Grid
  void set_u(Scalar*, int dims[3]) throw();

  /// Return a pointer to the right-hand side array associated with the Grid
  Scalar* get_f(int* nf0, int* nf1, int* nf2) throw();

  /// Set the right-hand side array associated with the Grid
  void set_f(Scalar*, int dims[3]) throw();

  /// Return a pointer to the specified Enzo data array associated with the Grid
  /// (does not check for non-NULL value)
  Scalar* get_E()    throw() { assert(E_);    return E_;    }; 
  Scalar* get_eta()  throw() { assert(eta_);  return eta_;  };
  Scalar* get_HI()   throw() { assert(HI_);   return HI_;   };
  Scalar* get_HeI()  throw() { assert(HeI_);  return HeI_;  };
  Scalar* get_HeII() throw() { assert(HeII_); return HeII_; };

  /// Set a pointer to the specified Enzo data array associated with the Grid
  /// (does not check for non-NULL value)
  void set_E(Scalar* E)       throw() { E_    = E;    };
  void set_eta(Scalar* eta)   throw() { eta_  = eta;  };
  void set_HI(Scalar* HI)     throw() { HI_   = HI;   };
  void set_HeI(Scalar* HeI)   throw() { HeI_  = HeI;  };
  void set_HeII(Scalar* HeII) throw() { HeII_ = HeII; };

  /// Get/Set routines for number of ghost zones in an Enzo grid (compared to current Grid)
  void get_Ghosts(int Ghosts[3][2]) throw() 
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<2; j++)
	Ghosts[i][j] = Ghosts_[i][j];
  };
  void set_Ghosts(int Ghosts[3][2]) throw() 
  { 
    for (int i=0; i<3; i++)
      for (int j=0; j<2; j++)
	Ghosts_[i][j] = Ghosts[i][j];
  };


  /// Deallocate storage for the u_ and f_ arrays
  void deallocate() throw();

  /// Deallocate storage for the u_ array
  void deallocate_u_() throw();

  /// Deallocate storage for the f_ array
  void deallocate_f_() throw();

  /// Allocate storage for the array of values associated with the Grid
  void allocate() throw();

  /// Allocate storage for the u_ array
  void allocate_u_(int offset) throw();

  /// Allocate storage for the f_ array
  void allocate_f_(int offset) throw();

  /// Input the Grid from the given string in compact format
  void input(std::string parms) throw();

  void input(int id, int id_parent, int ip, Scalar* xl, Scalar* xu, 
	     int* il, int* n) throw();

  //--------------------------------------------------------------------
  // Data access
  //--------------------------------------------------------------------

  /// Return the Grid's integer id
  int id() throw() { return id_; }; 

  /// Return the Grid's parent's id
  int id_parent() throw() { return id_parent_; }; 

  /// Set the Grid to be a child.  Should only be called once per child.
  void set_child(AMRsolve_Grid& child) throw() 
  { children0_[children0_.size()-1] = & child;
    children0_.push_back (0); };

  /// Return the ith child.
  AMRsolve_Grid& child(int i) { return *children0_.at(i); };

  /// Return the number of child Grids
  int num_children() { return children0_.size() - 1; };

  /// Set the Grid to be a neighbor.  Should only be called once per neighbor.
  void set_neighbor(AMRsolve_Grid& neighbor) throw() 
  { neighbors0_[neighbors0_.size()-1] = &neighbor;
    neighbors0_.push_back (0); };

  /// Determine the axis, face, and range of indices of zones adjacent 
  /// to a neighboring Grid.  Returns false if the neighbor is not
  /// actually a neighbor.  Indices are relative to the Grid.
  bool neighbor_shared_face(AMRsolve_Grid& neighbor, int& axis, int& face, 
			    int& il0, int& il1, int& iu0, int& iu1,
			    int iperiod[3], int& count) throw();

  /// Determine the axis, face, and range of indices of zones adjacent 
  /// to a Grid in the next-coarser level.  Returns false if
  /// the Grid is not actually adjacent.  Assumes the adjacent Grid is
  /// in the next-coarser level.   Indices are relative to the Grid.
  bool coarse_shared_face(AMRsolve_Grid& coarse, int& axis, int& face, 
			  int& il0, int& il1, int& iu0, int& iu1,
			  int iperiod[3], int& count) throw();

  /// Determine the "count"th axis (indexing from 0), face and
  /// corresponding range of coarse-grid indices of zones adjacent to
  /// the containing parent Grid, and increment "count".  Returns true
  /// if the returned values are valid, or false if there is no
  /// "count"th face.   Indices are relative to the parent.
  bool parent_shared_face(AMRsolve_Grid& parent, int& axis, int& face, 
			  int& il0, int& il1, int& iu0, int& iu1,
			  int& count) throw();

  /// Determine the "count"th axis (indexing from 0), face and
  /// corresponding range of coarse-grid indices of zones adjacent to
  /// the interior of the parent Grid, and increment "count".  Returns true
  /// if the returned values are valid, or false if there is no
  /// "count"th face.   Indices are relative to the Grid.
  bool parent_interior_face(AMRsolve_Grid& parent, int& axis, int& face, 
			    int& il0, int& il1, int& iu0, int& iu1,
			    int& count) throw();

  /// Return the ith neighbor
  AMRsolve_Grid& neighbor(int i) { return *neighbors0_.at(i); };

  /// Return the number of neighbors
  int num_neighbors() const { return neighbors0_.size() - 1; };

  /// Make g1 and g2 neighbors.  Should only be called once per pair of neighbors
  friend void assert_neighbors(AMRsolve_Grid& g1, AMRsolve_Grid& g2) throw()
  { g1.set_neighbor(g2); 
    g2.set_neighbor(g1); };

  /// Set level for this Grid.  0 = root
  void set_level(int level) throw() { level_ = level; };

  /// Return the level for this Grid
  int level() throw() { return level_; };

  /// Set the global lower index of the lower Grid unknown.  No error checking on i.
  void set_lower(int i0, int i1, int i2) throw()
  { il_[0]=i0; il_[1]=i1; il_[2]=i2; };

  /// Return the number of unknowns along the ith coordinate.  No error checking on i.
  int num_unknowns(int i) throw() { return n_[i]; };

  /// Return the total number of unknowns.
  int num_unknowns() throw() { return n_[0]*n_[1]*n_[2]; };

  /// Return the coordinates of the lower Grid vertex.  No error checking on i.
  void x_lower(Scalar &xl0, Scalar &xl1, Scalar &xl2) throw()
  { xl0 = xl_[0]; 
    xl1 = xl_[1]; 
    xl2 = xl_[2]; };

  /// Return the coordinates of the upper Grid vertex.  No error checking on i.
  void x_upper(Scalar &xu0, Scalar &xu1, Scalar &xu2) throw()
  { xu0 = xu_[0]; 
    xu1 = xu_[1]; 
    xu2 = xu_[2]; };

  /// Return the global lower index of the lower Grid unknown.  No error checking on i.
  int index_lower(int i) throw() { return il_[i]; };

  /// Return the global lower indices of the lower Grid unknown.  No error checking on il.
  void index_lower(int &il0, int &il1, int &il2) throw()
  { il0 = il_[0];
    il1 = il_[1];
    il2 = il_[2]; };

  /// Return the global upper index of the upper Grid unknown.  No error checking on i.
  int index_upper(int i) throw() { return il_[i] + n_[i] - 1; };

  /// Return the global upper indices of the upper Grid unknown.  No error checking on iu.
  void index_upper(int &iu0, int &iu1, int &iu2) throw()
  { iu0 = il_[0] + n_[0] - 1;
    iu1 = il_[1] + n_[1] - 1;
    iu2 = il_[2] + n_[2] - 1; };

  /// Return lower and upper global indices.  Upper indices are incremented by one.
  void indices(int ind[3][2]) throw()
  { ind[0][0] = il_[0];
    ind[1][0] = il_[1];
    ind[2][0] = il_[2];
    ind[0][1] = il_[0] + n_[0];
    ind[1][1] = il_[1] + n_[1];
    ind[2][1] = il_[2] + n_[2]; };

  /// Return lower and upper global indices.  Upper indices are not incremented by one.  
  /// Used primarily before Hypre calls
  void get_limits(int lower[3], int upper[3]) throw()
  { lower[0] = il_[0];
    lower[1] = il_[1];
    lower[2] = il_[2];
    upper[0] = il_[0] + n_[0] - 1;
    upper[1] = il_[1] + n_[1] - 1;
    upper[2] = il_[2] + n_[2] - 1; };

  /// Return the Grid size.  Used primarily before Hypre calls
  void get_size(int size[3]) throw()
  { size[0] = n_[0];
    size[1] = n_[1];
    size[2] = n_[2]; };
  
  /// Return the AMRsolve_Faces object for this Grid.  If not allocated yet,
  /// create a new AMRsolve_Faces object.
  AMRsolve_Faces& faces() throw()
  { return *(faces_ ? faces_ : faces_=new AMRsolve_Faces(n_)); };

  /// Return the mesh width along the given axis
  Scalar h(int i) throw() { return (xu_[i] - xl_[i]) / n_[i]; };

  /// Return the mesh width along the given axis
  void h(Scalar &h0, Scalar &h1, Scalar &h2) throw()
  { h0 = (xu_[0] - xl_[0]) / n_[0];
    h1 = (xu_[1] - xl_[1]) / n_[1];
    h2 = (xu_[2] - xl_[2]) / n_[2]; };    

  /// Return the Grid size along the given axis
  int n(int i) throw()  { return n_[i]; };

  /// Return the Grid size
  int n() throw() { return n_[0]*n_[1]*n_[2]; };    

  /// Return the center of the given zone
  void zone(int i0, int i1, int i2, Scalar &x0, Scalar &x1, Scalar &x2) throw()
  { x0 = xl_[0] + ((xu_[0] - xl_[0]) / n_[0]) * (i0 + 0.5); 
    x1 = xl_[1] + ((xu_[1] - xl_[1]) / n_[1]) * (i1 + 0.5); 
    x2 = xl_[2] + ((xu_[2] - xl_[2]) / n_[2]) * (i2 + 0.5); };

  /// Processor owner
  int ip() throw() { return ip_; };

  //--------------------------------------------------------------------
  // Query functions
  //--------------------------------------------------------------------

  /// Return true iff the two Grids are next to each other, or are the same
  /// Grid.  Assumes the Grids are in the same level.
  bool is_adjacent(AMRsolve_Grid& grid, Scalar period[3]) throw();

  /// Return true iff the Grid belongs to processor ip
  bool is_local() throw() { return mpi_.ip() == ip_; };

  /// Return true iff the Grid contains the other Grid (i.e. is an ancestor)
  bool contains(AMRsolve_Grid* grid) throw();
  
  //--------------------------------------------------------------------
  // Nonstencil / Matrix initialization functions
  //--------------------------------------------------------------------

  /// Return alias to the (i3[0],i3[1],i3[2])th counter
  int& counter(int i3[3])
  { 
    int i0=i3[0]-il_[0];
    int i1=i3[1]-il_[1];
    int i2=i3[2]-il_[2];
    assert(counters_);
    return counters_[index(i0,i1,i2,n_[0],n_[1],n_[2])]; 
  }

  /// Initialize the counters_ array to given value
  void init_counter(int value)
  {
    for (int i2=0; i2<n_[2]; i2++) {
      for (int i1=0; i1<n_[1]; i1++) {
	for (int i0=0; i0<n_[0]; i0++) {
	  counters_[index(i0,i1,i2,n_[0],n_[1],n_[2])] = value;
	}
      }
    }
  }

  //--------------------------------------------------------------------
  // STATIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Initialize static AMRsolve_Mpi object
  static void set_mpi(AMRsolve_Mpi& mpi) { mpi_ = mpi; };

  /// Set static AMRsolve_Domain object
  static void set_domain(AMRsolve_Domain& domain) { domain_ = domain; };

  static int index(int i0,int i1,int i2,int n0,int n1,int n2) 
  { return i0 + n0*(i1+n1*i2); }

};


//----------------------------------------------------------------------


class ItGridNeighbors
{

  /// @class    ItGridNeighbors
  /// @date     2007-05-02
  /// @brief    Class for iterating over all neighboring Grids of a AMRsolve_Grid
  ///
  /// Example usage:
  ///
  /// @verbatim
  ///     ItGridNeighbors it_neighbors (*grid);
  ///
  ///     while ((AMRsolve_Grid * neighbor = it_neighbors++)) {
  ///        neighbor->print();
  ///     }
  /// @endverbatim

private:

  unsigned int curr_;
  const AMRsolve_Grid* grid_;

public:

  ItGridNeighbors(AMRsolve_Grid& grid) throw () : curr_(0), grid_(&grid) {};

  ~ItGridNeighbors() throw() {};
  
  /// Iterate through all Grids in the AMRsolve_Grid.
  AMRsolve_Grid* operator++(int) { 
    if (curr_ == grid_->neighbors0_.size()) curr_ = 0;
    curr_ ++;
    return grid_->neighbors0_[curr_-1];
  }

};


//----------------------------------------------------------------------


class ItGridChildren
{

  /// @class    ItGridChildren
  /// @date     2007-05-02
  /// @brief    Class for iterating over all child Grids of a AMRsolve_Grid
  ///
  /// Example usage:
  ///
  /// @verbatim
  ///     ItGridChildren it_children (*grid);
  ///
  ///     while ((AMRsolve_Grid * child = it_children++)) {
  ///        child->print();
  ///     }
  /// @endverbatim

private:

  unsigned int curr_;
  const AMRsolve_Grid* grid_;

public:

  ItGridChildren(AMRsolve_Grid& grid) throw() : curr_(0), grid_(&grid) {};

  ItGridChildren(const ItGridChildren& it) : curr_(it.curr_), grid_(it.grid_) {};

  void operator=(const ItGridChildren& it) { curr_ = it.curr_;  grid_ = it.grid_; };

  ~ItGridChildren() throw() {};
  
  /// Iterate through all Grids in the AMRsolve_Grid.
  AMRsolve_Grid* operator++(int) { 
    if (curr_ == grid_->children0_.size()) curr_ = 0;
    curr_ ++;
    return grid_->children0_[curr_-1];
  }

};


#endif /* AMRSOLVE_GRID_H */
