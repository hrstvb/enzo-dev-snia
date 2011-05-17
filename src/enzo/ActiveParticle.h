/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/  modified2:  John Wise, Greg Bryan, Britton Smith, Cameron Hummels,
/              Matt Turk
/  date:       May, 2011 (converting from Star to ActiveParticle)
/
/  PURPOSE:
/
************************************************************************/
#ifndef __ACTIVE_PARTICLE_H
#define __ACTIVE_PARTICLE_H

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "StarBuffer.h"

struct ActiveParticleFormationData;
struct ActiveParticleFormationDataFlags;

class ActiveParticleType
{
    void static ConstructData(grid *_grid,
            ActiveParticleFormationDataFlags &flags,
            ActiveParticleFormationData &data);
    void static DestroyData(
            ActiveParticleFormationDataFlags &flags,
            ActiveParticleFormationData &data) { return; }
#ifdef ACTIVE_PARTICLE_IMPLEMENTED
  public:
    void static ConstructData(grid *_grid,
            ActiveParticleFormationDataFlags &flags,
            ActiveParticleFormationData &data);
    /* Several pure virtual functions */

    /* This should return the number of new star particles created, and should
     * create them. */
    ActiveParticleType();
    ActiveParticleType(grid *_grid, int _id, int _level);
    ActiveParticle(StarBuffer *buffer, int n);
    ActiveParticle(StarBuffer buffer) ;
    ~ActiveParticle();
    ActiveParticle* copy(void);

  int   ReturnID(void) { return Identifier; };
  double ReturnMass(void) { return Mass; };
  float ReturnBirthTime(void) { return BirthTime; };

  int   ReturnLevel(void) { return level; };
  void  ReduceLevel(void) { level--; };
  void  IncreaseLevel(void) { level++; };
  void  SetLevel(int i) { level = i; };
  void  SetGridID(int i) { GridID = i; };
  grid *ReturnCurrentGrid(void) { return CurrentGrid; };
  void  AssignCurrentGrid(grid *a) { this->CurrentGrid = a; };
  void  AddMass(double dM) { Mass += dM; };

  FLOAT *ReturnPosition(void) { return pos; }
  float *ReturnVelocity(void) { return vel; }
  void    ConvertMassToSolar(void);
  void    Merge(Star a);
  void    Merge(Star *a);
  bool    Mergable(Star a);
  bool  Mergable(Star *a);
  float Separation(Star a);
  float Separation(Star *a);
  float Separation2(Star a);
  float Separation2(Star *a);
  float RelativeVelocity2(Star a);
  float RelativeVelocity2(Star *a);
  void  UpdatePositionVelocity(void);
  void    DeleteCopyInGrid(void);
  int   DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[]);
  void    CopyToGrid(void);
  void  MirrorToParticle(void);
  virtual bool  IsARadiationSource(FLOAT Time) { return FALSE };
  int   DeleteParticle(LevelHierarchyEntry *LevelArray[]);
  int   DisableParticle(LevelHierarchyEntry *LevelArray[]);
  void  ActivateNewStar(FLOAT Time, float Timestep);

  int SphereContained(LevelHierarchyEntry *LevelArray[], int level, 
              float Radius);

  ActiveParticle* StarBufferToList(StarBuffer *buffer, int n);
  StarBuffer* StarListToBuffer(int n);
    
  protected:
    grid        *CurrentGrid;
    FLOAT     pos[MAX_DIMENSION];
    float         vel[MAX_DIMENSION];
    double            Mass;        // Msun
    float         BirthTime;
    PINT         Identifier;
    int         level;
    int         GridID;
    star_type     type;

    bool Active;

  private: /* Cannot be accessed by subclasses! */
  
#endif /* IMPLEMENTED */
  friend class grid;
  friend class ActiveParticleType_info;

};

struct ActiveParticleFormationData {
  /* This is where all the pointers that normally get passed into
   * formation routines gets placed. Things like fractional h2, dark
   * matter density, etc etc. Anything that's derived.  It's okay to
   * add to this.  */
  float *DarkMatterDensity;
  float *H2Fraction;
  float *CoolingTime;
  float *Temperature;
  float *TotalMetals;
  float DensityUnits;
  float LengthUnits;
  float TemperatureUnits;
  float TimeUnits;
  float VelocityUnits;
  double MassUnits;
  int DensNum;
  int Vel1Num;
  int Vel2Num;
  int Vel3Num;
  int MetalNum;
  int ColourNum;
  int level;
};

struct ActiveParticleFormationDataFlags {
  /* For every entry in the ActiveParticleFormationData struct, we
   * have a bool here. */
  bool DarkMatterDensity;
  bool H2Fraction;
  bool CoolingTime;
  bool Temperature;
  bool UnitConversions;
  bool DataFieldNumbers;
  bool MetalField;
};

const struct ActiveParticleFormationDataFlags flags_default = {
  false,    // DarkMatterDensity
  false,    // H2Fraction
  false,    // CoolingTime
  false,    // Temperature
  false,    // UnitConversions
  false,    // DataFieldNumbers
  false     // MetalField
};


//! maps the name of a plug-in to a pointer of the factory pattern
class ActiveParticleType_info;
typedef std::map<std::string, ActiveParticleType_info *> ActiveParticleMap;

ActiveParticleMap &get_active_particle_types();

void EnableActiveParticleType(char *active_particle_type_name);

class ActiveParticleType_info
{
    public:
       
       /* We will add more functions to this as necessary */
       ActiveParticleType_info(
           std::string this_name,
           int (*ffunc)(grid *thisgrid_orig),
           void (*dfunc)(ActiveParticleFormationDataFlags &flags)
           ){
        this->formation_function = ffunc;
        this->describe_data_flags = dfunc;
        get_active_particle_types()[this_name] = this;
       }

       static int count(){return get_active_particle_types().size();}

       int (*formation_function)(grid *thisgrid_orig);
       void (*describe_data_flags)(ActiveParticleFormationDataFlags &flags);
};

#endif

