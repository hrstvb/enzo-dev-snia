/***********************************************************************
/
/  CREATE CUSTOM MPI DATA TYPES
/
/  written by: John Wise
/  date:       February, 2011
/  modified1:  
/
***********************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "GroupPhotonList.h"
#include "Parallel.h"

namespace Parallel {

  int CreateMPITypes(void)
  {
#ifdef USE_MPI
    /* 1. MPI Buffer header */

    int i;
    mpi_header header;
    MPI_Datatype HeaderType;
    MPI_Datatype type[6] = {MPI_CHAR, IntDataType, IntDataType, IntDataType,
			    FLOATDataType, IntDataType};
    int blocklen[6] = {1, 2, 3, 3, 3, 3};
    MPI_Aint disp[6];
    MPI_Aint base;
    
    /* Compute the displacements inside the struct, then create header
       type */

    MPI_Address(&header, disp);
    MPI_Address(header.GridNum, disp+1);
    MPI_Address(header.GridOffset, disp+2);
    MPI_Address(header.GridDims, disp+3);
    MPI_Address(header.FArg, disp+4);
    MPI_Address(header.IArg, disp+5);
    base = disp[0];
    for (i = 0; i < 6; i++) disp[i] -= base;
    MPI_Type_struct(6, blocklen, disp, type, &HeaderType);
    MPI_Type_commit(&HeaderType);

    MPI_Header = HeaderType;

    /************************************************************************
      For most of the following types, we could construct a struct
      for heterogeneous systems, but for now, keep it as bytes.
    ************************************************************************/

    /* 2. Star object buffer. */

    MPI_Type_contiguous(sizeof(StarBuffer), MPI_BYTE, 
			&MPI_StarBuffer);
    MPI_Type_commit(&MPI_StarBuffer);

    /* 3. Particle entry for merging (CommunicationMergeStarParticle) */

    MPI_Type_contiguous(sizeof(ParticleEntry), MPI_BYTE,
			&MPI_ParticleEntry);
    MPI_Type_commit(&MPI_ParticleEntry);

    /* 4. Packed grid entry (CommunicationShareGrids) */

    MPI_Type_contiguous(sizeof(PackedGrid), MPI_BYTE,
			&MPI_PackedGrid);
    MPI_Type_commit(&MPI_PackedGrid);

    /* 5. Particle sharing list (CommunicationShareParticles) */

    MPI_Type_contiguous(sizeof(particle_data), MPI_BYTE,
			&MPI_ParticleShareList);
    MPI_Type_commit(&MPI_ParticleShareList);

    /* 6. Particle move list (CommunicationTransferParticles) */

    MPI_Type_contiguous(sizeof(ParticleMoveList), MPI_BYTE,
			&MPI_ParticleMoveList);
    MPI_Type_commit(&MPI_ParticleMoveList);

    /* 7. Star share list (CommunicationShareStars) */

    MPI_Type_contiguous(sizeof(star_data), MPI_BYTE,
			&MPI_StarShareList);
    MPI_Type_commit(&MPI_StarShareList);

    /* 8. Star move list (CommunicationTransferStars) */

    MPI_Type_contiguous(sizeof(StarMoveList), MPI_BYTE,
			&MPI_StarMoveList);
    MPI_Type_commit(&MPI_StarMoveList);

    /* 9. Photon move list (CommunicationTransferPhotons) */
    
    MPI_Type_contiguous(sizeof(GroupPhotonList), MPI_BYTE,
			&MPI_PhotonList);
    MPI_Type_commit(&MPI_PhotonList);

    /* 10. "two_int" structure with a grid number and processor number
       (DepositParticleMassFlaggingField) */

    MPI_Type_contiguous(sizeof(two_int), MPI_BYTE, &MPI_TwoInt);
    MPI_Type_commit(&MPI_TwoInt);

    /* 11. Photon entry (Grid_CommunicationSendPhotonPackages) */

    MPI_Type_contiguous(sizeof(PhotonBuffer), MPI_BYTE,
			&MPI_PhotonBuffer);
    MPI_Type_commit(&MPI_PhotonBuffer);

#endif /* USE_MPI */

    return SUCCESS;

  }

} // END NAMESPACE
