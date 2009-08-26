Eint32 enzo_main(Eint32 argc, char** argv);
int SetDefaultGlobalValues(TopGridData &MetaData);

#ifdef PARSE_GLOBAL_DATA
#include "svn_version.def"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "communication.h"
#include "CommunicationUtilities.h"
#ifdef TRANSFER
#include "PhotonCommunication.h"
#endif
#endif
