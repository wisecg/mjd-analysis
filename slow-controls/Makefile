# ${MGDODIR}/MJDB/Apps/SlowControls/Makefile
include $(MGDODIR)/buildTools/config.mk

# Give the list of applications, which must be the stems of cc files 
# with 'main'.  There can be more than one.
APPS = lnFills 
# The next two lines are important
SHLIB = 
ARCHIVE = 

# Include the correct flags 
# Modified by Clint to access MGDO and GAT classes.

INCLUDEFLAGS += -I.
INCLUDEFLAGS += -Wno-unused-parameter -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS) 
INCLUDEFLAGS += -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms -I$(MGDODIR)/Majorana -I$(MGDODIR)/Tabree
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/DCProcs -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/MGTEventProcessing #-I$(GATDIR)/buildtools 

LIBFLAGS += -L$(MGDODIR)/lib -lMGDOTabree -lMGDOMJDB $(ROOT_LIB_FLAGS) -lXMLParser
LIBFLAGS += -L$(MGDODIR)/lib -lMGDOBase -lMGDOMajorana -lMGDORoot -lMGDOTransforms -lMGDOTabree
LIBFLAGS += -v
LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATDCProcs -lGATMGOutputMCRunProcessing -lGATMGTEventProcessing #-lGATbuildtools 

include $(MGDODIR)/buildTools/BasicMakefile
# ${MGDODIR}/MJDB/Apps/SlowControls/Makefile