include $(MGDODIR)/buildTools/config.mk

# Give the list of applications, which must be the stems of cc files with 'main'.
APPS = spec-fit

# The next three lines are important
SHLIB =
ARCHIVE =
SOURCESSCRATCH = $(wildcard *.cc)
TAMDIR ?= $(ROOTSYS)
# Include the correct flags,
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms
INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS) -I$(TAMDIR)/inc -I$(TAMDIR)/include -I$(MGDODIR)/Tabree
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis -I$(GATDIR)/DCProcs -I$(GATDIR)/Apps
LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOMJDB -lMGDOTabree
LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATMGTEventProcessing -lGATMGOutputMCRunProcessing -lGATAnalysis -lGATMJDAnalysis -lGATDCProcs $(ROOT_LIB_FLAGS) -lSpectrum -lTreePlayer -lMinuit -L$(TAMDIR)/lib -lTAM

include $(MGDODIR)/buildTools/BasicMakefile

