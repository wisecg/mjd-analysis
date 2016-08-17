include $(MGDODIR)/buildTools/config.mk

# Give the application name, which must be the stem of cc file with 'main'.
APPS =  reportQuality

# The next two lines are important
SHLIB =
ARCHIVE =
TAMDIR ?= $(ROOTSYS)

# Set flags

INCLUDEFLAGS = -I. $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms

INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS) -I$(TAMDIR)/inc -I$(TAMDIR)/include -I$(MGDODIR)/Tabree -Wno-unused-parameter

INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis -I$(GATDIR)/DCProcs

LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOMJDB -lXMLParser -v -lMGDOTransforms -lMGDOTabree

LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATMGTEventProcessing -lGATMGOutputMCRunProcessing -lGATAnalysis -lGATMJDAnalysis -lGATDCProcs $(ROOT_LIB_FLAGS) -lSpectrum -lTreePlayer -L$(TAMDIR)/lib -lTAM

include $(MGDODIR)/buildTools/BasicMakefile

