include $(MGDODIR)/buildTools/config.mk

# Give the application name, which must be the stem of cc file with 'main'.  
APPS =  mjd_channel_quality
# The next two lines are important
SHLIB = 
ARCHIVE = 
TAMDIR ?= $(ROOTSYS)
# Include the correct flags, 
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms -I$(MGDODIR)/Tabree -I$(MGDODIR)/MJDB
INCLUDEFLAGS += -I$(MGDODIR)/Majorana $(ROOT_INCLUDE_FLAGS) -I$(TAMDIR)/inc -I$(TAMDIR)/include
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis
LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOTabree -lMGDOMJDB
LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATMGTEventProcessing -lGATMGOutputMCRunProcessing -lGATAnalysis -lGATMJDAnalysis $(ROOT_LIB_FLAGS) -lSpectrum -L$(TAMDIR)/lib -lTAM

include $(MGDODIR)/buildTools/BasicMakefile

