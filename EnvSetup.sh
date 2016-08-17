#!/bin/bash

# first two are for brew, they need to be in front.
export PATH="/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/Library/TeX/texbin"

alias pdsf="ssh -Y wisecg@pdsf.nersc.gov"
alias pdsf7="ssh -Y wisecg@pdsf7.nersc.gov"
alias pdsf8="ssh -Y wisecg@pdsf8.nersc.gov"
alias pdsf9="ssh -Y wisecg@pdsf9.nersc.gov"
alias pdsf10="ssh -Y wisecg@pdsf10.nersc.gov"
alias pdsf11="ssh -Y wisecg@pdsf11.nersc.gov"
alias ornl="ssh c2i@login1.ornl.gov"
alias vmini="ssh -Y wisecg@172.21.141.142"
alias rootmj="root -l ~/dev/rootlogon.C"
alias macspice="open ~/Applications/MacSpice/MacSpice.app/"
alias sub="sublime"
alias py2="python"
alias py3="python3"
alias tbr="rootmj -e \"TBrowser x\""
alias rsync="rsync -av --progress"
export PDSF=/global/homes/w/wisecg/
export DTN02=wisecg@dtn02.nersc.gov:/global/homes/w/wisecg
export PATH=${PATH}:/Users/wisecg/dev/vetoScan
source /Users/wisecg/ext/rw05/src/.radware.bashrc	# radware
export PATH=${PATH}:/Users/wisecg/ext/rw05/install/bin # radware
export PATH=${PATH}:/Users/wisecg/dev/channelSel/Cookies # mkcookie
export MJDDATADIR=/Users/wisecg/datasets/newModuleOne/

# ROOT 5 or ROOT 6 switch
# export ROOTSYS=/Users/wisecg/ext/ROOT/root_v5-34-00-patches
thisDir=$PWD
export ROOTSYS=/Users/wisecg/ext/ROOT/root_v6-06-00-patches
cd $ROOTSYS/bin
source thisroot.sh
cd $thisDir

# CLHEP
export CLHEP_BASE_DIR=/Users/wisecg/ext/CLHEP/2.3.1.0/CLHEP
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
export PATH=${PATH}:$CLHEP_BASE_DIR/bin
export DYLD_LIBRARY_PATH=$CLHEP_LIB_DIR:$DYLD_LIBRARY_PATH

# Python / PyRoot
export PYTHONDIR=/usr/lib/python2.7
export PYTHONPATH=${PYTHONPATH}:${HOME}/dev
export DYLD_LIBRARY_PATH=${ROOTSYS}/lib:${PYTHONDIR}/lib:${DYLD_LIBRARY_PATH}

# MGSW
export MGSWDIR=/Users/wisecg/ext
export TAMDIR=/Users/wisecg/ext/MGDO/tam
export MGDODIR=/Users/wisecg/ext/MGDO
export GATDIR=/Users/wisecg/ext/GAT
export ORDIR=/Users/wisecg/ext/OrcaRoot
export MJORDIR=/Users/wisecg/ext/MJOR
export PATH=${PATH}:$MGDODIR/bin:${GATDIR}/Apps:${GATDIR}/Scripts:${ORDIR}:${MJORDIR}:${MGDODIR}/lib
export DYLD_LIBRARY_PATH=${MGDODIR}/lib:${GATDIR}/lib:${ORDIR}/lib:${TAMDIR}/lib:${MJORDIR}:${ORDIR}/lib:${DYLD_LIBRARY_PATH}

# this is new in ROOT6
export ROOT_INCLUDE_PATH=$MGDODIR/Base:$MGDODIR/Gerda:$MGDODIR/GerdaTransforms:$MGDODIR/Majorana:$MGDODIR/MJDB:$MGDODIR/Root:$MGDODIR/Tabree:$MGDODIR/Tools:$MGDODIR/Transforms:$TAMDIR/inc:$MGDODIR/tam

