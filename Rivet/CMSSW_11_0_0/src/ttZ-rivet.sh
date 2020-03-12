#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc700
export RIVET_DATA_PATH=/afs/cern.ch/work/j/jthomasw/private/IHEP/LHCHXS-WG1/Rivet/CMSSW_11_0_0/src/GeneratorInterface/RivetInterface/data
export RIVET_REF_PATH=/afs/cern.ch/work/j/jthomasw/private/IHEP/LHCHXS-WG1/Rivet/CMSSW_11_0_0/src/GeneratorInterface/RivetInterface/data
export CMSSW_BASE=/afs/cern.ch/work/j/jthomasw/private/IHEP/LHCHXS-WG1/Rivet/CMSSW_11_0_0/src
pushd ${CMSSW_BASE}
echo ""
echo ${PWD}
ORIGINAL_DIR=${PWD}
eval `scram runtime -sh`
scram arch
popd
echo ${PWD}
which cmsRun
echo "input file: " $1
OFILE="TTZToLLNuNu.yoda"
echo "output file: " ${OFILE}
cmsRun /afs/cern.ch/work/j/jthomasw/private/IHEP/LHCHXS-WG1/Rivet/CMSSW_11_0_0/src/Rivet/HIGGSTOP/test/runRivetAnalyzer_TTH_TTZBCKG_cfg.py inputFiles="file:"$1
