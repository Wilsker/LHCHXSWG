# ttH Backgrounds Rivet

Package to run a study ttH backgrounds as a Rivet analysis using CMSSW framework.

## Author: Joshuha Thomas-Wilsker, Kirill Grevstov
### Institute: IHEP CAS, DESY

## Instructions
- Change to CMS directory and source CMSSW.
- Rivet analyser (python script) can be found here: CMSSW_XX/src/Rivet/TOP/runRivetAnalyzer_TTH_XXX.py
- This script is the CMSSW config that will load the necessary packages to run Rivet in CMSSW as well as setting up the input/output files and add the Rivet fragment to your process (see rivetAnalyzer.AnalysisName - can be found in CMSSW_10_6_0/src/Rivet/TOP/src/).
- Testing the code locally can be done simply using:
'''
cmsRun <CMSSW config>
'''
- One can then start editing scripts XXX-rivet.sub and XXX-rivet.sh which are used to run multiple jobs on condor.
- XXX-rivet.sub is the condor submission script that will send jobs to run XXX-rivet.sh for various input files.
- The submission script needs both the bash script to run the job and a list of files to run on e.g. XXX_fileList.txt.
- One can change the output names, condor queue etc as suits.
- Outputs, error and log files are stored in the appropriately named directories.
- To submit a job for each file in the file list simply run:
'''
condor_submit XXX-rivet.sub
'''
- Once jobs are done, the output yoda files will need merging. This can be done using yodamerge. Run the following command for more information:
'''
yodamerge -h
'''
- To print a list of files, separated by spaces in a single line, that can be used for merging try using:
'''
ls rivet_output/TTZ/* | tr '\n' ' '
'''
- Then use:
'''
yodamerge --output TTZToLLNuNuMerged.yoda <output from the above command>.
'''
- Once you have your merged .yoda file, you can create a root file from it:
'''
yoda2root TTZToLLNuNuMerged.yoda
'''
