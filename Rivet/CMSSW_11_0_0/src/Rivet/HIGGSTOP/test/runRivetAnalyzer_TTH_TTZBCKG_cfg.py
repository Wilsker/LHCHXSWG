import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("GEN")

options = VarParsing.VarParsing ('analysis')
options.inputFiles = '/store/mc/RunIIFall17MiniAODv2/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/90000/0016D8CA-2E32-E911-8820-246E96D10CBC.root'
options.maxEvents = -1
options.parseArguments()

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cff")
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    )

process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.rivetAnalyzer.HepMCCollection = cms.InputTag("genParticles2HepMC:unsmeared")
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2019_TTH_TTZBCKG')
process.rivetAnalyzer.OutputFile = cms.string("TTZToLLNuNu.yoda")

process.p = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.rivetAnalyzer)
