import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("GEN")

options = VarParsing.VarParsing ('analysis')
options.inputFiles = '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/00DFE1C9-BEAD-E811-A06B-0242AC130002.root'
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
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2019_TTH_TTWBCKG')
process.rivetAnalyzer.OutputFile = cms.string("TTWJetsToLNu.yoda")

process.p = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.rivetAnalyzer)
