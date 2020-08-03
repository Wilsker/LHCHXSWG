import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("GEN")

options = VarParsing.VarParsing ('analysis')
options.inputFiles = '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/A2184ECC-A9AF-E811-A9C6-FA163E042DE1.root'
#options.inputFiles = '/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/0EF562DB-C1AD-E811-9E06-E0071B7A46D0.root'
options.maxEvents = -1
options.parseArguments()

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cff")
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

#Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    )

process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.rivetAnalyzer.CrossSection = cms.double(0.196)
#process.rivetAnalyzer.CrossSection = cms.double(0.1960)#QCD NLO
#process.rivetAnalyzer.CrossSection = cms.double(0.01598)#NLO EWK Corrections only
process.rivetAnalyzer.HepMCCollection = cms.InputTag("genParticles2HepMC:unsmeared")
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMSATLAS_TTW_ttHBCKG')
process.rivetAnalyzer.OutputFile = cms.string("TTWJetsToLNu.yoda")
process.rivetAnalyzer.useLHEweights = cms.bool(True)
#process.rivetAnalyzer.LHEweightNumber = cms.int32(0)

process.p = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.rivetAnalyzer)
