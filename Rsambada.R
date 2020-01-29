##testing RSamBada for landscape genomics
library(R.SamBada)
library(installr)
library(devtools)
library(SNPRelate)
library(packcircles)

genoFile="new_char_Poly.ped"
genoFile

downloadSambada(directory=NULL)
prepareGeno(genoFile,outputFile="charr-mol.csv",FALSE,interactiveChecks=FALSE)

##sample IDs
idFile="charr-subset-id.csv"
readLines(idFile, n=20)

##longs, lats
locationFile= "charr-subset.csv"
readLines(locationFile, n=20)

##generate GDS file with SNPRelate
snpgdsPED2GDS("new_char_Poly.ped","new_char_Poly.map","new_char_Poly.gds")
gdsFile="new_char_Poly.gds"
gdsFile

##bioclim variables
envFile="charr-subset-env3.csv"
readLines(envFile, n=20) #View the first 20 lines of the file

##revised the prepareEnv function to include autosome.only=TRUE for SNPRelate
##load the preprocessing_edits script and run prepareEnv function (line 792-end)
##changed maxcorr to 0.7, mafThresh 0.01, missingnessThresh 0.05, ldThresh 1 (keep snps in LD)
prepareEnv(envFile=envFile, outputFile="charr-subset-env-export.csv", maxCorr=0.7, idName='short_name', genoFile=gdsFile, numPc=0.2, mafThresh=0.01, missingnessThresh=0.05, ldThresh=1, numPop=NULL, x='longitude', y='latitude', interactiveChecks=FALSE, locationProj=4326 )

envFile2="charr-subset-env-export.csv"
readLines(envFile2, n=20)

##from prepareGeno above
genoFile2="charr-mol.csv"

##run Sambada
sambadaParallel(genoFile=genoFile2, envFile=envFile2, idGeno='ID_indiv', idEnv='short_name', dimMax=2, cores=2, saveType='END ALL', populationVar='LAST', outputFile=file.path(tempdir(),'charr-subset-mol'))

prep = prepareOutput(sambadaname='charr-subset-mol', gdsFile = gdsFile, dimMax=2, popStr=TRUE, interactiveChecks=FALSE)

##marker information
Markers <- prep$sambadaOutput$Marker
write.csv(Markers,"markers.csv")

##these 6 variables are the uncorrelated set from the prepareEnv step
plotManhattan(prep, c('bio1','bio2','bio4','bio8','bio15','bio16'),chromo='all',valueName='pvalueG')

plotResultInteractive(preparedOutput=prep, varEnv='bio1', envFile=envFile2,species='drerio',pass=50000,x='longitude',y='latitude',gdsFile=gdsFile, IDCol='short_name', popStrCol='pop1')

##plot presence/absence of markers of interest across range
plotMap(envFile=envFile2, x='longitude', y='latitude', locationProj=4326,  popStrCol='pop1', gdsFile=gdsFile, marker='AX-182164053_GG', mapType='marker', varEnvName='bio1', simultaneous=FALSE)

##plot environmental variables across range
plotMap(envFile=envFile2, x='longitude', y='latitude', locationProj=4326,  popStrCol='pop1', gdsFile=gdsFile, marker='AX-182164053_GG', mapType='env', varEnvName='bio1', simultaneous=FALSE)

