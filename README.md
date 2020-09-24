# CMS DY - 2016 Data

## Source (Data and MC) files 

Source files can be downloaded at:

```
/afs/cern.ch/user/f/ftorresd/workSpace/public/forJamshid
```

## Setup CMSSW

```
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_10_1_9/ 
cmsenv
cd -
```

## Running the code

Data: `root -l -b -q NanoReaderData.cc++`

MC: `root -l -b -q NanoReaderMC.cc++`

All the produced outputs files will be saved at: `./outputs/`.
