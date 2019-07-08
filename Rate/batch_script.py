def makeBatchConfigFile( job_dir):

    config='#!/bin/sh'
    config+='\n'
    config+='#$ -S /bin/bash \n'
    config+='source /cvmfs/cms.cern.ch/cmsset_default.sh \n'
    config+='cd /u/user/moon/cmssw/CMSSW_8_0_26/src \n'
    config+='eval `scramv1 runtime -sh` \n'
    config+='cd ' + job_dir +  ' \n'
    config+='rm -rf test_C.so \n'
    config+='rm -rf test_C.d \n'
    config+='rm -rf result.log \n'
    config+='root -l -b < x_test.C \n'
    config+='rm -rf test_C.so \n'
    config+='rm -rf test_C.d \n'
    return config
