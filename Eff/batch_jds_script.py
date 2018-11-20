def makeBatchConfigjdsFile(job_dir):
    config ='executable = run.sh \n'
    config+='universe   = vanilla \n'
    config+='log        = log.log \n'
    config+='output     = out.out \n'
    config+='error      = err.err \n'
    config+='transfer_input_files = x_test.C \n'
    config+='should_transfer_files = YES \n'
    config+='when_to_transfer_output = ON_EXIT \n'
    config+='requirements = ( HasSingularity == true) \n'
    config+='+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el6:latest"  \n'
    config+='+SingularityBind = "/cvmfs, /cms, /share" \n'
    config+='queue'

    return config

