
#!/bin/bash

# Gaining access to the CRC modules
if [ -r /opt/crc/Modules/current/init/bash ]; then
    source /opt/crc/Modules/current/init/bash
fi

export PATH=$PATH:/bin:/usr/bin
export R_LIBS=/afs/crc.nd.edu/user/a/adraper2/condor/Rlibs

module load bio/R/3.4.0                  # Required modules

./runMCMCbatch.R $1 $2