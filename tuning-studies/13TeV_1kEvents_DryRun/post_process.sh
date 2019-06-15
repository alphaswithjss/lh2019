#!/bin/bash

prof2-runcombs ${PWD}/TuningSystem/tunepy8shower/ --pname=params.dat 0:1 -o runcomb.dat
prof2-ipol ${PWD}/TuningSystem/tunepy8shower/ --pname=params.dat --order=1
prof2-tune -d RivetRoutines/ --wfile=weightsFile ipol.dat -r ${PWD}/TuningSystem/tunepy8shower/ -o tune_DRY_RUN --debug --filter
prof2-envelopes ${PWD}/TuningSystem/tunepy8shower/ RivetRoutines/ --pname=params.dat
prof2-sens ipol.dat --grad
