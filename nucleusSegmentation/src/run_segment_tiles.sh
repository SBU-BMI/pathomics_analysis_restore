#!/bin/bash

export OMP_NUM_THREADS=1

mainSegmentFeatures -t tiles "$@"
