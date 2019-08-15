#!/bin/bash

python qmmm_spectra.py -f ../QMMM/frame_10000/job_1_static/checkpoint_iter_1.hdf5 -o static_spec.dat
python qmmm_spectra.py -f ../QMMM/frame_10000/job_0_vacuum/checkpoint_iter_1.hdf5 -o vacuum_spec.dat

gnuplot spectra.gp
