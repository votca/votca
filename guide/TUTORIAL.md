## VOTCA-XTP tutorial

This is the first iteration of the VOTCA-XTP tutorial. 

If you find bugs, weird behaviour or just incomprehensible instructions, please open an issue [here](https://github.com/votca/xtp-tutorials/issues)

We will try to fix it as fast as possible.

For other questions, we provide an email list:

[google group](https://groups.google.com/forum/?hl=de#!forum/votca-xtp)

and a slack channel:

[Slack](https://votca.slack.com/messages/C7XVBE9EG/?)

At the moment we are also working on a manual:

[Manual](http://doc.votca.org/xtp-manual.pdf)


### The tutorial is organised as follows:

The folder **GROMACS** contains **Methane** molecular dynamics data and a mapping file `system.xml`
using the `run.sh` script you can calculate all data necessary to run kmc simulations.
For some more expensive parts the parameters are calculated only for one pair.

You find exactly the same system in **KMC_Methane** with all data calculated. Here you can run kmc simulations. 
You find a `run.sh` in this folder as well to run kmc simulation.

The folder **LAMMPS** contains **Thiophene** molecular dynamics data and a mapping file `system.xml`
using the `run.sh` script you can calculate all data necessary to run kmc simulations.
For some more expensive parts the parameters are calculated only for one pair.

You find exactly the same system in **KMC_Thiophene** with all data calculated. Here you can run kmc simulations. 
You find a `run.sh` in this folder as well to run kmc simulation.

In the **LAMMPS** folders you will find a folder **plotting_scripts" you can run `.gp` files using `gnuplot`.
The `.py` scripts have a `-h` option, which tells you what they do.


In the folder **tools** are programs which are independent of state files:

**dftgwbse_CH4** runs a DFT and GWBSE on CH4 and then performs an espfit on it, outputting a `.mps` file
additionally a `.cube` file is produced for visualisation. It uses the internal DFT engine. 

**dftgwbse_CO_geoopt**  performs an excited state geometry optimisation using numerical derivatives, it is a bit more expensive







