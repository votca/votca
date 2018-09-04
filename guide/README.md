This is the first iteration of the VOTCA-XTP tutorial. 

If you find bugs, weird behaviour or just incomprehensible instructions, please open an issue [here](https://github.com/votca/xtp-tutorials/issues)

We will try to fix it as fast as possible.

For other questions, we provide an email list:

[google group](https://groups.google.com/forum/?hl=de#!forum/votca-xtp)

and a slack channel:

[Slack](https://votca.slack.com/messages/C7XVBE9EG/?)

At the moment we are also working on a manual:

[Manual](http://doc.votca.org/xtp-manual.pdf)


The tutorial is organised as follows:

the folder 'methane' contains a molecular dynamics data and a mapping file system.xml
using the 'run.sh' script you can calculate all data necessary to run kmc simulations.
For some more expensive parts the parameters are calculated only for one pair.
The KMC calculators are thus not tested

In the folder 'tools' are calculators which are independent of '.sql' files:

'dftgwbse_CH4' runs a DFT and GWBSE on CH4 and then performs an espfit on it, outputting a '.mps' file
additionally a '.cube' file is produced for visualisation. 

'dftgwbse_CO_geoopt'  performs an excited state geometry optimisation using numerical derivatives, it is a bit more expensive







