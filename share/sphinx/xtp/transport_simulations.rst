.. _sec_transport_simulations:

Transport Simulations
#####################

This chapter describes how to compute mobilities, times of flight and exciton diffusion lengths from a state file with rates, which calculator answers which question, and what to watch out for. The theory is in :ref:`the transport chapter <sec_transport_theory>`; the tutorial notebook ``xtp-tutorials/tools/transport_lattice/TRANSPORT.ipynb`` goes through all calculators on a small model system.

Which calculator
****************

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Question
     - Calculator
     - Method
   * - Mobility and diffusion tensor in the bulk, at a given carrier density
     - ``kmcmultiple``
     - KMC in a periodic box with a fixed number of carriers; exact single-carrier values from the master equation
   * - Time of flight through a slab, mobility from it
     - ``grw``
     - Exact mean hitting times (no sampling), one or two carriers
   * - The same with many carriers, pulses, the transient current or the distribution of transit times
     - ``kmctof``
     - KMC with source and sink regions
   * - Exciton lifetime, diffusion length, where excitons decay
     - ``kmclifetime``
     - KMC with site-dependent decay; exact single-exciton values
   * - A model system to test methods and settings
     - ``latticegen`` (``xtp_tools``)
     - Cubic lattice with Gaussian site energies, written as a state file
   * - A larger system with the statistics of a morphology (finite-size effects)
     - ``networkfit``, ``networkgen`` (``xtp_tools``)
     - Stochastic network model fitted to a state file, generated at any size

All four transport calculators run with ``xtp_run``, for example ``xtp_run -e grw -o grw.xml -f state.hdf5``; options can also be set on the command line with ``-c``, e.g. ``-c carriertype=hole field="1e7 0 0"``. They read the state file and never write to it.

Input
*****

The rates are Marcus rates :eq:`equ:theory:marcus` computed from the state file, so it must contain for the carrier type in question (``electron``, ``hole``, ``singlet`` or ``triplet``):

* site energies of the segments (internal energies from ``einternal`` plus electrostatic and polarization contributions, see the corresponding chapters),
* a neighbor list (``neighborlist``) with electronic couplings :math:`|J|^2` (e.g. from ``iqm`` or ``ipodcoupling``),
* reorganization energies (``einternal``) and, if wanted, outer-sphere contributions.

``eanalyze`` and ``ianalyze`` show the distributions of site energy differences and couplings and are worth a look before any transport calculation: the width of the density of states controls trapping, and pairs with tiny couplings hardly contribute.

Common options
==============

``carriertype``, ``temperature`` (K) and ``field`` (V/m, three components) define the rates. 1e7 V/m equals 100 kV/cm or 0.01 V/nm. Fields far above 1e8 V/m (1 MV/cm) are usually unphysical for organic semiconductors.

Segment types that should take no part in transport (e.g. an insulating matrix) are listed in ``ignoresegments`` (``kmcmultiple``, ``kmclifetime``) or ``tof.ignoretypes`` (``grw``, ``kmctof``): these segments get no hops.

Every calculator writes a short result file (``<calculator>_results.dat``) with ``key value unit`` lines in SI units (mobilities in cm\ :sup:`2`/Vs), meant for scripts and parameter scans, in addition to its log output.

Time of flight: ``grw`` and ``kmctof``
**************************************

Both use the same slab geometry, set in the ``tof`` block:

* ``axis``: the box vector (``a``, ``b`` or ``c``) across which the periodicity is cut. Transport runs along the normal of the face spanned by the other two box vectors; the sample stays periodic in those two directions.
* ``source.thickness`` and ``sink.thickness`` (nm): slabs at the two faces; and/or ``source.pattern``, ``sink.pattern``: segment type wildcards (both criteria must hold if both are given). Use patterns for explicit electrode or interface molecules.

The log reports the number of source and sink sites, the transport length :math:`L` and how many sources have no path to a sink (they are left out). The field should point along the normal so that it drives the carriers to the sink; a component along the periodic directions makes the rates violate detailed balance, which ``grw`` reports and handles with a non-symmetric solver.

``grw`` (graph random walk)
===========================

Solves for the mean hitting time of the sink without random numbers. ``numberofcarriers`` is 1 or 2; for two carriers ``observable`` selects ``first`` (time until the first carrier arrives) or ``transit`` (mean transit time per carrier). Results: the mean over a uniform injection on the source sites (``tau_mean``), its harmonic mean, the standard deviation of the hitting time, and the mobility :eq:`equ:theory:tofmobility`. ``sourcefile`` lists the time per source site or pair, ``occupationfile`` the expected time spent on every segment, which shows where the carriers are held up.

Solver settings rarely need changes. ``solver auto`` picks a direct factorization up to ``directlimit`` unknowns and an iterative solver above. The log reports the backward error; a warning means the solution is not reliable. Two carriers need about :math:`n^2/2` unknowns (``maxstates`` guards the memory): a few thousand sites are fine, :math:`10^5` are not.

``kmctof``
==========

Simulates ``runs`` independent experiments. In each, ``injection.pulses`` pulses of ``injection.carriers`` carriers each are placed on free source sites, ``injection.period`` apart, with a ``delta``, ``gaussian`` or ``rectangular`` time profile of ``injection.width``. Results:

* first arrival and mean transit time per carrier, with standard errors over the runs and the share of the slowest 1 % of runs in the sum,
* ``runfile``: one line per run,
* ``current.file``: transient current and collected fraction on a logarithmic time axis (``current.tmin``, ``current.tmax``, ``current.binsperdecade``).

Runs that reach ``maxtime`` or ``maxevents`` before all carriers arrived are *censored*: they are counted and left out of the means, which biases the times low. Raise the limits if any are reported. For one or two carriers, compare with ``grw`` on the same system: the two must agree within the standard errors.

Bulk transport: ``kmcmultiple``
*******************************

Places ``numberofcarriers`` carriers on random sites (``injectionpattern``) of the periodic box and follows them.

* ``runtime``: simulated time in seconds if at most 100, otherwise a number of KMC events (hops and attempts into occupied sites); ``outputtime`` is read in the same unit and sets the sampling of the trajectory (``trajectoryfile``, unwrapped positions of run 0) and of the ``timefile``.
* ``runs``: independent runs in parallel; with more than one run the errors are standard errors over the runs, otherwise they come from the time blocks.
* ``equilibration``: simulated time at the start of each run left out of all averages. The carriers start on random sites and relax into the low-energy tail of the density of states; the ``timefile`` shows the mean energy per carrier over time and tells how long this takes.
* ``blocks``: the minimum number of time blocks for the diffusion tensor. Blocks must be longer than the trapping times; if the diffusion tensor changes when ``blocks`` is reduced, the runs are too short.

Results: drift velocity, mobility :eq:`equ:theory:mobility` along the field, the diffusion tensor with its eigenvalues and, without field, the Einstein mobility. With ``exact`` (default) the drift velocity, mobility and diffusion tensor of a single carrier are also computed from the master equation; for one carrier they are the reference for the KMC, for several carriers the difference shows the effect of exclusion. ``occfile`` contains the occupation probability of every site, ``ratefile`` all rates.

The mobility from ``kmcmultiple`` is the stationary bulk mobility. A time-of-flight mobility from a thin slab can be considerably higher when the carriers cross the sample before they have relaxed into the tail of the density of states (dispersive transport); comparing both, and ToF results for different thicknesses, shows whether the sample is thick enough.

Excitons: ``kmclifetime``
*************************

Needs the lifetime of every site in ``lifetimefile``:

.. code-block:: xml

   <lifetimes>
     <site id="0">1.2e-9</site>
     <site id="1">0.9e-9</site>
     ...
   </lifetimes>

(seconds; ids are segment ids; segments of ignored types may be left out). ``numberofcarriers`` excitons are present at any time; whenever one decays, a new one is created on a random free site (``injectionpattern``). The simulation stops after ``numberofinsertions`` decays in total, distributed over ``runs``. ``carriertype`` is ``singlet`` by default.

Results: mean lifetime, hops per exciton, mean distance, diffusion length :math:`\sqrt{\langle|\mathbf{R}|^2\rangle}` and its components, and :math:`D=\langle|\mathbf{R}|^2\rangle/(6\langle\tau\rangle)`. With ``exact`` (default) lifetime and diffusion length of a single exciton are also computed from :eq:`equ:theory:lifetime`; ``numberofinsertions=0`` gives only these. ``trajectoryfile`` lists every decay, ``decayprobfile`` the local decay ratios :math:`k/(k+W)`, the fraction of the KMC decays and the exact decay probability per site, ``carrierenergy`` writes the mean energy of the excitons over time.

Output quantities
*****************

Times are in s, lengths in nm in the text files (m in the result files), energies in eV and mobilities in cm\ :sup:`2`/Vs unless stated otherwise. :math:`N_\mathrm{inj}` is the number of injected carriers, :math:`L` the transport length and :math:`F_n` the field component along the transport normal.

``grw``
=======

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Quantity
     - Definition
   * - ``tau_mean``
     - Mean hitting time :math:`\langle T\rangle = \frac{1}{N_s}\sum_s \tau_s` for uniform injection on the :math:`N_s` source sites (two carriers: source pairs). For ``observable=transit`` the mean transit time per carrier.
   * - ``tau_harmonic``
     - Harmonic mean :math:`N_s/\sum_s \tau_s^{-1}` over the injection states, as used in [Chen:2024]_. It weights fast sources more; the difference to ``tau_mean`` shows how much single sources matter.
   * - ``tau_std``
     - Standard deviation of the hitting time :math:`T` itself (not of its mean) for uniform injection, :math:`\sqrt{\langle T^2\rangle-\langle T\rangle^2}`, from the second moments. Not available for ``observable=transit``.
   * - ``mobility_mean``, ``mobility_harmonic``
     - :math:`L/(\tau |F_n|)` with ``tau_mean`` or ``tau_harmonic``.
   * - ``sourcefile``
     - One line per injection state: segment id, position along the normal and site energy of each carrier, then :math:`\tau_s` and the standard deviation of :math:`T` from this state.
   * - ``occupationfile``
     - Per segment: id, position along the normal, site energy, role (source, bulk, sink), expected time :math:`g_A` the segment is occupied before absorption (uniform injection), and the *fraction of the time of flight* :math:`g_A/(N_c\langle T\rangle)` with the number of carriers :math:`N_c`. The fractions add up to 1; sinks have 0.

``kmctof``
==========

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Quantity
     - Definition
   * - ``first_mean``, ``first_error``, ``first_std``
     - Time until the first carrier of a run reaches a sink: mean over the uncensored runs, its standard error and the standard deviation over runs.
   * - ``transit_mean``, ``transit_error``, ``transit_std``
     - Per run the mean of (absorption time - injection time) over its absorbed carriers; then mean, standard error and standard deviation over the uncensored runs. Equal to ``first`` for one carrier per run.
   * - ``transit_top_percent_share``
     - Sum of the per-run transit times of the slowest 1 % of the runs (at least one run) divided by the sum over all runs. For exponentially distributed transit times it is about 0.06 (:math:`(1+\ln 100)/100`); ``kmctof`` warns above 0.2, when a few runs dominate the mean.
   * - ``mobility_transit``
     - :math:`L/(\texttt{transit\_mean}\,|F_n|)`.
   * - ``censored``
     - Runs stopped by ``maxtime`` or ``maxevents`` before all carriers arrived; they are left out of the means.
   * - ``runfile``
     - Per run: run number, injected and absorbed carriers, censored (0/1), first arrival, mean transit time, hops, null events (attempts into occupied sites), time of the last event.
   * - ``current.file``: *current*
     - Transient current per injected carrier,

       .. math::

          I(t)=\frac{1}{L\,N_\mathrm{inj}\,\Delta t}\sum_{\text{hops in }[t,t+\Delta t)}\Delta x ,

       in 1/s, with :math:`\Delta x` the displacement of a hop along the transport normal (Ramo-Shockley). For :math:`N_\mathrm{inj}` injected carriers of charge :math:`q` the current in the external circuit is :math:`q\,N_\mathrm{inj}\,I(t)`. The integral :math:`\int I\,dt` is the mean displacement per injected carrier in units of :math:`L`; the parts before the first and after the last bin are given in the file header. It is close to 1 when all carriers cross the sample, but not exactly, since carriers start anywhere in the source slab and :math:`L` is measured from the mean source position to the first sink layer. Columns: lower and upper bin edge, geometric bin center, :math:`I`, collected fraction.
   * - ``current.file``: *collected fraction*
     - Number of carriers absorbed up to the upper edge of the bin, divided by :math:`N_\mathrm{inj}`.

``kmcmultiple``
===============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Quantity
     - Definition
   * - ``velocity_*``, ``velocity_error_*``
     - Drift velocity :math:`\mathbf{v}`: unwrapped displacement per carrier over the observation window (after ``equilibration``) divided by its length, averaged over runs; standard error over runs, or from the time blocks for one run.
   * - ``mobility``, ``mobility_error``
     - :math:`\mu=q\,\mathbf{v}\cdot\mathbf{F}/|\mathbf{F}|^2` with :math:`q=\pm1`: positive when the carriers drift along the force on them, for electrons and holes alike.
   * - ``D_**``, ``D_error_**``
     - Diffusion tensor :math:`D=\mathrm{Cov}(\Delta\mathbf{r}_\tau)/(2\tau)` from the displacements :math:`\Delta\mathbf{r}_\tau` of all carriers over time blocks of length :math:`\tau` (between ``blocks`` and twice as many blocks per run), with the drift removed; averaged over runs.
   * - ``einstein_mobility``
     - Without field: :math:`\mathrm{tr}(D)/(3k_BT)`.
   * - ``exact_*``
     - The same for one carrier from the master equation, :eq:`equ:theory:exactdrift` and :eq:`equ:theory:exactdiffusion`.
   * - ``occfile``
     - Per site the probability that it is occupied, i.e. the occupied time divided by the window, summed over carriers and averaged over runs; the sum over sites is the number of carriers.
   * - ``trajectoryfile``
     - Run 0: time, accepted hops so far, and the unwrapped position (nm) of every carrier, sampled every ``outputtime``.
   * - ``timefile``
     - Run 0 at the same times: hops, mean site energy of the carriers, mobility from the mean displacement since the first sample, mean displacement per carrier along the field, and the length of the mean displacement vector.

``kmclifetime``
===============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Quantity
     - Definition
   * - ``lifetime``, ``lifetime_error``
     - Mean time from creation to decay over all decays, and its standard error (assuming independent decays, exact for one carrier).
   * - ``mean_distance``, ``msd``, ``msd_*``
     - Mean of :math:`|\mathbf{R}|`, of :math:`|\mathbf{R}|^2` and of :math:`R_\alpha^2`, with :math:`\mathbf{R}` the unwrapped displacement from creation to decay.
   * - ``diffusion_length``
     - :math:`L_D=\sqrt{\langle|\mathbf{R}|^2\rangle}`; always at least ``mean_distance``.
   * - ``diffusion_coefficient``
     - :math:`\langle|\mathbf{R}|^2\rangle/(6\langle\tau\rangle)`; equal to :math:`D` for uniform decay rates in an isotropic medium.
   * - ``exact_*``
     - The same for one carrier from :eq:`equ:theory:lifetime`.
   * - ``trajectoryfile``
     - One line per decay: time of the decay, running number, carrier number, lifetime, hops, segment of the decay, displacement (nm), run.
   * - ``decayprobfile``
     - Per site: :math:`k/(k+W_\mathrm{out})` and :math:`k/(k+W_\mathrm{in})` with the sums of the outgoing and incoming hop rates (local measures of how likely an exciton decays there), the fraction of the KMC decays on the site, and the exact decay probability :math:`k_Ag_A` for one exciton.
   * - ``carrierenergy.outputfile``
     - Run 0, every ``outputsteps`` hops: time, hops, mean site energy of the excitons present, and its exponential average with weight ``alpha``.

Larger systems: stochastic network model
****************************************

A morphology of a few thousand molecules is often too small for transport: a box samples the occupied part of a Gaussian density of states of width :math:`\sigma` only if :math:`(\sigma/k_\mathrm{B}T)^2 \lesssim -5.7 + 1.05 \ln N` [Baumeier:2012]_, i.e. :math:`N\approx10^8` sites for :math:`\sigma = 3.7\,k_\mathrm{B}T` and more than :math:`10^{25}` for :math:`7.5\,k_\mathrm{B}T`. Below this size the mobility of a box depends on its deepest sites and changes by orders of magnitude from one box to the next. The stochastic network model generates hopping networks of any size with the statistics of the morphology. The model, the fit and the finite-size mathematics are described in :ref:`sec_stochastic_network_theory`; the tutorial ``xtp-tutorials/tools/stochastic_network/STOCHASTIC_NETWORK.ipynb`` goes through the workflow for MADN and bathocuproine (BCP).

Workflow
========

.. code-block:: bash

   # 1. fit the model to a state file with site energies, couplings and reorganization energies
   xtp_tools -e networkfit -c job_name=morph statefile=morph.hdf5 "carriers=h e" statsfile=morph_stats.dat
   # 2. generate a network of any size (and as many realizations as needed, seed=1,2,...)
   xtp_tools -e networkgen -c job_name=big parameters=morph_network.xml sites=64000 seed=1
   # 3. any transport calculator, e.g. the exact single-carrier mobility
   xtp_run -e kmcmultiple -c carriertype=hole field="1e7 0 0" runtime=1000 -f big.hdf5

``networkfit`` writes the parameter file (``parameterfile``, default ``<job_name>_network.xml``; lengths in nm, energies in eV, readable and editable), prints the fitted parameters, the traps it found and a validation table, and with ``statsfile`` writes all fitted functions for plots: :math:`D(r)`, :math:`H(r)`, :math:`g(r)` (morphology and model), :math:`f(r)`, pairs per site, :math:`\kappa(r)` with the fitted moving average, the fit quality for every :math:`l_c`, and the mean :math:`\log_{10}|J|^2` against distance. The fit of the point process takes about two minutes for 1000 sites; ``points.fit=false`` uses given values instead.

Important options of ``networkfit``:

* ``carriers``: carrier types to fit; the first two are correlated (site energies and coupling residuals). Carriers without energetic disorder in the state file get :math:`\sigma=0`.
* ``l_b``, ``l_c``: grids (first, last, step) of the short- and long-range moving averages. A warning is printed if the best value is at the end of a grid.
* ``traps.find`` (default true), ``traps.threshold``, ``traps.mass``, ``traps.maxsize``, ``traps.gap``: identification of traps from the slow modes of the master equation at ``temperature`` and zero field. Each trap is listed with its sites, energies, equilibrium weight, release and capture rates. Trap sites are left out of the bulk statistics and become templates. The dense eigenvalue problem limits this to ``traps.maxsites`` (5000) sites.
* ``validate`` (default 4): model boxes of the size of the morphology, compared with the morphology (pairs per site, nearest-neighbor distance, :math:`\sigma`, moments of :math:`\log_{10}|J|^2`, :math:`\kappa(r)`, energy correlation of the carriers, the exact zero-field diffusion coefficient of one carrier, traps found), as mean and standard deviation over the boxes. The transport of the morphology should lie within the spread of the boxes.

Important options of ``networkgen``: ``sites`` and ``box`` (by default the density and shape of the fitted box are kept), ``seed`` (the same network on every platform, up to the last digits of the math library), ``trapfactor`` (multiplies the trap densities; 0: no traps), ``energyfile`` (positions and site energies). It prints whether the network is large enough for the density of occupied states.

Recommendations
===============

* Check the validation table and the plots of the ``statsfile`` first. With about a thousand molecules, correlations of the energies over more than a third of the box are uncertain; the :math:`l_c` scan shows how well they are determined.
* Increase the size of the networks until the mobility stops changing, with several realizations per size; the box-to-box spread shrinks with size, the mean drifts as long as the network is too small.
* Do not average many small networks instead: plain means of their mobilities are biased, and a Boltzmann-weighted mean (:eq:`equ:theory:boltzmannmean`) only stands in for one larger network.
* One trap in a morphology fixes its density only to between 0.025 and 5.6 times the observed value (95 % Poisson interval): scan ``trapfactor``, including 0, to see how much the result depends on traps.
* The exact single-carrier values of ``kmcmultiple`` take seconds for :math:`10^4` sites and minutes for :math:`10^5`; beyond that use the kinetic Monte Carlo results.

Limits of the model: traps are inserted with the couplings of the bulk model, so traps that act through their couplings or reorganization energies are reproduced only through their energies; the energies of a site and its number of neighbors are independent in the model (in morphologies they are weakly correlated), so are couplings and energy differences; the mean of :math:`\log_{10}|J|^2` is linear in :math:`r`, which overestimates the couplings of the closest pairs slightly.

Checklist
*********

* Statistical errors assume independent samples. For ``kmctof`` they are over runs; check the share of the slowest runs as well. For ``kmcmultiple`` with one run they come from time blocks; prefer several runs.
* Compare KMC with the exact single-carrier values whenever they are available: they agree within a few standard errors for one carrier. Larger differences point to too short runs, too few blocks or a missing equilibration.
* Results do not depend on the number of threads (``-t``), only on ``seed`` and ``runs``.
* At high occupancy most KMC events are attempts into occupied sites; the fraction is printed. The results stay exact, but the simulation gets slower.
* ``latticegen`` writes model systems with known parameters (``dims``, ``spacing``, ``cutoff``, ``J2``, ``reorganization``, Gaussian site energies with ``energy.sigma``, uncorrelated or correlated). They are useful to test settings and to understand trends before running a real morphology.
