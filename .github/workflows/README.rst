Actions
=======

-  :code:`continuous-integration-workflow.yml` runs :code:`cmake`, :code:`make`, :code:`ctest` and :code:`make install` etc. for many different compiler,
   distribution, gromacs combinations and special builds for minimal, internal gromacs and coverage
     
     -  common setup parts done by `votca/actions/setup action <https://github.com/votca/votca/actions>`_
        it mainly generates the right cmake arguments (:code:`cmake_args`) for the combination of build parameters
     -  deploy of website and doxygen is happening in :code:`deploy.sh` (actual push only happens for master)
-  :code:`docker-build.yml`: does the build and deploy of the docker container from the `Dockerfile`
-  :code:`continuous-integration-workflow.yml` and :code:`docker-build.yml` are triggered for pushes to master and pull requests to
   master as well as on a schedule once a week to pull in updates from the votca/buildenv container
-  :code:`create-pr.yml`: will create a pull request when the :code:`update_(master|stable)_submodules` branch gets created (manually by us or automatically by
   the `forward action <https://github.com/votca/actions/tree/master/forward>`_ running a module)
-  :code:`create-stable-pr.yml`: will create a pull request from stable into master
-  :code:`release.yml`: can be triggered manually in the GitHub UI to create a release
-  :code:`stable_cron.yml`: workflow to run a CI and docker deploy for the :code:`stable` branch (workaround for this `bug <https://github.community/t/scheduled-builds-of-non-default-branch/16306>`_)
-  :code:`automerge.yml`: will automatically merge pull request that are marked with the :code:`automerge` label and fullfill all merge requirement, e.g. CI passed and review was made.
