Actions:

-    `continuous-integration-workflow.yml`: runs `cmake`, `make`, `ctest` and `make install` etc. for many different compiler,
     distro / gromacs combinations and special builds for minimal, internal gromacs and coverage
     -    common setup parts done by votca/actions/setup action that can be found here: https://github.com/votca/votca/actions 
          it mainly generates the right cmake arguments (`cmake_args`) for the combination of build parameters
     -    deploy of website and doxygen is happening in `deploy.sh` (for master only)
-    `docker-build.yml`: does the build and deploy of the docker container from the `Dockerfile`
-    `continuous-integration-workflow.yml` and `docker-build.yml` are triggered for pushes to master and pull requests to
     master as well as on a schedule once a week to pull in updates from the votca/buildenv container
-    `create-pr.yml`: will create a pull request when a `update_(master|stable)_submodules` branch is create (manually or by
     the forward action in a module)
