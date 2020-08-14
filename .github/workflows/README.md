Actions:

You have an automerge action which has not been explained. Does this allow you to merge to automatically when the ci completes successfully? What is the purpose of this action?
Also what is this line doing ` uses: "pascalgn/automerge-action@v0.9.0"`? Probably want to include a link https://github.com/pascalgn/automerge-action

-    `continuous-integration-workflow.yml`: runs `cmake`, `make`, `ctest` and `make install` etc. for many different compiler,
     distro / gromacs combinations and special builds for minimal, internal gromacs and coverage
     -    common setup parts done by votca/actions/setup action that can be found here: https://github.com/votca/votca/actions 
          it mainly generates the right cmake arguments (`cmake_args`) for the combination of build parameters
     -    deploy of website and doxygen is happening in `deploy.sh` (push happens for master only)

For the continuous workflow action is this scheduled to run every friday `  - cron:  '0 5 * * FRI'` if so it would probably be a good idea to have a badge displaying whether it is passing or not. 

Can you explain what the id setup is being used for, is this just a unique identifier so you can store ccache in seperate folders to have a fast build? 

`id: setup
 uses: votca/actions/setup@master`

-    `docker-build.yml`: does the build and deploy of the docker container from the `Dockerfile`
-    `continuous-integration-workflow.yml` and `docker-build.yml` are triggered for pushes to master and pull requests to
     master as well as on a schedule once a week to pull in updates from the votca/buildenv container


-    `create-pr.yml`: will create a pull request when a `update_(master|stable)_submodules` branch is created (manually or by
     the forward action in a module)

I don't understand what this is saying. This git action creates a pull request when a branch is created? When would you use this? Will this automatically update all the submodules and submit a pull request when the branch is created? Is this for use when a repo is updated, e.g. it will automatically try to update the votca/votca?  
