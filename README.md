# VOTCA's commonly used GitHub actions

These GitHub actions are used by multiple VOTCA modules.

Actions:

-   setup, initializes variables, like `cmake_args` to be used during the build.
-   forward, forward commit to master/stable in a module to votca/votca
-   changelog, create and commit a changelog entry
-   release\_changelog, create a md changelog file to be used with actions/create-release
-   copyright, update the copyright in all touched files
