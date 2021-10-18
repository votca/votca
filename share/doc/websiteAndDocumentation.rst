Website and Documentation
=========================

On this page you will find all the information you need to contribute to VOTCA's
website and documentation.

The website is build with Sphinx and uses a slightly modified version of the
Read the Docs theme. Parts of the documentation is automatically generated
based on the source code of VOTCA and jupyter notebooks, other parts are written
by hand. The markup language of Sphinx is reStructuredText, think of it as
Markdown on steroids. 

To start working on the website and documentation you will first need to be able
to build it, this is discussed first. Next, a brief overview of how to edit
existing pages and adding new ones is given. Finally, we discuss how to change the theme and publish the changes.

Building the documentation and website locally
----------------------------------------------

To be able to see and test changes you make to the website, you will need to
build it locally. You will first need to install and build VOTCA. Once you have
installed VOTCA, navigate to its build folder (if you followed the install
instructions, it will be called ``builddir``). In the build folder, run 

.. code-block:: bash

  make doc

This will probably output errors on the first few tries, read them carefully
since they will almost always be missing python packages. Find which
package/module is missing in the error message and simply install it with ``pip3
install <packagename>``.

A short list of packages you probably need to install manually:

- cma
- recommonmark
- nbsphinx
- sphinx_rtd_theme 

Once the errors are gone, the documentation is build and a static html website is
generated in the ``sphinx.html`` folder within the build folder. You can view
and test the website in any browser simply by navigating to the index.html file
within the ``sphinx.html`` folder.

Rebuilding the website after changes
------------------------------------

After you updated files you will need to rebuild the website to see your changes.

If you simply changed the contents of a .rst file you can rerun ``make doc`` to propagate the changes. If you refresh the browser, the changes will have taken
effect.

If you, however, changed the theme or added (deleted) a file to the website you
will need to reconfigure the website build. To make sure you start with a clean
slate, remove the two sphinx folders in the build folder (i.e. ``sphinx`` and
``sphinx.html``), next rerun the cmake configuration of votca and rerun ``make
doc``. 

Editing existing pages
----------------------

To edit existing pages, you only need limited knowledge of the reStructuredText
markup. You can find quick introductions `here <https://docutils.sourceforge.io/docs/user/rst/quickstart.html>`__ and `here <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`__. All pages are stored as
.rst files in the ``share/doc`` folder of the corresponding repo.

Writing a new page
------------------

To add a new page to the website, you will first need to decide where it belongs,
in the xtp, csg or votca repo. Once you have decided that, navigate to the
``share/doc`` folder within that repo. There will be either an ``index.rst``
file or something like ``XTP-MANUAL.rst``. This file is the main page of that repo on the website. If you look at the ``index.rst`` file in the ``votca`` repo you will
find all the subpages listed here, you need to make sure that your new page
can be found starting from this file. 

As an example, suppose we want to add a test page, called ``text.rst`` to the
developers section of the website. In ``votca/share/doc``, we would create the
new file. Next, we would update the file ``votca/share/doc/index.rst``, we need
to add our file to the table of contents of the site. Look for the table of
contents of the development section and add the test page. It will look
something like this

.. code-block:: text

  .. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Development

   DEVELOPERS_GUIDE.rst
   VOTCA_LANGUAGE_GUIDE.rst
   CODE_OF_CONDUCT.rst
   WEBSITE_AND_DOCUMENTATION.rst
   test.rst

After you have added content to the ``test.rst`` file, rebuild VOTCA and run
``make doc`` to rebuild the website.


Changing the theme
------------------

To change the theme, navigate to ``votca/share/doc/_static/css``, there you will
find the ``theme.css`` file which you can modify to change the look and feel of
the website.


Publishing your changes
-----------------------

All the changes so far are local, to get them published simply make a pull
request in the desired repo. Once the PR is accepted, the site will
automatically be updated.

More Advanced Stuff
-------------------

Except for the automatically generated pages, this website is a basic Sphinx
website, all the information you can find about Sphinx and reStructuredText on
the web is applicable here. So if you need something advanced, simply google it.
