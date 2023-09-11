############
Contributing
############

This document briefly describes how to contribute to
`mzSpecLib <https://github.com/hupo-psi/mzspeclib>`_.



Before you begin
################

If you have an idea for a feature, use case to add or an approach for a bugfix,
you are welcome to communicate it with the community by opening a
thread in `GitHub Issues <https://github.com/hupo-psi/mzspeclib/issues>`_.



Development setup
#################

Local install
*************

#. Setup Python 3, and preferably create a virtual environment.
#. Clone the `mzSpecLib repository <https://github.com/hupo-psi/mzspeclib>`_.
#. Use pip in editable mode to setup the development environment:

.. code-block:: sh

    pip install --editable .[test,docs]


Unit tests
**********

Run tests with ``pytest``:

.. code-block:: sh

    pytest ./tests


Documentation
*************

To work on the documentation and get a live preview, install the requirements
and run ``sphinx-autobuild``:

.. code-block:: sh

    pip install .[docs]
    sphinx-autobuild  ./docs/ ./docs/_build/

Then browse to http://localhost:8000 to watch the live preview.



How to contribute
#################

- Fork `mzSpecLib <https://github.com/hupo-psi/mzspeclib>`_ on GitHub to
  make your changes.
- Commit and push your changes to your
  `fork <https://help.github.com/articles/pushing-to-a-remote/>`_.
- Ensure that the tests and documentation (both Python docstrings and files in
  ``/docs/``) have been updated according to your changes. Python
  docstrings are formatted in the
  `numpydoc style <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
- Open a
  `pull request <https://help.github.com/articles/creating-a-pull-request/>`_
  with these changes. You pull request message ideally should include:

    - A description of why the changes should be made.
    - A description of the implementation of the changes.
    - A description of how to test the changes.

- The pull request should pass all the continuous integration tests which are
  automatically run by
  `GitHub Actions <https://github.com/hupo-psi/mzspeclib/actions>`_.
