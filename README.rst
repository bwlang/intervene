.. image:: https://raw.githubusercontent.com/asntech/intervene/master/docs/img/intervene_logo.png
   	:target: http://intervene.readthedocs.org
   	
Intervene
-----------

	a tool for intersection and visualization of multiple gene or genomic region sets



.. image:: https://img.shields.io/pypi/pyversions/intervene.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/intervene.svg
    :target: https://pypi.python.org/pypi/intervene

.. image:: https://anaconda.org/bioconda/intervene/badges/version.svg
	:target: https://anaconda.org/bioconda/intervene

.. image:: https://static.pepy.tech/badge/intervene
	:target: https://pypi.org/project/intervene/

.. image:: https://anaconda.org/bioconda/intervene/badges/downloads.svg
    :target: https://bioconda.github.io/recipes/intervene/README.html

.. image:: https://img.shields.io/github/issues/asntech/intervene.svg
	:target: https://github.com/asntech/intervene/issues

.. image:: https://img.shields.io/twitter/url/https/github.com/asntech/intervene.svg?style=social
	:target: https://twitter.com/intent/tweet?text=Intervene%20-%20a%20tool%20for%20intersection%20and%20visualization%20of%20multiple%20genomic%20region%20and%20gene%20sets%20https://github.com/asntech/intervene&url=%5Bobject%20Object%5D

Documentation
=============

**A detailed documentation is available in different formats:**  `HTML <http://intervene.readthedocs.org>`_ | `PDF <http://readthedocs.org/projects/intervene/downloads/pdf/latest/>`_ | `ePUB <http://readthedocs.org/projects/intervene/downloads/epub/latest/>`_

.. figure:: http://intervene.readthedocs.io/en/latest/_images/Intervene_sketch.png
   :width: 800px
   :align: left

Installation
============

Quick installation using Conda
------------------------------

.. code-block:: bash

	conda install -c bioconda intervene

This will install all the dependencies and you are ready to use Intervene.

Install using `pip`
-------------------
You can install Intervene from PyPi using pip.

.. code-block:: bash

	pip install intervene

Note: If you install using pip, make sure to install BEDTools and the R packages listed below.

Intervene requires the following Python modules and R packages:

	* Python (>= 3.10): https://www.python.org/
	* BedTools (latest): https://github.com/arq5x/bedtools2
	* pybedtools (>= 0.12.0): https://daler.github.io/pybedtools/
	* Pandas (>= 3.0): http://pandas.pydata.org/
	* Seaborn (>= 0.13): http://seaborn.pydata.org/
	* R (>= 3.0): https://www.r-project.org/
	* R packages: UpSetR (>= 1.4.0), corrplot

Install BEDTools
----------------
We are using pybedtools, which is a Python wrapper for BEDTools. BEDTools must be installed before using Intervene.

.. code-block:: bash

    conda install -c bioconda bedtools

Please read the instructions at https://github.com/arq5x/bedtools2 to install BEDTools, and make sure it is on your path.


Install required R packages
---------------------------

Intervene requires the R packages `UpSetR <https://cran.r-project.org/package=UpSetR>`_ and `corrplot <https://cran.r-project.org/package=corrplot>`_ for visualization, and `Cairo <https://cran.r-project.org/package=Cairo>`_ for high-quality vector and bitmap output.

.. code-block:: R

    install.packages(c("UpSetR", "corrplot", "Cairo"))

Install Intervene from source
=============================

.. code-block:: bash

    git clone https://github.com/asntech/intervene.git
    cd intervene
    pip install -e .

Development environment with pixi
----------------------------------

`pixi <https://pixi.sh>`_ manages all Python, R, and system dependencies in a reproducible environment.

.. code-block:: bash

    git clone https://github.com/asntech/intervene.git
    cd intervene
    pixi run dev-install

Useful development tasks:

.. code-block:: bash

    pixi run test        # run a basic smoke test
    pixi run lint        # run ruff linter
    pixi run typecheck   # run pyright type checker
    pixi run check       # lint + typecheck together

Releasing a new version
-----------------------

1. Update the version in ``intervene/__init__.py``
2. Commit and tag the release:

   .. code-block:: bash

       git commit -am "release vX.Y.Z"
       git tag vX.Y.Z

3. Build and publish to PyPI:

   .. code-block:: bash

       pixi run build     # builds sdist and wheel into dist/
       pixi run publish   # uploads to PyPI via twine

How to use Intervene
====================
Once you have installed Intervene, you can type:

.. code-block:: bash

	intervene --help

This will show the following help message.

.. code-block:: bash

	usage: intervene <subcommand> [options]
	    
	positional arguments <subcommand>:
	  {venn,upset,pairwise}
	                        List of subcommands
	    venn                Venn diagram of intersection of genomic regions or list sets (upto 6-way).
	    upset               UpSet diagram of intersection of genomic regions or list sets.
	    pairwise            Pairwise intersection and heatmap of N genomic region sets in <BED/GTF/GFF> format.

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit


to see the help for the three subcommands ``pairwise``, ``venn`` and ``upset`` type:

.. code-block:: bash
	
	intervene pairwise --help

	intervene venn --help

	intervene upset --help

Run Intervene on test data
--------------------------

To run Intervene using example data, use the following commands. To access the test data make sure you have ``sudo`` or ``root`` access.

.. code-block:: bash

	intervene pairwise --test

	intervene venn --test

	intervene upset --test

If you have installed Intervene locally from the source code, you may have problem to find test data. You can download the test data here https://github.com/asntech/intervene/tree/master/intervene/example_data and point to it using ``-i`` instead of ``--test``.

.. code-block:: bash

	./intervene/intervene venn -i intervene/example_data/ENCODE_hESC/*.bed       
  	./intervene/intervene upset -i intervene/example_data/ENCODE_hESC/*.bed      
  	./intervene/intervene pairwise -i intervene/example_data/dbSUPER_mm9/*.bed  

The above three test commands will generate the following three figures (a, b and c).

.. figure:: http://intervene.readthedocs.io/en/latest/_images/Intervene_plots.png
   :width: 800px
   :align: left

By default your results will stored in the current working directory with a folder named ``Intervene_results``. If you wish to save the results in a specific folder, you can type::

	intervene upset --test --output ~/path/to/your/folder

Interactive Shiny App
=====================
Intervene Shiny App is freely available at https://asntech.shinyapps.io/intervene or https://intervene.shinyapps.io/intervene

The source code for the Shiny app is available at https://github.com/asntech/intervene-shiny

Support
========
If you have questions, or found any bug in the program, please write to us at ``azez.khan[at]gmail.com``

Cite Us
=========
If you use Intervene please cite us: ``Khan A, Mathelier A. Intervene: a tool for intersection and visualization of multiple gene or genomic region sets. BMC Bioinformatics. 2017;18:287. doi: 10.1186/s12859-017-1708-7``

