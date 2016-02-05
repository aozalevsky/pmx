pmx
===

Introduction
------------

pmx (formerly pymacs) has started as a small bunch of classes to read
structure files such as pdb or gro and trajectory data in gromacs xtc
format. Over the years it has been extended towards a versatile (bio-)
molecular structure manipulation package with some additional
functionalities, e.g. gromacs file parsers and scripts for setup and
analysis of free energy calculations.

Citations:
----------

D. Seeliger and Bert L. de Groot, Biophys. J. 98(10):2309-2316 (2010)

V. Gapsys, S. Michielssens, D. Seeliger, B. L. de Groot. J. Comput.
Chem. 2014, DOI: 10.1002/jcc.23804

Purpose
-------

I mostly use pmx to write short scripts which perform some changes in
pdb files, e.g. changing atom or residue names, applying some geometric
transformations or doing some kind of analysis. The critical issue for
these things is usually not calculation time but straightforward
selection of some atoms/residues of interest, quick file parsing and
data visualization, which renders a well-organized data structure and
easy programming style more important than computation performance.
Hence, and ideal task for Python.

Installation
------------

Checkout the source code and run the usual python installation
``git clone https://github.com/dseeliger/pmx/ pmx cd pmx  sudo python setup.py install``

Software Requirements
---------------------

-  `numpy <http://numpy.scipy.org/>`_
-  `scipy <http://www.scipy.org/>`_
-  `matplotlib <http://matplotlib.org/>`_ ( for analysis scripts )

Getting Started
---------------

pmx stores structure data in Python classes. The "Model" class is the
uppermost class which contains severals lists of Atoms, Molecules and
Chains. The following script reads a pdb file, prints some atom
properties and writes the structure in gro format.

The figure above shows the most important data structure. A "Model"
instance contains list of chains, residues and atoms. A "Chain" instance
of residues and atoms and a "Molecule" instance of a list of atoms only.
Check the example scripts for how to navigate trough particular storage
classes.

pmx Modules
-----------

pmx currently contains the following modules. Click on the links below
for a (short) documentation

-  `atom <https://github.com/dseeliger/pmx/wiki/pmx_atom.md>`_
-  `molecule <https://github.com/dseeliger/pmx/wiki/pmx_molecule.md>`_
-  `chain <https://github.com/dseeliger/pmx/wiki/pmx_chain.md>`_
-  `model <https://github.com/dseeliger/pmx/wiki/model.md>`_
-  `atomselection <https://github.com/dseeliger/pmx/wiki/pmx_atomselection.md>`_
-  `options <https://github.com/dseeliger/pmx/wiki/pmx_options.md>`_
-  `library <https://github.com/dseeliger/pmx/wiki/pmx_library.md>`_
-  `ndx <https://github.com/dseeliger/pmx/wiki/pmx_ndx.md>`_
-  `geometry <https://github.com/dseeliger/pmx/wiki/pmx_geometry.md>`_
-  `parser <https://github.com/dseeliger/pmx/wiki/pmx_parser.md>`_
-  `xtc <https://github.com/dseeliger/pmx/wiki/pmx_xtc.md>`_
-  `builder <https://github.com/dseeliger/pmx/wiki/pmx_builder.md>`_
-  `forcefield <https://github.com/dseeliger/pmx/wiki/pmx_forcefield.md>`_
-  `rotamer <https://github.com/dseeliger/pmx/wiki/pmx_rotamer.md>`_
-  `ffparser <https://github.com/dseeliger/pmx/wiki/pmx_ffparser.md>`_

pmx Classes (the most important ones)
-------------------------------------

-  `Atom <https://github.com/dseeliger/pmx/wiki/pmx_atom.md>`_
-  `Molecule <https://github.com/dseeliger/pmx/wiki/pmx_molecule.md>`_
-  `Chain <https://github.com/dseeliger/pmx/wiki/pmx_chain.md>`_
-  `Model <https://github.com/dseeliger/pmx/wiki/model.md>`_
-  `Trajectory <https://github.com/dseeliger/pmx/wiki/pmx_trajectory.md>`_
-  `Topology <https://github.com/dseeliger/pmx/wiki/pmx_topology.md>`_
-  `IndexFile <https://github.com/dseeliger/pmx/wiki/pmx_indexfile.md>`_

Using pmx
---------

-  `Editing structure files, selecting atoms, residues, chain,
   etc. <https://github.com/dseeliger/pmx/wiki/editstruct.md>`_
-  `The pmx commandline parser, documenting programs and
   scripts <https://github.com/dseeliger/pmx/wiki/pmx_options.md>`_
-  `Building peptides and DNA strands from
   scratch <https://github.com/dseeliger/pmx/wiki/builder.md>`_
-  `Writing MD analysis tools, access xtc data from
   Python <https://github.com/dseeliger/pmx/wiki/trx.md>`_
-  `Modifying Gromacs
   topologies <https://github.com/dseeliger/pmx/wiki/topol.md>`_
-  `Generating Gromacs index
   files <https://github.com/dseeliger/pmx/wiki/pmx_indexfile.md>`_

