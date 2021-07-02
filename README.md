# Amyloid_fibril_width
Repository for David Boyer's AB Structure Project. Written by David Boyer and Nikos Mynhier. 

# Setting up an Environment
Dependencies (w/ versions) are listed in the environment.yml file. Install these dependencies using a conda environment. 

```
conda env create -f /path/to/gitclone/environment.yml
```

You will then want to activate the environment. 

```
conda activate amyloid_fibril_width
```

# Running the Workflow
In order to analyze fibril i and i+1 distances between symmetrically related atoms, calculate beta-sheet hydrogen bond lengths and off-axis tilts, and graph these results, please see the following instructions.

The input files needed are: (i) a helicized model of your amyloid fibril(s), (ii) a helicized model of your amyloid fibril(s) with hydrogens added (e.g., using phenix.reduce), and (iii) a helical parameters file listing the pdb id of the structure(s) and the helical twist and rise of the structure(s)

The two pdb files should be placed in different directories (e.g., pdb_noh/ and pdb_h/). You also need to make an output directory for where you would like the graphs to put written. 

In this example, we have:

pdb_noh/
pdb_h/
graphs/
helical_parameters.txt

In the first directory, we have the input pdb without hydrogens and the output csv file with interstrand distances and radii. In the second directory, we have the pdb file with hydrogens and the output files including the hydrogen bond lengths, off-axis tilts, and radii. In the third directory, there is the graphing output. To run the programs to produce the output files, specify in the macro_hbond.sh file the correct input directories, directory containing downloaded scripts, helical parameters, the graphing output directory, and the number of cpus to use for GNU's parallel program. 

The outputs can then be calculated by running the command: sh macro_hbond.sh

If you would like to add hydrogen bond length and planarity restraints during real space refinement with phenix.real_space_refine, you may use distrecord2phenix_withsidechains_planarity3.com on your structure. This will produce a .edits file that can be used during real space refinement of your structure using phenix.real_space_refine. distrecord2phenix_withsidechains_planarity3 relies on Fortran and hbplus (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/) to run.  

If you notice your structure does not follow helical symmetry, perhaps due to lack of helical symmetry employed during model refinement, you may use the program helical.com to helicize your structure. This program relies on pdbset from the CCP4 suite and the Fortran program beyond26chainsdiscont. In order to run this program successfully, you may need to remove all lines except for the ATOM and CRYST1 lines of the pdb file. In addition, it may be necessary for asymmetric fibrils with more than one protofilament to include all chains in the layer of the fibril to be symmetrized in one chain ID. At the moment, the program can only take one chain id as input.
