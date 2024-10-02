# converge_refmac
wrapper script for easy command-line refmac refinement with extensible options. 

By default, refinement is repeated until atoms stop moving. Easiest thing to do
is simply provide a *.pdb and *.mtz file on the command line and it will do the 
rest. Output file is always refmacout.pdb and refmacout.mtz. Best Rwork and Rfree
trials are automatically saved, as are the last 10 iterations. Log files sumarize
the results:
refmac_Rplot.txt : appended with the most recent R factors and geometry stats.
refmac_shifts.txt : appened with largest shifts since the last iteration
refmac_scales.txt : scale factors for model and bulk solvent
Partial structures for externally supplied bulk solvent are automatically detected
and incorporated. Ligands supplied as any number of *.cif files on the command line 
are automatically run through libcheck and incorporated. Other options include:
F000 estimation
anomalous refinement
pruning of atoms with too-high B factors
pruning of atoms that refine to zero occupancy
limit maximum shifts of xyz, B or occupancy with either clipping or overall scale-back
incremental nudge to occupancy for atoms with highest and lowest B factors.
Support for maximum run time for use in cluster environments.


