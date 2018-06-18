pdprop
======

Simulation of metal fatigue crack propagation with accounting for material memory effects.

The software can be used to simulate crack propagation in plates and pipes.
At this time the following geometries are incluced

    1. A full width surface crack in a plate 
    2. Elliptical surface crack in a plate
    3. Internal elliptical surface crack in pipe
    4. Center thru crack in a plate
    5. Edge thru crack in a plate
    6. Rod surface crack
    
Loading can be specified with a variable amplitude bending and membrane (axial) stress history.
The sequence effects of the stress history is modelled by rules that mimic the material's memory
of prior stress-strain behavior.  Estimation of the Stress Intensity follows the rules
set out by British Standard  BS7910.   For elliptical shaped cracks  a separate memory model
is run for each of the  "a" and "c"  crack lengths;  i.e. where "a" is the depth of the
crack,   and "c" is the 1/2 width of the crack on the surface.

  A example set of material data is provided:  A36 steel
with da/dn curves, crack initiation data, and a short stress history.  One of the three
problem folders also contains a sample output report generated from a simulation run.
The large output file created by a simulation is not, however, included in the folder.

Prerequisites: gfortran (a subset of gcc), gnuplot,  htmldoc,  bash, Linux

The code is contained in a compressed tar file.   It was developed on openSuse Linux.
There are some  bash scripts that help the user during problem setup, for final
report generation,  and for initial compile of all the software.  All of the runtime
code is written in  gfortran.   All code is GPL licensed.

A readme file in the tar file will explain how to setup the files and run the compile
script.  Each of the three crack type geometries is located in its own folder and a setup
script is provided to help the user with the problem definition.

It is expected that the user has an understanding of metal fatigue and the terms used in
BS7910.

Note:  Some example result files may be out of data due to software changes.

For reference on how the material memory logic works the reader is referred to:
Conle, A., T.R.Oxland, T.H.Topper, "Computer-Based Prediction of Cyclic Deformation
    and Fatigue Behavior," Low Cycle Fatigue ASTM STP 942, 1988, pp.1218-1236
    
An example of how both the crack initiation and the propagation software can be
applied to a fatigue life prediction case history is provided here:
  http://fde.uwaterloo.ca/Fde/Crackgrowth/Case5/case5_HeulerSeeger.html

A tutorial on fatigue testing and calculation can be found at
  http://fde.uwaterloo.ca/FatigueClass/classIndex.html


