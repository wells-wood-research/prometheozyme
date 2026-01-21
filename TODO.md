[] TODO add molecular orbital diagrams of amino acids' HOMOs and LUMOs to guide decisions on restraint angles (tutorial style)
[] TODO evaluate angle restraints
[] TODO evaluate reasonable backbone arrangement - what testable rules are there for reasonable arrangement?
[] TODO automatically write inputs to evaluate reactivity (e.g. with NEB) and generate protein (RFdiffusion)
[] TODO source of each individual final result is unknown because of renaming based on einter - write as REMARKs metadata information?
[] TODO add a slightly greasy polarisable medium (mimics protein environment) so lack of H-bonds isn't penalised so much
[] TODO clearer error messaging - the errors are there, but tend to be hidden by the code that comes next

[] TODO deal with this geom opt issue:
[] 
[]           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[]           !   SERIOUS PROBLEM WITH INTERNALS - ANGLE IS APPROACHING 180 OR 0 DEGREES   !
[]           !                       REBUILDING A NEW SET OF INTERNALS                    !
[]           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[] 
[] Making redundant internal coordinates   ... (2022 redundants).... [file orca_tools/Tool-Misc/qcredint.cpp, line 6411]: Error (ORCA_QCREDINT): For a linear bend you can only keep the initial value of your structure as constraint!!!
[] Options: use COPT (or coordsys cartesian in %geom)
[] Optimise per step
[] Maximum number of scan parameters is 3, cannot do 2 angles and 2 bonds
[] TODO ^ only optimise the restraints per this course, since host is frozen; keep previous geoms in keep but scan all new restrains; this way only the new restraints per host are limited to 3
[] TODO unfinished opts get copied - see last example, results should be 4 but there's 6
[] TODO implement multi-substrate docking (e.g. different conformers) properly
[x] TODO implement sampling conformations of substrate
[x] TODO implement custom output folder names
[x] TODO fix RMSD breaking at comparison between products from different ingredients (different number of atoms)
[] TODO add comment about fixHost being bad for any further than first course after init
[] TODO add warnings in logs (e.g. when fixHost is False past the first step)
[x] Don't fail on too small distances (unless specified)
[x] TODO is result of copt.xyz used in next steps?
[] Save at least ine example structure (per ingredient combination) in 'specials'
[] Grouping products by number of atoms is not enough for RMSD - different amino acids can have same number of atoms