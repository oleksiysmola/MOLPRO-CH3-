#!/bin/csh
#

cd $TMPDIR

module load molpro/2020.1/openmp

set point = $1
set endPoint = $2
echo $point 
while ($point < $endPoint)
    set line = `awk 'NR=='"$point" /scratch/scratch/zcaposm/CH3+/MOLPRO/Grids/MonteCarlo/GridMonteCarlo.txt`

    if ("$line" != "") then
        csh -f /home/zcaposm/Scratch/CH3+/MOLPRO/Grids/MonteCarlo/GenerateMolproScriptMonteCarlo.csh $line
    endif
@ point++
end