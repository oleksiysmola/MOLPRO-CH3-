#!/bin/csh
#
# Generation of the input file for MOLPRO
#

set pwd = `pwd`

set point = $1        
set directory = /scratch/scratch/zcaposm/CH3+/MOLPRO/Grids/1D/CH3+_1D_GRIDS
set fname = CH3+_1D_point_${point}

cat<<endb> ${fname}.inp
***, CH3+ Ground State Energy with CCSD(T)-F12 and cc-pVTZ-F12
memory,500,m;

geometry={angstrom
c 
h1 , 1, rch1 
h2 , 1, rch2, 2, ach12
h3 , 1, rch3, 2, ach23, 3, tau
}

rch1 = $2
rch2 = $3
rch3 = $4
ach12 = $5
ach23 = $6
tau = $7
Sa = $8
Sb = $9
rho = $10

! Set charge to +1
set,nelec=8
set,charge=+1

! Use the cc-pVTZ-F12 basis set
basis=cc-pVTZ-F12

hf

! Use explicitly correlated F12 methods
! First, MP2-F12 (useful for initial electronic energy)
{mp2-f12}

! If desired, perform CCSD(T)-F12 for more accurate results
{ccsd(t)-f12}

! Output the energy
xxx = "mmm"
point = ${point}
text ### CH3+
table,xxx,rch1,rch2,rch3,Sa,Sb,rho,energy,point
DIGITS, 0, 8, 8, 8, 8, 8, 8, 8, 4
save,CH3+_MP2_${point}.dat,new

--- End of Script ---
endb

module load molpro/2020.1/openmp
molpro ${fname}.inp
rm ${fname}.inp
cp ${fname}.out ${directory}
cp CH3+_MP2_${point}.dat ${directory}
