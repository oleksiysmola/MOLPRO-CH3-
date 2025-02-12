point::Int64 = 7462
jobSplit::Int64 = 9
numberOfJobs::Int64 = 1000

for i in 1:numberOfJobs
    run(`qsub -e CH3+_MonteCarlo_MP_VQZ_$(point)-$(point + jobSplit).e -o CH3+_MonteCarlo_MP_VQZ_$(point)-$(point + jobSplit).o -l h_rt="11:59:00" -l mem=20G -l tmpfs=100G RunMolproJobs.csh $(point) $(point + jobSplit)`)
    global point = point + jobSplit
    println(point)
end