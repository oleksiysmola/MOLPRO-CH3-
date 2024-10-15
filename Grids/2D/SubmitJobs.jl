point::Int64 = 3256
jobSplit::Int64 = 5
numberOfJobs::Int64 = 1

for i in 1:numberOfJobs
    run(`qsub -e CH3+_2D_MP_$(point)-$(point + jobSplit).e -o CH3+_2D_MP_$(point)-$(point + jobSplit).o -l h_rt="11:59:00" RunMolproJobs.csh $(point) $(point + jobSplit)`)
    global point = point + jobSplit
    println(point)
end