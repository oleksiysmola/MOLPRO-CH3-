using Printf

stretchesSpacing::Vector{Float64} = [0.0000, -0.025,0.025, 0.050, -0.0500, -0.1000, 0.1000, 0.12500, -0.125000, 0.1500, -0.1500, 0.200, 
   -0.200, 0.300, -0.300, 0.400, 0.500, 0.600, 0.700]
dihedralSpacing::Vector{Float64} = [0.0000, 1.000, -1.000, 2.500, -2.5000, 5.0000, -5.0000, 7.500, -7.500, 10.0000, -10.0000, 20.0000, 
    -20.0000, 30.0000, -30.0000, 50.0000, -50.0000, 60.0000]
inversionSpacing::Vector{Float64} = [0.0000, 1.000,  2.5000,  5.0000,  7.500, 10.0000, 15.000,  20.00, 
    30.000, 40.0000,  60.0000, 80.0000]

stretchesGrid::Int64 = size(stretchesSpacing)[1]
dihedralGrid::Int64 = size(dihedralSpacing)[1]
inversionGrid::Int64 = size(inversionSpacing)[1]
convertToRadians::Float64 = 2*pi/360

gridSize::Int64 = 9
SaEq::Float64 = 0.0
SbEq::Float64 = 0.0
rhoEq::Float64 = 90.0
rCHeq::Float64 = 1.08643631
aCHeq::Float64 = 120.00
tauEq::Float64 = 180.00
equilibriumGrid::Vector{Float64} = [rCHeq, rCHeq, rCHeq, aCHeq, aCHeq, tauEq, SaEq, SbEq, rhoEq]

function PrintGeometry(point::Int64, grid::Vector{Float64})
    @printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", point, grid[1], grid[2], grid[3], grid[4], grid[5], grid[6], grid[7], grid[8], grid[9])
end

point::Int64 = 101

for i in 1:3
    for k in i+1:3
        for j in 2:stretchesGrid
            displacementVector::Vector{Float64} = zeros(gridSize)
            displacementVector[i] = stretchesSpacing[j]
            for l in 2:stretchesGrid
                displacementVector[k] = stretchesSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
    for k in 7:8
        for j in 2:stretchesGrid
            displacementVector::Vector{Float64} = zeros(gridSize)
            displacementVector[i] = stretchesSpacing[j]
            for l in 2:dihedralGrid
                displacementVector[k] = dihedralSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
                d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
                grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
                grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
    for j in 2:stretchesGrid
        displacementVector::Vector{Float64} = zeros(gridSize)
        displacementVector[i] = stretchesSpacing[j]
        for l in 2:inversionGrid
            displacementVector[9] = inversionSpacing[l]
            grid::Vector{Float64} = equilibriumGrid + displacementVector
            d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
            d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
            grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
            grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
            grid[6] = acos(csc(grid[4]*convertToRadians)*csc(grid[5]*convertToRadians)*((sin(2*grid[9]*convertToRadians)*sin(d12/2)*sin(d13/2))^2 - sin(d12)*sin(d13)*sin(grid[9]*convertToRadians)^2))/convertToRadians
            PrintGeometry(point, grid)
            global point = point + 1
        end
    end
end
for i in 7:8
    for k in i+1:8
        for j in 2:dihedralGrid
            displacementVector::Vector{Float64} = zeros(gridSize)
            displacementVector[i] = dihedralSpacing[j]
            for l in 2:dihedralGrid
                displacementVector[k] = dihedralSpacing[l]
                grid::Vector{Float64} = equilibriumGrid + displacementVector
                d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
                d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
                grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
                grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
    for j in 2:dihedralGrid
        displacementVector::Vector{Float64} = zeros(gridSize)
        displacementVector[i] = dihedralSpacing[j]
        for l in 2:inversionGrid
            displacementVector[9] = dihedralSpacing[l]
            grid::Vector{Float64} = equilibriumGrid + displacementVector
            d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
            d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
            grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
            grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
            grid[6] = acos(csc(grid[4]*convertToRadians)*csc(grid[5]*convertToRadians)*((sin(2*grid[9]*convertToRadians)*sin(d12/2)*sin(d13/2))^2 - sin(d12)*sin(d13)*sin(grid[9]*convertToRadians)^2))/convertToRadians
            PrintGeometry(point, grid)
            global point = point + 1
        end
    end
end