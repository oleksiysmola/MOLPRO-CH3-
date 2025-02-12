using Printf
using Random


stretchesSpacing::Vector{Float64} = [0.0000, -0.025,0.025, 0.050, -0.0500, -0.1000, 0.1000, 0.12500, -0.125000, 0.1500, -0.1500, 0.200, 
   -0.200, 0.300, -0.300, 0.400, 0.500, 0.600]
dihedralSpacing::Vector{Float64} = [0.0000, 1.000, -1.000, 2.500, -2.5000, 5.0000, -5.0000, 7.500, -7.500, 10.0000, -10.0000, 20.0000, 
    -20.0000, 30.0000, -30.0000, 50.0000, -50.0000]
inversionSpacing::Vector{Float64} = [0.0000, 0.500, 1.000, 1.500,  2.000, 2.5000,  3.000, 3.5000, 4.000, 4.500, 5.0000, 6.000, 7.000, 7.500, 8.000, 9.00, 10.0000, 12.50, 15.000, 17.500,  20.00, 
    25.00, 30.000, 40.0000]

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

point::Int64 = 5412

numberOfGridPoints3D::Int64 = 5000
numberOfGridPoints4D::Int64 = 5000
numberOfGridPoints5D::Int64 = 5000
numberOfGridPoints6D::Int64 = 5000
degreesOfFreedom::Int64 = 6
equilibriumProbabilities::Vector{Float64} = [1/3, 1/3, 1/3, 0.5, 0.5, 0]

function GenerateMonteCarloGrid(dimension::Int64, numberOfGridPoints::Int64, maxGridRange::Float64)
    generatedPoint::Int64 = 1
    while generatedPoint <= numberOfGridPoints
        chosenCoordinates::Vector{Int64} = zeros(dimension)
        equilibriumProbabilitiesOfChosenCoordinates::Vector{Float64} = zeros(dimension)
        currentCoordinate::Int64 = 1
        while currentCoordinate <= dimension
            newCoordinate::Int64 = rand(1:degreesOfFreedom)
            if newCoordinate in chosenCoordinates
                continue
            else
                chosenCoordinates[currentCoordinate] = newCoordinate
                equilibriumProbabilitiesOfChosenCoordinates[currentCoordinate] = equilibriumProbabilities[newCoordinate]
                currentCoordinate += 1
            end
        end
        probabilitiesOfChosenCoordinates::Vector{Float64} = rand(dimension)
        if prod(abs.(probabilitiesOfChosenCoordinates .- equilibriumProbabilitiesOfChosenCoordinates)) > maxGridRange^dimension
            continue
        else
            grid::Vector{Float64} = equilibriumGrid .+ 0.0
            dihedralDisplacement::Vector{Float64} = [0.0, 0.0, 0.0]
                for j in 1:size(chosenCoordinates)[1]
                    if 1 <= chosenCoordinates[j] < 4
                        stretchDisplacement::Float64 = 0
                        if probabilitiesOfChosenCoordinates[j] < 1/3
                            # PDF p(x) = m1 x + c  (values given not normalised)
                            # c = 0.3
                            # m1 = 1
                            stretchDisplacement = 3*(-1 + sqrt(3)*sqrt(probabilitiesOfChosenCoordinates[j]))/10
                            grid[chosenCoordinates[j]] += stretchDisplacement
                        else
                            # PDF p(x) = m2 x + c  (values given not normalised)
                            # c = 0.3
                            # m2 = -1/2
                            stretchDisplacement = 3*(2-sqrt(6)*sqrt(1-probabilitiesOfChosenCoordinates[j]))/10
                            grid[chosenCoordinates[j]] += stretchDisplacement
                        end
                    elseif 4 <= chosenCoordinates[j] < 6
                        angleDisplacement::Float64 = real((25*(1-sqrt(3)*1im))/(-1+2*probabilitiesOfChosenCoordinates[j]+2*sqrt(Complex(-probabilitiesOfChosenCoordinates[j]+probabilitiesOfChosenCoordinates[j]^2)))^(1/3)+25*(1+sqrt(3)*1im)*(-1+2*probabilitiesOfChosenCoordinates[j]+2*sqrt(Complex(-probabilitiesOfChosenCoordinates[j]+probabilitiesOfChosenCoordinates[j]^2)))^(1/3))
                        dihedralDisplacement[chosenCoordinates[j] - 3] += angleDisplacement
                    else
                        inversionDisplacement::Float64 = 40*(1-sqrt(1-probabilitiesOfChosenCoordinates[j]))
                        dihedralDisplacement[3] += inversionDisplacement
                    end
                end
            grid[7:9] = dihedralDisplacement
            grid[9] += equilibriumGrid[9]
            d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
            d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
            grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
            grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
            if dihedralDisplacement[3] < 0.0000000001
                grid[6] = 180.000
            else
                grid[6] = acos(csc(grid[4]*convertToRadians)*csc(grid[5]*convertToRadians)*((sin(2*grid[9]*convertToRadians)*sin(d12/2)*sin(d13/2))^2 - sin(d12)*sin(d13)*sin(grid[9]*convertToRadians)^2))/convertToRadians
            end
            PrintGeometry(point, grid)
            generatedPoint += 1
            global point += 1
        end
    end
end

GenerateMonteCarloGrid(3, numberOfGridPoints3D, 0.20)
GenerateMonteCarloGrid(4, numberOfGridPoints4D, 0.15)
GenerateMonteCarloGrid(5, numberOfGridPoints5D, 0.10)
GenerateMonteCarloGrid(6, numberOfGridPoints6D, 0.05)

# for i in 1:3
#     for k in i+1:3
#         for j in 2:stretchesGrid
#             displacementVector::Vector{Float64} = zeros(gridSize)
#             displacementVector[i] = stretchesSpacing[j]
#             for l in 2:stretchesGrid
#                 displacementVector[k] = stretchesSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 PrintGeometry(point, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for k in 7:8
#         for j in 2:stretchesGrid
#             displacementVector::Vector{Float64} = zeros(gridSize)
#             displacementVector[i] = stretchesSpacing[j]
#             for l in 2:dihedralGrid
#                 displacementVector[k] = dihedralSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
#                 d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
#                 grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#                 grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#                 PrintGeometry(point, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for j in 2:stretchesGrid
#         displacementVector::Vector{Float64} = zeros(gridSize)
#         displacementVector[i] = stretchesSpacing[j]
#         for l in 2:inversionGrid
#             displacementVector[9] = inversionSpacing[l]
#             grid::Vector{Float64} = equilibriumGrid + displacementVector
#             d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
#             d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
#             grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#             grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#             grid[6] = acos(csc(grid[4]*convertToRadians)*csc(grid[5]*convertToRadians)*((sin(2*grid[9]*convertToRadians)*sin(d12/2)*sin(d13/2))^2 - sin(d12)*sin(d13)*sin(grid[9]*convertToRadians)^2))/convertToRadians
#             PrintGeometry(point, grid)
#             global point = point + 1
#         end
#     end
# end
# for i in 7:8
#     for k in i+1:8
#         for j in 2:dihedralGrid
#             displacementVector::Vector{Float64} = zeros(gridSize)
#             displacementVector[i] = dihedralSpacing[j]
#             for l in 2:dihedralGrid
#                 displacementVector[k] = dihedralSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
#                 d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
#                 grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#                 grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#                 PrintGeometry(point, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for j in 2:dihedralGrid
#         displacementVector::Vector{Float64} = zeros(gridSize)
#         displacementVector[i] = dihedralSpacing[j]
#         for l in 2:inversionGrid
#             displacementVector[9] = inversionSpacing[l]
#             grid::Vector{Float64} = equilibriumGrid + displacementVector
#             d12::Float64 = 4*pi/6 - (grid[7]/sqrt(6) - grid[8]/sqrt(2))*convertToRadians
#             d13::Float64 = 4*pi/6 - (grid[7]/sqrt(6) + grid[8]/sqrt(2))*convertToRadians
#             grid[4] = acos(cos(d12)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#             grid[5] = acos(cos(d13)*sin(grid[9]*convertToRadians)^2 + cos(grid[9]*convertToRadians)^2)/convertToRadians
#             grid[6] = acos(csc(grid[4]*convertToRadians)*csc(grid[5]*convertToRadians)*((sin(2*grid[9]*convertToRadians)*sin(d12/2)*sin(d13/2))^2 - sin(d12)*sin(d13)*sin(grid[9]*convertToRadians)^2))/convertToRadians
#             PrintGeometry(point, grid)
#             global point = point + 1
#         end
#     end
# end