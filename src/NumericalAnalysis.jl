module NumericalAnalysis
export romberg, int2d
export jacobi, jacobi!, eigenValues, eigenVectors

include("Integrate.jl")
using NumericalAnalysis.Integrate

include("Eigen.jl")
using NumericalAnalysis.Eigen

end