module SpatialIndexing
    using GeometryTypes: HyperCube, Vec3
    using LinearAlgebra: ⋅, norm
    sqrnorm(x) = x ⋅ x

    include("Octree.jl")

    export
        Octree,
        radius_search,
        radius_search_recursive
end
