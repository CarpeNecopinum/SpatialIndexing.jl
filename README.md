# SpatialIndexing.jl
Spatial Indexing Structures for Julia

## Octree

For now the only spatial indexing structure implemented.

Good for radius searches with large radii.

```
using GeometryTypes
using SpatialIndexing

xs = rand(Vec3f0, 1024^2)

tree = Octree(xs; leafsize = 16)

radius_search(oct, xs[1], 0.5f0)           
radius_search_recursive(oct, xs[1], 0.5f0)  # usually faster than the non-recursive variant, but allocates more
```
