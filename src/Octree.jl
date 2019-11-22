struct OctreeNode{T}
    center::Vec3{T}
    extent::T
    firstChild::Int32
    firstPoint::Int32
    lastPoint::Int32
end

@inline Base.length(node::OctreeNode) = node.lastPoint - node.firstPoint + 1
@inline with_child(node::OctreeNode{T}, child::Int) where {T} =
    OctreeNode{T}(node.center, node.extent, child, node.firstPoint, node.lastPoint)
@inline isleaf(node::OctreeNode) = node.firstChild == 0

mutable struct Octree{T}
    points::Vector{Vec3{T}}
    indices::Vector{Int32}
    cells::Vector{OctreeNode{T}}
end

leaves(oct::Octree) = filter(x->x.firstChild == 0, oct.cells)

function split_cell!(tree::Octree{T}, idx::Int) where {T}
    cell = tree.cells[idx] = with_child(tree.cells[idx], length(tree.cells) + 1)
    c = cell.center

    bins = ntuple(_->Int32[], 8)
    for i in cell.firstPoint:cell.lastPoint
        idx = tree.indices[i]
        delta = tree.points[idx] .< c
        bindex = (delta[1] | delta[2] << 1 | delta[3] << 2) + 1
        push!(bins[bindex], idx)
    end
    tree.indices[cell.firstPoint:cell.lastPoint] = vcat(bins...)

    current_index = cell.firstPoint
    ha = cell.extent / 2
    for i in eachindex(bins)
        code = i - 1
        nc = cell.center + ha * ((code .& Vec3(1,2,4) .== 0) * T(2) .- T(1))

        start_of_next = current_index + length(bins[i])
        push!(tree.cells, OctreeNode{T}(nc, ha, 0, current_index, start_of_next-1))
        current_index = start_of_next
    end
    nothing
end

function bounds_cube(pts)
    bbmin, bbmax = extrema(pts)
    center = (bbmin + bbmax) / 2
    width = maximum(center - bbmin)
    (center, width)
end

function Octree(data::AbstractVector{Vec3{T}}; leafsize::Int = 16) where {T}
    result = Octree{T}(copy(data), eachindex(data), OctreeNode[])
    bounds = bounds_cube(data)
    root = OctreeNode{T}(
        bounds[1],
        bounds[2],
        0,
        1,
        length(data))
    push!(result.cells, root)

    split_queue = Int[]
    if length(root) > leafsize
        push!(split_queue, 1) end

    while !isempty(split_queue)
        split_idx = popfirst!(split_queue)
        split_cell!(result, split_idx)

        splits = Int[]
        for i in eachindex(result.cells)[end-7:end]
            if length(result.cells[i]) > leafsize
                push!(splits, i) end end
        if length(splits) > 1
            append!(split_queue, splits)
        end
    end

    result.points .= data[result.indices]
    result
end

function radius_search_kernel!(idcs, tree, cell_idx, point, radius2)
    function in_searchsphere(cell)
        qprime = abs.(point - cell.center)
        sqrnorm(qprime .+ cell.extent) <= radius2
    end

    function overlaps_searchsphere(cell)
        qprime = abs.(point - cell.center)
        (minimum(qprime) <= cell.extent) || (sqrnorm(qprime .- cell.extent) <= radius2)
    end

    cell = tree.cells[cell_idx]
    if in_searchsphere(cell)
        append!(idcs, view(tree.indices, cell.firstPoint:cell.lastPoint))
    elseif isleaf(cell)
        for i in cell.firstPoint:cell.lastPoint
            if sqrnorm(point - tree.points[i]) <= radius2
                push!(idcs, tree.indices[i]) end end
    else
        for i in cell.firstChild:(cell.firstChild+7)
            if overlaps_searchsphere(tree.cells[i])
                radius_search_kernel!(idcs, tree, i, point, radius2) end end
    end
    idcs
end

function radius_search_recursive(tree::Octree{T}, point, radius) where {T}
    idcs = Int[]
    radius2 = radius * radius
    radius_search_kernel!(idcs, tree, 1, point, radius2)
end

function radius_search(tree::Octree{T}, point, radius) where {T}
    idcs = Int32[]
    search_queue = push!(sizehint!(Int32[], 32), 1)
    radius2 = radius * radius

    @inline function overlaps(cell, point, radius2)
        qprime = abs.(point - cell.center)
        (minimum(qprime) <= cell.extent) || (sqrnorm(qprime .- cell.extent) <= radius2)
    end
    stat = 0

    while !isempty(search_queue)
        cell = tree.cells[pop!(search_queue)]
        qprime = abs.(point - cell.center)
        if sqrnorm(qprime .+ cell.extent) <= radius2
            stat += length(cell.firstPoint:cell.lastPoint)
            append!(idcs, cell.firstPoint:cell.lastPoint)
        elseif isleaf(cell)
             for i in cell.firstPoint:cell.lastPoint
                if sqrnorm(point - tree.points[i]) <= radius2
                    push!(idcs, i) end end
        else
            for i in cell.firstChild:(cell.firstChild+7)
                if overlaps(tree.cells[i], point, radius2)
                    push!(search_queue, i) end end
        end
    end

    idcs .= @view tree.indices[idcs]
    idcs
end
