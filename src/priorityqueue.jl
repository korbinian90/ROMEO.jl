mutable struct PQueue{T}
    min::Int
    nbins::Int
    content::Vector{Vector{T}}
end

# initialize new queue
PQueue{T}(nbins) where T = PQueue(nbins + 1, nbins, [Vector{T}() for _ in 1:nbins])
PQueue(nbins, item::T) where T = push!(PQueue{T}(nbins), item, 1)

Base.isempty(q::PQueue) = q.min > q.nbins

function Base.push!(q::PQueue, item, weight)
    push!(q.content[weight], item)
    q.min = min(q.min, weight)
    return q
end

function Base.pop!(q::PQueue)
    elem = pop!(q.content[q.min])
    # increase smallestbin, if elem was last in bin
    while q.min ≤ q.nbins && isempty(q.content[q.min])
        q.min += 1
    end
    return elem
end
