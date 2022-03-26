module Chaos

using LinearAlgebra

export integer_partition, integer_partition_redundant, hopping_pair, average_spacing, state_evolution, reduced_density_matrix

"""
This function will give partitions of m into at most n integers. Please notice
that this function will give us redundant partition.
"""
function integer_partition_redundant(m::Int, n::Int)
    partitions = []

    if n == 2
        if m != 0
            for i::Int in 1:m
                append!(partitions, [[i, m-i]])
            end
        else
            append!(partitions, [[0,0]])
        end
    end

    if n > 2
        for i ::Int in 1:m
            for pair in integer_partition_redundant(m-i, n-1)
                to_be_appended = [pushfirst!(pair, i)]
                append!(partitions, to_be_appended)
            end
        end
    end
   return partitions
end

"""
Like the redundant version, this function will instead produce non-redundant
partitions.
"""
function integer_partition(m::Int, n::Int)::Vector{Vector{Int64}}
    partitions = []
    partitions_set = []
    partitions_redundant = integer_partition_redundant(m,n)
    
    for partition in partitions_redundant
        partition_set = Set(partition)
        append!(partitions_set, [partition_set])
    end

    partitions_set = Set(partitions_set)
    
    d = Dict() # this dict stores the occurance of each set
    for set in partitions_set
        get!(d, set, 0)
    end

    for partition in partitions_redundant
        set = Set(partition)
        if d[set] == 0
            append!(partitions, [partition])
            d[set] += 1
        else
            continue
        end
    end
    return partitions
end

"""
Finding the Bosonic operator pair within hopping Hamiltonian
"""
function hopping_pair(;m::Int64, Na::Int64, Hbs::Vector{Vector{Int64}})::Matrix{Float64}
    l = length(Hbs)
    H3s = zeros(Float64, l, l)
    for j::Int64 in 1:l
        if Hbs[j][m] < Na && Hbs[j][m+1] > 0
            Hb_copy = copy(Hbs[j])
            Hb_copy[m]::Int64 = Hbs[j][m] + 1
            Hb_copy[m+1]::Int64 = Hbs[j][m+1] - 1
            index::Int64 = findfirst(item -> item == Hb_copy, Hbs)
            H3s[j, index] += sqrt((Hbs[j][m]+1)*Hbs[j][m+1])
        end
    end
    return H3s
end

function average_spacing(differences::Vector{Float64}, i::Int64, v::Int64)
    drop_index = Int((v-1)/2)
    @assert i-drop_index > 0
    @assert i+drop_index <= length(differences)
    sum = 0
    for k in -drop_index:drop_index
        sum += differences[i+k]
    end
    return sum/v
end

"""
This function will give us the evolution state after time t.
Here init denotes the initial state, basis denotes the basis of the Hilbert space.
"""
function state_evolution(;H::Matrix{Float64},t::Float64,basis::Vector{Vector{Int64}},init::Vector{Int64})
    @assert t != 0.0
    ind = findall(x->x==init, basis)[1]
    eigenvals = eigvals(H)
    eigenvecs = eigvecs(H)
    len = length(eigenvals)
    ψ_t = zeros(Float64, length(basis),1)
    for i in 1:len
        vec = eigenvecs[:,i]
        val = eigenvals[i]
        ψ_t += vec[ind]*exp(-im*val*t)*vec
    end
    return vec(ψ_t)
end

"""
This function calculates the reduced_density_matrix by tracing out
subsystem B. Here Na is the number of atoms. Notice that this function
only works in 3-site system.
"""
function reduced_density_matrix(;state::Vector{T}, Na::Int64) where T <: Number
    N = length(state)
    
    rd_rank = 0
    int_rank = 0

    rank1 = Na+1
    rank2 = N/(Na+1)

    if rank1 > rank2
        int_rank = rank2
        rd_rank = rank1
    else
        int_rank = rank1
        rd_rank = rank2
    end

    rd_rank = Int(rd_rank)
    int_rank = Int(int_rank)

    d_mat = state * state' # constructing the density matrix.
    
    if typeof(state[1]) == Float64
        traces = Array{Float64}(undef, 1, rd_rank*rd_rank)
    else
        traces = Array{Complex}(undef, 1, rd_rank*rd_rank)
    end 
    #traces = Array{Float64}(undef, 1, rd_rank*rd_rank)
    #traces = Array{Number}(undef, 1, rd_rank*rd_ran)

    i = 1
    for j in 1:int_rank:N
        for k in 1:int_rank:N
            traces[i] = tr(d_mat[j:j+int_rank-1, k:k+int_rank-1])
            i += 1
        end
    end

    rd_mat = reshape(traces,(rd_rank, rd_rank))
    return transpose(rd_mat)
end


end
