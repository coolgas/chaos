module Chaos

using LinearAlgebra
using Combinatorics

export integer_partition, hopping_pair, average_spacing, state_evolution, reduced_density_matrix, reduced_density_matrix_3sites

"""
This function will give partitions of m into at most n integers.
"""
function integer_partition(;N::Int64, L::Int64)
    @assert typeof(N) == Int64 && typeof(L) == Int64
    @assert N != 0 && L != 0
    all_ps = Vector{Int64}[]
    equ_len_ps = Vector{Int64}[] # this list hold all the equal length partitions 
    for i in 1:L
        ps = partitions(N, i)
        for p in ps
            push!(all_ps, p)
        end         
    end
    
    for p in all_ps
        if length(p) == L
            push!(equ_len_ps, p)
        else
            p_mod = append!(p, zeros(Int64, 1, L-length(p)))
            push!(equ_len_ps, p_mod)
        end
    end
    
    Hbs = Vector{Int64}[]
    for p in equ_len_ps
        append!(Hbs, collect(multiset_permutations(p, L)))
    end
    
    Hbs_arranged = [Vector{Int64}[] for i in 1:N+1]
    indices = [i for i in 0:N]
    for vec in Hbs
        ind = findall(x->x==vec[1], indices)[1]
        push!(Hbs_arranged[ind], vec)
    end

    return vcat(Hbs_arranged...)
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
        ψ_t += conj(vec[ind])*exp(-im*val*t)*vec
    end
    return vec(ψ_t)
end

"""
This function will calculate the reduced_density_matrix by integrate 'int_rank' out
of the density matrix contructed from the state. Noticed that this function only
considers the bipartite system AB and only taking trace with respect to B.
"""
function reduced_density_matrix(;state::Array{T,1}, int_rank::Int) where T <: Number
    N = length(state)
    @assert N % 2 == 0
    @assert int_rank % 2 == 0
    @assert N % int_rank == 0

    rd_rank = Int(N / int_rank) # the output rank.
    d_mat = state * state' # constructing the density matrix.
    
    traces = Array{Float64}(undef, 1, rd_rank*rd_rank)
    i = 1
    for j in 1:int_rank:N
        for k in 1:int_rank:N
            traces[i] = tr(d_mat[j:j+int_rank-1, k:k+int_rank-1])
            i += 1
        end
    end

    rd_mat = reshape(traces,(rd_rank, rd_rank))
    return rd_mat'
end


end
