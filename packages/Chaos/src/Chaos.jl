module Chaos

using LinearAlgebra
using Combinatorics

export integer_partition, hopping_pair, hopping_pair_pbc, average_spacing, state_evolution, state_evolution_modified, state_evolution_random, reduced_density_matrix, reduced_density_matrix_modified, evolution_number

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
Finding the Bosonic operator pair within hopping Hamiltonian.
This is prepared for open boundary condition.
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

function hopping_pair_pbc(;m::Int64, Na::Int64, Hbs::Vector{Vector{Int64}})::Matrix{Float64}
    l = length(Hbs)
    ls = length(Hbs[1])
    H3s = zeros(Float64, l, l)
    if m < ls
        for j in 1:l
            if Hbs[j][m] < Na && Hbs[j][m+1] > 0
                Hb_copy = copy(Hbs[j])
                Hb_copy[m]::Int64 = Hbs[j][m] + 1
                Hb_copy[m+1]::Int64 = Hbs[j][m+1] - 1
                index::Int64 = findfirst(item -> item == Hb_copy, Hbs)
                H3s[j, index] += sqrt((Hbs[j][m]+1)*Hbs[j][m+1])
            end
        end

    elseif m == ls
        for j in 1:l
            if Hbs[j][m] < Na && Hbs[j][1] > 0
                Hb_copy = copy(Hbs[j])
                Hb_copy[m]::Int64 = Hbs[j][m] + 1
                Hb_copy[1]::Int64 = Hbs[j][1] - 1
                index::Int64 = findfirst(item -> item == Hb_copy, Hbs)
                H3s[j, index] += sqrt((Hbs[j][m]+1)*Hbs[j][1])
            end
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
    ??_t = zeros(Float64, length(basis),1)
    for i in 1:len
        vec = eigenvecs[:,i]
        val = eigenvals[i]
        ??_t += conj(vec[ind])*exp(-im*val*t)*vec
    end
    return vec(??_t)
end

function state_evolution_modified(;eigenvals::Vector{T},eigenvecs::Matrix{Float64},t::Float64,basis::Vector{Vector{Int64}},init::Vector{Int64}) where T <: Number
    @assert t != 0.0
    ind = findall(x->x==init, basis)[1]
    len = length(eigenvals)
    ??_t = zeros(Float64, length(basis), 1)
    for i in 1:len
        vec = eigenvecs[:,i]
        val = eigenvals[i]
        ??_t += conj(vec[ind])*exp(-im*val*t)*vec
    end
    return vec(??_t)
end

function state_evolution_random(;H::Matrix{Float64},t::Float64,basis::Vector{Vector{Int64}},init::Vector{T}) where T <: Number
    eigenvals = eigvals(H)
    eigenvecs = eigvecs(H)
    len = length(eigenvals)
    ??_t = zeros(Float64, length(basis), 1)
    for i in 1:len
        vec = eigenvecs[:,i]
        val = eigenvals[i]
        ??_i = vec'*init
        ??_t += ??_i*exp(-im*val*t)*vec
    end
    return vec(??_t)
end


# This function will integrate out the first "Rd" (stands for
# reduced dimension) sites of the system.
function reduced_density_matrix(;state::Vector{T}, Na::Int64, Ls::Int64, Rd::Int64) where T <: Number
    @assert Na != 0 && Rd != 0
    @assert Rd < Ls
    sys_bs = integer_partition(N=Na, L=Ls) # the set of all basis of the system 
    ind_lsts = Vector{Int64}[] # this set will be used to compute all non-zero elements of reduced density matrix
    
    # we first consider the case of zero particle
    init = vec(zeros(Int64, 1, Rd))
    ind_lst = findall(x->x[1:Rd]==init, sys_bs)
    push!(ind_lsts, ind_lst)
    
    # we then come to consider all the remaining cases
    for i in 1:Na
        hbs = integer_partition(N=i, L=Rd)
        for lst in hbs
            ind_lst = findall(x->x[1:Rd]==lst, sys_bs)
            push!(ind_lsts, ind_lst)
        end
    end
    
    # below we compute the reduced density matrix
    dmat = state*state'
    len = length(state)
    if typeof(state[1]) == Float64   
        rdmat = zeros(Float64, len, len)
    else
        rdmat = zeros(Complex, len, len)
    end
    
    for lst in ind_lsts
        for i in lst
            for j in lst
                rdmat[i,j] = dmat[i,j]
            end
        end
    end
    return rdmat
end

# Notice that this function only works in our paper i.e., triple-well system integrating out
# the first site.
function reduced_density_matrix_modified(;Na::Int64,state::Vector{T},Hbs::Vector{Vector{Int64}}) where T <: Number
    @assert Na > 0

    dmat = state*state' # the original density matrix
    len = length(state)
    rdmat = zeros(Complex, len, len) # initialize the reduced density matrix
    
    row_len = 0
    for num in 0:Na
        # these two list of indices is to locate
        # the corresponding value in dmat
        row_indices = findall(x->x[1]==num, Hbs)
        col_indices = findall(x->x[3]==num, Hbs)
        ind = row_indices[1]
        row_len = length(row_indices)
        dmat_values = []
        for i in row_indices
            for j in col_indices
                push!(dmat_values, dmat[i,j])
            end
        end
        dmat_values = reshape(dmat_values, row_len, row_len)
        rdmat[ind:ind+row_len-1,ind:ind+row_len-1] = dmat_values
    end
    return rdmat
end


"""
This function will give the evolution of number of atoms at i-th site
"""
function evolution_number(;H::Matrix{Float64}, init::Vector{Int64}, basis::Vector{Vector{Int64}}, site::Int64, ts::Vector{Float64})
    num_site = length(basis[1])
    @assert site >= 1 && site <= num_site
    
    vals = eigvals(H)
    vecs = eigvecs(H)
    
    m = findall(x->x==init, basis)[1]
    len_vals = length(vals)
    
    ns = Float64[]
    for t in ts 
        n = 0 # This will be the expectation value of number of the site
        for m??? in 1:len_vals
            n_i = basis[m???][site]
            prod_1 = 0
            prod_2 = 0
            for j in 1:len_vals
                prod_1 += vecs[m,j]*exp(1im*vals[j]*t)*conj(vecs[m???,j])
                prod_2 += vecs[m,j]*exp(-1im*vals[j]*t)*conj(vecs[m???,j])
            end
            n += prod_1*prod_2*n_i
        end
        append!(ns, n)
    end
    return ns
end
        




















end
