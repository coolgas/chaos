module Integer

export integer_partition, integer_partition_redundant, hopping_pair, average_spacing

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

end
