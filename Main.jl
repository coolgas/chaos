using Combinatorics
using LinearAlgebra
using LaTeXStrings
using GLMakie
using Chaos
using Statistics


Na = 60;
Ls = 3;

Hbs = integer_partition(N=Na, L=Ls);
Lb = length(Hbs)

# The hopping Hamiltonian
Hhop = let
	HLs = zeros(Float64, Lb, Lb);
	for i in 1:Ls-1
	    HLs += hopping_pair(m=i, Na=Na, Hbs=Hbs)
	end
	HLs + HLs'
end;

# The onsite interaction
Honsite = zeros(Float64, Lb, Lb);
for j in 1:Lb
    for k in 1:Lb
        if j == k
            sum = 0
            for l in 1:Ls
                sum += Hbs[j][l]^2 - Hbs[j][l]
            end
            Honsite[j,k] += sum
        end
    end
end

# The tilting Hamiltonian
function Htit(γ::Float64)
    Htit = zeros(Float64, Lb, Lb)
    for j in 1:Lb
        for k in 1:Lb
            if j == k
                for l in 1:Ls
                    Htit[j,k] += -(l-(Ls+1)/2)*γ*Hbs[j][l]
                end
            end
        end
    end
    return Htit
end

# Soft-core Interactions
Λ(delta::Int64, d::Float64, C6::Float64, R::Float64) = C6 / ((delta)^6*d^6 + R^6)

function Hsc(d::Float64, C6::Float64, R::Float64)
    Hsc = zeros(Float64, Lb, Lb)
    for j in 1:Lb
        for k in 1:Lb
            if j == k
                for l in 1:Ls
                    for m in 1:Ls
                        if l-m == 0
                            Hsc[j,k] += Λ(0,d,C6,R)*Hbs[j][l]*Hbs[j][m]
                        elseif abs(m-l) == 1
                            Hsc[j,k] += Λ(1,d,C6,R)*Hbs[j][l]*Hbs[j][m]
                        elseif abs(m-l) == 2
                            Hsc[j,k] += Λ(2,d,C6,R)*Hbs[j][l]*Hbs[j][m]
                        else
                            continue
                        end
                    end
                end
            end
        end
    end
    return Hsc
end

# The total Hamiltonian
function Htot(;γ::Float64, J::Float64, g::Float64, d::Float64, C6::Float64, R::Float64)
    Htot = Htit(γ) - J*Hhop + (g/2)*Honsite + (1/2)*Hsc(d,C6,R)
    return Htot
end

# In our case, we take average every 7 levels, we drop (7-1)/2 spacings at either end.
v = 5;

# This is defined for later convenience.
f(x) = exp(-x);
ff(x) = (π*x/2)*exp(-π*x^2/4);

# This is all the constants going to be used.
J = 1.0
g = -0.05
d = 1.5
C6 = 100.0
R = 3.0
gammas = 0:0.2:10
eigs = zeros(Float64, Lb, length(gammas));

# Here we try to plot γ = 0 case of level statistics
# to see if Makie actually work.
# I will not put too much energy on this.
eigs1 = eigvals(Htot(γ=0.0, J=J, g=g, d=d, C6=C6, R=R));
diffs1 = [eigs1[i+1]-eigs1[i] for i in 1:length(eigs1)-1];
index1 = Int((v-1)/2)
index2 = length(diffs1)-index1
diffs1_dropped = diffs1[index1:index2];
new_diffs1 = []
for i in 1:length(diffs1_dropped)
    drop_index = Int((v-1)/2)
    if i-drop_index>0 && i+drop_index<=length(diffs1_dropped)
        append!(new_diffs1, diffs1_dropped[i]/average_spacing(diffs1_dropped, i, v))
    end
end

fig = Figure()
ax = Axis(fig[1,1], xlabel="s", ylabel="P")
x = range(0, 7, length=50)
hist!(ax, new_diffs1, bins=:15, normalization=:pdf, xlabel="s")
lines!(ax, x, ff, color=:black, linestyle=:dash, linewidth=3.5)
lines!(ax, x, f, color=:red, linestyle=:dash, linewidth=3.5)
fig

# Everything set up, we are going to entanglement entropy
# This basically implements the general formula of entanglement entropy.
function EE(; γ::Float64, J::Float64, g::Float64, d::Float64, C6::Float64, R::Float64, Na::Int64)
    vec = eigvecs(Htot(γ=γ, J=J, g=g, d=d, C6=C6, R=R))[:,1]
    rdmat = reduced_density_matrix(state=vec, Na=Na, Ls=Ls, Rd=1)
    vals = svd(rdmat).S
    EE = -transpose(vals)*log.(vals)
    return EE
end

# This function find the evoling state
function evolEE(;H::Matrix{Float64},t::Float64,
basis::Vector{Vector{Int64}},init::Vector{Int})
    if t == 0.0
        EE = 0.0
    else
        vec = state_evolution(H=H,t=t,basis=basis,init=init)
        #rdmat = reduced_density_matrix(Na=Na, state=vec, Hbs=Hbs)
        rdmat = reduced_density_matrix(state=vec, Na=Na, Ls=Ls, Rd=1)
        vals = svd(rdmat).S
        #EE = -transpose(conj(vals).*vals)*log.(conj(vals).*(vals))
        #vals_squd = vals.*conj(vals)
        EE = -transpose(vals)*log.(vals)
    end
    return EE
end

# Calculating a bunch of Hamiltonians.
H0 = Htot(γ=0.0, J=J, g=g, d=d, C6=C6, R=R);
H1 = Htot(γ=1.0, J=J, g=g, d=d, C6=C6, R=R);
H2 = Htot(γ=2.0, J=J, g=g, d=d, C6=C6, R=R);
H3 = Htot(γ=3.0, J=J, g=g, d=d, C6=C6, R=R);
H4 = Htot(γ=4.0, J=J, g=g, d=d, C6=C6, R=R);
H5 = Htot(γ=5.0, J=J, g=g, d=d, C6=C6, R=R);
H6 = Htot(γ=6.0, J=J, g=g, d=d, C6=C6, R=R);
H7 = Htot(γ=7.0, J=J, g=g, d=d, C6=C6, R=R);
H8 = Htot(γ=8.0, J=J, g=g, d=d, C6=C6, R=R);
H9 = Htot(γ=9.0, J=J, g=g, d=d, C6=C6, R=R);
H10 = Htot(γ=10.0, J=J, g=g, d=d, C6=C6, R=R);

# Initializing a bunch of things.
ts = 0.0:0.05:6;
gammas_new = 0.0:1:8
ts1 = 0.0:0.2:5

# Initializing a bunch of things
aveEEs = zeros(Float64, 1, length(gammas_new))
std_ee = zeros(Float64, 1, length(gammas_new))
aveEEs_new = zeros(Float64, 1, length(gammas_new))
std_ee_new = zeros(Float64, 1, length(gammas_new))
aveEEs_newnew = zeros(Float64, 1, length(gammas_new))
std_ee_newnew = zeros(Float64, 1, length(gammas_new))

# Below we consider the case [20, 20, 20]
@time evolEE0s = [evolEE(H=H0,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts];
@time evolEE4s = [evolEE(H=H4,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts];
@time evolEE6s = [evolEE(H=H6,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts];
@time evolEE8s = [evolEE(H=H8,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts];
@time evolEE10s = [evolEE(H=H10,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts];

# Below we consider the case [0, 0, 60]
@time evolEE0s_new = [evolEE(H=H0,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts];
@time evolEE4s_new = [evolEE(H=H4,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts];
@time evolEE6s_new = [evolEE(H=H6,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts];
@time evolEE8s_new = [evolEE(H=H8,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts];
@time evolEE10s_new = [evolEE(H=H10,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts];

# Below we conside the case [10,40,10]
@time evolEE0s_newnew = [evolEE(H=H0,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts];
@time evolEE4s_newnew = [evolEE(H=H4,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts];
@time evolEE6s_newnew = [evolEE(H=H6,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts];
@time evolEE8s_newnew = [evolEE(H=H8,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts];
@time evolEE10s_newnew = [evolEE(H=H10,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts];


# [20 20 20] case of averages
@time for i in 1:length(gammas_new)
	gamma = gammas_new[i]
	H = Htot(γ=gamma, J=J, g=g, d=d, C6=C6, R=R)
	evolEEs = [evolEE(H=H,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts1]
	len = length(evolEEs)
	drop_rate = 0.2 # we drop the first 20% of data
	left_num = Int(floor(drop_rate*len))
	aveEEs[i] = mean(evolEEs[left_num:len])
    std_ee[i] = std(log.(evolEEs[left_num:len]))
	println(aveEEs)
end

# [0 0 60] case of averages
@time for i in 1:length(gammas_new)
	gamma = gammas_new[i]
	H = Htot(γ=gamma, J=J, g=g, d=d, C6=C6, R=R)
	evolEEs = [evolEE(H=H,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts1]
	len = length(evolEEs)
	drop_rate = 0.2 # we drop the first 20% of data
	left_num = Int(floor(drop_rate*len))
	aveEEs_new[i] = mean(evolEEs[left_num:len])
    std_ee_new[i] = std(log.(evolEEs[left_num:len]))
	println(aveEEs_new)
end

# [10 40 10] case of averages
@time for i in 1:length(gammas_new)
	gamma = gammas_new[i]
	H = Htot(γ=gamma, J=J, g=g, d=d, C6=C6, R=R)
	evolEEs = [evolEE(H=H,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts1]
	len = length(evolEEs)
	drop_rate = 0.2 # we drop the first 20% of data
	left_num = Int(floor(drop_rate*len))
	aveEEs_newnew[i] = mean(evolEEs[left_num:len])
    std_ee_newnew[i] = std(log.(evolEEs[left_num:len]))
	println(aveEEs_newnew)
end

gs = [0 4 6 8 10]

ee_dict1 = Dict(0=>evolEE0s, 4=>evolEE4s, 6=>evolEE6s, 8=>evolEE8s, 10=>evolEE10s)

ee_dict2 = Dict(0=>evolEE0s_new, 4=>evolEE4s_new, 6=>evolEE6s_new, 8=>evolEE8s_new, 10=>evolEE10s_new)

ee_dict3 = Dict(0=>evolEE0s_newnew, 4=>evolEE4s_newnew, 6=>evolEE6s_newnew, 8=>evolEE8s_newnew, 10=>evolEE10s_newnew)

for g in gs
    io = open("/Users/tianyiyan/Desktop/Chaos/data/ee_newnew/text$(g).txt", "w") do io
        for x in ee_dict3[g]
            println(io, x)
        end
    end
end


io = open("/Users/tianyiyan/Desktop/Chaos/data/aves/avee_newnew/aveEEs_newnew.txt", "w") do io
    for x in aveEEs_newnew
        println(io, x)
    end
end

std_ee_newnew

io = open("/Users/tianyiyan/Desktop/Chaos/data/aves/stdee_newnew/std_ee_newnew.txt", "w") do io
    for x in std_ee_newnew
        println(io, x)
    end
end

evolEE10s

io = open("/Users/tianyiyan/Desktop/Chaos/data/dense/ees/ee/text10.txt", "w") do io
    for x in evolEE10s
        println(io, x)
    end
end

