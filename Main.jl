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

fig = Figure();ax = Axis(fig[1,1], xlabel="s", ylabel="P")
x = range(0, 7, length=50)
hist!(ax, new_diffs1, bins=:15, normalization=:pdf, xlabel="s")
lines!(ax, x, ff, color=:black, linestyle=:dash, linewidth=3.5)
lines!(ax, x, f, color=:red, linestyle=:dash, linewidth=3.5)
fig


