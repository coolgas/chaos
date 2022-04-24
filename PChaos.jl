### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 407a8dea-26cf-4277-8157-827fae5b9b13
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()

    using Chaos
end

# ╔═╡ 2893af54-1f54-4fcd-ac54-2a918daea748
begin
	using Combinatorics
	using LinearAlgebra
	using LaTeXStrings
	using GLMakie
	using Statistics
end

# ╔═╡ 309117bd-3cdc-4a47-9831-aa7b6b5cef5f
Na = 60

# ╔═╡ f5d94d6f-f207-4247-81f6-e06fefe0228d
Ls = 3

# ╔═╡ 38195f3b-8078-44a2-b45f-3bb6f5339734
Hbs = let
	integer_partition(N=Na, L=Ls)
end


# ╔═╡ 59690aae-3b29-46b8-9f2a-2b482ffd27bb
Lb = let
	length(Hbs) # the dimension of the Hilbert space
end

# ╔═╡ 015762a8-bdb5-4070-aa59-3771754e6ac9
Lb_shouldbe = let
	binomial(Na+Ls-1, Na)
end

# ╔═╡ 540d4518-dd21-4021-8bd2-d4d3edba5330
# The hopping Hamiltonian
Hhop = let
	HLs = zeros(Float64, Lb, Lb);
	for i in 1:Ls-1
	    HLs += hopping_pair(m=i, Na=Na, Hbs=Hbs)
	end
	HLs + HLs'
end;

# ╔═╡ ff624ebf-b486-440f-a2c5-8cbb68c4440a
# The onsite interaction
begin
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
end

# ╔═╡ a21b8fb0-fc79-4cce-a4b7-6f7ca8eb7f2f
begin
	
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

end

# ╔═╡ 85436297-a22d-467f-9962-d870593f302c
begin
	# In our case, we take average every 7 levels, we drop (7-1)/2 spacings at either end.
	v = 5;
	
	# This is defined for later convenience.
	f(x) = exp(-x);
	ff(x) = (π*x/2)*exp(-π*x^2/4);
end

# ╔═╡ 55b4958b-fadb-4acc-b253-24509c35918d
# All the parameters we gonna use
begin 
	J = 1.0
	g = -0.05
	d = 1.5
	C6 = 100.0
	R = 3.0
	gammas = 0:0.2:10
	eigs = zeros(Float64, Lb, length(gammas));
end;

# ╔═╡ 46047bba-96df-48c3-b929-116b53f323a8
# Here we try to plot γ = 0 case of level statistics
# to see if Makie actually work.
# I will not put too much energy on this.
begin
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

end

# ╔═╡ 4ec20759-1705-47ec-bd82-aaaf0b079047
begin
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="s", ylabel="P")
	x = range(0, 7, length=50)
	hist!(ax, new_diffs1, bins=:15, normalization=:pdf, xlabel="s")
	lines!(ax, x, ff, color=:black, linestyle=:dash, linewidth=3.5)
	lines!(ax, x, f, color=:red, linestyle=:dash, linewidth=3.5)
	fig
end

# ╔═╡ fc1ef60f-6f95-4706-8d14-4df7420477f1
md"""
## Everything set up, below we come to the entanglement entropy.
"""

# ╔═╡ cf4c8537-f53d-4129-9f32-8dc192bf47f6
begin
	# This basically implements the general formula of entanglement entropy.
	function EE(; γ::Float64, J::Float64, g::Float64, d::Float64, C6::Float64, R::Float64, Na::Int64)
	    vec = eigvecs(Htot(γ=γ, J=J, g=g, d=d, C6=C6, R=R))[:,1]
	    rdmat = reduced_density_matrix(state=vec, Na=Na, Ls=Ls, Rd=1)
	    vals = svd(rdmat).S
	    EE = -transpose(vals)*log.(vals)
	    return EE
	end

	# I now only makes this function works given the initial state with equal
	# number of atoms in each site. But I potentially can extend it to arbitrary initial state.
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

	# This function works for 
	function evolEE_random(;H::Matrix{Float64},t::Float64,
	basis::Vector{Vector{Int64}},init::Vector{Float64})
	    
	    if t == 0.0
	        rdmat = reduced_density_matrix(state=init, Na=Na, Ls=Ls, Rd=1)
	        vals = svd(rdmat).S
	        EE = -transpose(vals)*log.(vals)
	    else
	        vec = state_evolution_random(H=H,t=t,basis=basis,init=init)
	        #rdmat = reduced_density_matrix(Na=Na, state=vec, Hbs=Hbs)
	        rdmat = reduced_density_matrix(state=vec, Na=Na, Ls=Ls, Rd=1)
	        vals = svd(rdmat).S
	        #EE = -transpose(conj(vals).*vals)*log.(conj(vals).*(vals))
	        #vals_squd = vals.*conj(vals)
	        EE = -transpose(vals)*log.(vals)
	    end
	    
	    return EE
	end
	
end
	

# ╔═╡ 21547d99-d365-41be-b047-3981f668a3e7
# We here prepare a bunch of Hamiltonians for later use.
begin
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
end;

# ╔═╡ e58bd42e-9000-4457-822e-20583c738a48
begin
	ts = 0.0:0.2:10;
	gammas_new = 0.0:1:8
	ts1 = 0.0:0.2:5
end

# ╔═╡ 36c5ed45-e785-435c-a6d9-78411552cd34
begin
	aveEEs = zeros(Float64, 1, length(gammas_new))
	std_ee = zeros(Float64, 1, length(gammas_new))
	aveEEs_new = zeros(Float64, 1, length(gammas_new))
	std_ee_new = zeros(Float64, 1, length(gammas_new))
	aveEEs_newnew = zeros(Float64, 1, length(gammas_new))
	std_ee_newnew = zeros(Float64, 1, length(gammas_new))
end

# ╔═╡ e0042278-087f-44dd-bb99-d9867c3d4d43
md"""
## Below we consider the case [20, 20, 20]
"""

# ╔═╡ e58536c2-78a4-4d4f-b76f-96085ca73641
@time evolEE0s = [evolEE(H=H0,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts];

# ╔═╡ 4c600027-7508-42fe-86e9-56a520a6bb3f
@time evolEE4s = [evolEE(H=H4,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts]

# ╔═╡ 95e2c945-8cdd-496d-b1d8-38268642a178
@time evolEE6s = [evolEE(H=H6,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts]

# ╔═╡ d0176929-6265-4179-aa3b-5e51624dcd91
@time evolEE8s = [evolEE(H=H8,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts]

# ╔═╡ 8c7366ec-6b50-418b-8043-0400ca8cbdc1
@time evolEE10s = [evolEE(H=H10,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts]

# ╔═╡ ee9ceb59-83b1-42cf-ace8-1286779144ce
md"""
## Below we consider the case [0, 0, 60]
"""

# ╔═╡ 1ca97531-5872-4967-9a10-fbc51453ed91
@time evolEE0s_new = [evolEE(H=H0,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts]

# ╔═╡ fe95646f-20a9-431b-9286-c125fc1ae998
@time evolEE4s_new = [evolEE(H=H4,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts]

# ╔═╡ 406ac126-9dfe-4c6a-8818-b0d471b2020c
@time evolEE6s_new = [evolEE(H=H6,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts]

# ╔═╡ df2d5c18-111f-4297-9ed1-2e897a8ce2c4
@time evolEE8s_new = [evolEE(H=H8,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts]

# ╔═╡ c194c089-87dc-410a-b209-2d0bb48ebfe4
@time evolEE10s_new = [evolEE(H=H10,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts]

# ╔═╡ 182552cf-c20d-4787-bdc4-fe5d0a40ce33
md"""
## Below we consider the case [10,40,10]
"""

# ╔═╡ cd7a1b1e-9c8e-47ee-aaef-24b2bb541227
@time evolEE0s_newnew = [evolEE(H=H0,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts]

# ╔═╡ 5d165176-e229-4337-af0b-748d5116138d
@time evolEE4s_newnew = [evolEE(H=H4,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts]

# ╔═╡ 546dca47-a0da-4521-84c4-1bbb90c179b9
@time evolEE6s_newnew = [evolEE(H=H6,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts]

# ╔═╡ b23b68eb-3b13-4db9-96dd-c5d390481311
@time evolEE8s_newnew = [evolEE(H=H8,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts]

# ╔═╡ c5d6d03e-5956-4e71-81f1-c5f6079c4671
@time evolEE10s_newnew = [evolEE(H=H10,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts]

# ╔═╡ adf1641d-25f6-4b7a-b5f6-6eeca8a0ca19
md"""
## [20 20 20] case of averages 
"""

# ╔═╡ 1d798ec4-8385-453f-a57f-db36fc421f1b
@time for i in 1:length(gammas_new)
	gamma = gammas_new[i]
	H = Htot(γ=gamma, J=J, g=g, d=d, C6=C6, R=R)
	evolEEs = [evolEE(H=H,t=t,basis=Hbs,init=vec([20 20 20])) for t in ts1]
	len = length(evolEEs)
	drop_rate = 0.2 # we drop the first 20% of data
	left_num = Int(floor(drop_rate*len))
	aveEEs[i] = mean(evolEEs[left_num:len])
	std_ee[i] = std(evolEEs[left_num:len])
	println(aveEEs)
end

# ╔═╡ a68426c6-13f3-4af7-abea-4978208c7596
md"""
## [0 0 60] case of averages
"""

# ╔═╡ db9a0449-b830-422c-a18a-9629b8b91cfc
@time for i in 1:length(gammas_new)
	gamma = gammas_new[i]
	H = Htot(γ=gamma, J=J, g=g, d=d, C6=C6, R=R)
	evolEEs = [evolEE(H=H,t=t,basis=Hbs,init=vec([0 0 60])) for t in ts1]
	len = length(evolEEs)
	drop_rate = 0.2 # we drop the first 20% of data
	left_num = Int(floor(drop_rate*len))
	aveEEs_new[i] = mean(evolEEs[left_num:len])
	std_ee_new[i] = std(evolEEs[left_num:len])
	display(aveEEs_new)
end


# ╔═╡ 794faca1-a799-4f83-b8e1-4c21ce1eec95
md"""
## [10 40 10] case of averages
"""

# ╔═╡ d4ae0b86-ce78-4c80-8002-a5a386fa0dcf
@time for i in 1:length(gammas_new)
	gamma = gammas_new[i]
	H = Htot(γ=gamma, J=J, g=g, d=d, C6=C6, R=R)
	evolEEs = [evolEE(H=H,t=t,basis=Hbs,init=vec([10 40 10])) for t in ts1]
	len = length(evolEEs)
	drop_rate = 0.2 # we drop the first 20% of data
	left_num = Int(floor(drop_rate*len))
	aveEEs_newnew[i] = mean(evolEEs[left_num:len])
	std_ee_newnew[i] = std(evolEEs[left_num:len])
	display(aveEEs_newnew)
end


# ╔═╡ 489534d2-e14c-4d29-9646-7217992b902b


# ╔═╡ Cell order:
# ╠═407a8dea-26cf-4277-8157-827fae5b9b13
# ╠═2893af54-1f54-4fcd-ac54-2a918daea748
# ╠═309117bd-3cdc-4a47-9831-aa7b6b5cef5f
# ╠═f5d94d6f-f207-4247-81f6-e06fefe0228d
# ╠═38195f3b-8078-44a2-b45f-3bb6f5339734
# ╠═59690aae-3b29-46b8-9f2a-2b482ffd27bb
# ╠═015762a8-bdb5-4070-aa59-3771754e6ac9
# ╠═540d4518-dd21-4021-8bd2-d4d3edba5330
# ╠═ff624ebf-b486-440f-a2c5-8cbb68c4440a
# ╠═a21b8fb0-fc79-4cce-a4b7-6f7ca8eb7f2f
# ╠═85436297-a22d-467f-9962-d870593f302c
# ╠═55b4958b-fadb-4acc-b253-24509c35918d
# ╠═46047bba-96df-48c3-b929-116b53f323a8
# ╠═4ec20759-1705-47ec-bd82-aaaf0b079047
# ╠═fc1ef60f-6f95-4706-8d14-4df7420477f1
# ╠═cf4c8537-f53d-4129-9f32-8dc192bf47f6
# ╠═21547d99-d365-41be-b047-3981f668a3e7
# ╠═36c5ed45-e785-435c-a6d9-78411552cd34
# ╠═e58bd42e-9000-4457-822e-20583c738a48
# ╠═e0042278-087f-44dd-bb99-d9867c3d4d43
# ╠═e58536c2-78a4-4d4f-b76f-96085ca73641
# ╠═4c600027-7508-42fe-86e9-56a520a6bb3f
# ╠═95e2c945-8cdd-496d-b1d8-38268642a178
# ╠═d0176929-6265-4179-aa3b-5e51624dcd91
# ╠═8c7366ec-6b50-418b-8043-0400ca8cbdc1
# ╠═ee9ceb59-83b1-42cf-ace8-1286779144ce
# ╠═1ca97531-5872-4967-9a10-fbc51453ed91
# ╠═fe95646f-20a9-431b-9286-c125fc1ae998
# ╠═406ac126-9dfe-4c6a-8818-b0d471b2020c
# ╠═df2d5c18-111f-4297-9ed1-2e897a8ce2c4
# ╠═c194c089-87dc-410a-b209-2d0bb48ebfe4
# ╠═182552cf-c20d-4787-bdc4-fe5d0a40ce33
# ╠═cd7a1b1e-9c8e-47ee-aaef-24b2bb541227
# ╠═5d165176-e229-4337-af0b-748d5116138d
# ╠═546dca47-a0da-4521-84c4-1bbb90c179b9
# ╠═b23b68eb-3b13-4db9-96dd-c5d390481311
# ╠═c5d6d03e-5956-4e71-81f1-c5f6079c4671
# ╠═adf1641d-25f6-4b7a-b5f6-6eeca8a0ca19
# ╠═1d798ec4-8385-453f-a57f-db36fc421f1b
# ╠═a68426c6-13f3-4af7-abea-4978208c7596
# ╠═db9a0449-b830-422c-a18a-9629b8b91cfc
# ╠═794faca1-a799-4f83-b8e1-4c21ce1eec95
# ╠═d4ae0b86-ce78-4c80-8002-a5a386fa0dcf
# ╠═489534d2-e14c-4d29-9646-7217992b902b
