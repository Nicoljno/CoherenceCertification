using EntanglementDetection
using Ket
using JuMP
using MosekTools
using LinearAlgebra
using MathOptInterface

"""
Definition of the quantum channels
- dephasing:             Λ(ρ)=v*ρ+(1-v)diag(ρ)
- depolarising:          Λ(ρ)=v*ρ+(1-v)I/d
- amplitude damping:     Λ(ρ)=K0ρK0†+K1ρK1†
"""
dephasing(ρ, v) = v*ρ + (1-v)*Diagonal(ρ)
depolarising(ρ, v) = v*ρ + (1-v)*I(size(ρ, 1))/size(ρ, 1)
function amplitude_damping(ρ, γ)
    E0 = [1 0; 0 sqrt(1 - γ)]
    E1 = [0 sqrt(γ); 0 0]
    return E0 * ρ * E0' + E1 * ρ * E1'
end

"""
Function emplojed to upper bound the joint measurability of binarization
- targets: the target states of the ensemble
- states: polytope approximation of the hypersphere
- μ: shrinking factor to define the outer approximation of the polytope
"""
function approx_jm_visibility_outer(targets, states, μ)
    d = size(targets[1], 1)
    m = length(targets)
    n = length(states)

    I_d = Matrix{ComplexF64}(I, d, d)

    rho = [Hermitian(states[l]*states[l]'/tr(states[l]*states[l]')) for l in 1:n]
    chi = [(1/μ) * rho[l] + (1 - 1/μ) * I_d / d for l in 1:n]

    model = Model(Mosek.Optimizer)
    set_silent(model)

    @variable(model, 1>= q[1:n] >= 0)
    @variable(model, 1>= r[1:n, 1:m] >= 0)
    @variable(model, 0 <= η <= 1)

    # sum_l q_l chi_l = I
    @constraint(model,
        sum(q[l] * chi[l] for l in 1:n) == I_d
    )
    
    # sum_l r_lx chi_l = η targets[x] + (1-η) I/d
    for x in 1:m
        @constraint(model,
            sum(r[l, x] * chi[l] for l in 1:n) ==
            depolarising(dephasing(targets[x], η), 0.99) #dephasing(t, 0.6)
            #η * targets[x] + (1 - η) * I_d / d
        )
    end

    for l in 1:n, x in 1:m
        @constraint(model, r[l, x] <= q[l])
    end

    @objective(model, Max, η)
    optimize!(model)

    return value(η), termination_status(model)
end

"""
Function emplojed to lower bound the joint measurability of binarization
- targets: the target states of the ensemble
- states: polytope approximation of the hypersphere
- vis_dep: depolarising channel visibility 
"""
function approx_jm_visibility(targets, states, vis_dep)
    d = size(targets[1], 1)
    m = length(targets)
    n = length(states)

    I_d = Matrix{ComplexF64}(I, d, d)
    rho = [Hermitian(states[l]*states[l]'/tr(states[l]*states[l]')) for l in 1:n]

    model = Model(Mosek.Optimizer)
    #set_silent(model)

    @variable(model, q[1:n] >= 0)
    @variable(model, r[1:n, 1:m] >= 0)
    @variable(model, 0 <= η <= 1)

    # sum_l q_l rho_l = I
    @constraint(model,
        sum(q[l] * rho[l] for l in 1:n) == I_d
    )

    # sum_l r_lx rho_l = Λ_η(targets[x]) 
    for x in 1:m
        @constraint(model,
            sum(r[l, x] * rho[l] for l in 1:n) ==
            depolarising(dephasing(targets[x], vis_dep), η)
            #η * targets[x] + (1 - η) * I_d / d
        )
    end

    # r_lx <= q_l ensures the second outcome is also valid
    for l in 1:n, x in 1:m
        @constraint(model, r[l, x] <= q[l])
    end

    @objective(model, Max, η)
    optimize!(model)

    return value(η), termination_status(model)
end

d=3
A = EntanglementDetection.CrossPolytopeSubdivision{Float64}(d, 3)
A_states = EntanglementDetection.CrossPolytopeSubdivision{Float64}(d, 3)

ket = Vector{ComplexF64}(undef, A.d)
states = [copy(EntanglementDetection._populate!(ket, A_states, i)) for i in 1:length(A_states)]

"""
targets=[]
X = shift(d, 1)
vecs = eigvecs(X)
for i=1:d
    push!(targets, ketbra(vecs[:,i]))
end
Z = clock(d, 1)
vecs = eigvecs(Z)
count = 1
for i=d+1:2*d
    push!(targets, ketbra(vecs[:,count]))
    global count = count+1
end

targets=[]
for i=1:200
    tmp = random_state_ket(d)
    push!(targets, tmp*tmp')
end
#"""
targets=[]
for i=1:length(states)
    push!(targets, Hermitian(states[i]*states[i]'/tr(states[i]*states[i]')))
end

res=[]
for n_ = 1:1
    local ηstar_in, status_in = approx_jm_visibility(targets, states, 1.0)
    println("status_in = ", status_in)
    push!(res, ηstar_in)
end
#ηstar_out, status_out = approx_jm_visibility_outer(targets, states, A.η)

#println("η*in = ", ηstar_in)
#println("η*out = ", ηstar_in/A.η)
#println("η*out = ", ηstar_out)
#println("status_out = ", status_out)

"""
idx = 0:20
col2 = 0.8 .+ 0.01 .* idx
col3 = Float64.(res)
col4 = col3 ./ A.η

M = hcat(idx, col2, col3, col4)
#"""
#writedlm("results_depol_damp.txt", M, '\t')

"""
the outer approximation can be obtained directly 
from inner result and shrinking factor
"""
(res, res/A_states.η)