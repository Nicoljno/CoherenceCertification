using LinearAlgebra
using Ket
using Nemo
using MAT
using NPZ


"""
Return the finite field F_q and a vector of its q elements, requires Nemo.  
"""

function paley_field(q::Int)
    K, _ = finite_field(q, "a")
    verts = collect(K)
    return K, verts
end

"""
Paley graph adjacency matrix P(q) for q prime, q ≡ 1 (mod 4)
Vertices are labeled 0, 1, ..., q-1
Matrix entry A[i,j] corresponds to vertices (i,j)
"""
function paley_adjacency(q::Int; res::Int = 1)

    if q % 4 != res
        error("q must satisfy q ≡ $(res) (mod 4).")
    end

    _, verts = paley_field(q)
    residues = Set(x^2 for x in verts if !iszero(x))
    delete!(residues, 0)

    A = Matrix{Int}(undef, q, q)

    for i in 1:q
        x = verts[i]
        for j in 1:q
            y = verts[j]
            A[i, j] = i != j && x - y in residues ?  1 : 0
        end
    end

    return A
end

function paley_seidel(q::Int; res::Int = 1)
    A = paley_adjacency(q, res=res)
    J = ones(Int, q, q)
    return res == 1 ? J - I - 2A : A-A'
end

"""
The conference matrix, sign indicates symmetrix or skew-symmetric
sign=1 → S is the adjacency matrix of a conference graph
sign=-1 → S is the adjacency matrix of a conference digraph
"""
function paley_conference(q::Int; res::Int = 1, sign::Int = 1)
    S=paley_seidel(q, res=res)
    C = zeros(Int, q+1, q+1)
    C[1,1] = 0
    C[2:end, 1] .= sign
    C[1, 2:end] .= 1
    C[2:end, 2:end] .= S
    return C
end

function paley_signature(q::Int)
    if q % 4 == 1
        C = paley_conference(q, res=1, sign=1)
        return ComplexF64.(C)
    elseif q % 4 == 3
        C = paley_conference(q, res=3, sign=-1)
        return 1im * ComplexF64.(C)
    else
        error("q must be an odd prime power.")
    end
end

function double_signature(S::AbstractMatrix, d::Int; ϵ::Int=1)
    n = size(S, 1)
    size(S, 2) == n || error("S must be square.")

    k = n - 2*d
    println(size(S, 1), " ", d)
    k in (-1, 0, 1) || error("Must have k∈(1,0,-1)")
    c = (n-2*d)*sqrt((n-1)/(d*(n-d)))
    β = - c + im*sqrt(1-c^2)
    I_n = Matrix{ComplexF64}(I, n, n)
    return [S S+β*I_n; S+β'*I_n -S]
end

"""
Gramm matrix for q prime power with 
q ≡ 1 mod 4 (obtaned from the Paley graph) 
q ≡ 3 mod 4 (obtaned from the Paley digraph)
"""
function paley_gram(q::Int; res::Int = 1, sign::Int = 1)
    if !(sign == 1 || sign == -1)
        error("sign must be +1 or -1")
    end
    C = paley_conference(q, res=res, sign=sign)
    return sign == 1 ? I + C / sqrt(q) : I + im * C / sqrt(q)
end

"""
Convert an ETF signature matrix S of a D × N ETF into its Gram matrix.
"""
function signature_to_gram(S::AbstractMatrix, D::Int)
    N = size(S, 1)
    size(S, 2) == N || error("S must be square.")
    α = sqrt((N - D) / (D * (N - 1)))
    I_N = Matrix{ComplexF64}(I, N, N)
    return Hermitian(I_N + α * ComplexF64.(S))
end

"""
Given a Gram matrix G of rank d, return a d × n matrix Φ
such that G = Φ' * Φ.

The columns of Φ are the ETF vectors.
"""
function gram_to_vectors(G; tol=1e-10)
    E = eigen(G)

    idx = findall(λ->real(λ)>(tol), E.values)
    λ = E.values[idx]
    V = E.vectors[:, idx]

    # G = V * Diagonal(λ) * V'
    # so Φ = Diagonal(sqrt.(λ)) * V' gives G = Φ'Φ
    Φ = Diagonal(sqrt.(λ)) * V'

    return Φ
end

function renes_strohmer_signature(q::Int)
    T = paley_seidel(q, res=3)
    dim = size(T, 1)
    G = (im.*T + 1/sqrt(dim)*ones(dim,dim) + sqrt(dim)*I)*sqrt(dim)/(dim+1)
    return (G-I)*sqrt(dim+1)
end

dims_graph=[2, 3, 4, 5, 6, 7, 9, 10, 12, 13]
vecs = Matrix{ComplexF64}[]
for dim in dims_graph
    local q=2*dim-1
    res=mod(q, 4)
    sign = res==1 ? 1 : -1
    local G=paley_gram(q, res=res, sign=sign)
    local GG=gram_to_vectors(G)
    push!(vecs, GG)
    if dim == 7
        local S_= paley_signature(dim)
        local S = double_signature(S_, 4)
        
        local G = signature_to_gram(S, 8)
        local GG = gram_to_vectors(G)
        push!(vecs, GG)
    end
    if dim == 10
        S = renes_strohmer_signature(11)
        S = double_signature(S, 6)
        G = signature_to_gram(S, 11)
        GG = gram_to_vectors(G)
        push!(vecs, GG)
    end


end

"""
Equiangular tight frames for dimensions ranging from 2 to 13
the results are saved as .mat or .npz file
"""

dims_graph=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

data_ = Dict{String, Matrix{ComplexF64}}()
for (i, GG) in enumerate(vecs)
    #data_["$(dims_graph[i])"] = GG
    data_["vec_$(dims_graph[i])"] = GG
end
#npzwrite("ETF_2d.npz", data_)
matwrite("ETF_2d.mat", data_)

overlaps=[]
for i=1:length(vecs)
    push!(overlaps, [abs(vecs[i][:,j]'*vecs[i][:,k]) for j=1:size(vecs[i], 2) for k=1:size(vecs[i], 2) if j!=k])
    println(all(abs.(overlaps[i].-overlaps[i][1]).<0.00001))
end

