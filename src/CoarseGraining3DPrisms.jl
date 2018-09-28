module CoarseGraining3DPrisms
using LinearAlgebra

greet() = print("Hello World!")


export numchop, delta, qn, qnfact, visqrt, trian, RacWig6j, Fsymb, Gsymb, Gsymb2


numchop(num::Real)=(abs(num) >= 1000*eps() ? num : zero(Real))

function numchop(num::Complex)
    numchop(imag(num)) == 0 ? numchop(real(num)) : Complex(numchop(real(num)), numchop(imag(num)))
end


const K = 2
const x = K+1
const y = K/2




#Define delta_{ijkl} -> coupling rules
function delta(i::Float64,j::Float64,k::Float64)
    sol = 0
    if i <=(j+k) && j<=(i+k) && k<=(i+j) && i+j+k <= K && 2*(i+j+k)%2 ==0
        sol = 1
    end
    return sol
end

#Define quantum numbers qn (this is real)
function qn(n::Float64)
    sol = (exp(pi*n*im/(K+2)) - exp(-pi*n*im/(K+2))) / (exp(pi*im/(K+2)) - exp(-pi*im/(K+2)))
    return real(sol)
end

#Define qn factorial
function qnfact(n::Float64)
    sol = 1
    for i in 1:n
        sol *= qn(i)
    end
    return sol
end

#Define square root of quantum dimension
function visqrt(i::Float64)
    sol = ((-1+0im)^i )*sqrt(qn(2*i+1))
    return sol
end

#Define triangle equality
function trian(i::Float64,j::Float64,k::Float64)
    sol = 0
    if delta(i,j,k) == 1
        sol = delta(i,j,k)*sqrt(qnfact(i+j-k)*qnfact(i-j+k)*qnfact(-i+j+k)/qnfact(i+j+k+1))
    end
    return sol
end

# Define Racah-Wigner six-j symbol
function RacWig6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    a = i+j+m; b = i+l+n; c = k+j+n; d = k+l+m;  e = i+j+k+l; f = i+k+m+n; g = j+l+m+n
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sumz = 0
        for z in max(a,b,c,d):min(e,f,g)
            sumz += (-1)^z *qnfact(z+1)/
                ((qnfact(e-z)*qnfact(f-z)*qnfact(g-z))* (qnfact(z-a)*qnfact(z-b)*qnfact(z-c)*qnfact(z-d)))
        end
        sol = trian(i,j,m)*trian(i,l,n)*trian(k,j,n)*trian(k,l,m)*sumz
    end
    return sol
end

#Define F-symbol
function Fsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = (-1+0im)^(i+j+k+l)*sqrt(qn(2*m+1)*qn(2*n+1)) * RacWig6j(i,j,m,k,l,n)
    end
    return sol
end

#Define G-symbol
function Gsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = Fsymb(i,j,m,k,l,n) /(visqrt(m)*visqrt(n))
    end
    return sol
end

#Define G-symbol 2
function Gsymb2(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = (-1+0im)^(i+j+k+l+m+n)*RacWig6j(i,j,m,k,l,n)
    end
    return sol
end


export prismA, prismB, blocksPrA, blocksPrB, svdA, svdB

#d1,d2,e1 is triangle we shall cut along
function prismA(ja1::Float64,jb1::Float64,jd1::Float64,je1::Float64,jf1::Float64,ja2::Float64,jb2::Float64,
		                 jd2::Float64,je2::Float64,jf2::Float64,jc::Float64,jg::Float64)
    sol = 0
    if delta(jd1,jd2,je1) != 0 && delta(jd1,jb1,jg) != 0 && delta(ja1,jd2,jg) != 0 && delta(ja1,jb1,je1) != 0 && delta(jd1,jc,jb2) != 0 && delta(jf1,jd2,jb2) != 0 && delta(jf1,jc,je1) != 0 && delta(jd1,ja2,jf2) != 0 && delta(jc,je2,jf2) != 0 && delta(jb2,je2,ja2) != 0
        dims = visqrt(ja1)*visqrt(jb1)*visqrt(jg)*visqrt(jf1)*visqrt(jc)*visqrt(jb2)*visqrt(je2)*visqrt(ja2)*visqrt(jf2)*visqrt(jd1)*visqrt(jd2)*visqrt(je1)
        sol =  dims*Gsymb(jd1,jd2,je1,ja1,jb1,jg) * Gsymb(jd1,jd2,je1,jf1,jc,jb2) * Gsymb(jd1,jc,jb2,je2,ja2,jf2)
    end
    return sol
end

function prismB(ja1::Float64,jb1::Float64,jd1p::Float64,je1::Float64,jf1::Float64,ja2::Float64,jb2::Float64,
		                 jd2p::Float64,je2::Float64,jf2::Float64,jc::Float64,jg::Float64)
    sol = 0
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y
        sol += numchop(Fsymb(ja1,jg,jd2,jb2,jf1,jd2p)*Fsymb(jb1,jg,jd1,ja2,jf2,jd1p)*prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
    end
    return sol
end


function blocksPrA(Prsm::Array{Float64,12},jd1::Float64,jd2::Float64,je1::Float64)
    mat = zeros(x^3,x^6)
    row = 1
    col = 1
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
        for jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y
            mat[row,col] =  Prsm[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)]
            col += 1
        end
        col = 1
        row += 1
    end
    return mat
end

function blocksPrB(Prsm::Array{Float64,12},jd2::Float64,jd1::Float64,je2::Float64)
    mat = zeros(x^3,x^6)
    row = 1
    col = 1
    for ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jg in 0.:0.5:y
        for jf2 in 0.:0.5:y, ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, jc in 0.:0.5:y
            mat[row,col] =  Prsm[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)]
            col += 1
        end
        col = 1
        row += 1
    end
    return mat
end







function svdA(prA::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    mat = blocksPrA(prA,ja,jb,jc)
    U, s, V = svd(mat)
    #V = ctranspose(V)
    return (U,s,V)
end

function svdB(prB::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e2
    mat = blocksPrB(prB,ja,jb,jc)
    U, s, V = svd(mat)
    return (U,s,V)
end



export UtensorA, UtensorB, VtensorA, VtensorB


function UtensorA(prA::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    UT = zeros(x,x,x)
    U = svdA(prA,ja,jb,jc)[1]
    s = svdA(prA,ja,jb,jc)[2]
    index = 1
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
        UT[Int(2*ja1+1),Int(2*jb1+1),Int(2*jg+1)] = numchop(U[index,1]*sqrt(s[1]))
        index += 1
    end
    return UT
end

function UtensorB(prB::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    UT = zeros(x,x,x)
    U = svdB(prB,ja,jb,jc)[1]
    s = svdB(prB,ja,jb,jc)[2]
    index = 1
    for ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jg in 0.:0.5:y
        UT[Int(2*ja2+1),Int(2*jb2+1),Int(2*jg+1)] = numchop(U[index,1]*sqrt(s[1]))
        index += 1
    end
    return UT
end

#VVtensors
function VtensorA(prA::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    VT = zeros(x,x,x,x,x,x)
    V = svdA(prA,ja,jb,jc)[3]
    s = svdA(prA,ja,jb,jc)[2]
    index = 1
    for jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y
        VT[Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1)] = numchop(V[index,1]*sqrt(s[1]))
        index += 1
    end
    return VT
end

function VtensorB(prB::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    VT = zeros(x,x,x,x,x,x)
    V = svdB(prB,ja,jb,jc)[3]
    s = svdB(prB,ja,jb,jc)[2]
    index = 1
    for jf2 in 0.:0.5:y, ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, jc in 0.:0.5:y
        VT[Int(2*jf2+1),Int(2*ja1+1),Int(2*jb1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*jc+1)] = numchop(V[index,1]*sqrt(s[1]))
        index += 1
    end
    return VT
end



export  fullUtensorA,  fullUtensorB, fullVtensorA, fullVtensorA

function fullUtensorA(prA::Array{Float64,12})
    fullUT = zeros(x,x,x,x,x,x)
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je1 in 0.:0.5:y
        fullUT[Int(2*jd1+1),Int(2*jd2+1),Int(2*je1+1),:,:,:] = UtensorA(prA,jd1,jd2,je1)
    end
    return fullUT
end

function fullVtensorA(prA::Array{Float64,12})
    fullVT = zeros(x,x,x,x,x,x,x,x,x)
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je1 in 0.:0.5:y
        fullVT[Int(2*jd1+1),Int(2*jd2+1),Int(2*je1+1),:,:,:,:,:,:] = VtensorA(prA,jd1,jd2,je1)
    end
    return fullVT
end

function fullUtensorB(prB::Array{Float64,12})
    fullUT = zeros(x,x,x,x,x,x)
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y
        fullUT[Int(2*jd1+1),Int(2*jd2+1),Int(2*je2+1),:,:,:] = UtensorB(prB,jd1,jd2,je2)
    end
    return fullUT
end

function fullVtensorB(prB::Array{Float64,12})
    fullVT = zeros(x,x,x,x,x,x,x,x,x)
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y
        fullVT[Int(2*jd1+1),Int(2*jd2+1),Int(2*je2+1),:,:,:,:,:,:] = VtensorB(prB,jd1,jd2,je2)
    end
    return fullVT
end


export prsmUAVB, prsmUBVA, prsmEff

function prsmUAVB(prA::Array{Float64,12},prB::Array{Float64,12})
    UA = fullUtensorA(prA)
    VB = fullVtensorB(prB)
    prsm = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
    for jb1 in 0.:0.5:y, jc1 in 0.:0.5:y, jd1 in 0.:0.5:y, jh1 in 0.:0.5:y, je11 in 0.:0.5:y, jf12 in 0.:0.5:y, ja2 in 0.:0.5:y, jf2 in 0.:0.5:y, jh2 in 0.:0.5:y, je21 in 0.:0.5:y, jg in 0.:0.5:y, jn in 0.:0.5:y
        ans = UA[Int(2*jb1+1),Int(2*jd1+1),Int(2*jg+1),Int(2*jn+1),Int(2*jh1+1),Int(2*je11+1)] *
                VB[Int(2*jd1+1),Int(2*ja2+1),Int(2*jf2+1),Int(2*je21+1),Int(2*jn+1),Int(2*jh2+1),Int(2*jf12+1),Int(2*je11+1),Int(2*jc1+1)]
        prsm[Int(2*jb1+1),Int(2*jc1+1),Int(2*jd1+1),Int(2*jh1+1),Int(2*je11+1),Int(2*jf12+1),Int(2*ja2+1),Int(2*jf2+1),Int(2*jh2+1),Int(2*je21+1),Int(2*jg+1),Int(2*jn+1)] = ans
    end
    return prsm
end



function prsmUBVA(prA::Array{Float64,12},prB::Array{Float64,12})
    UB = fullUtensorB(prB)
    VA = fullVtensorA(prA)
    prsm = Array{Float64}(x,x,x,x,x,x,x,x,x,x,x,x)
    for jb1 in 0.:0.5:y, jc1 in 0.:0.5:y, jd1 in 0.:0.5:y, jh1 in 0.:0.5:y, je11 in 0.:0.5:y, jf12 in 0.:0.5:y, ja2 in 0.:0.5:y, jf2 in 0.:0.5:y, jh2 in 0.:0.5:y, je21 in 0.:0.5:y, jg in 0.:0.5:y, jn in 0.:0.5:y
        ans = UB[Int(2*jb1+1),Int(2*jd1+1),Int(2*jg+1),Int(2*jn+1),Int(2*jh1+1),Int(2*je11+1)] *
                VA[Int(2*jd1+1),Int(2*ja2+1),Int(2*jf2+1),Int(2*je21+1),Int(2*jn+1),Int(2*jh2+1),Int(2*jf12+1),Int(2*je11+1),Int(2*jc1+1)]
        prsm[Int(2*jb1+1),Int(2*jc1+1),Int(2*jd1+1),Int(2*jh1+1),Int(2*je11+1),Int(2*jf12+1),Int(2*ja2+1),Int(2*jf2+1),Int(2*jh2+1),Int(2*je21+1),Int(2*jg+1),Int(2*jn+1)] = ans
    end
    return prsm
end


function prsmEff(prA::Array{Float64,12},prB::Array{Float64,12})
    pr1 = prsmUAVB(prA,prB)
    pr2 = prsmUBVA(prA,prB)
    prsm = zeros(x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x)
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jc1 in 0.:0.5:y, jd1 in 0.:0.5:y, je11 in 0.:0.5:y, je12 in 0.:0.5:y,jf1 in 0.:0.5:y, jh1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jc2 in 0.:0.5:y, jd2 in 0.:0.5:y, je21 in 0.:0.5:y,je22 in 0.:0.5:y, jf2 in 0.:0.5:y, jh2 in 0.:0.5:y, jf12 in 0.:0.5:y, jg in 0.:0.5:y
        sum = 0
        for jn in 0.:0.5:y
            ans1 = pr1[Int(2*jb1+1),Int(2*jc1+1),Int(2*jd1+1),Int(2*jh1+1),Int(2*je11+1),Int(2*jf12+1),Int(2*ja2+1),Int(2*jf2+1),Int(2*jh2+1),Int(2*je21+1),Int(2*jg+1),Int(2*jn+1)]
            ans2 = pr2[Int(2*jb2+1),Int(2*jc2+1),Int(2*jd2+1),Int(2*jh2+1),Int(2*je22+1),Int(2*jf12+1),Int(2*ja1+1),Int(2*jf1+1),Int(2*jh1+1),Int(2*je12+1),Int(2*jg+1),Int(2*jn+1)]
            sum+= (ans1*ans2)
        end
        prsm[Int(2*jb1+1),Int(2*jc1+1),Int(2*jd1+1),Int(2*je11+1),Int(2*ja2+1),Int(2*jf2+1),Int(2*je21+1),Int(2*jb2+1),Int(2*jc2+1),Int(2*jd2+1),Int(2*je12+1),Int(2*ja1+1),Int(2*jf1+1),Int(2*je22+1),Int(2*jh1+1),Int(2*jf12+1),Int(2*jh2+1),Int(2*jg+1)] = sum
    end
    return prsm
end



end # module
