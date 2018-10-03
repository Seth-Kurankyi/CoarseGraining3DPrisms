
using Test
using CoarseGraining3DPrisms
using LinearAlgebra

const K = 2
const x = K+1
const y = K/2

a = 0.0


@test delta(1.0,1.0,1.0) == 0

#@show PrA


PrA = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
    ans = numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
    PrA[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)] = ans
end

#PrA = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
ampsA = Float64[]
indxA = Array{Int64,1}[]
for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
    ans = numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
    if ans != 0 
        push!(ampsA,ans)
        push!(indxA, [Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)])
        #push!(lindx,)
    end
end



#PrB = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
#for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
#    ans = numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
   # PrB[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)] = ans
#end


#PrB = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
ampsB = Float64[]
indxB = Array{Int64,1}[]
for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
    ans = numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
    if ans != 0 
        push!(ampsB,ans)
        push!(indxB,[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)])
    end
end

#pos = findall(x-> x!=0, PrA )
#PrANZ = PrA[pos]
#posB = findall(x-> x!=0, PrB )
#PrANZ = PrB[posB]

#@show pos


function loppo(amps::Array{Float64,1},indx::Array{Array{Int64,1},1},jd1::Float64,jd2::Float64,je1::Float64)
    #posa = findall(x-> x!=0, prA )
    #PrANZ = PrA[pos]
    #a = Float64[] 
    #b = Array{Int64,1}[]
    c = []
    for i in 1:length(indx)# position of d1 = 3, d2= 8 and e1 = 4
        if indx[i][3] == Int(2*jd1+1) && indx[i][8] == Int(2*jd2+1) && indx[i][4] == Int(2*je1+1)
            #push!(a,amps[i])
            #push!(b,indx[i])
            push!(c,(indx[i],amps[i]))
        end
    end
    return c
end

function Nsmat(amps::Array{Float64,1},indx::Array{Array{Int64,1},1},jd1::Float64,jd2::Float64,je1::Float64)
    #lopV = loppo(amps,indx,jd1,jd2,je1)[1]
    lopC = loppo(amps,indx,jd1,jd2,je1)
    a = Array{Int64,1}[]
    #for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
    b = Array{Int64,1}[]
    for i in 1:length(lopC) # position of a1 = 1, b1 = 2 g = 12
        unique!(push!(a,[lopC[i][1][1],lopC[i][1][2],lopC[i][1][12]]))
        unique!(push!(b,[lopC[i][1][5],lopC[i][1][6],lopC[i][1][7],lopC[i][1][9],lopC[i][1][10],lopC[i][1][11]]))
    end
    dd = zeros(length(a),length(b))
    for i in 1:length(lopC)
        qa = findall(x-> x == [lopC[i][1][1], lopC[i][1][2], lopC[i][1][12]],a)
        qb = findall(x-> x == [lopC[i][1][5], lopC[i][1][6], lopC[i][1][7], lopC[i][1][9], lopC[i][1][10], lopC[i][1][11] ], b)
        dd[qa,qb] = [lopC[i][2]]
    end
    return a,b,dd
end

function blockReduceA(amps::Array{Float64,1},indx::Array{Array{Int64,1},1},jd1::Float64,jd2::Float64,je1::Float64)
    #a = Float64[] 
    #b = Array{Int64,1}[]
    ampInfo = [] # will contain indices and amplitudes
    for i in 1:length(indx)# position of d1 = 3, d2= 8 and e1 = 4
        if indx[i][3] == Int(2*jd1+1) && indx[i][8] == Int(2*jd2+1) && indx[i][4] == Int(2*je1+1)
            push!(ampInfo,(indx[i],amps[i]))
        end
    end
    nzcol = Array{Int64,1}[]
    #for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
    nzrow = Array{Int64,1}[]
    for i in 1:length(ampInfo) # position of a1 = 1, b1 = 2 g = 12
        unique!(push!(nzcol,[ampInfo[i][1][1],ampInfo[i][1][2],ampInfo[i][1][12]]))
        unique!(push!(nzrow,[ampInfo[i][1][5],ampInfo[i][1][6],ampInfo[i][1][7],ampInfo[i][1][9],ampInfo[i][1][10],ampInfo[i][1][11]]))
    end
    dd = zeros(length(nzcol),length(nzrow))
    if length(dd) != 0
        for i in 1:length(ampInfo)
            qa = findall(x-> x == [ampInfo[i][1][1], ampInfo[i][1][2], ampInfo[i][1][12]],nzcol)
            qb = findall(x-> x == [ampInfo[i][1][5], ampInfo[i][1][6], ampInfo[i][1][7], ampInfo[i][1][9], ampInfo[i][1][10], ampInfo[i][1][11] ], nzrow)
            dd[qa,qb] = [ampInfo[i][2]]
        end
        return nzcol,nzrow, dd
    else
        return 0,0,0
    end
end


#get the index for each of the shared triangular face.. return index values and amplitudes
# gives 30 elements for spins(0,0,0)
function loppo(amps::Array{Float64,1},indx::Array{Array{Int64,1},1},jd1::Float64,jd2::Float64,je1::Float64)
    #posa = findall(x-> x!=0, prA )
    #PrANZ = PrA[pos]
    a = Float64[] 
    b = Array{Int64,1}[]
    c = []
    for i in 1:length(indx)# position of d1 = 3, d2= 8 and e1 = 4
        if indx[i][3] == Int(2*jd1+1) && indx[i][8] == Int(2*jd2+1) && indx[i][4] == Int(2*je1+1)
            #push!(a,amps[i])
            #push!(b,indx[i])
            push!(c,(indx[i],amps[i]))
        end
    end
    return c
end

function Nsmat(amps::Array{Float64,1},indx::Array{Array{Int64,1},1},jd1::Float64,jd2::Float64,je1::Float64)
    #lopV = loppo(amps,indx,jd1,jd2,je1)[1]
    lopC = loppo(amps,indx,jd1,jd2,je1)
    a = Array{Int64,1}[]
    #for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
    b = Array{Int64,1}[]
    for i in 1:length(lopC) # position of a1 = 1, b1 = 2 g = 12
        unique!(push!(a,[lopC[i][1][1],lopC[i][1][2],lopC[i][1][12]]))
        unique!(push!(b,[lopC[i][1][5],lopC[i][1][6],lopC[i][1][7],lopC[i][1][9],lopC[i][1][10],lopC[i][1][11]]))
    end
    dd = zeros(length(a),length(b))
    for i in 1:length(lopC)
        qa = findall(x-> x == [lopC[i][1][1], lopC[i][1][2], lopC[i][1][12]],a)
        qb = findall(x-> x == [lopC[i][1][5], lopC[i][1][6], lopC[i][1][7], lopC[i][1][9], lopC[i][1][10], lopC[i][1][11] ], b)
        dd[qa,qb] = [lopC[i][2]]
    end
    return a,b,dd
end

function blockReduceA(amps::Array{Float64,1},indx::Array{CartesianIndex{12},1},jd1::Float64,jd2::Float64,je1::Float64)
    #a = Float64[] 
    #b = Array{Int64,1}[]
    ampInf = [] # will contain indices and amplitudes
    for i in 1:length(indx)# position of d1 = 3, d2= 8 and e1 = 4
        if indx[i][3] == Int(2*jd1+1) && indx[i][8] == Int(2*jd2+1) && indx[i][4] == Int(2*je1+1)
            push!(ampInf,(indx[i],amps[i]))
        end
    end
    # a and b are index and amplitudes
    c = []
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
        d = []
        for i in 1:length(b) # position of a1 = 1, b1 = 2 g = 12
            if b[i][1] == Int(2*ja1+1) && b[i][2] == Int(2*jb1+1) && b[i][12] == Int(2*jg+1)
                push!(d,a[i])
            end
        end
        if d != []
            push!(c,d)
        end
    end
    if length(c) !=0
        dd = zeros(length(c),length(c[1]))
        for i in 1:length(c)
            c[i] = convert(Array{Float64},c[i])
            dd[i,:] = c[i]
        end
        return dd
    else
        return 0
    end
end

function Nsmate(amps::Array{Float64,1},indx::Array{Array{Int64,1},1},jd1::Float64,jd2::Float64,je1::Float64)
    lopV = loppo(amps,indx,jd1,jd2,je1)[1]
    lopC = loppo(amps,indx,jd1,jd2,je1)[2]
    a = []
    #for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
    b = []
    for i in 1:length(lopC) # position of a1 = 1, b1 = 2 g = 12
        unique!(push!(a,(lopC[i][1],lopC[i][2],lopC[i][12])))
        unique!(push!(b,(lopC[i][5],lopC[i][6],lopC[i][7],lopC[i][9],lopC[i][10],lopC[i][11])))
    end
    dd = zeros(length(a),length(b))
    if length(dd) != 0
        if length(a)*length(b) == length(lopV)
            for i in 1:length(a), j in 1:length(b)
               dd[i,j] = lopV[(i-1)*length(b)+1] 
            end
            return dd
        elseif length(a) == length(b) == length(lopV)
            for i in 1:length(a)
                dd[i,i] = lopV[i]
            end
            return dd
        end
    else
        return 0
    end
    #return dd
end



#get reduced matrix
function smat(amps::Array{Float64,1},indx::Array{CartesianIndex{12},1},jd1::Float64,jd2::Float64,je1::Float64)
    lopV = loppo(amps,indx,jd1,jd2,je1)[1]
    lopC = loppo(amps,indx,jd1,jd2,je1)[2]
    b = []
    #for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
    #a = []
    for i in 1:length(lopC) # position of a1 = 1, b1 = 2 g = 12
        if lopC[i][1] == Int(2*ja1+1) && lopC[i][2] == Int(2*jb1+1) && lopC[i][12] == Int(2*jg+1)
            push!(a,lopV[i])
        end
    end
    if a != []
        push!(b,a)
    end
    #end
    if length(b) !=0
        dd = zeros(length(b),length(b[1]))
        for i in 1:length(b)
            b[i] = convert(Array{Float64},b[i])
            dd[i,:] = b[i]
        end
        return dd
    else
        return 0
    end
end
   

#Combine into one function
function blockReduceA(amps::Array{Float64,1},indx::Array{CartesianIndex{12},1},jd1::Float64,jd2::Float64,je1::Float64)
    #posa = findall(x-> x!=0, prA )
    #PrANZ = PrA[pos]
    a = [] 
    b = []
    for i in 1:length(indx)# position of d1 = 3, d2= 8 and e1 = 4
        if indx[i][3] == Int(2*jd1+1) && indx[i][8] == Int(2*jd2+1) && indx[i][4] == Int(2*je1+1)
            push!(a,amps[i])
            push!(b,indx[i])
        end
    end
    # a and b are index and amplitudes
    c = []
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
        d = []
        for i in 1:length(b) # position of a1 = 1, b1 = 2 g = 12
            if b[i][1] == Int(2*ja1+1) && b[i][2] == Int(2*jb1+1) && b[i][12] == Int(2*jg+1)
                push!(d,a[i])
            end
        end
        if d != []
            push!(c,d)
        end
    end
    if length(c) !=0
        dd = zeros(length(c),length(c[1]))
        for i in 1:length(c)
            c[i] = convert(Array{Float64},c[i])
            dd[i,:] = c[i]
        end
        return dd
    else
        return 0
    end
end


function svdR(prA::Array{Float64,12},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    mat = smat(prA,ja,jb,jc)
    U, s, V = svd(mat)
    #V = ctranspose(V)
    return (U,s,V)
end
    
function svdRA(amps::Array{Float64,1},indx::Array{CartesianIndex{12},1},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    mat = blockReduceA(amps,indx,ja,jb,jc)
    U, s, V = svd(mat)
    #V = ctranspose(V)
    return (U,s,V)
end



@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
    sing = svdRA(ampsA,indxA,ja,jb,jc)[2]
    if sing[1] != 0
        println(sing)
    end
end

#@time println(svdA(PrA,0.0,0.5,0.5)[2][1:6])
@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
    sing = svdA(PrA,ja,jb,jc)[2]
    if sing[1] != 0
        println(sing[1:3])
    end
end

#@time fullUtensorA(PrA)
#@time fullUtensorB(PrB)
#@time fullVtensorA(PrA)
#@time prsmUAVB(PrA,PrB)
#@time println(svdR(PrA,0.,0.5,0.5)[2])
#@time println(svdRA(PrA,0.,0.5,0.5)[2])


#@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
#    sing = svdRA(ampsB,indxB,ja,jb,jc)[2]
#    if sing[1] != 0
#        println(sing)
#    end
#end

#@time blockReduceA(ampsB,indxB,0.,0.5,1.)
#@time println(loppoF(PrA,0.,0.5,0.5))

CoarseGraining3DPrisms.greet()

