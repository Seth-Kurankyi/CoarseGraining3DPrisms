
using Test
using CoarseGraining3DPrisms
using LinearAlgebra

const K = 5
const x = K+1
const y = K/2

a = 0.0


#@test delta(1.0,1.0,1.0) == 1

#@show PrA


function dataPrA()
    PrA = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, 
                jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
        ans = numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
        PrA[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)] = ans
    end
    return PrA
end

function dataPrB2()
    PrB = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, 
                jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
        ans = numchop( prismB2(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
        PrB[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)] = ans
    end
    return PrB
end



#PrA = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
function dataA1()
    ampsInfo = []
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
        if numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg)) != 0
            ans = numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
            indx = [ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg]
            push!(ampsInfo,(indx,ans))
            #push!(lindx,)
        end
    end
    return ampsInfo
end

function dataB21()
    ampsInfo = []
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
        if numchop( prismB2(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg)) != 0
            ans = numchop( prismB2(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
            indx = [ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg]
            push!(ampsInfo,(indx,ans))
            #push!(lindx,)
        end
    end
    return ampsInfo
end

function blockReduceA1(ampJ::Array{Any,1},jd1::Float64,jd2::Float64,je1::Float64)
    #a = Float64[] 
    #b = Array{Int64,1}[]
    #ampInfo = [] # will contain indices and amplitudes
    # position of d1 = 3, d2= 8 and e1 = 4
    ampInfo = ampJ[findall(x-> x[1][3] == jd1 && x[1][8] == jd2 && x[1][4] == je1 , ampJ)]
    nzcol = Array{Float64,1}[]
    #for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jg in 0.:0.5:y
    nzrow = Array{Float64,1}[]
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


#PrB = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
#for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
#    ans = numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
   # PrB[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)] = ans
#end


#PrB = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
#ampsB = Float64[]
#indxB = Array{Int64,1}[]
#for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
#    if numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg)) != 0
#        ans = numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
#        push!(ampsB,ans)
#        push!(indxB,[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)])
#    end
#end

#pos = findall(x-> x!=0, PrA )
#PrANZ = PrA[pos]
#posB = findall(x-> x!=0, PrB )
#PrANZ = PrB[posB]

#@show pos

#get the index for each of the shared triangular face.. return index values and amplitudes
# gives 30 elements for spins(0,0,0)








    
function svdRA(ampInfo::Array{Any,1},ja::Float64,jb::Float64,jc::Float64)# parameters d1,d2,e1
    mat = blockReduceA1(ampInfo,ja,jb,jc)[3]
    U, s, V = svd(mat)
    #V = ctranspose(V)
    return (U,s,V)
end

#test1 = svdRA(ampsA,indxA,0.,0.,0.5)[3][2]
#@time @show test1

@time datA1 = dataA1();
@time datPA = dataPrA();

@time datB21 = dataB21();
@time datPB2 = dataPrB2();

#@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
#    sing = svdA(datPA,ja,jb,jc)[2]
#    if sing[1] != 0
#        println(sing[1:2])
#    end
#end

#@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
#    sing = svdRA(datA1,ja,jb,jc)[2]
#    if sing[1] != 0
#        println(sing[1:2])
#    end
#end


#@time println(svdA(PrA,0.0,0.5,0.5)[2][1:6])
#@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
#    sing = svdA(PrA,ja,jb,jc)[2]
#    if sing[1] != 0
#        println(sing[1:3])
#    end
#end

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

