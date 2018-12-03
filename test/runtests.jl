
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
function dataPrA()
    ampsInfo = []
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je1 in 0.:0.5:y
        if delta(jd1,jd2,je1) != 0 
            for jb1 in 0.:0.5:y, jg in 0.:0.5:y
                if delta(jd1,jb1,jg) != 0 
                    for ja1 in 0.:0.5:y
                        if delta(ja1,jd2,jg) != 0 && delta(ja1,jb1,je1) != 0
                            for jc in 0.:0.5:y, jb2 in 0.:0.5:y
                                if delta(jd1,jc,jb2) != 0
                                    for jf1 in 0.:0.5:y
                                        if delta(jf1,jd2,jb2) != 0 && delta(jf1,jc,je1) != 0
                                            for ja2 in 0.:0.5:y, jf2 in 0.:0.5:y
                                                if delta(jd1,ja2,jf2) != 0
                                                    for je2 in 0.:0.5:y
                                                        if delta(jc,je2,jf2) != 0 && delta(jb2,je2,ja2) != 0 #&& numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg)) != 0
                                                            #if numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg)) != 0
                                                            ans = numchop( prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
                                                            if ans != 0
                                                                indx = [ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg]
                                                                push!(ampsInfo,(indx,ans))
                                                                #push!(lindx,)
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
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
    mat = blockReduceA(ampInfo,ja,jb,jc)[3]
    U, s, V = svd(mat)
    #V = ctranspose(V)
    return (U,s,V)
end


#test1 = svdRA(ampsA,indxA,0.,0.,0.5)[3][2]
#@time @show test1

@time datA1 = dataPrA();
#@time datPA = dataPrA();

#@time datB21 = dataB21();
#@time datPB2 = dataPrB2();

#@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
#    sing = svdA(datPA,ja,jb,jc)[2]
#    if sing[1] != 0
#        println(sing[1:2])
#    end
#end

@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
    sing = svdRA(datA1,ja,jb,jc)[2]
    if sing[1] != 0
        println(sing[1:2])
    end
end





#@time @show blockReduceA1(datA1,0.,0.,0.)

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
#    sing = svdRA2(ampsB,indxB,ja,jb,jc)[2]
#    if sing[1] != 0
#        println(sing)
#    end
#end

#@time blockReduceA(ampsB,indxB,0.,0.5,1.)
#@time println(loppoF(PrA,0.,0.5,0.5))

CoarseGraining3DPrisms.greet()

