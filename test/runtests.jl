
using Test
using CoarseGraining3DPrisms
#using LinearAlgebra

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



PrB = zeros(x,x,x,x,x,x,x,x,x,x,x,x)
for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y, jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
    ans = numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
    PrB[Int(2*ja1+1),Int(2*jb1+1),Int(2*jd1+1),Int(2*je1+1),Int(2*jf1+1),Int(2*ja2+1),Int(2*jb2+1),Int(2*jd2+1),Int(2*je2+1),Int(2*jf2+1),Int(2*jc+1),Int(2*jg+1)] = ans
end


@time svdA(PrA,0.0,0.0,0.0)[2][1:3]
@time fullUtensorA(PrA)
@time fullUtensorB(PrB)
@time fullVtensorA(PrA)
@time prsmUAVB(PrA,PrB)


CoarseGraining3DPrisms.greet()

