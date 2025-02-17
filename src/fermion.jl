
function CalcSpinOcc(s::Int,Nsite::Int)
    return s & (2^Nsite-1), (s & ((2^Nsite-1) << Nsite)) >> Nsite
end

function CalcMagmom(s::Int,Nsite::Int)
    up,dn = CalcSpinOcc(s, Nsite)
    return (count1s(up) - count1s(dn))/2
end

function CalcNd(s::Int,Nsite::Int)
    up,dn = CalcSpinOcc(s, Nsite)
    return count1s(up & dn)
end

