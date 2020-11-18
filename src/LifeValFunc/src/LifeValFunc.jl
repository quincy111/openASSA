module LifeValFunc

using JuliaDB, DataFrames

using LifeValStructures

export Table_to_Policies, Table_Length

function unitPFunction(initVal, incr, rows)
    unitP = [initVal] .* (1 .+ incr) .^ (1 / 12) #unitP[1] is at end of first month
    for i in 2:rows
        append!(unitP, unitP[i-1] .* (1 .+ incr) .^ (1 / 12))
    end
    return(unitP)
end

function Table_to_Policies(policyTable::IndexedTable, productIndex, columnsIndex, PolicyType)
    #global PolicyType
    n = length(policyTable)
    # convert to policy and append to vector
    # change column order to correspond to fields                                                         of Policy
    orderedInput = JuliaDB.select(policyTable, columnsIndex)
    # create product identifier in vector
    product = JuliaDB.select(policyTable, productIndex)
    PolicyVector = Vector{LifeValStructures.Policy}(undef, n)
    # eventually create struct called policyVector
    for i in 1:n
        type = PolicyType[product[i]]
        PolicyVector[i] = type(orderedInput[i]...)
    end

    return PolicyVector
end

function Table_Length(policyTable::Array{IndexedTable}, runsettings)
    n = length(policyTable)
    maxlength = 0
    for i in 1:n
        if length(policyTable[i]) > maxlength
            maxlength = length(policyTable[i])
        end
    end
    reserveList = Array{Float64}(undef, maxlength, nrow(runsettings))
    return reserveList
end

end
