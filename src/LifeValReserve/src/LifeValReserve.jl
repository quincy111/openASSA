module LifeValReserve

using LifeValStructures

export Calc, Prem

# probability of surviving t years from age x given yearly survival probabilities
function tapxFunction(apx)
    tapx = [1.0]
    n = length(apx)
    for i in 2:n
        append!(tapx, tapx[i-1] .* apx[i])
    end
    return(tapx)
end

function escFunction(initVal, incr, rate)
    esc = [initVal]
    n = length(incr)
    for i in 2:n
        append!(esc, esc[i-1] .* (1 .+ rate .* incr[i]))
    end
    return(esc)
end

# reusable function to set expenses per policy when pricing, because P is unknown.
# common code for any product tyee
function perpolExp(policy, runtype, basisData)
    eBasis = basisData[7][runtype.runno, policy.spc].datfram

    exp_init = eBasis.INIT_EXP_PREM * (1 + runtype.exp_inc)

    #renewal expense also for increasing premiums need to recalc with prem=1
    j = 1:policy.projMax
    m = (runtype.valmonth .+ j) .% 12
    incr = map(m -> ifelse.(m == policy.esc_m, 1, 0), m)
    exp_ren_prem = eBasis.REN_EXP_PREM .* escFunction(1, incr, policy.prem_esc)

    return [exp_init, exp_ren_prem]
end

# reusable function to set basis for calculating reserves
# common code for any product tyee
function setBasis(policy, runtype, basisData, profit)

    # set mort basis
    mort = ifelse(runtype.morta==1, ifelse(policy.gender==0, :UXM_D, :UXF_D),
                            ifelse(policy.gender==0, :UXM_DS, :UXF_DS))

    mortBasis = ifelse(profit==0, basisData[1][runtype.runno, policy.spc].datfram[:, mort], #qxtables
                                   basisData[2][runtype.runno, policy.spc].datfram[:, mort]) #qxptables
    #qxData[runno].datfram[:, [:AGE, :mort]]
    # set wdl basis
    wdl = ifelse(runtype.wdlsa == 1, :uxw, :uxws)
    wdlBasis = ifelse(profit==0, basisData[3][runtype.runno, policy.spc].datfram[:, wdl],
                                basisData[4][runtype.runno, policy.spc].datfram[:, wdl])

    # set v basis
    v = ifelse(runtype.valR == 1, :vt, :vts)
    vBasis = ifelse(profit==0, basisData[5][runtype.runno, policy.spc].datfram[:, v],
                            basisData[6][runtype.runno, policy.spc].datfram[:, v])

    j = 1:policy.projMax
    age = policy.age .+ (j .- 1) /12
    term_wx = policy.termIF .+ j

    #read select mortality
    #for k in termwx:(termwx+policy.projMax)
    #    l += l
    #    if ceil.(Int, k/12) = 1

    #    elseif ceil.(Int, k/12) = 2

    #    else
    #        uxd[k] = mortBasis[floor.(Int, age)]
    #    end
    #end

    uxd = mortBasis[floor.(Int, age)]
    uxw = wdlBasis[ceil.(Int, term_wx/12)]
    aqx = 1 .- exp.((-uxd .-uxw)/12)
    aqxd = uxd ./(uxd + uxw) .*aqx
    awxd = uxw ./(uxd + uxw) .*aqx
    pop!(aqx)
    apx = append!([1.0], 1 .- aqx)
    tapx = tapxFunction(apx)
    vt = vBasis[j]

    #premium and SumAssured escalations
    m = (runtype.valmonth .+ j) .% 12
    incr = map(m -> ifelse.(m == policy.esc_m, 1, 0), m)
    prem = escFunction(policy.premium, incr, policy.prem_esc)
    sa =  escFunction(policy.SumAssured, incr, policy.sa_esc)

    #expenses
    eBasis = ifelse(profit==0, basisData[7][runtype.runno, policy.spc].datfram, basisData[8][runtype.runno, policy.spc].datfram)

    exp_init = (eBasis.INIT_EXP_AMT + eBasis.INIT_EXP_SA * policy.SumAssured +
                eBasis.INIT_EXP_PREM * policy.premium) *
                (1 + runtype.exp_inc)

    exp_ren = eBasis.REN_EXP_AMT * (1 + runtype.exp_inc)
    exp_ren = exp_ren[1]

    exp_ren_prem = eBasis.REN_EXP_PREM .* prem
    exp_ren_prem = exp_ren_prem[1]

    exp_claim = (eBasis.CLM_EXP_AMT .+ eBasis.CLM_EXP_SA .* sa) .*
                 (1 + runtype.exp_inc)

    #sur_claim = (eBasis[runtype.runno, policy.spc].datfram[:SUR_EXP_AMT]  .+
    #             eBasis[runtype.runno, policy.spc].datfram[:SUR_EXP_SA] .* sa) .*
    #             (1 + runtype.exp_inc)

    #mat_claim = (eBasis[runtype.runno, policy.spc].datfram[:MAT_EXP_AMT]  .+
    #              eBasis[runtype.runno, policy.spc].datfram[:MAT_EXP_SA] .* sa) .*
    #              (1 + runtype.exp_inc)

    exp_infl = (1 .+ eBasis.REN_EXP_INFL) * (1 + runtype.expinfl_inc) .- 1
    exp_infl = exp_infl[1]

    expenses = exp_ren_prem .+ (exp_ren/12) * (1 + exp_infl) .^((j .- 1)/12)

    return [aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa]
end

function setAlloc(policy, runno, prem, fmc_sens, basisData, profit)
    global chData

    alloc_P = basisData[9][runno, policy.spc].datfram.ALOC_P
    BOS = basisData[9][runno, policy.spc].datfram.BOS
    fmc = basisData[9][runno, policy.spc].datfram.FMC_PERC

    # set ch basis
    j = 1:policy.projMax
    term_group = policy.termIF .+ j
    rate_alloc = alloc_P[ceil.(Int, term_group/12)]
    rate_BOS = BOS[ceil.(Int, term_group/12)]
    rate_fmc = fmc[ceil.(Int, term_group/12)] * fmc_sens

    prem_alloc = prem .* rate_alloc .* rate_BOS

    return [prem_alloc, rate_fmc]
end

function calcUF(policy, prem_alloc, fmc, unitPrice)
    fundbef = Array{Float64}(undef, policy.projMax)
    fundbeg = Array{Float64}(undef, policy.projMax)
    fundmc = Array{Float64}(undef, policy.projMax)
    units = Array{Float64}(undef, policy.projMax)
    fv = Array{Float64}(undef, policy.projMax)

    units_purc = prem_alloc ./ unitPrice[1:policy.projMax]
    units_start = policy.units_1 #1 fund
    units[1] = (units_start + units_purc[1])
    fundbeg[1] = units[1] .* unitPrice[1]
    fundbef[1] = units[1] .* unitPrice[2]
    fundmc[1] = fundbef[1] * fmc[1]
    fv[1] = fundbef[1] - fundmc[1]

    for j in 2:policy.projMax
        units[j] = fv[j - 1] / unitPrice[j] + units_purc[j]
        fundbeg[j] = units[j] * unitPrice[j]
        fundbef[j] = units[j] * unitPrice[j + 1]
        fundmc[j] = fundbef[j] * fmc[j]
        fv[j] = fundbef[j] - fundmc[j]
    end

    return [units, fv, fundmc, fundbef, fundbeg]#, fundch]
end

function calcNUF(policy, netPrem, fundch, initExp, renExp, aqxd, apx, unitf, unitPrice)
    units_initexp = Array{Float64}(undef, policy.projMax)
    units_renexp = Array{Float64}(undef, policy.projMax)
    units = Array{Float64}(undef, policy.projMax)
    non_u_int = Array{Float64}(undef, policy.projMax)
    mort_rate = Array{Float64}(undef, policy.projMax)
    death_cost = Array{Float64}(undef, policy.projMax)
    cf = Array{Float64}(undef, policy.projMax)
    nonexp = Array{Float64}(undef, policy.projMax)
    reserve = Array{Float64}(undef, policy.projMax)
    fundmc = Array{Float64}(undef, policy.projMax)

    units_purc = netPrem ./ unitPrice[1:policy.projMax]
    units_start = policy.units_99 #NU FUND
    units_initexp = initExp ./ unitPrice[1]
    units_renexp = renExp ./ unitPrice[1:policy.projMax]

    units[1] = units_start + units_purc[1] - units_initexp[1]
    nonexp[1] = units_initexp[1] * unitPrice[1]
    non_u_int[1] = units[1] .* (unitPrice[2] .- unitPrice[1])
    fundmc[1] = fundch[1]
    mort_rate[1] = aqxd[1]
    death_cost[1] = max((policy.SumAssured - unitf[1]),0) * mort_rate[1]
    cf[1] = units[1] * unitPrice[1] + non_u_int[1] + fundmc[1] - death_cost[1]

    for j in 2:policy.projMax
        units[j] = units_purc[j] - units_renexp[j] # + fv[j - 1] / unitPrice[j] ?
        nonexp[j] = units_renexp[j] * unitPrice[j]
        non_u_int[j] = units[j] .* (unitPrice[j+1] .- unitPrice[j])
        fundmc[j] = fundch[j]
        mort_rate[j] = aqxd[j]
        death_cost[j] = max((policy.SumAssured - unitf[j]),0) * mort_rate[j]
        cf[j] = units[j] * unitPrice[j] + non_u_int[j] + fundmc[j] - death_cost[j]
    end

    #calc reserves
    if cf[policy.projMax]< 0
        reserve[policy.projMax] = -cf[policy.projMax] / unitPrice[policy.projMax] * unitPrice[policy.projMax - 1]
    end

    for j in (policy.projMax-1):-1:1
        if cf[j] < 0
            reserve[j]= (-cf[j] + apx[j]*reserve[j + 1]) / unitPrice[j + 1] * unitPrice[j]
        else
            reserve[j]= 0
        end
    end

    return [units, cf, non_u_int, death_cost, nonexp, reserve]#, fundch]
end

function Calc(policy::LifeValStructures.LifeAssurance, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy,runtype, basisData, 0)

    deathClaims = (sa .+ exp_claim) .* aqxd
    cf = deathClaims + expenses .- prem
    reserve = sum(cf .* tapx .* vt)
    return(reserve)
end

function Calc(policy::LifeValStructures.Annuity, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    deathClaims = (sa/12 .+ exp_claim) .* apx
    cf = deathClaims + expenses .- prem
    reserve = sum(cf .* tapx .* vt)
    return(reserve)
end

function Calc(policy::LifeValStructures.TermAssurance, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    deathClaims = (sa .+ exp_claim) .* aqxd
    cf = deathClaims + expenses .- prem
    reserve = sum(cf .* tapx .* vt)
    return(reserve)
end

function Calc(policy::LifeValStructures.EndowmentAssurance, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    deathClaims = (sa .+ exp_claim) .* aqxd
    deathClaims[policy.projMax]=deathClaims[policy.projMax]/aqxd[policy.projMax]
    cf = deathClaims + expenses .- prem
    reserve = sum(cf .* tapx .* vt)
    return(reserve)
end

function Calc(policy::LifeValStructures.EndowmentAssuranceUL, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    prem_alloc, fmc = setAlloc(policy, runtype.runno, prem, runtype.fmc_sens, basisData, 0)

    #currently only 1 fund, 100% allocation
    fmcp = (1 .+ fmc) .^ (1/12) .- 1 #fmc is annual, need to make it monthly

    #unit fund
    units, fundv, fundch, fundb, fundbeg = calcUF(policy, prem_alloc, fmcp, unitPrices[runtype.runno,1,1,:])
    #output=table((time = 1:policy.projMax, alloc = prem_alloc, fundvalbeg = fundbeg, fundvalbef = fundb, fmc = fundch, fundvalend = fundv); pkey = [:time])
    #writedlm("bu2.csv", output, ',')

    #non-unit fund
    non_units, noncf, nint, death_cost, nuexp, non_reserve = calcNUF(policy, (prem - prem_alloc), fundch, exp_init, expenses, aqxd, apx, fundv, unitPrices[runtype.runno,2,1,:])
    #output=table((time = 1:policy.projMax, alloc = (prem - prem_alloc), expenses = nuexp, interest = nint, fundmc = fundch, d_cost = death_cost, cf = noncf, v = non_reserve); pkey = [:time])
    #writedlm("nbu1.csv", output, ',')

    #profit test BASIS
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 1)
    prem_alloc, fmc = setAlloc(policy, runtype.runno, prem, runtype.fmc_sens, basisData, 1)
    #currently only 1 fund, 100% allocation
    fmcp = (1 .+ fmc) .^ (1/12) .- 1 #fmc is annual, need to make it monthly

    #unit fund
    units_p, fundv_p, fundch_p, fundb_p, fundbeg_p = calcUF(policy, prem_alloc, fmcp, unitPrices[runtype.runno,1,2,:])
    #output=table((time = 1:policy.projMax, alloc = prem_alloc, fundvalbeg = fundbeg_p, fundvalbef = fundb_p, fmc = fundch_p, fundvalend = fundv_p); pkey = [:time])
    #writedlm("bup2.csv", output, ',')

    non_units_p, noncf_p, nint_p, death_cost_p, nuexp_p, non_reserve_p = calcNUF(policy, (prem - prem_alloc), fundch_p, exp_init, expenses, aqxd, apx, fundv_p, unitPrices[runtype.runno,2,2,:])
    #output=table((time = 1:policy.projMax, alloc = (prem - prem_alloc), expenses = nuexp_p, interest = nint_p, fundmc = fundch_p, d_cost = death_cost_p, cf = noncf_p, v = non_reserve_p); pkey = [:time])
    #writedlm("nbup1.csv", output, ',')

    #TODO Profitvectors and signatures

    reserve = fundv[1]+non_reserve[2]  #2 because it is at start of year 1.
    return(reserve)
end

function Prem(policy::LifeValStructures.LifeAssurance, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy,runtype, basisData, 0)

    deathClaims = (sa .+ exp_claim) .* aqxd
    cf = deathClaims + expenses

    fixed, ren = perpolExp(policy,runtype, basisData)

    premium = (sum(cf .* tapx .* vt) .+ exp_init) ./ (sum((1 .- ren) .* tapx .* vt) .- fixed)
    return(premium[1])
end

function Prem(policy::LifeValStructures.Annuity, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    deathClaims = (sa/12 .+ exp_claim) .* apx
    cf = deathClaims + expenses

    fixed, ren = perpolExp(policy,runtype, basisData)

    premium = (sum(cf .* tapx .* vt) .+ exp_init) ./ (1 .- fixed) #lump sum
    return(premium[1])
end

function Prem(policy::LifeValStructures.TermAssurance, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    deathClaims = (sa .+ exp_claim) .* aqxd
    cf = deathClaims + expenses

    fixed, ren = perpolExp(policy,runtype, basisData)

    premium = (sum(cf .* tapx .* vt) .+ exp_init) ./ (sum((1 .- ren) .* tapx .* vt) .- fixed)
    return(premium[1])
end

function Prem(policy::LifeValStructures.EndowmentAssurance, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    deathClaims = (sa .+ exp_claim) .* aqxd
    deathClaims[policy.projMax]=deathClaims[policy.projMax]/aqxd[policy.projMax]
    cf = deathClaims + expenses

    fixed, ren = perpolExp(policy,runtype, basisData)

    premium = (sum(cf .* tapx .* vt) .+ exp_init) ./ (sum((1 .- ren) .* tapx .* vt) .- fixed)

    return(premium[1])
end

function Prem(policy::LifeValStructures.EndowmentAssuranceUL, runtype, basisData, unitPrices)
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 0)

    prem_alloc, fmc = setAlloc(policy, runtype.runno, prem, runtype.fmc_sens, basisData, 0)

    #currently only 1 fund, 100% allocation
    fmcp = (1 .+ fmc) .^ (1/12) .- 1 #fmc is annual, need to make it monthly

    #unit fund
    units, fundv, fundch, fundb, fundbeg = calcUF(policy, prem_alloc, fmcp, unitPrices[runtype.runno,1,1,:])
    #output=table((time = 1:policy.projMax, alloc = prem_alloc, fundvalbeg = fundbeg, fundvalbef = fundb, fmc = fundch, fundvalend = fundv); pkey = [:time])
    #writedlm("bu2.csv", output, ',')

    #non-unit fund
    non_units, noncf, nint, death_cost, nuexp, non_reserve = calcNUF(policy, (prem - prem_alloc), fundch, exp_init, expenses, aqxd, apx, fundv, unitPrices[runtype.runno,2,1,:])
    #output=table((time = 1:policy.projMax, alloc = (prem - prem_alloc), expenses = nuexp, interest = nint, fundmc = fundch, d_cost = death_cost, cf = noncf, v = non_reserve); pkey = [:time])
    #writedlm("nbu1.csv", output, ',')

    #profit test BASIS
    aqxd, tapx, vt, apx, expenses, exp_init, exp_claim, prem, sa = setBasis(policy, runtype, basisData, 1)
    prem_alloc, fmc = setAlloc(policy, runtype.runno, prem, runtype.fmc_sens, basisData, 1)
    #currently only 1 fund, 100% allocation
    fmcp = (1 .+ fmc) .^ (1/12) .- 1 #fmc is annual, need to make it monthly

    #unit fund
    units_p, fundv_p, fundch_p, fundb_p, fundbeg_p = calcUF(policy, prem_alloc, fmcp, unitPrices[runtype.runno,1,2,:])
    #output=table((time = 1:policy.projMax, alloc = prem_alloc, fundvalbeg = fundbeg_p, fundvalbef = fundb_p, fmc = fundch_p, fundvalend = fundv_p); pkey = [:time])
    #writedlm("bup2.csv", output, ',')

    non_units_p, noncf_p, nint_p, death_cost_p, nuexp_p, non_reserve_p = calcNUF(policy, (prem - prem_alloc), fundch_p, exp_init, expenses, aqxd, apx, fundv_p, unitPrices[runtype.runno,2,2,:])
    #output=table((time = 1:policy.projMax, alloc = (prem - prem_alloc), expenses = nuexp_p, interest = nint_p, fundmc = fundch_p, d_cost = death_cost_p, cf = noncf_p, v = non_reserve_p); pkey = [:time])
    #writedlm("nbup1.csv", output, ',')

    #TODO Profitvectors and signatures

    #TODO calc prem
    premium = fundv[1]+non_reserve[2]  #2 because it is at end of year 1.
    return(premium)
end

end
