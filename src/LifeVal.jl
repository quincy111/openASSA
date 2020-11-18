module LifeVal

using DrWatson
@quickactivate "LifeVal"

using JuliaDB, DelimitedFiles, DataFrames, XLSX, ProgressMeter #, Traceur
using Random, Distributions
#using IterableTables
using Pkg, Pkg.Artifacts
#using CUDA

using LifeValStructures
using LifeValData
using LifeValFunc
using LifeValReserve
using LifeValGenerate

#fooifier_path() = joinpath(artifact"fooifier", "bin", "fooifier" * (Sys.iswindows() ? ".exe" : ""))
  #println("Running the artifact")
  #res = read(`$(fooifier_path()) 5 10`, String)
  #println("The result of 2*5^2 - 10 == $res")

# --------------------------------------------------- #
#THE FOLLOWING SECTION SHOULD BE RUN ONCE FOR ANY CUSTOM INPUT
#dataFile = "ulData2.csv"
#policyData3 = loadtable(dataFile)
#save(policyData3, "ulData3")

function real_main(filePath, fileName)
    #@show ARGS
    #@show Base.PROGRAM_FILE
    #@show DEPOT_PATH
    #@show LOAD_PATH
    #@show pwd()
    #@show Base.active_project()
    #@show Threads.nthreads()
    #@show Sys.BINDIR
    #display(Base.loaded_modules)
    #println()

    #filePath=joinpath(dirname(@__FILE__), "test\\Basis2020127.xlsx")
    #filePath = joinpath(Pkg.dir("LifeVal"),"input") * "\\"
    #filePath = "c:\\julialifeval\\lifeval\\input\\"
    #import LifeVal; joinpath(dirname(pathof(LifeVal)), "..", paths...)
    #fileName = filePath * "Basis20201285v.xlsx"

    runsettings, spcode, names = LifeValData.Initialise(filePath, fileName)
    policyData, basisData, fData = LifeValData.Read_Data(runsettings, spcode, names, filePath)

    PolicyType = [LifeValStructures.LifeAssurance, LifeValStructures.TermAssurance, LifeValStructures.Annuity, LifeValStructures.EndowmentAssurance, LifeValStructures.LifeAssuranceUL, LifeValStructures.EndowmentAssuranceUL]
    inputOrder = (:id, :sex, :age, :Term, :Term_if, :Premium, :Sum_Assured, :spc,
        :prem_esc, :sa_esc, :esc_m, :units_1, :units_99, :alloc_1, :alloc_99)

    reserveList = LifeValFunc.Table_Length(policyData, runsettings)

    o = nrow(runsettings) #ProgressMeter
    q = Progress(o, 1) #ProgressMeter

    @time for runtype in eachrow(runsettings)
        # eventually create struct called policyList
        policyList = LifeValFunc.Table_to_Policies(policyData[runtype.runno], :prod, inputOrder, PolicyType)

        #CALC UNIT PRICES, this must be set per runno and only required up to max projTerm
        #ulList = filter!(policy.spcode > 49, policyList)
        #not sure how to remove this hardcoding!
        glob_proj_max = 421 #maximum(policyList[80001].projMax) + 1 #because permiums paid upfront, need 1 more at end
        unitpriceList = Array{Float64}(undef, nrow(runsettings), length(fData), 2, glob_proj_max)

        #for fund in fData[runtype.runno,:]
        #    unitpriceList[runtype.runno, fund.datfram.ID,1,:] = LifeValFunc.unitPFunction(fund.datfram.UP, (fund.datfram.V_IR - fund.datfram.V_IC), glob_proj_max)
        #    unitpriceList[runtype.runno, fund.datfram.ID,2,:] = LifeValFunc.unitPFunction(fund.datfram.UP / (1 + fund.datfram.V_IR .- fund.datfram.V_IC), (fund.datfram.V_IR .- fund.datfram.V_IC), glob_proj_max)
        #end
        #unitpriceList[runtype.runno, x.ID, 1, :] = map(x -> LifeValFunc.unitPFunction((x[1].datfram.UP, (x[1].datfram.V_IR .- x[1].datfram.V_IC), glob_proj_max)...), fData[runtype.runno,:])
        #unitpriceList[runtype.runno, x.ID, 2, :] = map(x -> LifeValFunc.unitPFunction((x[2].datfram.UP /(1+(x[2].datfram.V_IR .- x[2].datfram.V_IC)), (x[2].datfram.UP, (x[2].datfram.V_IR .- x[2].datfram.V_IC), glob_proj_max)...), fData[runtype.runno,:])

        for fund in fData[runtype.runno,:]
            unitpriceList[runtype.runno, fund.datfram.ID, 1, 1] = fund.datfram.UP #.* (1 .+ fund.datfram.V_IR .- fund.datfram.V_IC) .^ (1/12)
            unitpriceList[runtype.runno, fund.datfram.ID, 2, 1] = fund.datfram.UP
            for t in 2:glob_proj_max
                unitpriceList[runtype.runno, fund.datfram.ID, 1, t] = unitpriceList[runtype.runno, fund.datfram.ID, 1, (t-1)]  .* (1 .+ fund.datfram.V_IR .- fund.datfram.V_IC) .^ (1/12)
                unitpriceList[runtype.runno, fund.datfram.ID, 2, t] = unitpriceList[runtype.runno, fund.datfram.ID, 2, (t-1)]  .* (1 .+ fund.datfram.P_IR .- fund.datfram.P_IC) .^ (1/12)
            end
        end

        reserveList[:, runtype.runno] = map(x -> LifeValReserve.Calc(x, (runtype, basisData, unitpriceList)...), policyList)

        #outputtab = table(pol = policyList, res = reserveList)
        #writelm("out1.csv", outputtab, ",")
        #writedlm(filePath * "basis65-" * runtype.runno * ".csv", reserveList, ',')
        next!(q) #ProgressMeter
    end

    #combined = table(policyData,reserveList)

    writedlm(filePath * "basis65.csv", reserveList, ',')
end

function real_main_gpu()
    #filePath=joinpath(dirname(@__FILE__), "test\\Basis2020127.xlsx")
    #filePath = joinpath(Pkg.dir("LifeVal"),"input") * "\\"
    #filePath = "c:\\julialifeval\\lifeval\\input\\"
    #import LifeVal; joinpath(dirname(pathof(LifeVal)), "..", paths...)
    #fileName = filePath * "Basis20201281.xlsx"

    runsettings, spcode, names = LifeValData.Initialise(filePath, fileName)
    policyData, basisData, fData = LifeValData.Read_Data(runsettings, spcode, names, filePath)

    PolicyType = [LifeValStructures.LifeAssurance, LifeValStructures.TermAssurance, LifeValStructures.Annuity, LifeValStructures.EndowmentAssurance, LifeValStructures.LifeAssuranceUL, LifeValStructures.EndowmentAssuranceUL]
    inputOrder = (:id, :sex, :age, :Term, :Term_if, :Premium, :Sum_Assured, :spc,
        :prem_esc, :sa_esc, :esc_m, :units_1, :units_99, :alloc_1, :alloc_99)

    reserveList = LifeValFunc.Table_Length(policyData, runsettings)
    policyList1 = CuArray(policyData)

    o = nrow(runsettings) #ProgressMeter
    q = Progress(o, 1) #ProgressMeter

    @time for runtype in eachrow(runsettings)
        # eventually create struct called policyList
        policyList = LifeValFunc.Table_to_Policies(policyData[runtype.runno], :prod, inputOrder, PolicyType)
        policyList = CuArray(policyList)

        #CALC UNIT PRICES, this must be set per runno and only required up to max projTerm
        #ulList = filter!(policy.spcode > 49, policyList)
        #not sure how to remove this hardcoding!
        glob_proj_max = 421 #maximum(policyList[80001].projMax) + 1 #because permiums paid upfront, need 1 more at end
        unitpriceList = Array{Float64}(undef, nrow(runsettings), length(fData), 2, glob_proj_max)
        #for fund in fData[runtype.runno,:]
        #    unitpriceList[runtype.runno, fund.datfram.ID,:] = unitPFunction(fund.datfram.UP, (fund.datfram.V_IR .- fund.datfram.V_IC), glob_proj_max)
        #end

        for fund in fData[runtype.runno,:]
            unitpriceList[runtype.runno, fund.datfram.ID, 1, 1] = fund.datfram.UP #.* (1 .+ fund.datfram.V_IR .- fund.datfram.V_IC) .^ (1/12)
            unitpriceList[runtype.runno, fund.datfram.ID, 2, 1] = fund.datfram.UP
            for t in 2:glob_proj_max
                unitpriceList[runtype.runno, fund.datfram.ID, 1, t] = unitpriceList[runtype.runno, fund.datfram.ID, 1, (t-1)]  .* (1 .+ fund.datfram.V_IR .- fund.datfram.V_IC) .^ (1/12)
                unitpriceList[runtype.runno, fund.datfram.ID, 2, t] = unitpriceList[runtype.runno, fund.datfram.ID, 2, (t-1)]  .* (1 .+ fund.datfram.P_IR .- fund.datfram.P_IC) .^ (1/12)
            end
        end

        reserveList[:, runtype.runno] = map(x -> LifeValReserve.Calc(x, (runtype, basisData, unitpriceList)...), policyList)

        #outputtab = table(pol = policyList, res = reserveList)
        #writelm("out1.csv", outputtab, ",")
        #writedlm(filePath * "basis65-" * runtype.runno * ".csv", reserveList, ',')
        next!(q) #ProgressMeter
    end

    #combined = table(policyData,reserveList)

    writedlm(filePath * "basis65.csv", reserveList, ',')
end

function real_main_prem()

    #filePath = "c:\\julialifeval\\lifeval\\input\\"
    #import LifeVal; joinpath(dirname(pathof(LifeVal)), "..", paths...)
    #fileName = filePath * "Basis20201283 - test1.xlsx"
    #fileName = filePath * "Basis20201285p.xlsx"

    runsettings, spcode, names = LifeValData.Initialise(filePath, fileName)
    policyData, basisData, fData = LifeValData.Read_Data(runsettings, spcode, names, filePath)

    #for premiums we require age at inception
    #policyData.age = policyData.age .- policyData.Term_if / 12

    PolicyType = [LifeValStructures.LifeAssurance, LifeValStructures.TermAssurance, LifeValStructures.Annuity, LifeValStructures.EndowmentAssurance, LifeValStructures.LifeAssuranceUL, LifeValStructures.EndowmentAssuranceUL]
    inputOrder = (:id, :sex, :age, :Term, :Term_if, :Premium, :Sum_Assured, :spc,
        :prem_esc, :sa_esc, :esc_m, :units_1, :units_99, :alloc_1, :alloc_99)

    premList = LifeValFunc.Table_Length(policyData, runsettings)

    o = nrow(runsettings) #ProgressMeter
    q = Progress(o, 1) #ProgressMeter

    @time for runtype in eachrow(runsettings)
        # eventually create struct called policyList
        policyList = LifeValFunc.Table_to_Policies(policyData[runtype.runno], :prod, inputOrder, PolicyType)

        #CALC UNIT PRICES, this must be set per runno and only required up to max projTerm
        #ulList = filter!(policy.spcode > 49, policyList)
        #not sure how to remove this hardcoding!
        glob_proj_max = 1000 #maximum(policyList[80001].projMax) + 1 #because permiums paid upfront, need 1 more at end
        unitpriceList = Array{Float64}(undef, nrow(runsettings), length(fData), 2, glob_proj_max)

        for fund in fData[runtype.runno,:]
            unitpriceList[runtype.runno, fund.datfram.ID, 1, 1] = fund.datfram.UP #.* (1 .+ fund.datfram.V_IR .- fund.datfram.V_IC) .^ (1/12)
            unitpriceList[runtype.runno, fund.datfram.ID, 2, 1] = fund.datfram.UP
            for t in 2:glob_proj_max
                unitpriceList[runtype.runno, fund.datfram.ID, 1, t] = unitpriceList[runtype.runno, fund.datfram.ID, 1, (t-1)]  .* (1 .+ fund.datfram.V_IR .- fund.datfram.V_IC) .^ (1/12)
                unitpriceList[runtype.runno, fund.datfram.ID, 2, t] = unitpriceList[runtype.runno, fund.datfram.ID, 2, (t-1)]  .* (1 .+ fund.datfram.P_IR .- fund.datfram.P_IC) .^ (1/12)
            end
        end

        premList[:, runtype.runno] = map(x -> LifeValReserve.Prem(x, (runtype, basisData, unitpriceList)...), policyList)

        #p = Progress(length(policyList), barglyphs=BarGlyphs("[=> ]"))
        #progress_map(policyList, progress=p) do x
        #    premList[:, runtype.runno] .= LifeValReserve.Prem(x, (runtype, basisData, unitpriceList)...)
        #end


        #outputtab = table(pol = policyList, res = reserveList)
        #writelm("out1.csv", outputtab, ",")
        #writedlm(filePath * "basis65-" * runtype.runno * ".csv", reserveList, ',')
        next!(q) #ProgressMeter
    end

    writedlm(filePath * "prem67.csv", premList, ',')
end

#POLICYLIST test environment
#policyTable = policyData[1]
#productIndex = :prod
#columnsIndex = inputOrder

#test environment
#policyList = Table_to_Policies(policyData[2], :prod, inputOrder)
#policy=policyList[1]
#runtype=runsettings[1,:]

function real_main_gen()

    #filePath = "c:\\julialifeval\\lifeval\\input\\"
    #fileName = filePath * "rand_input2.xlsx"
    inputs = DataFrame(XLSX.readtable(fileName, "randominput")...)
    totpolicies = []
    totprems = []
    @time for poltype in eachrow(inputs)
        # eventually create struct called policyList
        poltypevec = Vector(poltype)
        policyList, premList =  LifeValGenerate.Gen(poltypevec)
        #premList =  LifeValGenerate.PremGen(poltypevec)
        for k in 1:length(policyList[:,1])
            polvect=policyList[k,:]
            push!(totpolicies,polvect)
        end
        for k in 1:length(premList[:,1])
            premvect=premList[k,:]
            push!(totprems,premvect)
        end
    end

    #policyTable=table(totpolicies, names = [:id,	:prod,	:age,	:SumAssured,	:Term,	:Term_if,	:Premium,	:sex,
    # :spc, :prem_esc, :sa_esc,	:esc_m,	:units_1,	:units_99,	:alloc_1,	:alloc_99,	:proj_max])

    writedlm(filePath * "inputs65.csv", totpolicies, ',')
    writedlm(filePath * "prems65.csv", totprems, ',')
    #dataFile = filePath * "inputs65.csv"
    #policyData3 = loadtable(dataFile)
    #save(policyTable, filePath * "ulData6")
    #save(totpolicies, filePath * "ulData6")
end

# ---------------------------------------------------------------------------- #

function julia_main()#(gpu::Int64)
    gpu = 0
    filePath = datadir("exp_raw") * "\\"
    #filePath = "C:/JuliaLifeVal/LifeVal/input/"
    inpu= "Basis2020126.xlsx"
    fileName = filePath * inpu

    try
        if gpu == 1
            real_main_gpu()
        elseif gpu == 2
            real_main_prem()
        elseif gpu == 10
            real_main_gen()
        else
            real_main(filePath, fileName)
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

# print("settings: 0 = val (default), 1 = gpu, 2=pricing, 3 = multicore, 10 = generate \n\n")
# print("What do you want to do ? \n\n")
# inpu = readline()
#num = parse(Int64, inpu)
# num=0

# print("Please copy path of input files ? \n\n")
# inpu = readline()
# inpu = inpu * "\\"
# filePath = inpu

# filePath = "C:/JuliaLifeVal/LifeVal/input/"

# print("Please copy name of runsettings (5 char minimum) ? \n\n")
# inpu = readline()
# if inpu[end-4:end] != ".xlsx"   #appends xlsx
#     inpu = inpu * ".xlsx"
# end
# fileName = filePath * inpu

# inpu= "Basis2020126.xlsx"
# fileName = filePath * inpu

# julia_main(num)

end # module

#LifeVal.julia_main(num)

#temp converting data
# dataFile = "c:\\JuliaLifeVal\\lifeval\\input\\prems65.csv"
# dataFile = "c:\\JuliaLifeVal\\lifeval\\input\\inputs65.csv"
# dataFile = "c:\\JuliaLifeVal\\lifeval\\input\\inputs-test1.csv"
# policyData3 = loadtable(dataFile)
# save(policyData3, "c:\\JuliaLifeVal\\lifeval\\input\\policyData65")
# save(policyData3, "c:\\JuliaLifeVal\\lifeval\\input\\policyData2020126")
