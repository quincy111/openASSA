module LifeVal

#using Distributed

#@everywhere begin
    using DrWatson
    @quickactivate "LifeVal"
#end

using JuliaDB, DelimitedFiles, DataFrames, XLSX, ProgressMeter #, Traceur
using Random, Distributions #not technically required, but added here so that all modules use same project

using Pkg, Pkg.Artifacts
#using CUDA, Distributed    #wip

using LifeValStructures
using LifeValData
using LifeValFunc
using LifeValReserve
using LifeValGenerate

function real_main()
    runsettings, spcode, names, outfile = LifeValData.Initialise(filePath, fileName)
    policyData, basisData, fData = LifeValData.Read_Data(runsettings, spcode, names, filePath)

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

    writedlm(filePath * outfile, reserveList, ',')
end

#THE GPU ENVIRONMENT IS NOT ACTIVATED AT THE MOMENT DUE TO CuArray ERROR (bits??)
function real_main_gpu()

    runsettings, spcode, names, outfile = LifeValData.Initialise(filePath, fileName)
    policyData, basisData, fData = LifeValData.Read_Data(runsettings, spcode, names, filePath)

    reserveList = LifeValFunc.Table_Length(policyData, runsettings)
    policyList1 = CuArray(policyData)  #this is wip

    o = nrow(runsettings) #ProgressMeter
    q = Progress(o, 1) #ProgressMeter

    @time for runtype in eachrow(runsettings)
        # eventually create struct called policyList
        policyList = LifeValFunc.Table_to_Policies(policyData[runtype.runno], :prod, inputOrder, PolicyType)
        policyList = CuArray(policyList)  #this is an ERROR???

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

    writedlm(filePath * outfile, reserveList, ',')
end

#MAKE SURE TO USE THE PREMS OUTPUT FOR PREMIUMS!
function real_main_prem()
    runsettings, spcode, names, outfile = LifeValData.Initialise(filePath, fileName)
    policyData, basisData, fData = LifeValData.Read_Data(runsettings, spcode, names, filePath)

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

        #INSTEAD OF THE MAP ABOVE, THE FOLLOWING CAN ALSO BE USED, 2X SLOWER...
        #p = Progress(length(policyList), barglyphs=BarGlyphs("[=> ]"))
        #progress_map(policyList, progress=p) do x
        #    premList[:, runtype.runno] .= LifeValReserve.Prem(x, (runtype, basisData, unitpriceList)...)
        #end

        #outputtab = table(pol = policyList, res = reserveList)
        #writelm("out1.csv", outputtab, ",")
        #writedlm(filePath * "basis65-" * runtype.runno * ".csv", reserveList, ',')
        next!(q) #ProgressMeter
    end

    writedlm(filePath * outfile, premList, ',')
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
    inputs = DataFrame(XLSX.readtable(fileName, "randominput")...)
    totpolicies = []
    totprems = []
    @time for poltype in eachrow(inputs)
        poltypevec = Vector(poltype)
        policyList, premList =  LifeValGenerate.Gen(poltypevec)
        for k in 1:length(policyList[:,1])
            polvect=policyList[k,:]
            push!(totpolicies,polvect)
        end
        for k in 1:length(premList[:,1])
            premvect=premList[k,:]
            push!(totprems,premvect)
        end
    end

    writedlm(filePath * "inputs69.csv", totpolicies, ',')
    writedlm(filePath * "prems69.csv", totprems, ',')
end

# ---------------------------------------------------------------------------- #

function julia_main(gpu::Int64)
    try
        if gpu == 1 #using CUDA
            real_main_gpu()
        elseif gpu == 3  #using Distributed (NOT ENABLED YET)
            real_main_cpu()
        elseif gpu == 4  #e.g. using CUDA and Distributed (NOT ENABLED YET)
            real_main_cpu()
        elseif gpu == 10  #generating policies
            real_main_gen()
        elseif gpu == 20
            real_main_prem()
        else   #using Threads
            real_main()
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

#THE FOLLOWING IS TESTED IN PROD ENVIRONMENT, NOT RUN IN TEST ENVIRONMENT
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

#fooifier_path() = joinpath(artifact"fooifier", "bin", "fooifier" * (Sys.iswindows() ? ".exe" : ""))

#FIXED VALUES FOR ALL RUNS / CURRENT SETTINGS FOR ALL
PolicyType = [LifeValStructures.LifeAssurance, LifeValStructures.TermAssurance, LifeValStructures.Annuity, LifeValStructures.EndowmentAssurance, LifeValStructures.LifeAssuranceUL, LifeValStructures.EndowmentAssuranceUL]
inputOrder = (:id, :sex, :age, :Term, :Term_if, :Premium, :Sum_Assured, :spc,
    :prem_esc, :sa_esc, :esc_m, :units_1, :units_99, :alloc_1, :alloc_99)

filePath = datadir("exp_raw") * "\\"   #current fixed input/output directory

print("Do you want to repeat last run settings (y/n) ? \n\n")
inpu = readline()
if inpu == "y" || inpu == "Y" #note not working for GEN currently
    fileName = filePath * "lastrun.xlsx"
    print("settings: 0 = val (default), 1 = gpu, 2-9 = multicore/+, 10 = generate, 20=pricing \n\n")
    print("What do you want to do ? \n\n")
    inpu = readline()
    num = parse(Int64, inpu)
else
    print("settings: 0 = val (default), 1 = gpu, 2-9 = multicore/+, 10 = generate, 20=pricing \n\n")
    print("What do you want to do ? \n\n")
    inpu = readline()
    num = parse(Int64, inpu)

    # print("Please copy path of input files ? \n\n") #ENABLE FOR CUSTOM DIR
    # inpu = readline()
    # filePath = inpu * "\\"

    print("Please copy name of runsettings (5 char minimum) ? \n\n")
    inpu = readline()
    if inpu[end-4:end] != ".xlsx"   #appends xlsx
        inpu = inpu * ".xlsx"
    end
    fileName = filePath * inpu
end

julia_main(num)

end # module

#LifeVal.julia_main(num)

#temp converting data (again, just lazy to go and find it each time)
 # filePath = datadir("exp_raw") * "\\"
 # dataFile = filePath * "inputsA1.csv"
 # policyDataT = loadtable(dataFile)
 # JuliaDB.save(policyDataT, filePath * "policyDataA")
