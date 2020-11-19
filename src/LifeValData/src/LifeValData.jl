module LifeValData

using XLSX, DataFrames, JuliaDB, DelimitedFiles

using LifeValStructures

export Initialise, Read_Data

# set some global constants that need to be worked into functions at some point
function Initialise(filePath, fileName)
    main = DataFrame(XLSX.readtable(fileName, "MAIN")...)
    runsettings = DataFrame(XLSX.readtable(fileName, main.RS[1])...)
    spcode = DataFrame(XLSX.readtable(fileName, main.SP[1])...)

        #writes these as lastrun SETTINGS
    XLSX.writetable(filePath * "lastrun.xlsx", MAIN =( DataFrames._columns(main), DataFrames._names(main) ),
        TEMP1 =( DataFrames._columns(runsettings), DataFrames._names(runsettings) ),
        TEMP2 =( DataFrames._columns(spcode), DataFrames._names(spcode) ), overwrite=true)

    XLSX.openxlsx(filePath * "lastrun.xlsx", mode="rw") do xf
        sheet = xf[2]
        XLSX.rename!(sheet, main.RS[1])
        sheet = xf[3]
        XLSX.rename!(sheet, main.SP[1])
    end

    #read bases - this should happen before reading policies for e.g. impact of unit pricing on policies
    qxname = Vector{String}()
    wxname = Vector{String}()
    vtname = Vector{String}()
    expname = Vector{String}()
    chname = Vector{String}()
    fname = Vector{String}()
    fbname = Vector{String}()

    #qxLookup = unique(spcode, :qxbasis)
    for runtype in eachrow(runsettings)
        qxn = filePath * runtype.qxname
        wxn = filePath * runtype.wxname
        vtn = filePath * runtype.vtname
        expn = filePath * runtype.expname
        chn = filePath * runtype.chname
        fn = filePath * runtype.fname
        fbn = runtype.fbasis

        push!(qxname, qxn)
        push!(wxname, wxn)
        push!(vtname, vtn)
        push!(expname, expn)
        push!(chname, chn)
        push!(fname, fn)
        push!(fbname, fbn)
    end

    names = DataFrame([qxname, wxname, vtname, expname, chname, fname, fbname])
    rename!(names, [:qxname,:wxname,:vtname,:expname,:chname,:fname,:fbname])

    outpn = main.outv[1]

    return [runsettings, spcode, names, outpn]
end

function Read_Data(runsettings, spcode, names, filePath)
    #read policyData - []
    policyData = Vector{IndexedTable}(undef, nrow(runsettings))
    fB = DataFrame(XLSX.readtable(names.fname[1], names.fbname[1])...) #temporary, assume all f has same #s
    fData = Array{LifeValStructures.fTables}(undef, nrow(runsettings), nrow(fB))
    for runtype in eachrow(runsettings)
        policyData[runtype.runno] = load(filePath * runtype.policyfile)

        #fund data
        if runtype.fbasis != ""
            fB = DataFrame(XLSX.readtable(names.fname[runtype.runno], runtype.fbasis)...)
            for fbItem in eachrow(fB)
                fData[runtype.runno, fbItem.ID]=LifeValStructures.fTables(runtype.runno, fbItem.ID, DataFrame(fbItem))
            end
        end
    end

    #find all (unique) spcodes
    spcodes = []
    for runno in eachrow(spcode)
        spcod = runno.spc
        #if spcod !in spcodes
            push!(spcodes, spcod)
        #end
    end
    spcodes = unique(spcodes)

    qxData = Array{LifeValStructures.qxTables}(undef, nrow(runsettings), 100)
    qxpData = Array{LifeValStructures.qxTables}(undef, nrow(runsettings), 100)
    wxData = Array{LifeValStructures.wxTables}(undef, nrow(runsettings), 100)
    wxpData = Array{LifeValStructures.wxTables}(undef, nrow(runsettings), 100)
    vtData = Array{LifeValStructures.vtTables}(undef, nrow(runsettings), 100)
    vtpData = Array{LifeValStructures.vtTables}(undef, nrow(runsettings), 100)
    expData = Array{LifeValStructures.expTables}(undef, nrow(runsettings), 100)
    exppData = Array{LifeValStructures.expTables}(undef, nrow(runsettings), 100)
    chData = Array{LifeValStructures.chTables}(undef, nrow(runsettings), 100)

    for runno in eachrow(spcode)
        runtemp = runno.runno
        spcod = runno.spc

        qxB = DataFrame(XLSX.readtable(names.qxname[runtemp], runno.qxbasis)...)
        qxData[runtemp,spcod]=LifeValStructures.qxTables(runtemp,spcod,qxB)

        if runno.qxbasisp != ""
            qxpB = DataFrame(XLSX.readtable(names.qxname[runtemp], runno.qxbasisp)...)
            qxpData[runtemp,spcod]=LifeValStructures.qxTables(runtemp,spcod,qxpB)
        end

        wxB = DataFrame(XLSX.readtable(names.wxname[runtemp], runno.wxbasis)...)
        wxData[runtemp,spcod]=LifeValStructures.wxTables(runtemp,spcod,wxB)

        if runno.wxbasisp != ""
            wxpB = DataFrame(XLSX.readtable(names.wxname[runtemp], runno.wxbasisp)...)
            wxpData[runtemp,spcod]=LifeValStructures.wxTables(runtemp,spcod,wxpB)
        end

        vtB = DataFrame(XLSX.readtable(names.vtname[runtemp], runno.vtbasis)...)
        vtData[runtemp,spcod]=LifeValStructures.vtTables(runtemp,spcod,vtB)

        if runno.vtbasisp != ""
            vtpB = DataFrame(XLSX.readtable(names.vtname[runtemp], runno.vtbasisp)...)
            vtpData[runtemp,spcod]=LifeValStructures.vtTables(runtemp,spcod,vtpB)
        end

        expB = DataFrame(XLSX.readtable(names.expname[runtemp], runno.ebasis)...)
        expItem = filter(row -> row.SPCODE == runno.spc, expB) #don't need every spcode on spcode level
        expData[runtemp,spcod]=LifeValStructures.expTables(runtemp,spcod,expItem)

        if runno.ebasisp != ""
            exppB = DataFrame(XLSX.readtable(names.expname[runtemp], runno.ebasisp)...)
            exppItem = filter(row -> row.SPCODE == runno.spc, exppB) #don't need every spcode on spcode level
            exppData[runtemp,spcod]=LifeValStructures.expTables(runtemp,spcod,exppItem)
        end

        if runno.chbasis != ""
            chB = DataFrame(XLSX.readtable(names.chname[runtemp], runno.chbasis)...)
            chData[runtemp,spcod]=LifeValStructures.chTables(runtemp,spcod,chB)
        end
    end

    basisData = [qxData, qxpData, wxData, wxpData, vtData, vtpData, expData, exppData, chData]
    #names!(basisData, [:qxData,:qxpData,:wxData,:wxpData,:vtData,:vtpData,:expData,:exppData,:chData,:fData])

    return [policyData, basisData, fData]
end

end
