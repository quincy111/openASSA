module LifeValStructures

using DataFrames

export qxTables, wxTables, vtTables, expTables, chTables, fTables

export Policy, LifeAssurance, LifeAssuranceUL, TermAssurance, Annuity, EndowmentAssurance, EndowmentAssuranceUL

const MAX_AGE = 108

abstract type typeBasis end

struct qxTables <: typeBasis
    #id::string
    runno::Int64
    spcode::Int64
    datfram::DataFrame
    agemax_m::Int64
    agemax_f::Int64
    function qxTables(runno, spcode, datfram)
        #
        agemax_m=117
        for n in eachrow(datfram)
            if  n.QXM == 1
                agemax_m = n.AGE - 1
                break;
            end
        end
        #
        agemax_f=109
        for n in eachrow(datfram)
            if  n.QXF == 1
                agemax_f = n.AGE - 1
                break;
            end
        end
        return(new(runno, spcode, datfram, agemax_m, agemax_f))
    end
end

struct wxTables <: typeBasis
    #id::string
    runno::Int64
    spcode::Int64
    datfram::DataFrame
    function wxTables(runno, spcode, datfram)
        return(new(runno, spcode, datfram))
    end
end

struct vtTables <: typeBasis
    #id::string
    runno::Int64
    spcode::Int64
    datfram::DataFrame
    function vtTables(runno, spcode, datfram)
        return(new(runno, spcode, datfram))
    end
end

struct expTables <: typeBasis
    #id::string
    runno::Int64
    spcode::Int64
    datfram::DataFrame  #this should be a vector!
    function expTables(runno, spcode, datfram)
        return(new(runno, spcode, datfram))
    end
end

struct chTables <: typeBasis
    #id::string
    runno::Int64
    spcode::Int64
    datfram::DataFrame  #this should be a vector!
    function chTables(runno, spcode, datfram)
        return(new(runno, spcode, datfram))
    end
end

struct fTables <: typeBasis
    #id::string
    runno::Int64
    fundid::Int64
    datfram::DataFrame #this should be a vector!
    function fTables(runno, fundid, datfram)
        return(new(runno, fundid, datfram))
    end
end

# create Policy type
abstract type Policy end

# Life Assurance
# ---------------------------------------------------------------------------- #
struct LifeAssurance <: Policy
    id::Int64
    gender::Bool
    age::Int64
    term::Int64
    termIF::Int64
    premium::Float64
    SumAssured::Float64
    spc::Int64
    prem_esc::Float64
    sa_esc::Float64
    esc_m::Int64
    units_1::Float64
    units_99::Float64
    alloc_1::Float64
    alloc_99::Float64
    projMax::Int64


    function LifeAssurance(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99)
        # check for negative premium
        if premium < 0
            error("NegativePremium Error")
        end

        if age>MAX_AGE
            error("MaxAge Error")
        end

        global MAX_AGE
        projMax = 12(MAX_AGE-age)
        return(new(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99, projMax))
    end
end

struct LifeAssuranceUL <: Policy
    id::Int64
    gender::Bool
    age::Int64
    term::Int64
    termIF::Int64
    premium::Float64
    SumAssured::Float64
    spc::Int64
    prem_esc::Float64
    sa_esc::Float64
    esc_m::Int64
    units_1::Float64
    units_99::Float64
    alloc_1::Float64
    alloc_99::Float64
    projMax::Int64

    function LifeAssuranceUL(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99)
        # check for negative premium
        if premium < 0
            error("NegativePremium Error")
        end

        if age>MAX_AGE
            error("MaxAge Error")
        end

        global MAX_AGE
        projMax = 12(MAX_AGE-age)
        return(new(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99, projMax))
    end
end

struct Annuity <: Policy
    id::Int64
    gender::Bool
    age::Int64
    term::Int64
    termIF::Int64
    premium::Float64
    SumAssured::Float64
    spc::Int64
    prem_esc::Float64
    sa_esc::Float64
    esc_m::Int64
    units_1::Float64
    units_99::Float64
    alloc_1::Float64
    alloc_99::Float64
    projMax::Int64


    function Annuity(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99)
        # check for negative premium
        if premium != 0
            #error("Premium Error")
            #just set it to 0
            premium = 0
        end

        if prem_esc != 0
            #error("Premium Error")
            #just set it to 0
            prem_esc = 0
        end

        if age>MAX_AGE
            error("MaxAge Error")
        end

        global MAX_AGE
        projMax = 12(MAX_AGE-age)
        return(new(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99, projMax))
    end
end

struct TermAssurance <: Policy
    id::Int64
    gender::Bool
    age::Int64
    term::Int64
    termIF::Int64
    premium::Float64
    SumAssured::Float64
    spc::Int64
    prem_esc::Float64
    sa_esc::Float64
    esc_m::Int64
    units_1::Float64
    units_99::Float64
    alloc_1::Float64
    alloc_99::Float64
    projMax::Int64

    function TermAssurance(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99)
        # check for negative premium
        if premium < 0
            error("NegativePremium Error")
        end

        if age>MAX_AGE
            error("MaxAge Error")
        end

        global MAX_AGE
        projMax = term-termIF
        return(new(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99, projMax))
    end
end

struct EndowmentAssurance <: Policy
    id::Int64
    gender::Bool
    age::Int64
    term::Int64
    termIF::Int64
    premium::Float64
    SumAssured::Float64
    spc::Int64
    prem_esc::Float64
    sa_esc::Float64
    esc_m::Int64
    units_1::Float64
    units_99::Float64
    alloc_1::Float64
    alloc_99::Float64
    projMax::Int64


    function EndowmentAssurance(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99)
        # check for negative premium
        if premium < 0
            error("NegativePremium Error")
        end

        if age>MAX_AGE
            error("MaxAge Error")
        end

        global MAX_AGE
        projMax = term-termIF
        return(new(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, alloc_1, alloc_99, projMax))
    end
end

struct EndowmentAssuranceUL <: Policy
    id::Int64
    gender::Int64
    age::Int64
    term::Int64
    termIF::Int64
    premium::Float64
    SumAssured::Float64
    spc::Int64
    prem_esc::Float64
    sa_esc::Float64
    esc_m::Int64
    #units::Vector
    #allocPrem::Vector
    units_1::Float64
    units_99::Float64
    allocPrem_1::Float64
    allocPrem_99::Float64
    projMax::Int64


    function EndowmentAssuranceUL(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, allocPrem_1, allocPrem_99)
        # check for negative premium
        if premium < 0
            error("NegativePremium Error")
        end

        if age>MAX_AGE
            error("MaxAge Error")
        end

        global MAX_AGE
        projMax = term-termIF
        return(new(id, gender, age, term, termIF, premium, SumAssured, spc, prem_esc, sa_esc, esc_m, units_1, units_99, allocPrem_1, allocPrem_99, projMax))
    end
end

end # module
