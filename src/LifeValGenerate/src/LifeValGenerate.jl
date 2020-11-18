module LifeValGenerate

using Random, Distributions

export Gen

Random.seed!(2020)

function gen_normal(meanv, sdv, minv, maxv, fact)
    if meanv==99
        val = 0
    else
        x = rand(Normal(meanv,sdv),1)
        y=x[1]
        val = round(y, digits=0)
        if val < minv
            val=minv
        elseif val > maxv
            val=maxv
        end
    end
    return (val * fact)
end

function Gen(v::Vector)
    policies = Array{Int64}(undef, v[4], 17)
    prems = Array{Int64}(undef, v[4], 17)

    for n in 1:v[4]
        age = gen_normal(v[6],v[7],v[8],v[9],1)
        sa = gen_normal(v[11],v[12],v[13],v[14],1000)
        term =gen_normal(v[16],v[17],v[18],v[19],1)

        if v[19] != 99  #term_if, ceiling if term applicable
            term_if =gen_normal(v[20],v[21],v[22],term-1,1) #with ceiling
        else
            term_if =gen_normal(v[20],v[21],v[22],99,1) #no ceiling
        end

        #alloc to policies in force
        policies[n,1]=n   #id
        policies[n,2]=v[1] #prod
        policies[n,3]=age #age
        policies[n,4]=sa #sumassured
        policies[n,5]=term #term
        policies[n,6]=term_if
        #current age...
        policies[n,3]=policies[n,3]+policies[n,6] #age is current age for reserve calc
        #terms are in months!
        policies[n,5]=term*12  #term converted to months
        policies[n,6]=term_if*12  #term_if converted to months
        policies[n,7]=0
        policies[n,8]=rand(0:1)  #gender, 0 males otherwise 1 females
        policies[n,9]=v[2] #spcode
        policies[n,10]=0
        policies[n,11]=0
        policies[n,12]=rand(1:12) #prem-inc month
        policies[n,13]=0
        policies[n,14]=0
        policies[n,15]=0
        policies[n,16]=0
        policies[n,17]=0

        prems[n,1]=n   #id
        prems[n,2]=v[1] #prod
        prems[n,3]=age #age
        prems[n,4]=sa #sumassured
        prems[n,5]=term*12 #term
        prems[n,6]=0
        prems[n,7]=0
        prems[n,8]=rand(0:1)  #gender, 0 males otherwise 1 females
        prems[n,9]=v[2] #spcode
        prems[n,10]=0
        prems[n,11]=0
        prems[n,12]=rand(1:12) #prem-inc month
        prems[n,13]=0
        prems[n,14]=0
        prems[n,15]=0
        prems[n,16]=0
        prems[n,17]=0
    end

    return [policies, prems]
end

end # module

#d = Normal(50,10)

#x = rand(d, 100)

#quantile.(Normal(), [0.5, 0.95])
#pdf.(Normal(),[0.5, 0.95])
#cdf.(Normal(), [0.5, 0.95])

#Binomial(p) # Discrete univariate
#Cauchy(u, b)  # Continuous univariate
#Multinomial(n, p) # Discrete multivariate
#Wishart(nu, S)  # Continuous matrix-variate

#fit(Normal, x)

#test environment
#Pkg.add("DataFrames")
#Pkg.add("XLSX")
#Pkg.add("Random")
#Pkg.add("Distributions")
#using DataFrames, XLSX
#using Random, Distributions

#filePath = "c:\\julialifeval\\lifeval\\input\\"
#fileName = filePath * "rand_input2.xlsx"
#inputs = DataFrame(XLSX.readtable(fileName, "randominput")...)
#totprems = []
#totprems = []
#poltype=inputs[1,:]

#for poltype in eachrow(inputs)
    # eventually create struct called policyList
    #poltypevec = Vector(poltype)
    #policyList, premList =  Gen(poltypevec)
    #premList =  LifeValGenerate.PremGen(poltypevec)
    #for k in 1:length(policyList[:,1])
    #    polvect=policyList[k,:]
    #    push!(totprems,polvect)
    #end
    #for k in 1:length(premList[:,1])
    #    premvect=premList[k,:]
    #    push!(totprems,premvect)
    #end
#end

#writedlm(filePath * "inputs65.csv", totprems, ',')
#writedlm(filePath * "prems65.csv", totprems, ',')
