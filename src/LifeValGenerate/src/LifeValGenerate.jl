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

function gen_beta(alphav, betav, minv, maxv, fact)
    x = rand(Beta(alphav,betav),1)
    y=x[1]*fact
    val = round(y, digits=-3)
    if val < minv
        val=minv
    elseif val > maxv
        val=maxv
    end
    return (val)
end

#the Gen function below generates 2 sets of output, the IF book and the at inception book
#code is crudely written!  our appologies, to be revised but not top of list!
function Gen(v::Vector)
    policies = Array{Int64}(undef, v[4], 17)
    prems = Array{Int64}(undef, v[4], 17)

    for n in 1:v[4]
        #generate from Normal and Beta distributions!
        age = gen_normal(v[6],v[7],v[8],v[9],1)
        sa = gen_beta(v[11],v[12],v[13],v[14],v[23])
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
