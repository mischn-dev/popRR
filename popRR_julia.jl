# perform the snp calling using a variant caller that can create allelic depth values for each locus (bcftools version 1.10.2)
bcftools mpileup -Ou -q 20 -Q 25 -a FORMAT/AD  -f $REF $mydir/*filtered.sorted.bam | bcftools call -vm -Ov > $mydir/variants.vcf
# extract the alleic depth information from the raw vcf file -- only 2 alleles per locus allowed
vcftools --vcf /mnt/d/GeneticMapProject/4-SimReads/SNPcalling/SNPcalling_.vcf --out  /mnt/d/GeneticMapProject/4-SimReads/SNPcalling/outfile --max-alleles 2 --min-alleles 2 --minQ 40 --extract-FORMAT-info AD


using DataFrames, Statistics, CSV, StatsBase


function popRR(input, windowsize=10, Popsize=300)
                # load the file to the julia env - one genotype per file expected 
                snp = CSV.read(input, DataFrame; delim="\t")
                WS = windowsize * 1000000 # calculate the actual window size 
                # seperate the libraries in dicts
                Libs = Dict()
                # store the meta data 
                MPO = DataFrame(Chr = String[], POS = Int64[], Readcount = Int64[], Ref = Float64[], Alt = Float64[], DistReal = Float64[], LDposReal2 = Float64[], RR = Float64[], sample = String[])
                nam = names(snp) # get the names to loop over 

                ## 1. run the genetic map estimation 
                for I in 3:length(nam)

                                                    # the user might give either a vector of population sizes or a single value - if all samples have the same population size
                                                    if typeof(Popsize) == Int64
                                                                                popsize = Popsize
                                                    elseif  typeof(Popsize) == Vector{Int64}
                                                                                popsize = Popsize[I-2]
                                                    end

                                                    println("Calling RR for sample $(nam[I])")
                                                        inter = snp[.!(occursin.("0,0", snp[!,I])),[1,2,I]]
                                                        inter[!,:ReadcountRef] .= 0
                                                        inter[!,:ReadcountAlt] .= 0

                                                        # extract the Read count for each allele
                                                        Threads.@threads for rf in 1:size(inter,1)
                                                                                m = parse.(Int64, split(inter[rf,3], ","))
                                                                                inter[rf, 4:5] = m # push to ref and alt column
                                                        end

                                                        ## remove the raw column 
                                                        select!(inter, Not(3))

                                                        # calculate read count + relative allele frequency
                                                        inter[!,:Readcount] .= inter.ReadcountRef + inter.ReadcountAlt
                                                        inter.Ref .= inter.ReadcountRef ./ inter.Readcount
                                                        inter.Alt .= inter.ReadcountAlt ./ inter.Readcount

                                                        # remove the recount counts single columns
                                                        select!(inter, Not(3,4))

                                                        # remove rows with missing data 
                                                        inter = inter[.!(isnan.(inter.Ref)),:]

                                                        # calculate the genetic map, based on the AF and distance
                                                        inter[!,:DistReal] .= 0
                                                        inter[!,:AFvar] .= 0.0
                                                        for i in 2:size(inter,1)
                                                                            inter[i,:DistReal] = inter.POS[i] - inter.POS[i-1]
                                                                            inter[i,:AFvar] = abs(inter.Alt[i] - inter.Alt[i-1])
                                                           end

                                                        inter.DistReal = Float64.(inter.DistReal)
                                                        inter[inter.DistReal .<= 1,:DistReal] .= 1.1
                                                        # normalize the ditsance between the markers
                                                        inter.DistReal = log10.(inter.DistReal)
                                                        # split the chromosomes and calculate the D' for the markers
                                                        inter[!,:Dr] .= 0.0
                                                        inter[!,:LDposReal] .= 0.0

                                                        # split by chromosome
                                                        mp = Dict()
                                                                        # get the name of the first column to change it 
                                                        NAME = names(inter)[1]
                                                        rename!(inter, NAME => "Chr")
                                                        
                                                        for i in unique(inter.Chr); mp[i] = inter[inter.Chr .== i,:]; end

                                                        # calculate the genetic map from AF
                                                        for i in collect(keys(mp))
                                                                                mp[i][!,:Dr] .= (mp[i][!,:AFvar] .* (1 ./ mp[i][!,:DistReal]) )
                                                                                mp[i][1,:Dr] = 0.0 # first position does not have any AF variation, therefore set it to 0

                                                                            for j in 2:size(mp[i],1); mp[i][j,:LDposReal] = mp[i].Dr[j] + mp[i].LDposReal[j-1]; end # calculate the genetic position
                                                            end

                                                        # rejoin the dict in one DF
                                                        mpo = DataFrame(mp[collect(keys(mp))[1]][1,:])
                                                        for i in collect(keys(mp)); mpo = vcat(mpo, mp[i]); end
                                                        delete!(mpo, 1)

                                                        # sort the file by position - make sure the second column is named "POS"
                                                        NAME = names(inter)[2]
                                                        rename!(inter, NAME => "POS")

                                                        sort!(mpo, [:Chr, :POS])

                                                        # recalculate the map extention by adding the nlr regression

                                                        Marker = nrow(mpo) # marker count - we check the extend of the variant file and use this value 
                                                        marker = log2(Marker)
                                                        individuals = log2(popsize)

                                                        mpo.LDposReal2 = mpo.LDposReal .* (7958.9212 * exp(-0.5401 * marker) *  exp(0.3491 * individuals) + (691.0495 / sqrt(Marker)))

                                                        ## hypothesis - the sequencing introduces an error, that results in overestimation of genetic map length. For the first sample tested, it was 45 times too long. This could be assosiated with the read depth and the Allele frequency
                                                        # the overestimation depends on various factors, which are partly unknown.
                                                        # The user has the option to adjust the genetic map length to a commonly expected value - like exemplarily 300 cM 
                                                        corrector = 300 / maximum(mpo.LDposReal2) # needs to be saved somewhere

                                                        # remove intermediate columns
                                                        select!(mpo, Not(7,8,9))

                                                        ## calculate the cM/BP ratio 
                                                        mpo[!,:RR] .= 0.0
                                                        for k in 2:nrow(mpo);  
                                                                                            mpo.RR[k] = mpo.LDposReal2[k] - mpo.LDposReal2[k-1] # calcuate the cM part 
                                                                                            mpo[k,:DistReal] = mpo.POS[k] - mpo.POS[k-1] # calulcate the distance once more
                                                                                        end
                                                        # the distance in the beginning of the chromosome is false - set to 0 
                                                            mpo[mpo.DistReal .< 0,:DistReal] .= 0
                                                            mpo[mpo.RR .< 0, :RR] .= 0

                                                        mpo.RR .= mpo.RR ./ mpo.DistReal

                                                        # multiply with 1 mio to get the cM/MB 
                                                        mpo.RR *= 1000000

                                                        
                                                        ### generate the cM/MB values in windows
                                                        er = Dict(); for i in unique(mpo.Chr);  er[Symbol(i)] = mpo[mpo.Chr.==i,:]; end

                                                        stepper = WS
                                                        
                                                        for opm in collect(keys(er))

                                                                        i = Symbol(opm)
                                                                        asx = collect(range(1, maximum(er[i].POS), step=stepper)) # the range to serach in # to overlap, each step will be added 500mb to the bottom and 500 mb to the top end
                                                                        asx = DataFrame(hcat(asx, 1:length(asx)), :auto)
                                                                        asx[!,:start] .= asx.x1 .- stepper/2
                                                                        asx[!,:END] .= asx.x1 .+ stepper/2
                                                        
                                                                        er[i][!,:Group] .= ""
                                                                        er[Symbol(string("Groups_", String(i)))] = DataFrame(Group = Int64[], pos = Int64[], RRwindow = Float64[], Markercount = Int64[])                                                       
                                                        
                                                                        for itt in 1:nrow(asx)-1
                                                                                                                        run = string("m", asx[itt,:x2], "_")
                                                                                                                        start = asx[itt,:start]
                                                                                                                        steP = asx[itt+1, :END]
                                                                                                                        er[i][.&(er[i].POS .< steP, er[i].POS .> start),:Group] .= string.(er[i][.&(er[i].POS .< steP, er[i].POS .> start),:Group], run)
                                                                                                        end
                                                        
                                                                        for m in 1:nrow(asx)-1
                                                                                                                                    if any(occursin.("m$(m)_", er[i].Group))
                                                                                                                                    # do this for each replicate seperatly
                                                                                                                                    u = median(er[i][.&(occursin.("m$(m)_", er[i].Group)),8])
                                                                                                                                    POS = round(mean(er[i][.&(occursin.("m$(m)_", er[i].Group)), 2]))
                                                                                                                                    MK = nrow(er[i][.&(occursin.("m$(m)_", er[i].Group)),:])
                                                                                                                                    push!(er[Symbol(string("Groups_", String(i)))], [m,POS, u, MK])
                                                                                                                                    end
                                                                        end
                                                        end
                                                        
                                                        if typeof(inter.Chr[1]) == String7 ||  typeof(inter.Chr[1]) == String
                                                                                        sink = DataFrame(Group = Int64[], pos = Int64[], RRwindow = Float64[], Markercount = Int64[], Chr = String[])
                                                        elseif typeof(inter.Chr[1]) == Int64 
                                                                                        sink = DataFrame(Group = Int64[], pos = Int64[], RRwindow = Float64[], Markercount = Int64[], Chr = Int64[])
                                                        else
                                                            println("Your chromosomes seem to have an unexpected type. Please user either Strings (chr1) or integer (1) values")
                                                        end

                                                        for i in collect(keys(er))[occursin.("Groups", String.(collect(keys(er))))]
                                                            er[i][!,:Chr] .= String.(split(String(i),"_")[2])
                                                            append!(sink, er[i])
                                                        end
                                                        
                                                        sort!(sink, [:Chr, :pos])

                                                        filename = nam[I]

                                                        Libs[Symbol(filename)] = sink

                                                        # save meta files of mpo to a DF
                                                        mpo[!,:sample] .= filename 
                                                        # some inputs might give the chromosome as integer - adjust this, so that MPO can be created
                                                        if isempty(MPO) && typeof(mpo.Chr[1]) == Int64;
                                                                                                        MPO.Chr = Int64[]
                                                        end
                                                        # add the metadata of sample I to file                 
                                                        append!(MPO, mpo) 

                end

                ## merge the populations into a single file by the Group column 
                # generate a column to merge with 
                for samples in collect(keys(Libs))
                                                    Libs[samples][!,:Match] .= string.(Libs[samples].Chr, "_", Libs[samples].Group)
                end   
                # basic info to merge with 
                BASE1 = Libs[collect(keys(Libs))[1]][!,[6, 1, 5, 2]]
                BASE2 = Libs[collect(keys(Libs))[1]][!,[6, 1, 5, 2]]
                
                # merge the samples in one file
                for samples in collect(keys(Libs))
                                                    BASE1 = outerjoin(BASE1, Libs[samples][!,[6,3]], on=:Match)
                                                    BASE2 = outerjoin(BASE2, Libs[samples][!,[6,4]], on=:Match)
                                                    rename!(BASE1, "RRwindow" => String(samples))
                                                    rename!(BASE2, "Markercount" => String(samples))
                end   

                sort!(BASE1, [:Chr, :pos])
                sort!(BASE2, [:Chr, :pos])

                # push RR and Markercount into a new dict 
                out = Dict()
                out[:RecRate] = BASE1
                out[:Markercount] = BASE2 
                out[:SingleSNPinfo] = MPO
                
                return out
    end

# end of function


### example call style 
## input values 
 input = "D:/GeneticMapProject/4-SimReads/SNPcalling/outfile.AD.FORMAT"
 windowsize = 10 # in MB 
 Popsize = [300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50] # population size 

 result = popRR(input, windowsize, Popsize)

result = @time popRR(input, 50, Popsize)
