# load arguments from external source (RScript)
args <- commandArgs(trailingOnly = TRUE)

input = args[1] # path to file
winSize = args[2] # genomic window size to work on
Popsize = args[3] # path to text file with info

## load packages 
# data.table
if (any(installed.packages()[,1] == "data.table")) { 
                                                    library(data.table)
  
                                              }else{
                                                    install.packages("data.table")
                                                    library(data.table)
                                                    }

# stringr
if (any(installed.packages()[,1] == "data.table")) { 
                                                    library(stringr)
  
                                                  }else{
                                                    install.packages("stringr")
                                                    library(stringr)
                                                  }

### example call style 
## input values 
#input = "D:/GeneticMapProject/4-SimReads/SNPcalling/outfile.AD.FORMAT"
#windowsize = 50 # in MB 
#Popsize = c(300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50, 300, 100, 50) # population size 

# result = popRR(input, windowsize, Popsize)


popRR = function(input, windowsize=10, Popsize){
                                                
                                                  
                                                # load the file to env    
                                                snp = fread(input)
                                                WS = windowsize * 1000000 # calculate the actual window size 
                                                # seperate the libraries in dicts
                                                Libs = list()
                                                # store the meta data 
                                                MPO = data.frame(Chr = character(), POS = integer(), Readcount = integer(), Ref = numeric(), Alt = numeric(), DistReal = numeric(), LDposReal2 = numeric(), RR = numeric(), sample = character())
                                                nam = colnames(snp) # get the names to loop over 
                                                
                                                # run the genetic map and RR estimation for each sample seperatly
                        for (I in 3:length(nam)){
                                                
                                                # Popsize might be either a single integer value or a vector of integers - depending on user input 
                                                if (length(Popsize) > 1){
                                                                        popsize = Popsize[I-2]
                                                                 }else {
                                                                        popsize = Popsize[1]
                                                                 }
                                                
                                                print(paste("Calling RR for sample", nam[I]))  
                                                # convert data table to data frame
                                                snp = as.data.frame(snp)
                                                
                                                # remove missing calls for the sample 
                                                inter = snp[snp[,I] != "0,0",c(1,2,I)]
                                                inter$ReadcountRef = as.integer(str_split_fixed(inter[,3], ",", 2)[,1])
                                                inter$ReadcountAlt = as.integer(str_split_fixed(inter[,3], ",", 2)[,2])
                                                
                                                # remove column 3 
                                                inter = inter[,-3]
                                                
                                                # calculate read count + relative allele frequency
                                                inter$Readcount = inter$ReadcountRef + inter$ReadcountAlt
                                                inter$Ref = inter$ReadcountRef / inter$Readcount
                                                inter$Alt = inter$ReadcountAlt / inter$Readcount  
                                                
                                                # remove read count single columns
                                                inter = inter[,-c(3,4)]
                                                
                                                # rename the first column
                                                colnames(inter)[1] = "Chr"
                                                
                                                # split into list entries, each chromosome separated from each other
                                                mp = list()
                                                for (chr in unique(inter$Chr)){mp[[chr]] = inter[inter$Chr == chr,]}

                                  for (chr in names(mp)){
                                                                      
                                                                      # calculate the genetic map, based on the AF and distance
                                                                      mp[[chr]]$DistReal = 0
                                                                      mp[[chr]]$AFvar = 0 
                                                                      
                                                                      for (SNPs in 2:nrow(mp[[chr]])){
                                                                                                      mp[[chr]][SNPs,"DistReal"] = mp[[chr]]$POS[SNPs] - mp[[chr]]$POS[SNPs-1]
                                                                                                      mp[[chr]][SNPs,"AFvar"] = abs(mp[[chr]]$Alt[SNPs] - mp[[chr]]$Alt[SNPs-1])
                                                                      } # loop Distreal, afvar
                                                                      
                                                                      # transform the Distance to logDist 
                                                                      mp[[chr]][mp[[chr]]$DistReal <= 1, "DistReal"] = 1.1
                                                                      mp[[chr]]$logDist = log10(mp[[chr]]$DistReal)
                                                                      
                                                                      # pepare some additional columns
                                                                      mp[[chr]]$Dr = 0
                                                                      mp[[chr]]$LDposReal = 0
                                                                      
                                                                      # calculate the genetic map from AF 
                                                                      mp[[chr]]$Dr = mp[[chr]]$AFvar * (1 / mp[[chr]]$logDist)
                                                                      
                                                                      # sum up Dr along the chromosome to generate a "genetic map"
                                                                      for (gm in 2:nrow(mp[[chr]])){
                                                                                                    mp[[chr]][gm,"LDposReal"] = mp[[chr]]$Dr[gm] + mp[[chr]]$LDposReal[gm-1]
                                                                                                  }# genetic map generation loop (gm)
                                                
                                                      } # loop chr
                              
                                                
                                   # rejoin the mp list entries in one DF
                                  mpo = mp[[unique(inter$Chr)[1]]][1,]
                                  for (chr in unique(inter$Chr)){mpo = rbind(mpo, mp[[chr]])}
                                  
                                  # delete first row of mpo
                                  mpo = mpo[-1,]
                                                                  
                                  # recalculate the map extention by adding the nlr regression
                                  Marker = nrow(mpo) # marker count - we check the extend of the variant file and use this value 
                                  marker = log2(Marker)
                                  individuals = log2(popsize)
                                                                    
                                  mpo$LDposReal2 = mpo$LDposReal * (7958.9212 * exp(-0.5401 * marker) *  exp(0.3491 * individuals) + (691.0495 / sqrt(Marker)))
                                                                    
                                                                    
                                  ## calculate the cM/BP ratio 
                                  mpo$RR = 0.0
                                  for (k in 2:nrow(mpo)){
                                                              mpo$RR[k] = mpo$LDposReal2[k] - mpo$LDposReal2[k-1] # calcuate the cM part 
                                                                mpo[k,"DistReal"] = mpo$POS[k] - mpo$POS[k-1] # calulcate the distance once more
                                                        }
                                  
                                  # the distance in the beginning of the chromosome is false - set to 0 
                                  mpo[mpo$DistReal < 0,"DistReal"] = 0
                                  mpo[mpo$RR < 0, "RR"] = 0
                                                                    
                                  mpo$RR = mpo$RR / mpo$DistReal
                                                                    
                                  # multiply with 1 mio to get the cM/MB 
                                  mpo$RR = mpo$RR * 1000000
                                                                    
                                                    
                                  ### generate the cM/MB values in windows
                                  # split chr to seperate df
                                   er = list() 
                                   for (i in unique(mpo$Chr)){er[[i]] = mpo[mpo$Chr==i,]}
                                                                    
                                   stepper = WS # define the fwindow size to estimate the median RR for

             for (opm in names(er)){
        
                              i = opm
                              # generate the windows
                              asx = seq(1, max(er[[i]]$POS), by = stepper)
                              asx = as.data.frame(asx)
                              asx$x2 = c(1:length(asx$asx))
                              asx$start = asx$asx - stepper/2
                              asx$END = asx$asx + stepper/2
                              
                              # create a file where to store the window RR results in
                              er[[i]]$Group = ""
                              yxc = paste0("Groups_", i)
                              er[[yxc]] = data.frame(Group = integer(), pos = integer(), RRwindow = numeric(), Markercount = numeric())
                              
                              # add the group tag to each SNP
                              for (itt in c(1:nrow(asx)-1)){
                                                            run = paste0("m", asx[itt,"x2"], "_")
                                                            start = asx[itt, "start"]
                                                            steP = asx[itt+1, "END"]
                                                            er[[i]][er[[i]]$POS < steP & er[[i]]$POS > start, "Group"] = paste0(er[[i]][er[[i]]$POS < steP & er[[i]]$POS > start, "Group"], run)
                                } # itt loop - set the groups
          
                              for (m in 1:nrow(asx)-1){
                                
                                                        if (any(grep(paste0("m", m, "_"), er[[i]]$Group))){
                                                                                                            u = median(er[[i]][grepl(paste0("m", m, "_"), er[[i]]$Group),12])
                                                                                                            POS = round(mean(er[[i]][grepl(paste0("m", m, "_"), er[[i]]$Group),2]))
                                                                                                            MK = nrow(er[[i]][grepl(paste0("m", m, "_"), er[[i]]$Group),])
                                                                                                            # append the info to the new df storing the recombination rate in a genomic window
                                                                                                            er[[yxc]][nrow(er[[yxc]])+1,] = c(m,POS,u,MK)
                                                              } # if group present
                                
                                
                                            }# loop m 
                                        
                                      
                                      
                        } # opm loop 
                              
                              # save the results of all chromosomes in a single df
                              if (typeof(inter$Chr[1]) == "character"){
                                                                      sink = data.frame(Group = integer(), pos = integer(), RRwindow = numeric(), Markercount = numeric(), Chr = character())
                              }else if (typeof(inter$Chr[1]) == "double"){
                                                                      sink = data.frame(Group = integer(), pos = integer(), RRwindow = numeric(), Markercount = numeric(), Chr = integer())
                                
                              }else{print("Your chromosomes seem to have an unexpected type. Please user either Strings (chr1) or integer (1) values")}
                              
                            
                            ####       
                              for (i in  names(er)[grepl("Groups_", names(er))]){
                                                                                  er[[i]]$Chr = str_split_fixed(i, "_", 2)[,2]
                                                                                  sink = rbind(sink, er[[i]])
                                } # loop i - concat RR windows for chromosomes 
                                      
                              sink = sink[order(sink$Chr, sink$pos),]
                              
                              filename = nam[I] # sample name to store under
                              Libs[[filename]] = sink  # store the df in a list - transferd later to a new df
                              
                              # save intermediate files
                              mpo$Sample = filename
                              
                              # some inputs might give the chromosome as integer - adjust this, so that MPO can be created
                              if (nrow(MPO) == 0 && typeof(mpo$Chr[1]) == "integer"){MPO$Chr = integer()}
                              # add the metadata of sample I to file                 
                              MPO = rbind(MPO, mpo) 
                              
             } # loop I
                                
                              ## merge the populations into a single file by the Group column 
                              # generate a column to merge with 
                              for (samples in names(Libs)){
                                                                   Libs[[samples]]$Match = paste0(Libs[[samples]]$Chr, "_", Libs[[samples]]$Group)
                              }
                              # basic info to merge with 
                              BASE1 = Libs[[names(Libs)[1]]][,c(6,1,5,2)]
                              BASE2 = Libs[[names(Libs)[1]]][,c(6,1,5,2)]
                              
                              # merge the samples in one file
                              for (samples in names(Libs)){
                                                            
                                                          BASE1 = merge(BASE1, Libs[[samples]][,c(6,3)], by = "Match", all = T) 
                                                          BASE2 = merge(BASE2, Libs[[samples]][,c(6,4)], by = "Match", all = T) 
                                                          
                                                          colnames(BASE1)[colnames(BASE1) == "RRwindow"] = samples 
                                                          colnames(BASE2)[colnames(BASE2) == "RRwindow"] = samples 
                                                          
                                } # samples loop - merge to one file
                              
                              BASE1 = BASE1[order(BASE1$Chr, BASE1$pos),]
                              BASE2 = BASE2[order(BASE2$Chr, BASE2$pos),]
                              
                              # remove the match column
                              BASE1 = BASE1[,-1]
                              BASE2 = BASE2[,-1]
                              
                              # remove the distreal column
                              MPO = MPO[,-c(6,7,8,9,10)]
                              # rename the columns of MPO 
                              colnames(MPO)[6] = "GeneticMapPosition"
                              # round values 
                              MPO$Ref = round(MPO$Ref, digits=4)
                              MPO$Alt = round(MPO$Alt, digits=4)
                              MPO$GeneticMapPosition = round(MPO$GeneticMapPosition, digits=6)
                              MPO$RR = round(MPO$RR, digits=6)
                              
                              out = list()
                              out[["Recrate"]] = BASE1
                              out[["Markercount"]] = BASE2 
                              out[["SingleSNPinfo"]] = MPO  
                              
                              return(out)
} # function popRR

result = popRR(input, windowsize = winSize, Popsize = Popsize)


#write result to file 
loc = getwd()
w1 = paste(loc, "Recrate.txt", sep = "/")
w2 = paste(loc, "Markercount.txt", sep = "/")
w3 = paste(loc, "SingleSNPinfo.txt", sep = "/")

write.table(result$Recrate, w1, quote = F, row.names = F)
write.table(result$Markercount, w2, quote = F, row.names = F)
write.table(result$SingleSNPinfo, w3, quote = F, row.names = F)


