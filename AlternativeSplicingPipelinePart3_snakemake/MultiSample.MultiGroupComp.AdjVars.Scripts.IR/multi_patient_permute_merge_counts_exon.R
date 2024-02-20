## calculate log odds ratio for all clusters with >=2 junctions per cluster and total reads in genotyped cells >= 5
## permute by shuffling genotype within each patient 
## calculate one odds ratio by pseudobulking across all patients 
## adds in regularization factor of 1e-5 to all junctions 
## if cluster contains a junction that is NA, makes the values = 0 otherwise leads to NA later on 

library(Matrix)
library(tidyverse)
library(matrixStats)

args = commandArgs(TRUE)
path.to.split = args[1]
path.to.cell.meta = args[2]
comp.groups.column = args[3] #e.g. genotype labels
cell.groups.column = args[4] #e.g. column with cell type labels 
nperm = args[5]
sample.names = args[6]
output.dir = args[7]
output.file = args[8]
group1_name = args[9]
group2_name = args[10]

print(args)

nperm = as.numeric(nperm)

path.to.three.matrix = paste(path.to.split, "/counts_files", sep = "")
path.to.three.data = paste(path.to.split, "/data_tables", sep = "")

sample.names = unlist(strsplit(sample.names, split = ","))
print(sample.names)

#read in matrix as named list of matrix
setwd(path.to.three.matrix)
file.paths = list.files(path.to.three.matrix)
suffix <- lapply(strsplit(file.paths[1], split= "mtx"), "[",2)
files <- list()

for (sample in sample.names ){
        matrix.path <- paste0(sample,".mtx", suffix)
        files <- c(files , matrix.path)
}

files <- unlist(files)
three.mtx.list = lapply(files, function(x) read.table(file=x, header = TRUE))
print(files)
print(sample.names)
names(three.mtx.list) = sample.names
print(names(three.mtx.list))



message("Matrix loaded")

#load cell metadata information (should be save the version from the previous script that performs the same action?)
##Then we won't have to include these arguments again??

cell.meta.list = list()
setwd(path.to.cell.meta)
i = 1
for (file in list.files(path.to.cell.meta)){
  
  cell.meta = as.data.frame(read.table(file, stringsAsFactors = F))
  cell.meta = cell.meta[!is.na(cell.meta[,cell.groups.column]),]
  cell.meta = cell.meta[,c(cell.groups.column, comp.groups.column)]
  print("Removing Pattern on cell barcode")
  rownames(cell.meta) = unlist(lapply(strsplit(rownames(cell.meta), split="_"), "[",1))
  ##You will have the same cell barcodes in the colnames of the three prime and five prime matircies
  sum(rownames(cell.meta) %in% colnames(three.mtx.list[[i]]))
  cell.meta = cell.meta[rownames(cell.meta) %in% colnames(three.mtx.list[[i]]),]
  ##keeps only the cell belonging to the 2 groups you're comparing
  cell.meta.sub = cell.meta[!is.na(cell.meta[,comp.groups.column]),]
  
  cell.meta.list[[sample.names[i]]] = cell.meta.sub
  i = i+1
}

message("Cell meta data loaded")

setwd(path.to.three.data)
files = list.files(path.to.three.data)
##We created 2 files but we are only reading one ????? remove this 
three.data.comb = read.csv(files[1])



message("All Data Loaded")


for (sample in sample.names){
  three.mtx.list[[sample]]$sample = sample
}




three.prime.clusters = as.character(unique(three.data.comb$five_prime_ID))


three.obs.ratio = list()

for(cluster in three.prime.clusters){
  subset = as.data.frame(three.data.comb[which(three.data.comb$five_prime_ID == cluster),])
  subset$group1.all = sum(subset$obs.group1)
  subset$group2.all = sum(subset$obs.group2)
  subset$group1.remain = subset$group1.all - subset$obs.group1
  subset$group2.remain = subset$group2.all - subset$obs.group2
  for (junc in subset$exon_coordinates){
    three.obs.ratio[junc] = log(((subset[which(subset$exon_coordinates == junc),"obs.group1"] + 1e-5)/(subset[which(subset$exon_coordinates == junc),"group1.remain"] + 1e-5))/(((subset[which(subset$exon_coordinates == junc),"obs.group2"]+1e-5)/(subset[which(subset$exon_coordinates == junc),"group2.remain"]+1e-5))))
  }
}

three.obs.ratio.num = as.numeric(three.obs.ratio)
names(three.obs.ratio.num) = names(three.obs.ratio)


three.data.comb$alt_three_prime_exon_coordinates = paste(three.data.comb$exon_coordinates, three.data.comb$five_prime_ID, sep = ":")

#create two data frames with final output
three.sample.output = data.frame(three.obs.logOR.ratio = three.obs.ratio.num, exon_coordinates = names(three.obs.ratio))
three.sample.output = left_join(three.sample.output, three.data.comb, by = "exon_coordinates")



message("Observed difference calculated")

#initialize to calculate whether or not obs OR > shf OR
temp.three = rep(0,length(three.obs.ratio.num)) # Create an empty vector to update in each iteration

#initiate shuffled data frame, do this only once
three.shf.data.list = list()
#five.shf.data.list = list()
for (sample in sample.names){
  three.shf.data.list[[sample]] = data.frame(exon_coordinates = three.data.comb$exon_coordinates, five_prime_ID=three.data.comb$five_prime_ID)
}


for(x in 0:nperm){

  set.seed(x)

  #cycle through each sample and shuffle genotypes before calculating OR
  #use only filtered data frames based on previous filtering done
  
  for (sample in sample.names) {
    
    #Looks at all cells in metadata table when shuffling annotations. 
    cell.annotation_vec = as.character(cell.meta.list[[sample]][,comp.groups.column])
    names(cell.annotation_vec) = rownames(cell.meta.list[[sample]])
    orig.names = names(cell.annotation_vec)
    shf.annotation = sample(cell.annotation_vec, size = length(cell.annotation_vec), replace = F)
    names(shf.annotation) = orig.names

    group2 = names(shf.annotation)[shf.annotation == group2_name]
    group1 = names(shf.annotation)[shf.annotation == group1_name]

    three.shf.data = three.shf.data.list[[sample]]

    ##Only grabbing counts from the per patient count matrix
    three.shf.data$shf.group2 = rowSums(three.mtx.list[[sample]][,colnames(three.mtx.list[[sample]]) %in% group2])
    three.shf.data$shf.group1 = rowSums(three.mtx.list[[sample]][,colnames(three.mtx.list[[sample]]) %in% group1])
    three.shf.data$total.reads.per.junction = three.shf.data$shf.group2+three.shf.data$shf.group1
    three.shf.data.list[[sample]] = three.shf.data
    
  }
  
  for (sample in sample.names){
    three.shf.data.list[[sample]]$sample = sample
    #five.shf.data.list[[sample]]$sample = sample
  }
  
  ##Binding the results from the per patient shuffled counts so that we can collapse them 
  three.shf.data = bind_rows(three.shf.data.list)
  
  three.shf.data.comb = three.shf.data %>% group_by(exon_coordinates, five_prime_ID) %>% summarise(shf.group2 = sum(shf.group2), shf.group1 = sum(shf.group1), total.reads.per.junction = sum(total.reads.per.junction))
  
  three.shf.data.comb[is.na(three.shf.data.comb)] = 0
  
  three.shf.ratio = list()
  
  #log OR for each junction after shuffling genotypes within each sample
  for(clust in three.prime.clusters){
    subset = as.data.frame(three.shf.data.comb[which(three.shf.data.comb$five_prime_ID == clust),])
    subset$group1.all = sum(subset$shf.group1)
    subset$group2.all = sum(subset$shf.group2)
    subset$group1.remain = subset$group1.all - subset$shf.group1
    subset$group2.remain = subset$group2.all - subset$shf.group2
    for (junc in subset$exon_coordinates){
      three.shf.ratio[junc] = log(((subset[which(subset$exon_coordinates == junc),"shf.group1"] + 1e-5)/(subset[which(subset$exon_coordinates == junc),"group1.remain"] + 1e-5))/(((subset[which(subset$exon_coordinates == junc),"shf.group2"]+1e-5)/(subset[which(subset$exon_coordinates == junc),"group2.remain"]+1e-5))))
    }
  }

  three.shf.ratio.num = as.numeric(three.shf.ratio)
  names(three.shf.ratio.num) = names(three.shf.ratio)
  
  
  #create output data frame with results for each sample
  three.shf.sample.output = data.frame(three.shf.logOR.ratio = three.shf.ratio.num, alt_three_prime_exon_coordinates = names(three.shf.ratio))
  


  three.shf.diff = abs(three.obs.ratio.num) > abs(three.shf.ratio.num)

  temp.three = temp.three + three.shf.diff

  if(x == nperm) cat("Permuted differences calculated")
}

pvals.three = 1 - temp.three/(nperm + 1)

message("Creating final data frames")
final.three = data.frame(pvalue = pvals.three, three.obs.logOR.ratio = three.obs.ratio.num, exon_coordinates = names(three.obs.ratio.num))
final.three = left_join(final.three,three.sample.output)

three.cluster.cov = final.three %>% group_by(five_prime_ID) %>% summarise(three.group1.cluster.cov = sum(obs.group1), three.group2.cluster.cov = sum(obs.group2))
final.three = left_join(final.three, three.cluster.cov, by= "five_prime_ID")


message("Done with permutations!")

message("writing output")
setwd(output.dir)
three.filename = paste("./alt_three_prime/", output.file, ".csv", sep = "")
write.csv(final.three, file = three.filename, quote = FALSE, row.names = FALSE)

message("Done!!")

