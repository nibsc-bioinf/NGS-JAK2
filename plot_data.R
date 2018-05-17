library(ggplot2)
library(reshape2)
library(egg)
library(data.table)
library(plyr)

# empty results table
results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(results) <- c('projid', 'sample', 'type', 'rep', 'value')

# project numbers
projects <- c("26", "93")

# read in data from TSV file
all <- fread("/data/output/tsv_files/results.tsv",showProgress = T)
#rename columns
names(all) <- c("project", "percentT", "run_no", "read_id", "orientation", "start", "hc_start", "sc_start", "sc_end", "hc_end", 
                "end", "allele", "ins_count", "ins_dist", "del_count", "del_dist", "mut_count", "snp_count")


# function for calculating T percent, T count, G count, and read count in sample
pGT <- function(reads){
  percentT<-c(sum(reads$allele=="T"))/nrow(reads)
  t<-c(sum(reads$allele=="T"))
  g<-c(sum(reads$allele=="G"))
  c(percentT,nrow(reads),t,g)
}

# function for getting results for each sample and project
get_results <- function(i, samp, projid){
  # Currently adjusts distance of mutation site from read end. This can be adjusted as required.
  b <- abs(all$sc_start)>=i & abs(all$sc_end)>=i & abs(all$percentT==samp) & abs(all$project==projid)
  # performs the pGT function for all reads, then for forward reads and then for reverse reads.
  data <- c(pGT(all[b]), pGT(all[b &  abs(all$orientation=='F')]), pGT(all[b & abs(all$orientation!='F')]))
  samp_results <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(results) <- c('projid', 'sample', 'type', 'rep', 'value')
  # loops through results and added result to a dataframe
  for (j in 1:12){
    result <- data.frame(projid, samp, j, i, data[j])
    colnames(result) <- c("projid", 'sample', 'type', 'rep', 'value')
    samp_results <- rbind(samp_results, result)
  }
  print(samp_results)
  samp_results
}

# loops through all data and gets results which are appended to the results table
for (proj in projects){
  if (proj == "26"){
    samps <- c("0", "0-1", "0-6", "1", "10", "90", "100")
  }
  if (proj == "93"){
    samps <- c("0", "0-1", "0-4", "0-7", "1", "10", "30", "90", "100")
  }
  for (samp in samps){
    print(proj)
    print(samp)
    # 100 loops which increased the minimum allowed distance of the mutation site to the read end
    for (i in 1:100){
      print(i)
      new_results <- get_results(i, samp, proj)
      results <- rbind(results, new_results)
    }
  }
}


# creating subsets of the data based on the information recorded, type numbers renamed to read orientation
percentT <- subset(results, type == 1 | type == 5 | type == 9)
percentT$type <- mapvalues(percentT$type, from = c(1,5,9), to = (c("All","Forward","Reverse")))
nReadsAll <- subset(results, type == 2 | type == 6 | type == 10)
nReadsAll$type <- mapvalues(nReadsAll$type, from = c(2,6,10), to = (c("All","Forward","Reverse")))
nReadsT <- subset(results, type == 3 | type == 7 | type == 11)
nReadsT$type <- mapvalues(nReadsT$type, from = c(3,7,11), to = (c("All","Forward","Reverse")))
nReadsG <- subset(results, type == 4 | type == 8 | type == 12)
nReadsG$type <- mapvalues(nReadsG$type, from = c(4,8,12), to = (c("All","Forward","Reverse")))

# plots produced (needs refactoring - each plot currently produced individually and then stacked together using ggarrange)
p1 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "0",], aes(x=rep, y = value * 100, colour = type, group = type)) + 
  geom_point(size=1) +
  geom_hline(yintercept=0, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank()) +
  labs(title = "0%", y = "Percent T allele")

p2 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "0-1",], aes(x=rep, y = value * 100, colour = type, group = type)) + 
  geom_point(size=1) +
  geom_hline(yintercept=0.03, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title = "0.1%")

p3 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "0-6",], aes(x=rep, y = value * 100, colour = type, group = type)) + 
  geom_point(size=1) +
  geom_hline(yintercept=0.6, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title = "0.6%")

p4 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "1",], aes(x=rep, y = value * 100, colour = type, group = type)) +
  geom_point(size=1) +
  geom_hline(yintercept=1, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title = "1%")

p5 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "10",], aes(x=rep, y = value * 100, colour = type, group = type)) +
  geom_point(size=1) +
  geom_hline(yintercept=10.8, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title = "10%")

p6 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "90",], aes(x=rep, y = value * 100, colour = type, group = type)) +
  geom_point(size=1) +
  geom_hline(yintercept=89.5, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title = "90%")

p7 <- ggplot(data=percentT[percentT$projid=="26" & percentT$sample == "100",], aes(x=rep, y = value * 100, colour = type, group = type)) +
  geom_point(size=1) +
  geom_hline(yintercept=100, linetype="dashed") + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(title = "100%")

p8 <- ggplot(data=nReadsT[nReadsT$projid=="26" & nReadsT$sample == "0",], aes(x=rep, y = value, colour = type, group = type)) +
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank()) +
  labs(y = "Minor Allele Count (T or G)")

p9 <- ggplot(data=nReadsT[nReadsT$projid=="26" & nReadsT$sample == "0-1",], aes(x=rep, y = value, colour = type, group = type)) +
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

p10 <- ggplot(data=nReadsT[nReadsT$projid=="26" & nReadsT$sample == "0-6",], aes(x=rep, y = value, colour = type, group = type)) +
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p11 <- ggplot(data=nReadsT[nReadsT$projid=="26" & nReadsT$sample == "1",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p12 <- ggplot(data=nReadsT[nReadsT$projid=="26" & nReadsT$sample == "10",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p13 <- ggplot(data=nReadsG[nReadsG$projid=="26" & nReadsG$sample == "90",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p14 <- ggplot(data=nReadsG[nReadsG$projid=="26" & nReadsG$sample == "100",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

p15 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "0",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank()) +
  labs(y = "Total read count")

p16 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "0-1",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p17 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "0-6",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p18 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "1",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.y=element_blank()) +
  labs(x = "Minimum distance from end of read")

p19 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "10",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p20 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "90",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

p21 <- ggplot(data=nReadsAll[nReadsAll$projid=="26" & nReadsAll$sample == "100",], aes(x=rep, y = value, colour = type, group = type)) + 
  geom_point(size=1) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# aligns and produces plots
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, ncol = 7)

### Run up to this point to produce plots. Code below is additional and can be used instead to build results also split by run.

#######################################################################################################################
# does all combinations for all runs for all samples.. very slow.


results <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(results) <- c('projid', 'sample', 'run', 'type', 'rep', 'value')

get_results <- function(i, samp, projid, run_no){
  b <- abs(all$sc_start)>=i & abs(all$sc_end)>=i & abs(all$percentT==samp) & abs(all$project==projid) & abs(all$run_no==run_no)
  data <- c(pGT(all[b]), pGT(all[b &  abs(all$orientation=='F')]), pGT(all[b & abs(all$orientation!='F')]))
  samp_results <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(results) <- c('projid', 'sample', 'run', 'type', 'rep', 'value')
  for (j in 1:12){
    result <- data.frame(projid, samp, run_no, j, i, data[j])
    colnames(result) <- c("projid", 'sample', 'run', 'type', 'rep', 'value')
    samp_results <- rbind(samp_results, result)
  }
  print(samp_results)
  samp_results
}


for (proj in projects){
  if (proj == "26"){
    samps <- c("0", "0-1", "0-6", "1", "10", "90", "100")
    runs <- 5
  }
  if (proj == "93"){
    samps <- c("0", "0-1", "0-4", "0-7", "1", "10", "90", "100")
  }
  for (samp in samps){
    print(proj)
    print(samp)
    if (proj == "93" & (samp == "0" | samp == "100")){
      runs <- 16
    }
    if (proj == "93" & (samp == "30" | samp == "90")){
      runs <- 4
    }
    if (proj == "93" & (samp == "0-1" | samp == "0-4" | samp == "0-7" | samp == "1" | samp == "10")){
      runs <- 12
    }
    for (run_no in 1:runs){
      for (i in 1:100){
        print(i)
        new_results <- get_results(i, samp, proj, run_no)
        results <- rbind(results, new_results)
      }
    } 
  }
}


