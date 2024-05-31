dir1 = "~/LC716/bam.sc/"      #### Path to bam/sam file ######
dir2 = "~/LC716/chimeric.sc.pairs.mapq40/"     #### Path to save all intermediate files ######
mm10.chrom.size = "~/mm10.chrom.sizes"   ######## Path to mm10 chrom sizes file ########
mm10.blacklist.region = "~/Droplet.FIREscore/mm10.blacklist.bed" ####### Path to mm10 Black list regions #########

files = list.files(path=dir1 ,pattern=".int.sam", all.files=FALSE, full.names=FALSE)

for (file in files){
  modified_string <- sub("\\.int.sam$", "", file)
  command <- paste0("pairtools parse2 -c ",mm10.chrom.size," --assembly mm10 --add-pair-index --expand --drop-sam --min-mapq 40 ",dir1, file," -o ",dir2, modified_string,".pairs")
  print(command)
  system(command)
  
}


files = list.files(path=dir2 ,pattern=".pairs", all.files=FALSE, full.names=FALSE)

for (file in files){
  
  modified_string <- sub("\\.pairs$", "", file)
  
  command <- paste0("pairtools sort -o ",dir2,modified_string,"_sorted.pairs ",dir2, modified_string,".pairs")
  print(command)
  system(command)
  
  command <- paste0("pairtools dedup --mark-dups --output-stats ",dir2, modified_string,".dedup.stats -o ",dir2, modified_string,"_dedup.pairs ",dir2, modified_string,"_sorted.pairs")
  print(command)
  system(command)
  
  command <- paste0("pairtools select '(pair_type==\"UU\") and ((chrom1!=chrom2) or ((chrom1==chrom2) and (abs(pos1 - pos2) > 10000)))' ",dir2, modified_string,"_dedup.pairs > ",dir2, modified_string,".Noselfligation.UU.pairs")
  print(command)
  system(command)
  
}


files = list.files(path=dir2, pattern=".Noselfligation.UU.pairs", all.files=FALSE, full.names=FALSE)
for (file in files){
  modified_string <- sub("\\.Noselfligation.UU.pairs$", "", file)
  
  a<-read.table(paste0(dir2, file), head=F)
  a<-as.data.frame(a)
  a$ID1 <- 0
  a$ID2 <- 0
  
  x<-read.table(mm10.blacklist.region, head=F)
  x<-as.data.frame(x)
  colnames(x)<-c('chr','start','end')
  nrow(x) # 164 blacklist
  sum(x$end-x$start)/1e3 # 81.46Kb in total
  x<-x[order(x$chr, x$start),]
  
  for(i in 1:nrow(x))
  {
    index1 <- I(a[,2]==x$chr[i] & a[,3]>=x$start[i] & a[,3]<=x$end[i])
    index2 <- I(a[,4]==x$chr[i] & a[,5]>=x$start[i] & a[,5]<=x$end[i])
    if(sum(index1)>=1)
    {
      a[index1==1, ]$ID1 <- 1
    }
    if(sum(index2)>=1)
    {
      a[index2==1, ]$ID2 <- 1
    }
  }
  
  table(a$ID1, a$ID2)
  
  # remove reads overlapping with mm10 blacklist region
  y <- a[ a$ID1==0 & a$ID2==0 ,]
  dim(y) # 92,452 contacts left
  sum(y[,2]!=y[,4]) # 38,203 inter-chr contacts
  sum(y[,2]==y[,4]) # 54,249 intra-chr contacts >10Kb
  
  write.table(y[,1:10], file=paste0(dir2, modified_string,".UU.NoSelfLigation.NoBlackList.pairs"), row.names=F, col.names=F, sep='\t', quote=F)
  
}


files = list.files(path=dir2, pattern=".UU.NoSelfLigation.NoBlackList.pairs", all.files=FALSE, full.names=FALSE)
for (file in files){
  
  modified_string <- sub("\\.UU.NoSelfLigation.NoBlackList.pairs$", "", file)
  y<-read.table(paste0(dir2,file),head=F)
  dim(y)
  
  nrow(y) # 57,983 rows
  nrow(unique(y[,-1])) # 57,983
  
  PE <- table(y[,1])
  length(PE) # 41,216 PE reads
  
  length(PE[PE==1]) # 32,349 PE reads only support pair-wise contacts
  length(PE[PE>=2]) #  8,867 PE reads support multi-way contacts (>=3)
  
  u <- y[y[,1] %in% names(PE[PE==1]),] # pair-wise contacts
  v <- y[y[,1] %in% names(PE[PE>=2]),] # multi-way contacts
  
  # summarize # of unique 10Kb bin for each PE read supporting multi-way contacts
  rec <- as.data.frame(matrix(0, nrow=length(PE[PE>=2]), ncol=2))
  colnames(rec)<-c('PE', 'Num.10Kb.Bin')
  rec$PE <- names(PE[PE>=2])
  
  for(i in 1:nrow(rec))
  {
    if(i%%1000==0)
    {
      print(c(i, i/nrow(rec))) 
    }
    
    w <- v[ v[,1]==rec$PE[i] ,]  
    
    w.l <- w[,2:3]
    w.r <- w[,4:5]
    colnames(w.l)<-c('chr','pos')
    colnames(w.r)<-c('chr','pos')
    anchor <- unique(rbind(w.l, w.r))
    anchor <- anchor[order(anchor$chr, anchor$pos),]
    anchor$Bin10Kb <- floor(anchor$pos/1e4)*1e4
    anchor <- unique(anchor[,-2])
    rec$Num.10Kb.Bin[i] <- nrow(anchor)
  }
  
  table(rec$Num.10Kb.Bin)
  
  nrow(rec) # 8,867 PE reads supporting multi-way contacts (>=3)
  out <- rec[ rec$Num.10Kb.Bin>=3 ,]
  nrow(out) # 5,770 PE reads connecting 3 unique 10Kb bins
  
  table(rec[,2])
  # 3,097 PE reads connect 2 unique 10Kb bins (need to remove)
  # 5,017 PE reads connect 3 unique 10Kb bins 
  #   651 PE reads connect 4 unique 10Kb bins
  #    92 PE reads connect 5 unique 10Kb bins
  #     9 PE reads connect 6 unique 10Kb bins
  #     1 PE read  connect 7 unique 10Kb bins
  
  z<-y[y[,1] %in% out$PE,]
  nrow(z)               # 19,245 contacts
  length(unique(z[,1])) #  5,770 PE reads, each PE read connects >=3 unique 10Kb bins
  
  write.table(z, file=paste0(dir2, modified_string,'.Noselfligation.NoBlackList.Multi10KbBin.UU.pairs'),row.names=F, col.names=F, sep='\t', quote=F)
}
