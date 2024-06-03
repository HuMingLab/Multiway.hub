dir<-'./example.data/chimeric.sc.pairs.mapq40/'   #### Path to directory where all intermediate result files were saved in step1 (dir2 in step1)

cellnames = "./example.data/LC716_MAPQ40_cell_names_Shreya_light.txt"     ###### Path to example data cellnames file
color.loc = './example.data/color_20_YangOrder.txt'                       ###### Path to Color assigned to each celltype 

mm10.blacklist = '10KbBin_mm10_blacklist.txt'     ##### Path to mm10 Blacklist region
annotation = 'F_GC_M_MboI_10Kb_el.mm10.txt'       ##### Path to mm10 10Kb MboI annotation file

metadata = '052124_LC716_MAPQ40_metadata_MultiContact.txt'      ###### Filename to save multiway Contacts metadata file
final.result = 'hub_20celltypes_MAPQ40.txt'                     ###### Filename to save final hub file


a <- read.table(cellnames, head=F)
a<-as.data.frame(a)
colnames(a)<-c('ID','CellType')
dim(a)

a$NumContact <- 0
a$NumMultiContact <- 0
a$Num10KbBin.with.MultiContact <- 0

for(i in 1:nrow(a))
{
  if(i%%100==0)
  {
    print(c(i, i/nrow(a))) 
  }
  
  u<-read.table(paste(dir, a$ID[i], '.', a$CellType[i],'.UU.NoSelfLigation.NoBlackList.pairs', sep=''), head=F)  
  
  z<-read.table(paste(dir, a$ID[i], '.', a$CellType[i], '.Noselfligation.NoBlackList.Multi10KbBin.UU.pairs', sep=''), head=F)
  
  a$NumContact[i] <- nrow(u)
  a$NumMultiContact[i] <- nrow(z) 
  
  z.left <- z[,2:3]
  z.right <- z[,4:5]
  colnames(z.left)<-c('chr','pos')
  colnames(z.right)<-c('chr','pos')
  z.out <- unique( rbind(z.left, z.right) )
  dim(z.out)
  z.out <- z.out[order(z.out$chr, z.out$pos),]
  
  z.out$Bin <- floor(z.out$pos/1e4)*1e4
  final<-unique(z.out[,-2])
  dim(final) # 16,357 unique 10Kb bins
  
  a$Num10KbBin.with.MultiContact[i]<- nrow(final)
  
  write.table(final, file=paste(dir,'Bin.with.MultiContact.', a$ID[i], '.txt', sep=''),
              row.names=F, col.names=T, sep='\t', quote=F)
}
write.table(a, file=metadata,row.names=F, col.names=T, sep='\t', quote=F)

##########################################################################################################################################

options(scipen=999)


ann    <- read.table(annotation, sep='')
ann <- as.data.frame(ann)
colnames(ann)<- c('chr','start','end','F','GC','M')
ann$ID <- paste(ann$chr, ann$start, sep='_')

z<-read.table(metadata, head=T)
p<-read.csv(color.loc, sep='\t', head=T)

CELL20 <- as.data.frame( matrix(0, nrow=20, ncol=2) )
colnames(CELL20)<-c('CellType','CellNum')
CELL20[,1]<-p[,1]
for(i in 1:20)
{
  CELL20[i,2] <- sum(z$CellType == CELL20[i,1])  
}

for(ID in 1:20)
{
  CELL.name <- CELL20[ID, 1]
  CELL.num  <- CELL20[ID, 2]
  
  z.cell <- z[z$CellType == CELL.name,]
  rec <- as.data.frame( matrix(0, nrow=nrow(ann), ncol=CELL.num) )
  colnames(rec)<- z.cell$ID
  
  for(i in 1:nrow(z.cell))
  {
    x <- read.table(paste(dir,"/", 'Bin.with.MultiContact.',z.cell$ID[i],'.txt', sep=''), head=T)  
    x$ID <- paste(x$chr, x$Bin, sep='_')
    rec[,i]<- as.numeric( I(ann$ID %in% x$ID) )
  }
  
  rec.sum <- apply(rec, 1, sum)
  out <- cbind(ann[,1:6], rec.sum)
  write.table(out, file=paste(dir,"/",'freq_',CELL.name,'.txt', sep=''),row.names=F, col.names=T, sep='\t', quote=F)
}


#####################################################################################################################################

options(scipen=999)
f<-function(KEY)
{
  a0<-read.table(paste(dir,"/",'freq_',KEY,'.txt', sep=''), head=T)
  #dim(a0) # 272,566
  
  # keep autosomes
  a1 <- a0[ a0$chr %in% paste('chr', 1:19, sep='')  ,] 
  #dim(a1) # 246,285
  
  # keep 10Kb bins with mappability >0.8
  a <- a1[a1$M>0.8,]
  #dim(a) # 220,634
  
  a$ID <- paste(a$chr, a$start, sep='_')
  
  BLACK <- read.table(mm10.blacklist,head=T)
  
  #sum(a$ID %in% BLACK$ID)
  
  a<-a[!a$ID %in% BLACK$ID, ]
  #dim(a) # 220,580 
  
  #barplot( table(a$rec.sum) )
  
  # trim top 1%
  v <- sort(a$rec.sum)
  v.cut <- v[ round(length(v)*0.99) ]
  v.trim <- v[ v<=v.cut ]
  
  #barplot(table(v.trim))
  #print(c(mean(v.trim), sd(v.trim)))
  
  v.mean <- mean(v.trim)
  v.sd <- sd(v.trim)
  
  a$Zscore <- round((a$rec.sum - v.mean)/v.sd,4)
  
  out <- a[ a$Zscore > 1.96 ,]
  #dim(out)
  
  print(c(KEY, nrow(out), nrow(out)/nrow(a)))
  
  write.table(out, file=paste(dir,"/",'hub_', KEY,'.txt', sep=''),
              row.names=F, col.names=T, sep='\t', quote=F)
}
ann <- read.csv(color.loc, sep='\t', head=T)
for(i in 1:20)
{
  f(ann[i,1])  
}


#######################################################################################################################

ann <- read.csv(color.loc, sep='\t', head=T)

x<-read.table(annotation, head=F)
dim(x)
x<-as.data.frame(x)
colnames(x)<-c('chr','start','end','F','GC','M')
x <- x[ x$chr %in% paste('chr', 1:19, sep='') & x$M>0.8 ,] 
dim(x) # 220,634

x$ID <- paste(x$chr, x$start, sep='_')
BLACK <- read.table(mm10.blacklist, head=T)
sum(x$ID %in% BLACK$ID)

x<-x[!x$ID %in% BLACK$ID, ]
dim(x) # 220,580

options(scipen=999)

rec<-as.data.frame(matrix(0,nrow=nrow(x),ncol=20))
colnames(rec)<-ann[,1]

for(i in 1:20)
{
  p<-read.table(paste(dir,"/",'hub_',ann[i,1],'.txt', sep=''),head=T)   
  p$ID <- paste(p$chr, p$start, sep='_')
  rec[,i] <- as.numeric(I(x$ID %in% p$ID))
}

final<-cbind(x, rec)
write.table(final, file=final.result,
            row.names=F, col.names=T, sep='\t', quote=F)


