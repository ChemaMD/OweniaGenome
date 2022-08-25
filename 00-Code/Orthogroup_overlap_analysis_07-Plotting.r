require(tidyverse)
require(pheatmap)
require(RColorBrewer)


## 1. Owenia fusiformis vs. Capitella teleta ALL

elts=read.table('01-ctel_ofus_phyper_ALL.tsv',h=T)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% arrange(paj) 

elts[elts$pv>0.05,]
elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l
elts.l

elts.m <- elts.m[order(rownames(elts.m)),]
elts.m <- elts.m[,order(colnames(elts.m))]
elts.m <- cbind(elts.m[-c(2:4)],elts.m[c(2:4)])
elts.m <- rbind(elts.m[-c(2:4),],elts.m[c(2:4),])

table(-log10(elts.m)<0.01)

elts.final <- -log10(elts.m)
write.table(elts.final, "07-Ofus2Ctel_ALL_log10.txt", 
            sep ='\t', quote = F, row.names = T)


## 2. Owenia fusiformis vs. Capitella teleta TFs

elts=read.table('04-ctel_ofus_phyper_TF.tsv',h=T)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% arrange(paj) 

elts[elts$pv>0.05,]
elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l
elts.l

elts.m <- elts.m[order(rownames(elts.m)),]
elts.m <- elts.m[,order(colnames(elts.m))]
elts.m <- cbind(elts.m[-c(2:4)],elts.m[c(2:4)])
elts.m <- rbind(elts.m[-c(2:4),],elts.m[c(2:4),])

table(-log10(elts.m)<0.01)

elts.final <- -log10(elts.m)
write.table(elts.final, "10-Ofus2Ctel_TF_log10.txt", 
            sep ='\t', quote = F, row.names = T)



## 3. Owenia fusiformis vs. Dimorphilus gyrociliatus ALL

elts=read.table('02-dgyr_ofus_phyper_ALL.tsv',h=T)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% arrange(paj) 

elts[elts$pv>0.05,]
elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l
elts.l

elts.m <- elts.m[order(rownames(elts.m)),]
elts.m <- elts.m[,order(colnames(elts.m))]
elts.m <- rbind(elts.m[-c(2:4),],elts.m[c(2:4),])

table(-log10(elts.m)<0.01)

elts.final <- -log10(elts.m)
write.table(elts.final, "08-Ofus2Dgyr_ALL_log10.txt", 
            sep ='\t', quote = F, row.names = T)




## 4. Owenia fusiformis vs. Dimorphilus gyrociliatus TF

elts=read.table('05-dgyr_ofus_phyper_TF.tsv',h=T)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% arrange(paj) 

elts[elts$pv>0.05,]
elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l
elts.l

elts.m <- elts.m[order(rownames(elts.m)),]
elts.m <- elts.m[,order(colnames(elts.m))]
elts.m <- rbind(elts.m[-c(2:4),],elts.m[c(2:4),])

table(-log10(elts.m)<0.01)

elts.final <- -log10(elts.m)
write.table(elts.final, "11-Ofus2Dgyr_TF_log10.txt", 
            sep ='\t', quote = F, row.names = T)





## 5. Capitella teleta vs. Dimorphilus gyrociliatus ALL

elts=read.table('03-ctel_dgyr_phyper_ALL.tsv',h=T)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% arrange(paj) 

elts[elts$pv>0.05,]
elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l
elts.l

elts.m <- elts.m[order(rownames(elts.m)),]
elts.m <- elts.m[,order(colnames(elts.m))]
elts.m <- cbind(elts.m[-c(2:4)],elts.m[c(2:4)])

table(-log10(elts.m)<0.01)

elts.final <- -log10(elts.m)
write.table(elts.final, "09-Ctel2Dgyr_ALL_log10.txt", 
            sep ='\t', quote = F, row.names = T)



## 6. Capitella teleta vs. Dimorphilus gyrociliatus TF

elts=read.table('06-ctel_dgyr_phyper_TF.tsv',h=T)

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% arrange(paj) 

elts[elts$pv>0.05,]
elts %>% select(cl1,cl2,pv) %>% pivot_wider(names_from = cl1,id_cols=cl2,values_from=pv)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=paj)%>%
  column_to_rownames(var='cl2') -> elts.m

elts %>% mutate(paj=p.adjust(pv,method='BH')) %>% select(cl1,cl2,paj) %>%
  group_by(cl1) %>% mutate(best=max(-log10(paj)),mv=ifelse(-log10(paj)==best,'*','')) %>%
  pivot_wider(names_from = cl1,id_cols=cl2,values_from=mv) %>%
  column_to_rownames(var='cl2') -> elts.l
elts.l

elts.m <- elts.m[order(rownames(elts.m)),]
elts.m <- elts.m[,order(colnames(elts.m))]
elts.m <- rbind(elts.m[-c(2:4),],elts.m[c(2:4),])

table(-log10(elts.m)<0.01)

elts.final <- -log10(elts.m)
write.table(elts.final, "12-Ctel2Dgyr_TF_log10.txt", 
            sep ='\t', quote = F, row.names = T)



## 7. Plot normalised quadrants

data <- read.table("14-Summary_input_data.txt", sep ='\t', header = T, row.names = 1)

heatmap_color <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(300)
heatmap_color <- heatmap_color[150:300]
heatmap_color[1] <- "#FFFFFF"

pheatmap(data,
         col = heatmap_color,
         breaks = seq(0, 2, length.out = 100), 
         cluster_cols = FALSE,
         cluster_rows = FALSE, 
         display_numbers = TRUE,
         border_color = NA, 
         cellwidth = 30, 
         cellheight = 30)


