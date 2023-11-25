# packages
library(tidyverse)
library(edgeR)
library(readxl)

# load all data 
load(file='raw_data/smRa.Rdata')

# 
miRa <- smRa[["miRa"]]

miRa1 <-miRa %>% select(Mature_ID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                        ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(miRa1)) {
  out <- sum(miRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

miRa1_num <- cbind(miRa1, num) %>% 
  as_tibble() 

sum(miRa1_num$num)/(nrow(miRa1)*84)
z <- miRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)

miRa1 %>% ggplot(aes(num))+geom_bar()

mirna <- miRa1_num %>% filter(num<42) %>% select(-num)


for (i in 2:43) {
  mirna1 <- column_to_rownames(mirna, var ="Mature_ID") %>% 
    select(colnames(mirna)[[i]], colnames(mirna)[[i+42]])
  dge <- DGEList(counts=mirna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  mirna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue, normEXP$FDR)
  colnames(mirna_fc)[1] <- "Mature_ID"
  colnames(mirna_fc)[2] <- str_c(str_sub(colnames(mirna)[[i]], 1, 8), "_FC")
  colnames(mirna_fc)[3] <- str_c(str_sub(colnames(mirna)[[i]], 1, 8), "_pvalue")
  colnames(mirna_fc)[4] <- str_c(str_sub(colnames(mirna)[[i]], 1, 8), "_FDR")
  mirna <- merge(mirna, mirna_fc, by="Mature_ID")
}

# select , 508 x 42
miR <-  mirna %>%  select(Mature_ID, contains("_FC"), contains("_pvalue")) %>% as_tibble() 


###########
piRa <- smRa[["piRa"]]

piRa1 <-piRa %>% select(SeqID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                        ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(piRa1)) {
  out <- sum(piRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

piRa1_num <- cbind(piRa1, num) %>% 
  as_tibble() 

sum(piRa1_num$num)/(nrow(piRa1)*84)
z <- piRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
piRa1 %>% ggplot(aes(num))+geom_bar()

pirna <- piRa1_num %>% filter(num<42) %>% select(-num)

for (i in 2:43) {
  pirna1 <- column_to_rownames(pirna, var ="SeqID") %>% 
    select(colnames(pirna)[[i]], colnames(pirna)[[i+42]])
  dge <- DGEList(counts=pirna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  pirna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue , normEXP$FDR)
  colnames(pirna_fc)[1] <- "SeqID"
  colnames(pirna_fc)[2] <- str_c(str_sub(colnames(pirna)[[i]], 1, 8), "_FC")
  colnames(pirna_fc)[3] <- str_c(str_sub(colnames(pirna)[[i]], 1, 8), "_pvalue")
  colnames(pirna_fc)[4] <- str_c(str_sub(colnames(pirna)[[i]], 1, 8), "_FDR")
  pirna <- merge(pirna, pirna_fc, by="SeqID")
}

# select 353x 42
piR <-  pirna %>%  select(SeqID, contains("_FC"), contains("_pvalue")) %>% as_tibble() 

##########
snoRa <- smRa[["snoRa"]]

snoRa1 <-snoRa %>% select(SeqID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                          ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(snoRa1)) {
  out <- sum(snoRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

snoRa1_num <- cbind(snoRa1, num) %>% 
  as_tibble() 

sum(snoRa1_num$num)/(nrow(snoRa1)*84)
z <- snoRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
snoRa1 %>% ggplot(aes(num))+geom_bar()

snorna <- snoRa1_num %>% filter(num<42) %>% select(-num)

for (i in 2:43) {
  snorna1 <- column_to_rownames(snorna, var ="SeqID") %>% 
    select(colnames(snorna)[[i]], colnames(snorna)[[i+42]])
  dge <- DGEList(counts=snorna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  snorna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue , normEXP$FDR)
  colnames(snorna_fc)[1] <- "SeqID"
  colnames(snorna_fc)[2] <- str_c(str_sub(colnames(snorna)[[i]], 1, 8), "_FC")
  colnames(snorna_fc)[3] <- str_c(str_sub(colnames(snorna)[[i]], 1, 8), "_pvalue")
  colnames(snorna_fc)[4] <- str_c(str_sub(colnames(pirna)[[i]], 1, 8), "_FDR")
  snorna <- merge(snorna, snorna_fc, by="SeqID")
}

# select 629x 41
snoR <-  snorna %>%  select(SeqID, contains("_FC"), contains("_pvalue")) %>% as_tibble() 

#######
snRa <- smRa[["snRa"]]

snRa1 <-snRa %>% select(SeqID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                        ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(snRa1)) {
  out <- sum(snRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

snRa1_num <- cbind(snRa1, num) %>% 
  as_tibble() 

sum(snRa1_num$num)/(nrow(snRa1)*84)
z <- snRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
snRa1 %>% ggplot(aes(num))+geom_bar()

snRna <- snRa1_num %>% filter(num<42) %>% select(-num)

for (i in 2:43) {
  snRna1 <- column_to_rownames(snRna, var ="SeqID") %>% 
    select(colnames(snRna)[[i]], colnames(snRna)[[i+42]])
  dge <- DGEList(counts=snRna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  snRna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue, normEXP$FDR)
  colnames(snRna_fc)[1] <- "SeqID"
  colnames(snRna_fc)[2] <- str_c(str_sub(colnames(snRna)[[i]], 1, 8), "_FC")
  colnames(snRna_fc)[3] <- str_c(str_sub(colnames(snRna)[[i]], 1, 8), "_pvalue")
  colnames(snRna_fc)[4] <- str_c(str_sub(colnames(snRna)[[i]], 1, 8), "_FDR")
  snRna <- merge(snRna, snRna_fc, by="SeqID")
}

# select 2000x 42
snR <-  snRna %>%  select(SeqID, contains("_FC"), contains("_pvalue")) %>% as_tibble() 


#######
tRa <- smRa[["tRa"]]

tRa1 <-tRa %>% select(SeqID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                      ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(tRa1)) {
  out <- sum(tRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

tRa1_num <- cbind(tRa1, num) %>% 
  as_tibble() 

sum(tRa1_num$num)/(nrow(tRa1)*84)
z <- tRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
tRa1 %>% ggplot(aes(num))+geom_bar()

tRna <- tRa1_num %>% filter(num<42) %>% select(-num)

for (i in 2:43) {
  tRna1 <- column_to_rownames(tRna, var ="SeqID") %>% 
    select(colnames(tRna)[[i]], colnames(tRna)[[i+42]])
  dge <- DGEList(counts=tRna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  tRna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue, normEXP$FDR)
  colnames(tRna_fc)[1] <- "SeqID"
  colnames(tRna_fc)[2] <- str_c(str_sub(colnames(tRna)[[i]], 1, 8), "_FC")
  colnames(tRna_fc)[3] <- str_c(str_sub(colnames(tRna)[[i]], 1, 8), "_pvalue")
  colnames(tRna_fc)[4] <- str_c(str_sub(colnames(tRna)[[i]], 1, 8), "_FDR")
  tRna <- merge(tRna, tRna_fc, by="SeqID")
}

# select 2523x 41
tR <-  tRna %>%  select(SeqID, contains("_FC"), contains("_pvalue")) %>% as_tibble() 


#######
yRa <- smRa[["yRa"]]

yRa1 <-yRa %>% select(SeqID, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                      ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count")) 

num <- c() 
for (i in 1:nrow(yRa1)) {
  out <- sum(yRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

yRa1_num <- cbind(yRa1, num) %>% 
  as_tibble() 

sum(yRa1_num$num)/(nrow(yRa1)*84)
z <- yRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
yRa1 %>% ggplot(aes(num))+geom_bar()

yRna <- yRa1_num %>% filter(num<42) %>% select(-num)

for (i in 2:43) {
  yRna1 <- column_to_rownames(yRna, var ="SeqID") %>% 
    select(colnames(yRna)[[i]], colnames(yRna)[[i+42]])
  dge <- DGEList(counts=yRna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=10000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  yRna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue, normEXP$FDR)
  colnames(yRna_fc)[1] <- "SeqID"
  colnames(yRna_fc)[2] <- str_c(str_sub(colnames(yRna)[[i]], 1, 8), "_FC")
  colnames(yRna_fc)[3] <- str_c(str_sub(colnames(yRna)[[i]], 1, 8), "_pvalue")
  colnames(yRna_fc)[4] <- str_c(str_sub(colnames(yRna)[[i]], 1, 8), "_FDR")
  yRna <- merge(yRna, yRna_fc, by="SeqID")
}

# select 700x 42
yR <-  yRna %>%  select(SeqID, contains("_FC"), contains("_pvalue")) %>% as_tibble() 


########## mR
mRa <- smRa[["mRa"]]

mRa1 <- mRa %>% select(Gene_Symbol, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                       ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count"))


num <- c() 
for (i in 1:nrow(mRa1)) {
  out <- sum(mRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

mRa1_num <- cbind(mRa1, num) %>% 
  as_tibble() 

sum(mRa1_num$num)/(nrow(mRa1)*84)
z <- mRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
mRa1 %>% ggplot(aes(num))+geom_bar()

mrna <- mRa1_num %>% filter(num<42) %>% select(-num)


for (i in 2:43) {
  mrna1 <- column_to_rownames(mrna, var ="Gene_Symbol") %>% 
    select(colnames(mrna)[[i]], colnames(mrna)[[i+42]])
  dge <- DGEList(counts=mrna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=30000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  mrna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue, normEXP$FDR)
  colnames(mrna_fc)[1] <- "Gene_Symbol"
  colnames(mrna_fc)[2] <- str_c(str_sub(colnames(mrna)[[i]], 1, 8), "_FC")
  colnames(mrna_fc)[3] <- str_c(str_sub(colnames(mrna)[[i]], 1, 8), "_pvalue")
  colnames(mrna_fc)[4] <- str_c(str_sub(colnames(mrna)[[i]], 1, 8), "_FDR")
  mrna <- merge(mrna, mrna_fc, by="Gene_Symbol")
}

### select 14908x 42
mR <-  mrna %>%  select(Gene_Symbol, contains("_FC"), contains("_pvalue")) %>% as_tibble() 



########## lncR
lncRa <- smRa[["lncRa"]]

lncRa1 <- lncRa %>% select(Gene_Symbol, ends_with("01_Read_Count"), ends_with("Pre_Read_Count"), ends_with("-1_Read_Count"),
                           ends_with("02_Read_Count"), ends_with("After_Read_Count"), ends_with("-2_Read_Count"))


num <- c() 
for (i in 1:nrow(lncRa1)) {
  out <- sum(lncRa1[i,]==0)
  num <- rbind(num, out)
} # num - NA number

lncRa1_num <- cbind(lncRa1, num) %>% 
  as_tibble() 

sum(lncRa1_num$num)/(nrow(lncRa1)*84)
z <- lncRa1_num %>% filter(num<42) 
sum(z$num)/(nrow(z)*84)
lncRa1 %>% ggplot(aes(num))+geom_bar()

lncRna <- lncRa1_num %>% filter(num<42) %>% select(-num)


for (i in 2:43) {
  lncRna1 <- column_to_rownames(lncRna, var ="Gene_Symbol") %>% 
    select(colnames(lncRna)[[i]], colnames(lncRna)[[i+42]])
  dge <- DGEList(counts=lncRna1, group=c(1,2))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group, dge$samples)
  rownames(design) <- rownames(dge$samples)
  dge_dispersion <- estimateGLMCommonDisp(dge)
  dge_tagwise <- estimateGLMTagwiseDisp(dge_dispersion)
  plotBCV(dge_tagwise)
  fit <- glmFit(dge_tagwise, dispersion = dge_tagwise$tagwise.dispersion)
  lrt <- glmLRT(fit)
  norm <- round(cpm(dge_tagwise, normalized.lib.sizes = T))
  normEXP <- merge(norm, topTags(lrt, n=30000), by.x="row.names", by.y="row.names")
  names(normEXP)[1]
  head(normEXP)
  lncRna_fc <- data.frame(normEXP$Row.names, normEXP$logFC, normEXP$PValue, normEXP$FDR)
  colnames(lncRna_fc)[1] <- "Gene_Symbol"
  colnames(lncRna_fc)[2] <- str_c(str_sub(colnames(lncRna)[[i]], 1, 8), "_FC")
  colnames(lncRna_fc)[3] <- str_c(str_sub(colnames(lncRna)[[i]], 1, 8), "_pvalue")
  colnames(lncRna_fc)[4] <- str_c(str_sub(colnames(lncRna)[[i]], 1, 8), "_FDR")
  lncRna <- merge(lncRna, lncRna_fc, by="Gene_Symbol")
}

# 3151 x 41
lncR <-  lncRna %>%  select(Gene_Symbol, contains("_FC"), contains("_pvalue")) %>% as_tibble() 

rt_fc <- list(miR, piR, mR, lncR, snoR, snR, tR, yR)
names(rt_fc) <- c("miR", "piR", "mR", "lncR", "snoR", "snR", "tR", "yR")
save(rt_fc, file="process_data/rt_fc1.Rdata")


rt_fc <- list(mirna, pirna, mrna, lncRna, snorna, snRna, tRna, yRna)
names(rt_fc) <- c("miR", "piR", "mR", "lncR", "snoR", "snR", "tR", "yR")
save(rt_fc, file="process_data/rt_fc2.Rdata")
