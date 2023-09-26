library(readxl)
library(tidyverse)
library(lubridate)
library(edgeR)
library(lubridate)
library(Hmisc) # rcorr
library(leaps) # rugsubset
library(xlsx)

### loading fc miR: 508 (4), piR:353(2), snoR: 629 (3), snR: 2000 (0),tR:2523(1), yR:700 (1),lncR: 3151 (11), mR:14908 (33) -24772 (55) # at least 1 read count
load(file='process_data/rt_fc1.Rdata')
miR <- rt_fc[["miR"]]; piR <- rt_fc[["piR"]]; snoR <- rt_fc[["snoR"]]
mR <- rt_fc[["mR"]]; lncR <- rt_fc[["lncR"]]; snR <- rt_fc[["snR"]]
tR <- rt_fc[["tR"]]; yR <- rt_fc[["yR"]]

Cx <- read_excel("process_data/group.xlsx", sheet = "Cx_e") 

cx <- Cx %>% mutate(ID=substr(SID, 1,8), recur_m=ifelse(recur1>1,1,0)) %>% select(ID, css, group1, Hb0:min_PLT,  Hb1:Mo3) 


########### mR
mR <- mR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="Gene_Symbol") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID

mR1 <- mR %>% 
    inner_join(cx, by="ID") 
mRa <- mR1 %>% column_to_rownames(var="ID")
x<-select_if (mRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)

mR_cor <- result %>% as_tibble() %>% mutate(Gene_Symbol = colnames(result)) %>% select(Gene_Symbol, css:Mo3) %>% as.data.frame()
mR_cor <- mR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 

t <- mR_cor %>% select(-(Gene_Symbol:group1), -(Hb2:Mo2), -(Hb3:Mo3))

### hematological correlation evlautation (total)
t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                                     PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                                     ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                                     ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                                     Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                                     NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                                     PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                                     LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))
# number of relevant parameters (sum)
for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
    }
mR_cor$sum  <- rowSums(t[,1:21])
mR_cor$total <- rowSums(t[,22:29])

# at least MR relevant 3types of hemaological paramteres 
mR_sel <- mR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))

# fc_selection
mR_fc <- mR1 %>% select(ID, css, group1, mR_sel$Gene_Symbol)

## mR_sel parmater details
mR_sel1 <- t %>% mutate(Gene_Symbol=mR_cor$Gene_Symbol)  %>% select(-(Hb0:LMR1)) %>%
    inner_join(mR_sel,by="Gene_Symbol")  %>% select(Gene_Symbol, Hb:LMR, total, css:sum) %>% 
    rename("ID"="Gene_Symbol") %>% mutate(type=rep("mR", nrow(mR_sel)))



################# miR
miR <- miR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="Mature_ID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
miR1 <- miR %>% 
    inner_join(cx, by="ID") 
miRa <- miR1 %>% column_to_rownames(var="ID")
x<-select_if (miRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
miR_cor <- result %>% as_tibble() %>% mutate(Mature_ID = colnames(result)) %>% select(Mature_ID, css:Mo3) %>% as.data.frame() 
##
miR_cor <- miR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 


t <- miR_cor %>% select(-(Mature_ID:group1), -(Hb2:Mo2), -(Hb3:Mo3))


t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                               PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                               ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                               ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                               Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                               NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                               PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                               LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))

for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
miR_cor$sum  <- rowSums(t[,1:21])
miR_cor$total <- rowSums(t[,22:29])


miR_sel <- miR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))

miR_fc <- miR1 %>% select(ID, miR_sel$Mature_ID) 

miR_sel1 <- t %>% mutate(Mature_ID=miR_cor$Mature_ID) %>% select(-(Hb0:LMR1)) %>% 
    inner_join(miR_sel, by="Mature_ID") %>% select(Mature_ID, Hb:LMR, total, css:sum) %>% rename("ID"="Mature_ID") %>% mutate(type=rep("miR", nrow(miR_sel)))

#piR######3

piR <- piR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="SeqID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
piR1 <- piR %>% 
    inner_join(cx, by="ID") 
piRa <- piR1 %>% column_to_rownames(var="ID")
x<-select_if (piRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
piR_cor <- result %>% as_tibble() %>% mutate(SeqID = colnames(result)) %>% select(SeqID, css:Mo3) %>% as.data.frame()
piR_cor <- piR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 


t <- piR_cor %>% select(-(SeqID:group1), -(Hb2:Mo2), -(Hb3:Mo3))

t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                               PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                               ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                               ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                               Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                               NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                               PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                               LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))


for (i in 1:22) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
piR_cor$sum  <- rowSums(t[,1:21])
piR_cor$total <- rowSums(t[,22:29])

piR_sel <- piR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))

piR_fc <- piR1 %>% select(ID, piR_sel$SeqID) 

piR_sel1 <- t %>% mutate(SeqID=piR_cor$SeqID) %>% select(-(Hb0:LMR1)) %>% 
    inner_join(piR_sel, by="SeqID") %>% select(SeqID, Hb:LMR, total, css:sum) %>% rename("ID"="SeqID") %>% mutate(type=rep("piR", nrow(piR_sel)))

######## snoR
snoR <- snoR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="SeqID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
snoR1 <- snoR %>% inner_join(cx, by="ID") 
snoRa <- snoR1 %>% column_to_rownames(var="ID")
x<-select_if (snoRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
snoR_cor <- result %>% as_tibble() %>% mutate(SeqID = colnames(result)) %>% select(SeqID, css:Mo3) %>% as.data.frame()
snoR_cor <- snoR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 

t <- snoR_cor %>% select(-(SeqID:group1), -(Hb2:Mo2), -(Hb3:Mo3))
t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                                PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                                ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                                ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                                Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                                NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                                PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                                LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))


for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
snoR_cor$sum  <- rowSums(t[,1:21])
snoR_cor$total <- rowSums(t[,22:29])

# seperately saved to rna_cor.Rdata -snoR_cor

snoR_sel <- snoR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))

snoR_fc <- snoR1 %>% select(ID, snoR_sel$SeqID) 

snoR_sel1 <- t %>% mutate(SeqID=snoR_cor$SeqID) %>% select(-(Hb0:LMR1)) %>% 
    inner_join(snoR_sel, by="SeqID") %>% select(SeqID, Hb:LMR, total, css:sum) %>% rename("ID"="SeqID") %>% mutate(type=rep("snoR", nrow(snoR_sel)))


########### snR

snR <- snR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="SeqID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
snR1 <- snR %>% 
    inner_join(cx, by="ID") 
snRa <- snR1 %>% column_to_rownames(var="ID")
x<-select_if (snRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
snR_cor <- result %>% as_tibble() %>% mutate(SeqID = colnames(result)) %>% select(SeqID, css:Mo3)  %>% as.data.frame()
snR_cor <- snR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 

t <- snR_cor %>% select(-(SeqID:group1), -(Hb2:Mo2), -(Hb3:Mo3))
t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                               PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                               ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                               ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                               Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                               NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                               PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                               LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))


for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
snR_cor$sum  <- rowSums(t[,1:21])
snR_cor$total <- rowSums(t[,22:29])

# seperately saved to rna_cor.Rdata -snR_cor

snR_sel <- snR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))
snR_fc <- snR1 %>% select(ID, snR_sel$SeqID)  

snR_sel1 <- t %>% mutate(SeqID=snR_cor$SeqID) %>% select(-(Hb0:LMR1)) %>% 
    inner_join(snR_sel, by="SeqID") %>% select(SeqID, Hb:LMR, total, css:sum) %>% rename("ID"="SeqID") %>%  mutate(type=rep("snR", nrow(snR_sel)))

###
tR <- tR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="SeqID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
tR1 <- tR %>% 
    inner_join(cx, by="ID") 
tRa <- tR1 %>% column_to_rownames(var="ID")
x<-select_if (tRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
tR_cor <- result %>% as_tibble() %>% mutate(SeqID = colnames(result)) %>% select(SeqID, css:Mo3) %>% as.data.frame()
tR_cor <- tR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 

t <- tR_cor %>% select(-(SeqID:group1), -(Hb2:Mo2), -(Hb3:Mo3))

t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                              PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                              ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                              ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                              Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                              NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                              PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                              LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))


for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
tR_cor$sum  <- rowSums(t[,1:21])
tR_cor$total <- rowSums(t[,22:29])

tR_sel <- tR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))
tR_fc <- tR1 %>% select(ID, tR_sel$SeqID) 
tR_sel1 <- t %>% mutate(SeqID=tR_cor$SeqID) %>% select(-(Hb0:LMR1)) %>% 
    inner_join(tR_sel, by="SeqID") %>% select(SeqID, Hb:LMR, total, css:sum) %>% rename("ID"="SeqID") %>%  mutate(type=rep("tR", nrow(tR_sel)))

########## yR
yR <- yR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="SeqID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
yR1 <- yR %>% 
    inner_join(cx, by="ID") 
yRa <- yR1 %>% column_to_rownames(var="ID")
x<-select_if (yRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
yR_cor <- result %>% as_tibble() %>% mutate(SeqID = colnames(result)) %>% select(SeqID, css:Mo3) %>% as.data.frame()
yR_cor <- yR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 

t <- yR_cor %>% select(-(SeqID:group1), -(Hb2:Mo2), -(Hb3:Mo3))
t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                              PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                              ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                              ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                              Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                              NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                              PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                              LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))


for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
yR_cor$sum  <- rowSums(t[,1:21])
yR_cor$total <- rowSums(t[,22:29])
yR_sel <- yR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))
yR_fc <- yR1 %>% select(ID, yR_sel$SeqID)  
yR_sel1 <- t %>% mutate(SeqID=yR_cor$SeqID) %>% select(-(Hb0:LMR1)) %>% 
    inner_join(yR_sel, by="SeqID") %>% select(SeqID, Hb:LMR, total, css:sum) %>% rename("ID"="SeqID") %>%  mutate(type=rep("yR", nrow(yR_sel)))

### lncR
lncR <- lncR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="Gene_Symbol") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
lncR1 <- lncR %>% 
    inner_join(cx, by="ID") 
lncRa <- lncR1 %>% column_to_rownames(var="ID")
x<-select_if (lncRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
lncR_cor <- result %>% as_tibble() %>% mutate(Gene_Symbol = colnames(result)) %>% select(Gene_Symbol, css:Mo3) %>% as.data.frame()
lncR_cor <- lncR_cor[-((nrow(result)-32):nrow(result)),] %>% as_tibble() 

t <- lncR_cor %>% select(-(Gene_Symbol:group1), -(Hb2:Mo2), -(Hb3:Mo3))
t <- t %>% mutate(Hb=ifelse(Hb0==0 & min_Hb==0 & Hb1==0, 0, 1), 
                                PLT=ifelse(PLT0==0 & min_PLT==0 & PLT1==0, 0, 1),
                                ANC=ifelse(ANC0==0 & min_ANC==0 & ANC1==0, 0, 1),
                                ALC=ifelse(ALC0==0 & min_ALC==0 & ALC1==0, 0, 1),
                                Mo=ifelse(Mo0==0 & min_Mo==0 & Mo1==0, 0, 1),
                                NLR=ifelse(NLR0==0 & NLR1==0, 0, 1),
                                PLR=ifelse(PLR0==0 & PLR1==0, 0, 1),
                                LMR=ifelse(LMR0==0 & LMR1==0, 0, 1))


for (i in 1:21) {
    t[,i]=ifelse(t[,i]==0,0,1)
}
lncR_cor$sum  <- rowSums(t[,1:21])
lncR_cor$total <- rowSums(t[,22:29])
lncR_sel <- lncR_cor %>% filter((total>2 & abs(css)>0 & abs(group1)>0)|(abs(css)>0.5 & abs(group1)>0)) %>% arrange(desc(sum))
lncR_fc <- lncR1 %>% select(ID, lncR_sel$Gene_Symbol) 
lncR_sel1 <- t %>% mutate(Gene_Symbol=lncR_cor$Gene_Symbol)  %>% select(-(Hb0:LMR1)) %>%
    inner_join(lncR_sel,by="Gene_Symbol")  %>% select(Gene_Symbol, Hb:LMR, total, css:sum) %>% rename("ID"="Gene_Symbol") %>%  mutate(type=rep("lncR", nrow(lncR_sel)))


## rna_cor 
rna_cor <- list(mR_cor, lncR_cor, miR_cor, piR_cor, snoR_cor, snR_cor, tR_cor, yR_cor)
names(rna_cor) <- c("mR_cor", "lncR_cor", "miR_cor", "piR_cor","snoR_cor", "snR_cor", "tR_cor", "yR_cor")
save(rna_cor, file="process_data/rna_cor.Rdata")


## network
all_fc <- mR_fc %>% inner_join(lncR_fc, by="ID") %>% 
    inner_join(miR_fc, by="ID") %>% 
    inner_join(piR_fc, by="ID") %>% 
    inner_join(snoR_fc, by="ID") %>%
    inner_join(snR_fc, by="ID") %>% 
    inner_join(tR_fc, by="ID")%>% 
    inner_join(yR_fc, by="ID")

rf <- Cx %>% mutate(ID1=ID, ID=substr(SID, 1,8), recur_m=ifelse(recur1>1,1,0)) %>% select(ID, ID1) %>% inner_join(all_fc, by="ID")  %>% mutate(ID=ID1) %>% select(-ID1)
write.xlsx(rf, "process_data/reference.xlsx", sheetName = "network_fc", append = T)

library(igraph)
c <- all_fc %>% select(-ID,-css,-group1)

##
c <- data.frame(t(c))
c[is.na(c)] <- 0
rownames(c)

g <- graph.adjacency(
    as.matrix(as.dist(cor(t(c), method="pearson"))),
    mode="undirected",
    weighted=TRUE,
    diag=FALSE
)
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
E(g)[which(E(g)$weight<0)]$color <- "darkblue"
E(g)[which(E(g)$weight>0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight<=0.5)])
g <- delete.vertices(g, degree(g)==0)
V(g)$name <- V(g)$name
V(g)$shape <- "circle"
V(g)$color <- "skyblue"
V(g)$vertex.frame.color <- "white"
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(c, 1, mean)) + 1.0) * 5
edgeweights <- E(g)$weight * 4
mst <- mst(g, algorithm="prim")

par(mfrow=c(1,1))


mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1


plot(
    mst,
    layout=layout_as_tree(g),
    edge.curved=FALSE,
    vertex.size=vSizes,
    vertex.label.dist=0,
    vertex.label.color="black",
    asp=FALSE,
    vertex.label.cex=0.8,
    edge.width=edgeweights,
    edge.arrow.mode=F,
    edge.arrow.size=1,
    main="")

central <- data.frame(degree(mst))
central$ID <- rownames(central)
central$closeness_centrality <- closeness(mst)
central$betweeness_centrality <- betweenness(mst)
central$eigen_centrality <-  eigen_centrality(mst)$vector
central$community <- mst.communities$membership
central$modularity <- mst.communities$modularity
central <- central %>% arrange(community)  


## hematological parameter correlation number
all_sel1 <- mR_sel1 %>% bind_rows(lncR_sel1, miR_sel1, piR_sel1, snoR_sel1, tR_sel1, yR_sel1)


central <- central %>% inner_join(all_sel1, by="ID")

central %>% as_tibble() %>% arrange(desc(abs(group1))) %>% select(ID, css, group1, community, type, total) %>% as.data.frame()

library(xlsx)
write.xlsx(all_sel1, "process_data/group.xlsx", sheetName = "all_rna1", append = T)
write.xlsx(central, "process_data/group.xlsx", sheetName = "network_css_05", append = T)
write.xlsx(central, "process_data/group.xlsx", sheetName = "network_css_06", append = T)


library(leaps)

c1 <- mR %>% inner_join(miR, by="ID") %>%
    inner_join(piR, by="ID")%>%
    inner_join(snoR, by="ID") %>%
    inner_join(snR, by="ID") %>%
    inner_join(lncR, by="ID") %>%
    inner_join(yR, by="ID") %>% 
    inner_join(cx, by="ID") %>% select(css, central$ID)


c1[is.na(c1)] <- 0


g <-  regsubsets(css ~ .,
                 data = c1,
                 nbest = 1,       # 1 best model for each number of predictors
                 nvmax = 10,    # NULL for no limit on number of variables
                 force.in = NULL, force.out = NULL,
                 method = "exhaustive")
out <- summary(g)
out$adjr2

### figure2D, Subgroup selection for stage

g$xnames <- ifelse(str_sub(g$xnames,1,4)=="`hsa", str_sub(g$xnames,6,(nchar(g$xnames)-1)), g$xnames)
g$xnames <- ifelse(str_sub(g$xnames,1,3)=="URS", paste("URS",str_sub(g$xnames,4,nchar(g$xnames)),sep="\n"), g$xnames)
g$xnames <- ifelse(str_sub(g$xnames,1,3)=="LOC", paste("LOC",str_sub(g$xnames,4,nchar(g$xnames)),sep="\n"), g$xnames)

plot(g, scale = "adjr2")


####
