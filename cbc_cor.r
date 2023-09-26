library(readxl)
library(tidyverse)
library(reshape)
library(ggvenn)
library(ggpubr)
library(reshape)
library(ggVennDiagram)
library(xlsx)
load(file="process_data/rna_cor.Rdata")

mR_cor <- rna_cor[["mR_cor"]] %>% dplyr::rename("ID"="Gene_Symbol") %>% mutate(type="mR")
lncR_cor <- rna_cor[["lncR_cor"]] %>% dplyr::rename("ID"="Gene_Symbol") %>% mutate(type="lncR")
miR_cor <- rna_cor[["miR_cor"]] %>% dplyr::rename("ID"="Mature_ID") %>% mutate(type="miR")
piR_cor <- rna_cor[["piR_cor"]] %>% dplyr::rename("ID"="SeqID") %>% mutate(type="piR")
snoR_cor <- rna_cor[["snoR_cor"]] %>% dplyr::rename("ID"="SeqID") %>% mutate(type="snoR")
snR_cor <- rna_cor[["snR_cor"]] %>% dplyr::rename("ID"="SeqID") %>% mutate(type="snR")
tR_cor <- rna_cor[["tR_cor"]] %>% dplyr::rename("ID"="SeqID") %>% mutate(type="tR")
yR_cor <- rna_cor[["yR_cor"]] %>% dplyr::rename("ID"="SeqID") %>% mutate(type="yR")
all_cor <- mR_cor %>% bind_rows(lncR_cor, miR_cor, piR_cor, snoR_cor, snR_cor, tR_cor, yR_cor)

m <- c("Hb0", "Hb1", "Hb2", "min_Hb", "Hb3", "PLT0", "PLT1", "PLT2", "min_PLT", "PLT3",
"ANC0", "ANC1", "ANC2", "min_ANC", "ANC3", "Mo0", "Mo1", "Mo2", "min_Mo", "Mo3",
"ALC0", "ALC1", "ALC2", "min_ALC", "ALC3")

cbc <- c()
for (i in seq(m)){
c <- all_cor %>% filter(abs(get(m[i]))>0) %>% mutate(cbc=m[i]) 
cbc <- rbind(cbc, c)
}


central <- read_excel("process_data/group.xlsx", sheet = "network_css_05") 
c <- central %>% select(ID) %>% inner_join(cbc, by="ID") 
c1 <- central %>% select(ID) %>% inner_join(cbc, by="ID") %>% select(ID, type, cbc) %>%
group_by(ID) %>% summarise(n=n()) 
c<- c %>% inner_join(c1, by="ID") %>% select(ID, cbc, n) %>% cast(ID~cbc) %>% as_tibble()

c[is.na(c)] <- 0 
c[,2:24] <- ifelse(c[,2:24]>0,1,0)

c1 <-c %>% mutate(ALC=ifelse(ALC0==1|ALC1==1|ALC2==1|min_ALC==1,1,0),
ANC=ifelse(ANC0==1|ANC1==1|ANC3==1|min_ANC==1,1,0),
Hb=ifelse(Hb0==1|Hb1==1|Hb2==1|Hb3==1|min_Hb==1,1,0),
Mo=ifelse(Mo0==1|Mo1==1|Mo2==1|Mo3==1|min_Mo==1,1,0),
PLT=ifelse(PLT0==1|PLT1==1|PLT2==1|PLT3==1|min_PLT==1,1,0)) %>% 
select(ID, Hb, PLT, ANC, Mo, ALC)

c1[,2:6]<- ifelse(c1[,2:6]==1, c1$ID, 0)

x <- list(ALC=c1$ALC, ANC=c1$ANC, Hb=c1$Hb, Mo=c1$Mo, PLT=c1$PLT)
network <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle("Genes related to Network")+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))


c1[,2:6]<- ifelse(c1[,2:6]==0, 0, 1)
c2 <- c1 %>% mutate(total=rowSums(c1[,2:6]))
colSums(c1[,2:6]) %>% as_tibble() %>% 
mutate(pct=round(value/sum(value)*100,1), ID=c("Hb", "PLT", "ANC", "Mo", "ALC"))
c2 %>% arrange(desc(total)) 


all <- cbc %>%  select(ID, type, cbc) %>% group_by(ID) %>% summarise(n=n()) 
a <- cbc %>%  inner_join(all, by="ID") %>% select(ID, cbc, n) %>% cast(ID~cbc) %>% as_tibble()
a[is.na(a)] <- 0 
a[,2:26] <- ifelse(a[,2:26]>0,1,0)

#write.xlsx(a, "cbc_cor.xlsx", sheetName = "all", append = T)

##CBC all
all1 <- a %>% mutate(ALC=ifelse(ALC0==1|ALC1==1|ALC2==1|ALC3==1|min_ALC==1,1,0),
ANC=ifelse(ANC0==1|ANC1==1|ANC2==1|ANC3==1|min_ANC==1,1,0),
Hb=ifelse(Hb0==1|Hb1==1|Hb2==1|Hb3==1|min_Hb==1,1,0),
Mo=ifelse(Mo0==1|Mo1==1|Mo2==1|Mo3==1|min_Mo==1,1,0),
PLT=ifelse(PLT0==1|PLT1==1|PLT2==1|PLT3==1|min_PLT==1,1,0)) %>% 
select(ID, Hb, PLT, ANC, Mo, ALC)
all1[,2:6]<- ifelse(all1[,2:6]>0, all1$ID, 0)
x <- list(ALC=all1$ALC, ANC=all1$ANC, Hb=all1$Hb, Mo=all1$Mo, PLT=all1$PLT)
all <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle("Genes related to at least 1 \n among  CBC0, CBC1, CBC2, CBC3, and min_CBC")+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))
all1[,2:6]<- ifelse(all1[,2:6]==0, 0, 1)
colSums(all1[,2:6]) %>% as_tibble() %>% 
mutate(pct=round(value/sum(value)*100,1), ID=c("Hb", "PLT", "ANC", "Mo", "ALC"))

ggarrange(network, all, labels = c("A", "B"), nrow=1, ncol=2)


# CBC0

z1 <- cbc %>% filter(cbc=="Hb0"|cbc=="PLT0"|cbc=="ANC0"|cbc=="Mo0"|cbc=="ALC0")
z <-cbc %>% filter(cbc=="Hb0"|cbc=="PLT0"|cbc=="ANC0"|cbc=="Mo0"|cbc=="ALC0") %>% group_by(ID) %>% 
summarise(n=n()) 
z <-z %>% inner_join(z1, by="ID") %>% select(ID, cbc, n) %>% 
cast(ID~cbc) %>% as_tibble() 
z[is.na(z)]<-0
z[,2:6] <- ifelse(z[,2:6]>0,z$ID, 0)
x <- list(ALC0=z$ALC0, ANC0=z$ANC0, Hb0=z$Hb0, Mo0=z$Mo0, PLT0=z$PLT0)
ven0 <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle("Genes related to CBC0")+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))

#CBC1
z1 <- cbc %>% filter(cbc=="Hb1"|cbc=="PLT1"|cbc=="ANC1"|cbc=="Mo1"|cbc=="ALC1")
z <-cbc %>% filter(cbc=="Hb1"|cbc=="PLT1"|cbc=="ANC1"|cbc=="Mo1"|cbc=="ALC1") %>% group_by(ID) %>% 
    summarise(n=n()) 
library(reshape)
z <-z %>% inner_join(z1, by="ID") %>% select(ID, cbc, n) %>% 
    cast(ID~cbc) %>% as.tibble() 
z[is.na(z)]<-0
z[,2:6] <- ifelse(z[,2:6]>0,z$ID,0)
x <- list(ALC1=z$ALC1, ANC1=z$ANC1, Hb1=z$Hb1, Mo1=z$Mo1, PLT1=z$PLT1)
ven1 <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle(expression(bolditalic(paste("Genes related to CBC1(",sqrt(CBC0%*%min_CBC),")"))))+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))

#CBC2
z1 <- cbc %>% filter(cbc=="Hb2"|cbc=="PLT2"|cbc=="ANC2"|cbc=="Mo2"|cbc=="ALC2")
z <-cbc %>% filter(cbc=="Hb2"|cbc=="PLT2"|cbc=="ANC2"|cbc=="Mo2"|cbc=="ALC2") %>% group_by(ID) %>% 
    summarise(n=n()) 
z <-z %>% inner_join(z1, by="ID") %>% select(ID, cbc, n) %>% 
    cast(ID~cbc) %>% as.tibble() 
z[is.na(z)]<-0
z[,2:6] <- ifelse(z[,2:6]>0,z$ID,0)
x <- list(ALC2=z$ALC2, ANC2=z$ANC2, Hb2=z$Hb2, Mo2=z$Mo2, PLT2=z$PLT2)
ven2 <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle("Genes related to CBC2")+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))

#CBC3
z1 <- cbc %>% filter(cbc=="Hb3"|cbc=="PLT3"|cbc=="ANC3"|cbc=="Mo3"|cbc=="ALC3")
z <-cbc %>% filter(cbc=="Hb3"|cbc=="PLT3"|cbc=="ANC3"|cbc=="Mo3"|cbc=="ALC3") %>% group_by(ID) %>% 
    summarise(n=n()) 
library(reshape)
z <-z %>% inner_join(z1, by="ID") %>% select(ID, cbc, n) %>% 
    cast(ID~cbc) %>% as.tibble() 
z[is.na(z)]<-0
z[,2:6] <- ifelse(z[,2:6]>0,z$ID,0)
x <- list(ALC3=z$ALC3, ANC3=z$ANC3, Hb3=z$Hb3, Mo3=z$Mo3, PLT3=z$PLT3)
ven3 <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle(expression(bolditalic(paste("Genes related to CBC3 (log ", frac(CBC2,CBC0),")"))))+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))


#min_CBC
z1 <- cbc %>% filter(cbc=="min_Hb"|cbc=="min_PLT"|cbc=="min_ANC"|cbc=="min_Mo"|cbc=="min_ALC")
z <-cbc %>% filter(cbc=="min_Hb"|cbc=="min_PLT"|cbc=="min_ANC"|cbc=="min_Mo"|cbc=="min_ALC") %>% group_by(ID) %>% 
    summarise(n=n()) 
z <-z %>% inner_join(z1, by="ID") %>% select(ID, cbc, n) %>% 
    cast(ID~cbc) %>% as.tibble() 
z[is.na(z)]<-0
z[,2:6] <- ifelse(z[,2:6]>0,z$ID,0)
x <- list(min_ALC=z$min_ALC, min_ANC=z$min_ANC, min_Hb=z$min_Hb, min_Mo=z$min_Mo, min_PLT=z$min_PLT)
ven4 <- ggVennDiagram(x, label_alpha = 0.6)+scale_fill_gradient(low="white",high = "red")+ggtitle("Genes related to min_CBC")+theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"))


ggarrange(all, ven0,ven1,ven2,ven3,ven4, labels = c("A","B","C","D","E","F"), nrow=3, ncol=2, font.label = list(size = 35, color = "black"))



