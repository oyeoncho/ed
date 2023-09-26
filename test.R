library(readxl)
library(tidyverse)
library(reshape)
library(ggvenn)
library(ggpubr)
library(reshape)
library(ggVennDiagram)
library(xlsx)
library(reshape2)
library(ggpubr)

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


#CBC0~CBC3 or CBC1 ~CBC3 

#122
PLT <- cbc %>% filter(PLT0 !=0 & PLT3!=0) %>% select(ID, css, PLT0, PLT3) %>% distinct %>% filter((PLT0 <0 &PLT3>0) | (PLT0 >0 &PLT3<0)) %>%
    ggplot(aes(PLT0, PLT3)) + geom_point(size=3)+ xlab("PLT0-RNA correlation")+ylab("PLT3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("PLT0-PLT3 (122 RNAs)")

#8

PLT1<- cbc %>% filter(PLT0==0 & PLT1 !=0 & PLT3!=0) %>% select(ID, css, PLT1, PLT3) %>% distinct %>% filter((PLT1 <0 &PLT3>0) | (PLT1 >0 &PLT3<0)) %>% 
    ggplot(aes(PLT1, PLT3)) + geom_point(size=3)+ xlab("PLT1-RNA correlation")+ylab("PLT3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("PLT1-PLT3 (8 RNAs)")+
    lims(x=c(-0.4,0.4), y=c(-0.4,0.4))

#238
ANC <- cbc %>% filter(ANC0 !=0 & ANC3!=0) %>% select(ID, css, ANC0, ANC3) %>% distinct %>% filter((ANC0 <0 &ANC3>0) | (ANC0 >0 &ANC3<0)) %>%
    ggplot(aes(ANC0, ANC3)) + geom_point(size=3)+ xlab("ANC0-RNA correlation")+ylab("ANC3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("ANC0-ANC3 (238 RNAs)")

#24
ANC1 <- cbc %>% filter(ANC0==0 & ANC1 !=0 & ANC3!=0) %>% select(ID, css, ANC1, ANC3) %>% distinct %>% filter((ANC1 <1 &ANC3>0) | (ANC1 >0 &ANC3<0)) %>% 
    ggplot(aes(ANC1, ANC3)) + geom_point(size=3)+ xlab("ANC1-RNA correlation")+ylab("ANC3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("ANC1-ANC3 (24 RNAs)")


#2067
Hb <- cbc %>% filter(Hb0 !=0 & Hb3!=0) %>% select(ID, css, Hb0, Hb3) %>% distinct %>%  filter((Hb0 <0 &Hb3>0) | (Hb0 >0 &Hb3<0)) %>%
    ggplot(aes(Hb0, Hb3)) + geom_point(size=3)+ xlab("Hb0-RNA correlation")+ylab("Hb3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("Hb0-Hb3 (2067 RNAs)")

#20
Hb1 <- cbc %>% filter(Hb0==0 & Hb1 !=0 & Hb3!=0)  %>% select(ID, css, Hb1, Hb3) %>% distinct %>% filter((Hb1 <0 &Hb3>0) | (Hb1 >0 &Hb3<0)) %>%
    ggplot(aes(Hb1, Hb3)) + geom_point(size=3)+ xlab("Hb1-RNA correlation")+ylab("Hb3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("Hb1-Hb3 (20 RNAs)")

#145
ALC <- cbc %>% filter(ALC0 !=0 & ALC3!=0) %>% select(ID, css, ALC0, ALC3) %>% distinct %>% filter((ALC0 <0 &ALC3>0) | (ALC0 >0 &ALC3<0)) %>%
    ggplot(aes(ALC0, ALC3)) + geom_point(size=3)+ xlab("ALC0-RNA correlation")+ylab("ALC3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("ALC0-ALC3 (145 RNAs)")

#0
ALC1<- cbc %>% filter(ALC0 ==0 & ALC1!=0 & ALC3!=0) %>% select(ID, css, ALC1, ALC3) %>% distinct %>% filter((ALC1 <0 &ALC3>0) | (ALC1 >0 &ALC3<0)) %>%
    ggplot(aes(ALC1, ALC3)) + geom_point(size=3)+ xlab("ALC1-RNA correlation")+ylab("ALC3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("ALC1-ALC3 (0 RNAs)")

#185
Mo <- cbc %>% filter(Mo0 !=0 & Mo3!=0) %>% select(ID, css, Mo0, Mo3) %>% distinct %>% filter((Mo0 <0 &Mo3>0) | (Mo0 >0 &Mo3<0)) %>%
    ggplot(aes(Mo0, Mo3)) + geom_point(size=3)+ xlab("Mo0-RNA correlation")+ylab("Mo3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("Mo0-Mo3 (185 RNAs)")

#16
Mo1 <- cbc %>% filter(Mo0==0 & Mo1 !=0 & Mo3!=0) %>% select(ID, css, Mo1, Mo3) %>% distinct %>% filter((Mo1 <0 &Mo3>0) | (Mo1 >0 &Mo3<0)) %>%
    ggplot(aes(Mo1, Mo3)) + geom_point(size=3)+ xlab("Mo1-RNA correlation")+ylab("Mo3-RNA correlation")+
    theme(axis.text = element_text(size=15,face='bold'), axis.title = element_text(size=20,face='bold'), title=element_text(size=20, face="bold"))+ggtitle("Mo1-Mo3 (16 RNAs)")
###

ggarrange(PLT, PLT1,ANC,ANC1, Hb, Hb1, ALC, ALC1, Mo, Mo1,  labels=c("A","B","C","D","E","F","G","H", "I","J"), nrow=5, ncol=2, font.label = list(size = 30, color = "black"))

### in css +, CBC0(1)-CBC3
PLT <- cbc %>% filter(PLT0 !=0 & PLT3!=0 & css !=0) %>% select(ID, css, PLT0, PLT3) %>% distinct %>% filter( (css>0 & PLT0>0)|(css<0&PLT0<0)) %>% 
    melt(id=c("ID", "css")) %>%   ggplot(aes(value, css, color=variable)) +geom_point(size=5)+
    ggtitle("RNAs associated with PLT0, PLT3, and survival (6 RNAs)")+ xlab("RNA-PLT0 or RNA-PLT3 correlation")+ ylab("RNA-survival correlation")+
    theme(axis.text = element_text(size=12,face='bold'), axis.title = element_text(size=12,face='bold'), title=element_text(size=12, face="bold"), 
          legend.text = element_text(size=12, face='bold'), legend.title = element_blank(), legend.position = "top")

#6
PLT1 <- cbc %>% filter(PLT0 !=0 & PLT3!=0 & css !=0) %>% select(ID, css, PLT0, PLT3, type) %>% distinct %>% filter( (css>0 & PLT0>0)|(css<0&PLT0<0)) %>% arrange(css) %>%
mutate(cbc="PLT") %>% dplyr::rename("CBC0"="PLT0", "CBC3"="PLT3") 

ANC <- cbc %>% filter(ANC0 !=0 & ANC3!=0 & css !=0) %>% select(ID, css, ANC0, ANC3) %>% distinct %>% filter( (css>0 & ANC0>0)|(css<0&ANC0<0)) %>% 
    melt(id=c("ID", "css")) %>%  ggplot(aes(value, css, color=variable)) +geom_point(size=5)+
    ggtitle("RNAs associated with ANC0, ANC3, and survival (8 RNAs)")+ xlab("RNA-ANC0 or RNA-ANC3 correlation")+ ylab("RNA-survival correlation")+
    theme(axis.text = element_text(size=12,face='bold'), axis.title = element_text(size=12,face='bold'), title=element_text(size=12, face="bold"), 
          legend.text = element_text(size=12, face='bold'), legend.title = element_blank(), legend.position = "top")

#8
ANC1 <- cbc %>% filter(ANC0 !=0 & ANC3!=0 & css !=0) %>% select(ID, css, ANC0, ANC3, type) %>% distinct %>% filter( (css>0 & ANC0>0)|(css<0&ANC0<0)) %>% arrange(css) %>%
mutate(cbc="ANC") %>% dplyr::rename("CBC0"="ANC0", "CBC3"="ANC3") 

#2
ANC2<-cbc %>% filter(ANC0==0 & ANC1 !=0 & ANC3!=0 & css !=0) %>% select(ID, css, ANC1, ANC3, type) %>% distinct %>% filter( (css>0 & ANC1>0)|(css<0&ANC1<0)) %>% arrange(css) %>%
mutate(cbc="ANC") %>% dplyr::rename("CBC1"="ANC1", "CBC3"="ANC3") 

ANCa <- cbc %>% filter(ANC0==0 & ANC1 !=0 & ANC3!=0 & css !=0) %>% select(ID, css, ANC1, ANC3) %>% distinct %>% filter((css>0 & ANC1>0)|(css<0&ANC1<0))%>%
    melt(id=c("ID", "css")) %>%  ggplot(aes(value, css, color=variable)) +geom_point(size=5)+
    ggtitle("RNAs associated with ANC1, ANC3, and survival (2 RNAs)")+ xlab("RNA-ANC1 or RNA-ANC3 correlation")+ ylab("RNA-survival correlation")+
    theme(axis.text = element_text(size=12,face='bold'), axis.title = element_text(size=12,face='bold'), title=element_text(size=12, face="bold"), 
          legend.text = element_text(size=12, face='bold'), legend.title = element_blank(), legend.position = "top")+lims(x=c(-0.4,0.4), y=c(-0.4,0.4))


Hb <- cbc %>% filter(Hb0 !=0 & Hb3!=0 & css !=0) %>% select(ID, css, Hb0, Hb3) %>% distinct %>% 
     melt(id=c("ID", "css")) %>%  ggplot(aes(value, css, color=variable)) +geom_point(size=5)+
    ggtitle("RNAs associated with Hb0, Hb3, and survival (58 RNAs)")+ xlab("RNA-Hb0 or RNA-Hb3 correlation")+ ylab("RNA-survival correlation")+
    theme(axis.text = element_text(size=12,face='bold'), axis.title = element_text(size=12,face='bold'), title=element_text(size=12, face="bold"), 
          legend.text = element_text(size=12, face='bold'), legend.title = element_blank(), legend.position = "top")
#58
Hb1 <-cbc %>% filter(Hb0 !=0 & Hb3!=0 & css !=0) %>% select(ID, css, Hb0, Hb3, type) %>% distinct %>% arrange(css) %>%
mutate(cbc="Hb") %>% dplyr::rename("CBC0"="Hb0", "CBC3"="Hb3") 

Hba <- cbc %>% filter(Hb0==0 & Hb1 !=0 & Hb3!=0 & css !=0) %>% select(ID, css, Hb1, Hb3) %>% distinct %>% 
    melt(id=c("ID", "css")) %>%  ggplot(aes(value, css, color=variable)) +geom_point(size=5)+
    ggtitle("RNAs associated with Hb1, Hb3, and survival (2 RNAs)")+ xlab("RNA-Hb1 or RNA-Hb3 correlation")+ ylab("RNA-survival correlation")+
    theme(axis.text = element_text(size=12,face='bold'), axis.title = element_text(size=12,face='bold'), title=element_text(size=12, face="bold"), 
          legend.text = element_text(size=12, face='bold'), legend.title = element_blank(), legend.position = "top")+lims(x=c(-0.4,0.4), y=c(-0.4,0.4))

#2
Hb2<-cbc %>% filter(Hb0==0 & Hb1 !=0 & Hb3!=0 & css !=0) %>% select(ID, css, Hb1, Hb3, type) %>% distinct %>% arrange(css) %>%
mutate(cbc="Hb") %>% dplyr::rename("CBC1"="Hb1", "CBC3"="Hb3") 

Mo <- cbc %>% filter(Mo0 !=0 & Mo3!=0 & css !=0) %>% select(ID, css, Mo0, Mo3) %>% distinct %>%  
    melt(id=c("ID", "css")) %>%  ggplot(aes(value, css, color=variable)) +geom_point(size=5)+
    ggtitle("RNAs associated with Mo0, Mo3, and survival (5 RNAs)")+ xlab("RNA-Mo0 or RNA-Mo3 correlation")+ ylab("RNA-survival correlation")+
    theme(axis.text = element_text(size=12,face='bold'), axis.title = element_text(size=12,face='bold'), title=element_text(size=12, face="bold"), 
          legend.text = element_text(size=12, face='bold'), legend.title = element_blank(), legend.position = "top")
#5
Mo1 <-cbc %>% filter(Mo0 !=0 & Mo3!=0 & css !=0) %>% select(ID, css, Mo0, Mo3, type) %>% distinct %>% arrange(css) %>%
mutate(cbc="Mo") %>% dplyr::rename("CBC0"="Mo0", "CBC3"="Mo3") 

#0 Mo2
ggarrange(PLT, ANC, ANCa, Hb,Hba, Mo,  labels=c("A","B","C","D","E","F"), nrow=2, ncol=3, font.label = list(size = 30, color = "black"))

#81
cor <- PLT1 %>% bind_rows(ANC1, ANC2, Hb1, Hb2, Mo1)
c <- cor %>% distinct(ID)
cor %>% arrange(desc(abs(css))) %>% as.data.frame() %>% arrange(ID)

## cbc0-cbc3-css
cor %>% filter(type=="mR") %>% distinct(ID) #58
cor %>% filter(type=="miR") %>% distinct(ID) #7
cor %>% filter(type=="lncR") %>% distinct(ID) #4
cor %>% filter(type=="tR") %>% distinct(ID) #3
cor %>% filter(type=="snR") %>% distinct(ID) #6
cor %>% filter(type=="piR") %>% distinct(ID)#1

###
write.xlsx(cor, "process_data/reference.xlsx", sheetName = "cbc0_cbc3_css", append = T)

#### cbc0(1)-cbc3-css on network
central <- read_excel("process_data/group.xlsx", sheet = "network_css_05") 
central %>% inner_join(c, by="ID") %>% select(ID) 
central %>% filter(abs(Hb3)>0 | abs(PLT3)>0 | abs(Mo3)>0 |abs(ANC3)>0) %>% select(ID)

m <- central %>% select(ID) %>% inner_join(cbc, by="ID")
m %>% filter(PLT3 !=0 ) %>% select(ID, css,PLT0, PLT1, min_PLT, PLT2, PLT3) %>% distinct()
m %>% filter(ANC3 !=0 ) %>% select(ID, css,ANC0, ANC1, min_ANC, ANC2, ANC3) %>% distinct()
m %>% filter(Hb3 !=0 ) %>% select(ID, css, Hb0, Hb1, min_Hb, Hb2, Hb3) %>% distinct()
m %>% filter(Mo3 !=0 ) %>% select(ID, css, Mo0, Mo1, min_Mo, Mo2, Mo3) %>% distinct()


##########


z <- read_excel("process_data/reference.xlsx", sheet = "network_fc") %>% dplyr::rename("No"="...1")
Cx <- read_excel("process_data/group.xlsx", sheet = "Cx_e") 
z <- Cx %>% select(ID, pathology, Age, stage, fu_date, recur_date, recur1, Hb0:Mo3) %>%  inner_join(z, by="ID")

summary(lm(css~`hsa-miR-574-3p`+LINC01003+ACOT9,data=z))


library(survival)
library(survminer)
library(ggpubr)


q1 <- z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% arrange(desc(test)) %>%
    select(ID, Age, pathology, stage,css, fu_date,  recur1,  recur_date, `hsa-miR-574-3p`, LINC01003, ACOT9, test) %>%
    ggplot(aes(`hsa-miR-574-3p`, test))+ geom_point(size=3)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("miR-574-3p (log2FC)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"))+
    geom_abline(intercept= -0.3179, slope=1.2951, color='blue', size = 1.5)+
    annotate(geom="text", x=0.5, y=5, label=expression(R^2==0.1327), size=7, color="blue")

q2 <- z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% arrange(desc(test)) %>%
    select(ID, Age, pathology, stage,css, fu_date,  recur1,  recur_date, `hsa-miR-574-3p`, LINC01003, ACOT9, test) %>%
    ggplot(aes(LINC01003, test))+ geom_point(size=3)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("LINC01003 (log2FC)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"))+
    geom_abline(intercept= -0.2429, slope=-0.8889, color='blue', size = 1.5)+
    annotate(geom="text", x=0.5, y=5, label=expression(R^2==0.5865), size=7, color="blue")

q3 <- z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% arrange(desc(test)) %>%
    select(ID, Age, pathology, stage,css, fu_date,  recur1,  recur_date, `hsa-miR-574-3p`, LINC01003, ACOT9, test) %>%
    ggplot(aes(ACOT9, test))+ geom_point(size=3)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("ACOT9 (log2FC)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"))+
    geom_abline(intercept= -0.1570, slope=-0.8609, color='blue', size = 1.5)+
    annotate(geom="text", x=0.5, y=5, label=expression(R^2==0.2046), size=7, color="blue")

q4 <- z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% arrange(desc(test)) %>%
    select(ID, Age, pathology, stage,css, fu_date,  recur1,  recur_date, `hsa-miR-574-3p`, LINC01003, ACOT9, test) %>%
    ggplot(aes(LINC01003+ACOT9, test))+ geom_point(size=3)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("ACOT9+LINC01003 (log2FC)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"))+
    geom_abline(intercept= -0.008864, slope=-1.026713, color='blue', size = 1.5)+
    annotate(geom="text", x=0.5, y=5, label=expression(R^2==0.9215), size=7, color="blue")


q5 <- z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% arrange(desc(test)) %>%
    select(ID, Age, pathology, stage,css, fu_date,  recur1,  recur_date, `hsa-miR-574-3p`, LINC01003, ACOT9, test) %>%
    ggplot(aes(`hsa-miR-574-3p`-LINC01003, test))+ geom_point(size=3)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("miR-574-3p-LINC01003 (log2FC)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"))+
    geom_abline(intercept= -0.21933, slope=0.95203, color='blue', size = 1.5)+
    annotate(geom="text", x=0.5, y=5, label=expression(R^2==0.7257), size=7, color="blue")

q6 <-z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% arrange(desc(test)) %>%
    select(ID, Age, pathology, stage,css, fu_date,  recur1,  recur_date, `hsa-miR-574-3p`, LINC01003, ACOT9, test) %>%
    ggplot(aes(`hsa-miR-574-3p`-ACOT9, test))+ geom_point(size=3)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("miR-574-3p-ACOT9 (log2FC)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"))+
    geom_abline(intercept= -0.1549, slope=0.8049, color='blue', size = 1.5)+
    annotate(geom="text", x=0.5, y=5, label=expression(R^2==0.2738), size=7, color="blue")


ggarrange(q1, q2, q3, q4, q5, q6, labels=c("A","B", "C", "D","E","F"), nrow=3, ncol=2, font.label = list(size = 35, color = "black"))

z1 <- z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, 
                   test_m= ifelse(test>=0, 1, 0))

fit=survfit(Surv(fu_date, css==1)~test_m, data=z1)
summary(fit, time=30)
summary(z1$fu_date)

z %>% mutate(test=`hsa-miR-574-3p`-LINC01003-ACOT9, survival=ifelse(css==1,"death","survival")) %>% 
    ggboxplot(x="survival",y="test", palette = "jco", add="jitter", color="survival")+stat_compare_means(label="p.signif", label.x=1.5, label.y=10, cex=10)+
    ylab("miR-574-3p-LINC01003-ACOT9 (log2FC)")+xlab("")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")+ coord_cartesian(ylim=c(-8,10))+
    scale_y_continuous(breaks=seq(-10,12,2))


ggsurvplot(fit, xlab="Months", ylab="Disease specific survival rate", pval=T, pval.size=8,fun="pct", conf.int=F, risk.table = T, size=1, linetype = c("dashed","solid"),
           palette = c("blue","red"), legend= "none", legend.title="miR-574-3p-LINC01003\n-ACOT9(log2FC)", legend.lab=c("< 0",  "\u2265 0"),
           break.time.by=6, xlim=c(0,48), risk.table.height=0.3, ylim=c(0,100), ggtheme = theme_classic2(base_size=20, base_family = "Arial"),
           font.family = "Arial" )


fit=survfit(Surv(fu_date, css==1)~1, data=z1)
summary(fit, time=30)


