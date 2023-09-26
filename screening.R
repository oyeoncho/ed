library(readxl)
library(tidyverse)
library(lubridate)
library(ggpubr)
library(xlsx)

## 460
Cx <- read_excel("raw_data/Cx2023.xlsx", sheet = "all") %>% 
    mutate(start_date= as.Date(start_date, origin = "1900-01-01"), end_date= as.Date(end_date, origin = "1900-01-01"))

## Cohort 2
Cx_e <- Cx %>% filter(is.na(SID)==F)  %>%
    mutate(stage1=ifelse(stage=="2A"|stage=="2B"|stage=="3A"|stage=="3B"|stage=="3C1", 0,1), 
           RT_field=ifelse(`RT field`=="L5",0,1))

# stage 2A-3C1 vs 3C2-4A (1B/4B remove) & remove cohort2
Cx <- Cx %>% filter(is.na(SID)==T)
Cx <- Cx %>% filter(stage!="1B" & stage!="4B" & EQD2> 50) %>% 
    mutate(stage1=ifelse(stage=="2A"|stage=="2B"|stage=="3A"|stage=="3B"|stage=="3C1", 0,1))

################# group1 : progression within 12 months and dead within 15 months after end of treatment 
Cx %>% filter(css==1) %>% mutate(z=ifelse(recur_date <=12, "Progression (<= 12 months)","Progression (> 12 months)")) %>% ggplot(aes(z, fu_date))+geom_boxplot()+
    scale_y_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,120,6))+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")+
    ylab("Death date (month)")+xlab("")

Cx %>% filter(css==1) %>% ggplot(aes(recur_date, fu_date))+geom_point(size=3)+
    scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,120,6))+
    scale_y_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,120,6))+
    xlab("Progression date (month)")+ylab("Death date (month)")+coord_fixed()+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")+
    geom_hline(aes(yintercept=30), color="blue", linetype="dashed", linewidth=1)+
    geom_hline(aes(yintercept=15), color="red", linetype="dashed", linewidth=1)

#################group2 : progression but survive at least more than 30 months after end of treatment 
Cx %>% filter(css==1) %>% ggplot(aes(fu_date))+geom_density()+
    scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,100,6))+
    geom_vline(aes(xintercept=30), color="blue", linetype="dashed", linewidth=1)+
    geom_vline(aes(xintercept=15), color="red", linetype="dashed", linewidth=1)+
    xlab("Death date (month)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")

################### group 3 : no progression after at least more than 60months after end of treatment
Cx %>% filter(recur1>0) %>% ggplot(aes(recur_date))+geom_density()+
    scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,100,6))+
    geom_vline(aes(xintercept=60), color="black", linetype="dashed", linewidth=1)+
    geom_vline(aes(xintercept=12), color="black", linetype="solid", linewidth=1)+
    xlab("Progression date (month)")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")

#  definition, death =< 15mo, recur but survival > 30 mo, > 60 mo fu recur x
Cx1 <- Cx %>% filter(is.na(SID)) %>% mutate(group=ifelse(recur1==0 & recur_date >60, 3, ifelse(css==1 & fu_date <= 15, 1, ifelse(recur1>0 & fu_date >30, 2, NA))), 
                     RT_field=ifelse(`RT field`=="L5",0,1)) %>% filter(is.na(group)==F) 

library(moonBook)
# n=174
Cx1$z <- ifelse(Cx1$recur_date<=12,1,0)
mytable(group ~ Age+pathology+stage1+RT_field+duration+EQD2+recur_date+recur1+fu_date+css, data=Cx1)

#cohort 2
mytable(~ Age+pathology+stage1+RT_field+duration+EQD2+recur_date+recur1+fu_date+css, data=Cx_e)

## PSM, optimal 
library(MatchIt)
# 1 vs 3 (E vs no recur)
Cx2 <- Cx1 %>% filter(group != 2) %>% mutate(group=ifelse(group==3, 0, 1), path = ifelse(pathology=="sqcc",0,1))
match_fit <- matchit(group ~ path+stage1+RT_field, method="optimal",distance = "glm", data=Cx2)
matched_dat <- match.data(match_fit)
mytable(group ~ Age+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration, data=matched_dat,method=3,catMethod=0)
a11 <- matched_dat 
summary(match_fit)
plot(match_fit, type = "jitter", interactive = FALSE)


# 1 vs 2 (E vs recrur)
Cx2 <- Cx1 %>% filter(group != 3) %>% mutate(group=ifelse(group==2, 0, 1), path = ifelse(pathology=="sqcc",0,1))
match_fit <- matchit(group ~ path+stage1+RT_field, method="optimal",distance = "glm", data=Cx2)
matched_dat <- match.data(match_fit)
mytable(group ~ Age+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration, data=matched_dat,method=3, catMethod=0)
a22 <- matched_dat 
plot(match_fit, type = "jitter", interactive = FALSE)


# 2 vs 3 (recur vs no recur)
Cx2 <- Cx1 %>% filter(group != 1) %>% mutate(group=ifelse(group==3, 0, 1), path = ifelse(pathology=="sqcc",0,1))
match_fit <- matchit(group ~ path+stage1+RT_field, method="optimal",distance = "glm", data=Cx2)
matched_dat <- match.data(match_fit)
mytable(group ~ Age+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration, data=matched_dat,method=3, catMethod=0)
a33 <- matched_dat 
plot(match_fit, type = "jitter", interactive = FALSE)

### all date used in marching
all <- a11 %>% bind_rows(a22) %>% bind_rows(a33) %>% distinct(ID, .keep_all = T)  


# CBC data, PSM comparison
CBC_a <- read_excel("raw_data/CBC.xlsx", sheet = "CBC_a") # CBC for cohort 1
CBC_e <- read_excel("raw_data/CBC.xlsx", sheet = "CBC_e") # CBC for cohort 2

# pre CBC define (match)
c0 <- c()
for (i in 1:nrow(all)){ 
    m0 <- CBC_a %>% filter(ID==all$ID[i]) %>% filter(diff<1 & is.na(NLR)==F) %>% arrange(desc(diff)) %>% slice(1) %>% 
        select(ID, Hb0=Hb, PLT0=PLT, ANC0=ANC, ALC0=ALC, Mo0=Mo, NLR0=NLR, PLR0=PLR, LMR0=LMR)
    c0 <- rbind(c0, m0)
}

x <- c0 %>% right_join(all, by="ID") %>% filter(is.na(Hb0)) %>% select(ID)

for (i in 1:nrow(x)) {
    m0 <- CBC_a %>% filter(ID==x$ID[i]) %>% filter(diff<5 & is.na(NLR)==F) %>% arrange(desc(diff)) %>% slice(1) %>% 
        select(ID, Hb0=Hb, PLT0=PLT, ANC0=ANC, ALC0=ALC, Mo0=Mo, NLR0=NLR, PLR0=PLR, LMR0=LMR)
    c0 <- rbind(c0, m0)
}

# CBC2 define
c4 <- c()
for (i in 1:nrow(all)) {
    m4 <- CBC_a %>% filter(ID==all$ID[i]) %>% filter(diff>1 & is.na(ANC)==F) %>% mutate(f=abs(14-diff)) %>% arrange(f) %>% slice(1) %>%
           select(ID, Hb2=Hb, PLT2=PLT, ANC2=ANC, ALC2=ALC, Mo2=Mo)
    c4 <- rbind(c4, m4)
}

CBC_a %>% ggplot(aes(diff, PLT))+geom_point()+lims(x=c(0,100))

# CBC defined during RT (match)
c1 <- c();m <-c()
for (i in 1:nrow(all)) {
    m1 <- CBC_a %>% filter(ID==all$ID[i]) %>% filter(diff < 50 & diff>0 & is.na(NLR)==F) %>% arrange(diff) %>% 
        mutate(min_ANC=min(ANC), min_ALC=min(ALC), min_Mo=min(Mo)) %>% slice(1) %>% select(ID, min_ANC:min_Mo)
    m2 <- CBC_a %>% filter(ID==all$ID[i]) %>% filter(diff < 50 &diff>0 & is.na(Hb)==F) %>% arrange(diff) %>% 
        mutate(min_Hb=min(Hb)) %>% slice(1) %>%  select(min_Hb)
    m3 <- CBC_a %>% filter(ID==all$ID[i]) %>% filter(diff < 50 & diff>0 & is.na(PLT)==F) %>% arrange(diff) %>% 
        mutate(min_PLT=min(PLT)) %>% slice(1) %>%  select(min_PLT)
    m <- cbind(m1, m2,m3)
    c1 <- rbind(c1,m)
}


# pre CBC , CBC2 define (exosome)
e0 <- c() ; m <- c()
for (i in 1:nrow(Cx_e)){ 
    m0 <- CBC_e %>% filter(ID==Cx_e$ID[i]) %>% filter(diff<1 & is.na(NLR)==F) %>% arrange(desc(diff)) %>% slice(1) %>% 
        select(ID, Hb0=Hb, PLT0=PLT, ANC0=ANC, ALC0=ALC, Mo0=Mo, NLR0=NLR, PLR0=PLR, LMR0=LMR)
    m4 <- CBC_e %>% filter(ID==Cx_e$ID[i]) %>% mutate(f=abs(14-diff)) %>% arrange(f) %>% slice(1) %>%
           select(Hb2=Hb, PLT2=PLT, ANC2=ANC, ALC2=ALC, Mo2=Mo)
    m <- cbind(m0,m4)
    e0 <- rbind(e0,m)
}



# CBC defined during RT (exosome)
e1 <- c(); m <- c()
for (i in 1:nrow(Cx_e)) {
    m1 <- CBC_e %>% filter(ID==Cx_e$ID[i]) %>% filter(diff < 50 & diff>0 & is.na(NLR)==F) %>% arrange(diff) %>% 
        mutate(min_ANC=min(ANC), min_ALC=min(ALC), min_Mo=min(Mo)) %>% slice(1) %>% select(ID, min_ANC:min_Mo)
    m2 <- CBC_e %>% filter(ID==Cx_e$ID[i]) %>% filter(diff < 50 &diff>0 & is.na(Hb)==F) %>% arrange(diff) %>% 
        mutate(min_Hb=min(Hb)) %>% slice(1) %>%  select(min_Hb)
    m3 <- CBC_e %>% filter(ID==Cx_e$ID[i]) %>% filter(diff < 50 & diff>0 & is.na(PLT)==F) %>% arrange(diff) %>% 
        mutate(min_PLT=min(PLT)) %>% slice(1) %>%  select(min_PLT)
    m <- cbind(m1, m2,m3)
    e1 <- rbind(e1,m)
}

# join CBC (cohort1)
all1 <- all %>% left_join(c0, by="ID") %>%  select(ID, Hb0:LMR0)%>% 
    inner_join(c1, by="ID") %>% inner_join(c4, by="ID")

# join CBC (cohort2)
all_e <- Cx_e %>% inner_join(e0, by="ID") %>% select(ID, Hb0:Mo2) %>% inner_join(e1, by="ID")



### coreelation CBC0,CBC1, CBC2, min_CBC ~CBC3
library(corrplot)
library(colorRamp2)
library(Hmisc) 
cor <- all  %>% inner_join(all1, by="ID") %>%   mutate(Hb1=(Hb0*min_Hb)^0.5, ANC1=(ANC0*min_ANC)^0.5, PLT1=(PLT0*min_PLT)^0.5, ALC1=(ALC0*min_ALC)^0.5, Mo1=(Mo0*min_Mo)^0.5, 
                                                NLR1=ANC1/ALC1, PLR1=PLT1/ALC1, LMR1=ALC1/Mo1,
                                            Hb3=log(Hb2/Hb0), PLT3=log(PLT2/PLT0), ALC3=log(ALC2/ALC0), ANC3=log(ANC2/ANC0), Mo3=log(Mo2/Mo0), css1=ifelse(css==1 & fu_date<=15,1,0)) %>%
    select(css1, ANC0, PLT0, Hb0, ALC0, Mo0,  ANC1, PLT1, Hb1, ALC1, Mo1, 
    min_ANC, min_PLT, min_Hb, min_ALC, min_Mo, ANC3, PLT3, Hb3, ALC3, Mo3)
x<-select_if (cor, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
col1 <- colorRampPalette(c("blue", "white", "red"))
corrplot(result, type="upper", addrect = 2, col=col1(100))
summary(lm(Hb3~Hb0, data=cor))

q <- cor %>% mutate(ED=ifelse(css1==1, "yes", "no")) %>% ggplot(aes(Hb0, Hb3, color=ED))+geom_point(size=3)+
    geom_abline(intercept= 0.426151, slope=-0.041347, color='black', linewidth = 1.5)+
    annotate("text",10,0.6,label=expression(paste(R^2,"=0.2676")), size=10)+theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "top")+
    ylab("Hb3")+xlab("Hb0 (g/dL)")+ 
    scale_y_continuous(breaks=seq(-0.8,1,0.1))+
    scale_x_continuous(breaks=seq(5,15,1))+
    ggtitle("Correlation between Hb0 and HB3 in 120 patients")+
    geom_hline(yintercept = median(cor$Hb3), linewidth=1.5, lty=2)+
    annotate("text",6,0, label="Median Hb3",size=6)


w <- cor %>% mutate(ED=ifelse(css1==1, "yes", "no")) %>% 
    ggboxplot(x="ED",y="Hb3", palette = "jco", add="jitter", color="ED")+stat_compare_means(label="p.signif", label.x=1.5, label.y=0.4, cex=10)+
    ylab("Hb3")+xlab("Early death")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")+
    scale_y_continuous(breaks=seq(-0.5,0.8,0.2))


### table 1 PSM CBC comparision

#  ED (G1) vs no progression (G3)
data <- a11 %>% inner_join(all1, by="ID") %>% 
    mutate(Hb1=(Hb0*min_Hb)^0.5, ANC1=(ANC0*min_ANC)^0.5, PLT1=(PLT0*min_PLT)^0.5, ALC1=(ALC0*min_ALC)^0.5, Mo1=(Mo0*min_Mo)^0.5, 
           NLR1=ANC1/ALC1, PLR1=PLT1/ALC1, LMR1=ALC1/Mo1,
           Hb3=log(Hb2/Hb0), PLT3=log(PLT2/PLT0), ALC3=log(ALC2/ALC0), ANC3=log(ANC2/ANC0), Mo3=log(Mo2/Mo0)) %>% as.data.frame()

com1 <- mytable(group ~ Age+pathology+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration+
Hb0+min_Hb+Hb1+Hb2+Hb3+PLT0+min_PLT+PLT1+PLT2+PLT3+Mo0+min_Mo+Mo1+Mo2+Mo3+ANC0+min_ANC+ANC1+ANC2+ANC3+ALC0+min_ALC+ALC1+ALC2+ALC3+
NLR0+NLR1+PLR0+PLR1+LMR0+LMR1+alpha, data=data ,method=3,catMethod=0, digits=2)

com1

## ED (G1) vs progression (G2)
data <- a22 %>% inner_join(all1, by="ID") %>% 
    mutate(Hb1=(Hb0*min_Hb)^0.5, ANC1=(ANC0*min_ANC)^0.5, PLT1=(PLT0*min_PLT)^0.5, ALC1=(ALC0*min_ALC)^0.5, Mo1=(Mo0*min_Mo)^0.5, 
           NLR1=ANC1/ALC1, PLR1=PLT1/ALC1, LMR1=ALC1/Mo1,
           Hb3=log(Hb2/Hb0), PLT3=log(PLT2/PLT0), ALC3=log(ALC2/ALC0), ANC3=log(ANC2/ANC0), Mo3=log(Mo2/Mo0)) %>% as.data.frame()

com2 <- mytable(group ~ Age+pathology+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration+
Hb0+min_Hb+Hb1+Hb2+Hb3+PLT0+min_PLT+PLT1+PLT2+PLT3+Mo0+min_Mo+Mo1+Mo2+Mo3+ANC0+min_ANC+ANC1+ANC2+ANC3+ALC0+min_ALC+ALC1+ALC2+ALC3+
NLR0+NLR1+PLR0+PLR1+LMR0+LMR1+alpha, data=data ,method=3,catMethod=0, digits=2)

com2

## progression (G2) vs non-progression (G3)
data <- a33 %>% inner_join(all1, by="ID") %>% 
    mutate(Hb1=(Hb0*min_Hb)^0.5, ANC1=(ANC0*min_ANC)^0.5, PLT1=(PLT0*min_PLT)^0.5, ALC1=(ALC0*min_ALC)^0.5, Mo1=(Mo0*min_Mo)^0.5, 
           NLR1=ANC1/ALC1, PLR1=PLT1/ALC1, LMR1=ALC1/Mo1,
           Hb3=log(Hb2/Hb0), PLT3=log(PLT2/PLT0), ALC3=log(ALC2/ALC0), ANC3=log(ANC2/ANC0), Mo3=log(Mo2/Mo0)) %>% as.data.frame()
com3 <- mytable(group ~ Age+pathology+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration+
Hb0+min_Hb+Hb1+Hb2+Hb3+PLT0+min_PLT+PLT1+PLT2+PLT3+Mo0+min_Mo+Mo1+Mo2+Mo3+ANC0+min_ANC+ANC1+ANC2+ANC3+ALC0+min_ALC+ALC1+ALC2+ALC3+
NLR0+NLR1+PLR0+PLR1+LMR0+LMR1+alpha, data=data ,method=3,catMethod=0, digits=2)

com3

## cohort 2 PSM comarison

Cx_e1 <- Cx_e %>% inner_join(all_e, by="ID")  %>% 
    mutate(Hb1=(Hb0*min_Hb)^0.5,  ANC1=(ANC0*min_ANC)^0.5, PLT1=(PLT0*min_PLT)^0.5, ALC1=(ALC0*min_ALC)^0.5, Mo1=(Mo0*min_Mo)^0.5, 
           NLR1=ANC1/ALC1, PLR1=PLT1/ALC1, LMR1=ALC1/Mo1,
           Hb3=log(Hb2/Hb0), PLT3=log(PLT2/PLT0), ALC3=log(ALC2/ALC0), ANC3=log(ANC2/ANC0), Mo3=log(Mo2/Mo0),
           group1=ifelse(recur1>1 & recur_date<=12,1,0), path=ifelse(pathology=="sqcc",0,1))


# patients who did not ensure either ED or non-ED 
z <- Cx_e1 %>% filter(ID!="141eaa12" & ID!="d9af1daa")

match_fit <- matchit(css ~ Age+path+stage1, method="optimal",distance = "glm",  data=z, ratio=2)
matched_dat <- match.data(match_fit)
ee <- mytable(css ~ Age+pathology+path+stage1+RT_field+SCC0+SCC2+Cyfra0+Cyfra2+EQD2+duration+
            Hb0+min_Hb+Hb1+Hb2+Hb3+PLT0+min_PLT+PLT1+PLT2+PLT3+Mo0+min_Mo+Mo1+Mo2+Mo3+ANC0+min_ANC+ANC1+ANC2+ANC3+ALC0+min_ALC+ALC1+ALC2+ALC3+
            NLR0+NLR1+PLR0+PLR1+LMR0+LMR1+alpha, data=matched_dat ,method=3,catMethod=0, digits=2)
ee

# corplot CBC for cohort2
cor <- Cx_e1 %>% select(css,ANC0, PLT0, Hb0, ALC0, Mo0,  ANC1, PLT1, Hb1, ALC1, Mo1, min_ANC, min_PLT, min_Hb, min_ALC, min_Mo, ANC3, PLT3, Hb3, ALC3, Mo3)
x<-select_if (cor, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)
col1 <- colorRampPalette(c("blue", "white", "red"))
corrplot(result, type="upper", addrect = 2, col=col1(100))

e <-  cor %>% mutate(ED=ifelse(css==1, "yes", "no")) %>% ggplot(aes(Hb0, Hb3, color=ED))+geom_point(size=3)+
    geom_abline(intercept= 0.48428, slope=-0.0457, color='black', linewidth = 1.5)+
    annotate("text",12,0.2,label=expression(paste(R^2,"=0.3423")), size=10)+theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "top")+
    ylab("Hb3")+xlab("Hb0 (g/dL)")+
    scale_y_continuous(breaks=seq(-0.8,1,0.1))+
    scale_x_continuous(breaks=seq(5,15,1))+
    ggtitle("Correlation between Hb0 and HB3 in 42 patients (exosome)")+
    geom_hline(yintercept = median(cor$Hb3), linewidth=1.5, lty=2)+
    annotate("text",9,-0.05, label="Median Hb3",size=6)

r <- cor %>% mutate(ED=ifelse(css==1, "yes", "no")) %>% 
    ggboxplot(x="ED",y="Hb3", palette = "jco", add="jitter", color="ED")+stat_compare_means(label="p.signif", label.x=1.5, label.y=0.15, cex=10)+
    ylab("Hb3")+xlab("Early death")+
    theme(axis.title =element_text(face="bold", size=20), axis.text = element_text(size=20, face="bold"), legend.position = "NA")+
    scale_y_continuous(breaks=seq(-0.5,0.5,0.2))

ggarrange(q,e,w,r, labels=c("A","B","C","D"), nrow=2, ncol=2, font.label = list(size = 30, color = "black"))

write.xlsx(Cx_e1, "process_data/group.xlsx", sheetName = "Cx_e", append = T)


