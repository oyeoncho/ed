library(readxl)
library(tidyverse)
library(lubridate)
library(edgeR)
library(lubridate)
library(Hmisc) # rcorr
library(leaps) # rugsubset
library(xlsx)
library(ggVennDiagram)
library(ggpubr)

### loading fc 50%
load(file='process_data/rt_fc1.Rdata')
miR <- rt_fc[["miR"]]; mR <- rt_fc[["mR"]]; lncR <- rt_fc[["lncR"]]

mR <- mR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="Gene_Symbol") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
miR <- miR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="Mature_ID") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID
lncR <- lncR %>% select(-contains("pvalue")) %>% 
    column_to_rownames(var="Gene_Symbol") %>% t %>% 
    as.data.frame() %>% rownames_to_column(var="ID") %>% 
    as_tibble() %>% mutate(ID=substr(ID, 1,8))  # key ID

test <- miR %>% inner_join(lncR, by="ID") %>% select(ID, `hsa-miR-574-3p`, LINC01003) %>%
inner_join(mR,by="ID")
mRa <- test %>% column_to_rownames(var="ID")
x<-select_if (mRa, is.numeric)
rx <- rcorr (as.matrix (x))$r
result <-ifelse (rcorr (as.matrix (x))$P<0.05, rx, 0)

a0 <- result %>% as_tibble() %>% mutate(ID = colnames(result)) %>% select(ID, ACOT9, `hsa-miR-574-3p`, LINC01003) %>%
filter(abs(ACOT9)> 0 & ID!="ACOT9") %>% arrange(desc(abs(ACOT9))) 
aa0 <- a0 %>% slice(1:100)

a1 <- result %>% as_tibble() %>% mutate(ID = colnames(result)) %>% select(ID, ACOT9, `hsa-miR-574-3p`, LINC01003) %>%
filter(abs(`hsa-miR-574-3p`)> 0 & ID!="hsa-miR-574-3p") %>% arrange(desc(abs(`hsa-miR-574-3p`))) 
aa1 <- a1 %>% slice(1:100)

a2 <- result %>% as_tibble() %>% mutate(ID = colnames(result)) %>% select(ID, ACOT9, `hsa-miR-574-3p`, LINC01003) %>%
filter(abs(LINC01003)> 0 & ID!="LINC01003") %>% arrange(desc(abs(LINC01003))) 

aa2 <- a2 %>% slice(1:100)

x =list(ACOT9=a0$ID, `hsa-miR-574-3p`=a1$ID, LINC01003=a2$ID)
g1 <- ggVennDiagram(x, label_alpha = 0.6, label_size = 8)+scale_fill_gradient(low="white",high = "red")+ggtitle("219 mRNAs significantly correlated with\nall three RNAs")+theme(
    plot.title = element_text(color="red", size=25, face="bold.italic"))

x =list(ACOT9=aa0$ID, `hsa-miR-574-3p`=aa1$ID, LINC01003=aa2$ID)
g2 <-ggVennDiagram(x, label_alpha = 0.6, label_size = 8)+scale_fill_gradient(low="white",high = "red")+ggtitle("100 mRNAs highly correlated with\n each of 3 RNAs")+theme(
    plot.title = element_text(color="red", size=25, face="bold.italic"))

ggarrange(g1, g2, nrow=1, ncol=2)

# lncr & mr & mir
c1 <- result %>% as_tibble() %>% mutate(ID = colnames(result)) %>% select(ID, ACOT9, `hsa-miR-574-3p`, LINC01003) %>%
filter(abs(LINC01003) > 0 & abs(`hsa-miR-574-3p`) > 0 & abs(ACOT9)>0) %>% select(ID)

write.xlsx(aa0, "process_data/function.xlsx", sheetName = "mR", append = T)
write.xlsx(aa1, "process_data/function.xlsx", sheetName = "miR", append = T)
write.xlsx(aa2, "process_data/function.xlsx", sheetName = "lncR", append = T)
write.xlsx(c1, "process_data/function.xlsx", sheetName = "mR_lncR_miR", append = T)
