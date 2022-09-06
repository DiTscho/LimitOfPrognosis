library("tidyverse")
library("survival")
library("survminer")
source("00_0.R")

df = readRDS("../database/platform.rds")


### tmp ####

df$expression_data[[3]] = filter_and_annotate_mb(df$expression_data[[3]])
df$expression_data[[3]][1:5,1:5]


n = 3 #METABRIC
E = df$expression_data[[n]]
C = df$clinical_data[[n]]



MB = cbind.data.frame(C$tumor_size, 
                      C$Lymph.nodes.examined.positive, 
                      E$ESR1, E$PGR, 
                      C$Overall.Survival..Months./12, 
                      C$Patient.s.Vital.Status)
colnames(MB) = c("tumorsize", "lymph_nodes", "ESR1", "PGR", "survival_time", "status")

MB = MB %>%
  filter(status == "Died of Disease" | status == "Living",
         is.na(lymph_nodes) == FALSE) %>% 
  mutate(Q = PGR/ESR1,
         status = ifelse(status == "Living" | survival_time >= 10, 0, 1),
         node_status = case_when(
           lymph_nodes == 0 ~ 1,
           lymph_nodes <= 3 ~ 2,
           lymph_nodes >= 4 ~ 3
         ))

m = coxph(Surv(survival_time, status) ~ tumorsize + lymph_nodes + PGR, data = MB)
m



