install.packages("ggsurvfit")
install.packages("tidycmprsk")
install.packages("gtsummary")

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(tidycmprsk)
library(survminer)

# Load dataset 
lung <- lung %>% mutate(status = recode(status, "1" = 0, "2" = 1))
lung
s1 <- survfit(Surv(lung$time, lung$status) ~ 1, data = lung)

survfit(Surv(time, status) ~ 1, data = lung) %>% ggsurvfit() +
  labs(x = "Days", 
       y = "Overall survival probability") +
  add_confidence_interval()

summary(survfit(Surv(time, status) ~ 1, data = lung), times = 365.25)


# Load dataset
lung <- as_tibble(lung)

# Survival object
s <- Surv(lung$time, lung$status)
class(s)

# Fit a survival curve
survfit(s ~ 1)
sfit <- survfit(Surv(time, status) ~ 1, data = lung)
summary(sfit)
ggsurvplot(sfit)

# Kaplan-Meier Plots
sfit <- survfit(Surv(time, status) ~ sex, data = lung)
summary(sfit)
plot(sfit)
ggsurvplot(sfit, conf.int = T, pval = T, risk.table = T, legend.labs = c("Male", "Female"),
           legend.title = "Sex", palette = c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for Lung Cancer Survival", 
           risk.table.height = 0.15)

# Colon cancer data
colon <- as_tibble(colon)
# Select death data
colondeath <- filter(colon, etype==2) 
# Male versus female
sfit <- survfit(Surv(time, status) ~ sex, data = colondeath)
summary(sfit, times = seq(0, 2000, 500))
ggsurvplot(sfit, conf.int = T, pval = T, risk.table = T, legend.labs = c("Male", "Female"), 
           legend.title = "Sex", palette = c("dodgerblue2", "orchid2"),
           title = "Kaplan-Meier Curve for Colon Cancer Survival", 
           risk.table.height = 0.15)

# Differentiation level
sfit <- survfit(Surv(time, status) ~ differ, data = colondeath)
ggsurvplot(sfit, conf.int = F, pval = T, risk.table = T, 
           legend.labs = c("differ=1", "differ=2", "differ=3"), 
           legend.title = "Differ", palette = c("red", "green", "blue"), 
           title = "Kaplan-Meier Curve for Colon Cancer Survival", 
           risk.table.height = 0.15)
# Lymph nodes
sfit <- survfit(Surv(time, status) ~ node4, data = colondeath)
ggsurvplot(sfit, conf.int = F, pval = T, risk.table = T, 
           legend.labs = c("node4=0", "node4=1"), 
           legend.title = "Node", palette = c("dodgerblue2", "orchid2"), 
           title = "Kaplan-Meier Curve for Colon Cancer Survival", 
           risk.table.height = 0.15)

# Cox regression
ggsurvplot(survfit(Surv(time, status) ~ nodes, data=colondeath))
sfit <- coxph(Surv(time, status) ~ sex + age + ph.ecog + ph.karno + pat.karno + meal.cal + wt.loss, 
              data=lung)

sfit <- survfit(Surv(time, status) ~ rx, data = colondeath)
ggsurvplot(sfit, conf.int = F, pval = T, risk.table = T)
coxfit <- coxph(Surv(time, status) ~ sex + age + nodes + rx, data = colondeath)
coxfit


