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

lung <- lung %>% mutate(status = recode(status, "1" = 0, "2" = 1))
lung
s1 <- survfit(Surv(lung$time, lung$status) ~ 1, data = lung)

survfit(Surv(time, status) ~ 1, data = lung) %>% ggsurvfit() +
  labs(x = "Days", 
       y = "Overall survival probability") +
  add_confidence_interval()

summary(survfit(Surv(time, status) ~ 1, data = lung), times = 365.25)
