# Compiled Files of Charcoal Process Point Model using Big Woods Lily Lake
# Aidan Draper

rm(list=ls())

library(readxl)
require(rbacon)

setwd("~/Documents/Junior_Year/DISC_REU/DISC_bayesian_model/")

lily.lake <- read_excel("bigwoods.xls", sheet = 14, col_names = TRUE)

lake.dat <- data.frame(age=rep(NA,length(lily.lake$Date)))

lake.dat$age <- with(lily.lake, round(Date) - 2018) # calculate YBP
lake.dat$sed.rate <- with(lily.lake, Depth/5) # NOTE: do not have age of sediment core so we left it fixed
lake.dat$influx <- lily.lake$`Char Flux`
lake.dat$count <- lily.lake$Count

# NOTE: we cannot compute any of these without volume or sediment core form time
#lake.dat$char <- with(lily.lake, )
#lake.dat$offset



