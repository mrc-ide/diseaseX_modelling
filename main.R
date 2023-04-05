
library(rstudioapi)
library(squire)
library(nimue)
library(squire.page) 
library(data.table)
library(ggplot2)
library(parallel)
library(tidyverse)
library(countrycode)
library(furrr)
library(zoo)
library(tibble)
library(dplyr)
library(purrr)
library(furrr)
library(tidyr)
library(wesanderson)
library(patchwork)

setwd(dirname(getActiveDocumentContext()$path)) 

# colours

#cols1 <- c("#08519c", "#440154ff", "#9ecae1","#993366","#e6b3cc")
cols1 <- c("#08519c", "#440154ff", "#9ecae1","#993366","#333333")

cols2 <- c( "#993366","#440154ff", "#e6b3cc")
cols3 <- c("#8c2d04", "#fe9929", "#fed976")

lightgreen <- "#bbdf27"
teal <- "#25848e"
darkpurple <- "#440154"