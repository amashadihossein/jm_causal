library(simcausal)
library(ggplot2)
library(ggdag)
require(visNetwork, quietly = TRUE)
library(patchwork)
library(dplyr)
library(survival)
library(survminer)
library(stringr)
# library(JM)
# library(nlme)
library(conflicted)

# pkg conflict preference
#------------------------
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("area", "patchwork")

source("./R/util.R")
y_spline <- yspline_get()