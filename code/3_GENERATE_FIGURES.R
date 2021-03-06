# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
# created by Ben Williams (Ben.Williams@alaska.gov);Nov 3, 2016; 2017-7-4
# Changes made by Sara Miller (Sara.Miller@alaska.gov); October 2020

LowerB <- 20000 #lower bound of recommended escapement goal range
UpperB <- 75000 #upper bound of recommended escapement goal range
SMSY <- 55863.05 #Lambert W version of SMSY from file stats.csv
lnalpha.c <-2.32010222657675 
lnalpha <- 1.7449322109841
beta <- 0.0000139956245496603 #let's try and grab these last 4 parameter outputs from the "stats.csv" output

#load----
library(plyr)
library(tidyverse)
library(reshape2)
library(grid)
library(gsl)
library(scales)
library(backports)
library(FField)
library(gg.gap)

devtools::install_github("commfish/fngr")
library(fngr)
library(extrafont)
loadfonts(device="win")
windowsFonts(Times=windowsFont("Times New Roman"))
theme_set(theme_report())

source('code/functions.r')

if(!dir.exists(file.path("output", "rjags_Full_Basecase", "processed"))){dir.create(file.path("output", "rjags_Full_Basecase", "processed"))}

# data----
spawnrecruitdat <- read.csv("data/Coghill_Sock.csv") #Load Data File (make sure this file is updated)
coda <- read.csv("output/rjags_Full_Basecase/coda.csv") 

# data clean----
# profile spawnrecruitdat
coda %>% 
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
         MSY.c = R.msy.c-S.msy.c, 
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

# analysis----
# create function for probability profiles and figures
profile(i=10, z=500, xa.start=0, xa.end=10000,lnalpha.c, beta) #can change i,z, xa.start, xa.end
#profile(i=10, z=500, xa.start=0, xa.end=1000,lnalpha.c, beta) #can change i,z, xa.start, xa.end
QM <- read.csv("output/rjags_Full_BaseCase/processed/QM.csv")
CI <- read.csv("output/rjags_Full_BaseCase/processed/CI.csv")
num <- nrow(QM)
QM %>%
  dplyr::select(c(escapement)) -> x
coda %>%
  dplyr::select(c(lnalpha.c, beta)) %>%
  filter(row_number()==1:50) %>%
  slice(rep(1:n(), each = num)) %>%
  cbind(., x) %>% #lnalpha.c, beta, and S 
  mutate(recruitment = escapement*exp(lnalpha.c-beta*escapement),
         variable = rep(1:50,each=num),
         year = "") -> dataset

QM %>%
  dplyr::select(c(escapement)) %>%
  mutate (lnalpha.c = lnalpha.c,
          beta = beta,
          recruitment = escapement*exp(lnalpha.c-beta*escapement),
          variable = 51,
          year = "")%>%
  rbind(., dataset) -> dataset

spawnrecruitdat %>%
  mutate(escapement = spawn,
         lnalpha.c = NA,
         beta = NA,
         recruitment =recruit,
         variable = 52) %>%
  dplyr::select(c(escapement, lnalpha.c, beta, recruitment, variable, year)) %>%
  rbind(., dataset) %>%
  arrange(variable, escapement) -> dataset

ggplot(data=dataset, aes(x=escapement, y=recruitment, group=variable))+
  geom_line(data=subset(dataset,dataset$variable<52),linetype="solid", size=0.5, color="grey80")+
  scale_y_continuous(labels = comma,breaks = seq(0, 1500000, 250000), limits = c(0, 1500000))+
  scale_x_continuous(labels= comma,breaks = seq(0, 450000, 100000), limits = c(0, 450000))+
  ylab("Recruits (R)")+xlab("Spawners (S)")+
  geom_line(data=dataset, aes(x=escapement, y=escapement, group=1),linetype="solid", size=1)+#replacement line
  geom_vline(xintercept = 20000,linetype = "solid", color = "red", size=1)+
  geom_vline(xintercept = 75000,linetype = "solid", color = "red", size=1)+
  geom_line(data=subset(dataset,variable==51),colour = "black", lty=2, size=2)+
  geom_point(data=subset(dataset,variable==52),colour = "black", pch=16, size=1)+
  geom_text(size=3, data=dataset, aes(x=escapement, y=recruitment, group=52, label=year,family="Times", 
                                      hjust = -0.1, vjust= -0.4)) +
  theme_classic(base_size = 14)
ggsave("output/rjags_Full_BaseCase/processed/horsetail.png", dpi = 500, height = 6, width = 8, units = "in")

#Alternative horsetail plot
dataset %>%
  filter (variable %in% c(52)) %>%
  mutate(escapement1 = escapement) %>%
  dplyr::select(-c(lnalpha.c, beta, escapement)) %>%
  mutate(escapement = 'NA',
         Median = 'NA',
         q95 = 'NA',
         q90 ='NA',
         q10 ='NA',
         q5 = 'NA') -> dataset
CI %>%
  mutate(year = 'NA',
         variable = 51,
         recruitment = 'NA',
         escapement1 ='NA')%>%
  rbind(., dataset) %>%
  mutate_if(is.character, as.numeric) -> dataset1

x.fact <- 100/max(dataset1$escapement) 
y.fact <- 100/max(dataset1$recruitment1)
coords <- FFieldPtRep(coords = cbind(dataset1$escapement1 * x.fact, dataset1$recruitment * y.fact), rep.fact = 40)
x.t <- coords$x/x.fact 
y.t <- coords$y/y.fact 
# https://stackoverflow.com/questions/61564313/remove-axis-titles-in-gg-gap-plots
ggplot(data=dataset1, aes(x=escapement, y=Median, group=variable)) + 
  geom_line(size=1, lty=2, group=51) +
  geom_ribbon(aes(ymin = q5, ymax = q95, group=51), alpha=.08) +
  geom_ribbon(aes(ymin = q10, ymax = q90, group=51), alpha=.08) +
  xlab('Spawners (S)') +
  ylab('Recruits (R)') +
  scale_y_continuous(labels = comma, breaks = seq(0, 1300000, 100000), limits = c(0, 1300000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 225000, 25000), limits = c(0, 225000)) +
  geom_line(aes(x=escapement, y=escapement, group=51),linetype="solid", size=1) +
  geom_point(data=dataset1, aes(x=x.t, y=y.t, group=52),pch=16, size=1) +
  geom_point(data=dataset1, aes(x=escapement1, y=recruitment, group=52),pch=16, size=1) +
  geom_vline(xintercept = 20000,linetype = "solid", color = "red", size=1)+
  geom_vline(xintercept = 75000,linetype = "solid", color = "red", size=1)+
  geom_vline(xintercept = 60000,linetype = "dotted", color = "black", size=1)+
  geom_text(size=3, data=dataset1, aes(x=escapement1, y=recruitment, group=52, label=year,family="Times", 
                                       hjust = -0.1, vjust= -0.4)) +
  theme_report(base_size = 14)
#gg <- gg.gap(plot=p,segments=c(800000,1200000),tick_width=(100000), ylim=c(0,1300000))
ggsave("output/rjags_Full_BaseCase/processed/horsetail2.png", dpi = 500, height = 6, width = 8, units = "in")

