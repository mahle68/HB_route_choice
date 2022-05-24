#honey buzzard ssf in inla

library(INLA)
library(tidyverse)
library(ggregplot)
#library(survival)
#library(jtools) #to plot coeffs of clogit
library(magrittr)

#open data
load("/home/enourani/Desktop/Hester_HB/HB_annotated_data.RData")

#add variables for ID and year
all_data <- annotated_data %>% 
  mutate(migration = as.numeric(ifelse(migration == "Adult", 5, migration)), #create a numeric variable by giving adults age 5
         id1 = factor(id),
         id2 = factor(id),
         yr1 = factor(year),
         yr2 = factor(year))

all_data <- all_data %>% 
  drop_na(tail) %>%
  mutate_at(c("tail", "BLH", "migration" , "vertical_pressure"),
            list(z = ~(scale(.)))) %>% 
  rename(stratum = step_id_)

#
### use clogit to see
library(survival)
f <- case_ ~  + 
  BLH_z * category +
  tail_z * category + 
  vertical_pressure_z * category +
  strata(stratum)

f2 <- case_ ~  + 
  BLH_z * tail_z + 
  strata(stratum)

f2 <- case_ ~  + 
  vertical_pressure_z * tail_z + 
  strata(stratum)

ssf <- clogit(f, data = all_data)

#separate models with no interaction terms
ssf_all <- clogit(f2, data = all_data)

d <- all_data[all_data$category == "Adult",]
ssf_a <-  clogit(f2, data = d)

int <- all_data[all_data$category == "Intermediate",]
ssf_i <-  clogit(f2, data = int)

j <- all_data[all_data$category == "First year",]
ssf_j <-  clogit(f2, data = j)


plot_summs(list(ssf_a, ssf_i, ssf_j))

#include category as a continous variable!!

####
#model formula
formula <- case_ ~ -1 + tail_z * BLH_z * scale(as.numeric(category)) +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+

#use blh
formula_n <- case_ ~ -1 + tail_z * BLH_z * migration_z + #instead of three levels, migration has five
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+

#use vertical velocity
formula_n_vv <- case_ ~ -1 + tail_z * vertical_pressure_z * migration_z  +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+

formula_i <- case_ ~ -1 + tail_z * BLH_z * migration_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(id1, BLH_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(id2, tail_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) #+
#  f(yr1, BLH_z, model = "iid",
#    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
#  f(yr2, tail_z,  model = "iid",
#    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

formula_v <- case_ ~ -1 + tail_z * vertical_pressure_z * migration_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(id1, vertical_pressure_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(id2, tail_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) #+
#  f(yr1, BLH_z, model = "iid",
#    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
#  f(yr2, tail_z,  model = "iid",
#    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))



n <- 500

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(case_ = NA,
         category = sample(unique(all_data$category), n, replace = T),
          migration_z = sample(seq(min(all_data$all_data), max(all_data$all_data), length.out = 10), n, replace = T),
         BLH_z = sample(seq(min(all_data$BLH_z), max(all_data$BLH_z), length.out = 10), n, replace = T), #regular intervals, so we can make a raster later on
         tail_z = sample(seq(min(all_data$tail_z), max(all_data$tail_z), length.out = 10), n, replace = T)) %>% 
  full_join(all_data)

all_data <- readRDS("/home/enourani/Desktop/Hester_HB/inla_input.rds")

mean.beta <- 0
prec.beta <- 1e-4 

#Model for predictions
(b <- Sys.time())
M_pred_c <- inla(formula, family = "Poisson", 
                control.fixed = list(
                  mean = mean.beta,
                  prec = list(default = prec.beta)),
                data = new_data, 
                num.threads = 10,
                control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
                control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b #500 missing values, no random effects: 17 mins

Efxplot(M_pred) +
  xlim(-.2,0.1)


#model without missing values
(b <- Sys.time())
M_v <- inla(formula_v, family = "Poisson", 
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 data = all_data, 
                 num.threads = 10,
                 control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
                 control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b #with ind id as random effect: 5 min in base r; with vertical velocity and ind ID, 8 min

Efxplot(list(M_n,M_n_3)) 

Efxplot(list(M_n,M_n_vv)) 

#i expect blh to have performed better. test this:
M_n$mlik[1] < M_n_vv$mlik[1] 

mean(M_n$cpo$cpo[1]) < mean(M_n_vv$cpo$cpo[1]) 

saveRDS(M_n, file = "/home/enourani/Desktop/Hester_HB/M_n_5yrs_norndm.rds") #this is better than the 3 year version... has lower cpo and mlik too
saveRDS(M_n_vv, file = "/home/enourani/Desktop/Hester_HB/M_n_5yrs_norndm_vv.rds")


mean(M_v$cpo$cpo) #model with vertical v and random effect for individual id
# 0.9799288

M_v$mlik
# integration -38533.82
#Gaussian -38534.20

mean(M_i$cpo$cpo) #model with blh and random effect for individual id
#0.9799288

M_i$mlik
# integration : -38518.27 
#Gaussian : -38517.73 (this one is reported)


######### ind variation plots ####################################################

#original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)

#extract info for making the plots later on 

tab_blh <- data.frame(ID = as.factor(M_i$summary.random$id1$ID),
                     mean = M_i$summary.random$id1$mean,
                     IClower = M_i$summary.random$id1[, 4],
                     ICupper = M_i$summary.random$id1[, 6])


tab_wspt <- data.frame(ID = as.factor(M_i$summary.random$id2$ID),
                       mean = M_i$summary.random$id2$mean,
                       IClower = M_i$summary.random$id2[, 4],
                       ICupper = M_i$summary.random$id2[, 6])

saveRDS(tab_vv, file = "/home/enourani/Desktop/Hester_HB/tab_vv_M_v.rds")
saveRDS(tab_wspt, file = "/home/enourani/Desktop/Hester_HB/tab_wspt_M_v.rds")


saveRDS(tab_blh, file = "/home/enourani/Desktop/Hester_HB/tab_blh_M_i.rds")
saveRDS(tab_wspt, file = "/home/enourani/Desktop/Hester_HB/tab_wspt_M_i.rds")


#plot
X11(width = 9, height = 9)

par(mfrow = c(1,2), bty="n",
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 2, 0.5, 1)
)


plot(0, bty = "l", labels = FALSE, tck = 0, xlim = c(-0.2,0.15), ylim = c(0.5,28.5), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_vv$mean, as.numeric(tab_vv$ID) - 0.2, col = "darkgoldenrod2", pch = 19, cex = 1.3)
arrows(tab_vv$IClower, as.numeric(tab_vv$ID) - 0.2,
       tab_vv$ICupper, as.numeric(tab_vv$ID) - 0.2,
       col = "darkgoldenrod2", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line

points(tab_wspt$mean, as.numeric(tab_wspt$ID) , col = "cornflowerblue", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) ,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) ,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 2) 

axis(side= 1, at = c(-2,-1,0,1,2), labels = c(-2,-1,0,1,2), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4), 
     labels =  tab_dt$ID, 
     tick = T ,col = NA, col.ticks = 1, 
     tck = -.015 , 
     las = 2) 

#add legend
legend(x = 1.25, y = 4.5, legend = c( "Wind support var", "Wind support", expression(italic(paste(Delta,"T")))), 
       col = c("pink1", "cornflowerblue","darkgoldenrod1"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.9)


### extract info for coeff plots

# posterior means of coefficients
graph <- as.data.frame(summary(M_i)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)


#saveRDS(graph, file = "/home/enourani/Desktop/Hester_HB/graph_M_v.rds")
saveRDS(graph, file = "/home/enourani/Desktop/Hester_HB/graph_M_i.rds")



##########include movement kernel
f_mv <- case_ ~ -1 + log(sl_) + cos(ta_) + tail_z * BLH_z * migration_z + #instead of three levels, migration has five
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+

f_mv_i <- case_ ~ -1 + log(sl_) * cos(ta_) * migration_z + #instead of three levels, migration has five
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+

f_mv_ii <- case_ ~ -1 + scale(log(sl_)) * tail_z * BLH_z * migration_z + #instead of three levels, migration has five
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+


(b <- Sys.time())
M_mv_ii <- inla(f_mv_ii, family = "Poisson", 
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = all_data, 
               num.threads = 10,
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b #no random effect: 57 sec; random effect of individual id


Efxplot(list(M_mv_i, M_mv_ii))


M_mv_i$mlik[1] < M_mv_ii$mlik[1] 
