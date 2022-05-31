#Honey buzzard ssf in inla
#Elham Nourani, PhD.
#MPI of Animal Behavior. Konstanz, DE. May 23. 2022

library(INLA)
library(tidyverse)
library(ggregplot)
library(magrittr)
library(rNDOW) #for ggTraj

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

#### STEP 1: model building ---------------------------------------------------------------------

#model formula
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


# build the model (I did this using base R, not rstudio to save RAM)

all_data <- readRDS("/home/enourani/Desktop/Hester_HB/inla_input.rds")

mean.beta <- 0
prec.beta <- 1e-4 


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



#extract CPO and Marginal likelihood values
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


#### STEP 2: extract info from model output for plotting ---------------------------------------------------------------------

#original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)

#extract info individual variation plots. done separately for each model
tab_blh <- data.frame(ID = as.factor(M_i$summary.random$id1$ID),
                     mean = M_i$summary.random$id1$mean,
                     IClower = M_i$summary.random$id1[, 4],
                     ICupper = M_i$summary.random$id1[, 6])


tab_wspt <- data.frame(ID = as.factor(M_i$summary.random$id2$ID),
                       mean = M_i$summary.random$id2$mean,
                       IClower = M_i$summary.random$id2[, 4],
                       ICupper = M_i$summary.random$id2[, 6])

#saveRDS(tab_vv, file = "/home/enourani/Desktop/Hester_HB/tab_vv_M_v.rds")
#saveRDS(tab_wspt, file = "/home/enourani/Desktop/Hester_HB/tab_wspt_M_v.rds")


#saveRDS(tab_blh, file = "/home/enourani/Desktop/Hester_HB/tab_blh_M_i.rds")
#saveRDS(tab_wspt, file = "/home/enourani/Desktop/Hester_HB/tab_wspt_M_i.rds")



# extract info for coefficient plots

# posterior means of coefficients
graph <- as.data.frame(summary(M_i)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)


#saveRDS(graph, file = "/home/enourani/Desktop/Hester_HB/graph_M_v.rds")
#saveRDS(graph, file = "/home/enourani/Desktop/Hester_HB/graph_M_i.rds")



#### STEP 3: coefficient plots ---------------------------------------------------------------------

gr_blh <- readRDS("/home/enourani/Desktop/Hester_HB/graph_M_i.rds")
gr_vv <- readRDS("/home/enourani/Desktop/Hester_HB/graph_M_v.rds")


#X11(width = 8.4, height = 3)
png("/home/enourani/Desktop/Hester_HB/coeff_plot.png", units = "in", width = 8.4, height = 3, res = 300)

par(mfrow = c(1,2),
    cex = 0.7,
    oma = c(0,14,0,0),
    mar = c(3, 1, 3, 1),
    bty = "l",
    #font = 3,
    font.axis = 3
)


for(graph in list(gr_blh,gr_vv)){
  

VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

min <- min(graph$Lower,na.rm = T)
max <- max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-0.05,0.14), ylim = c(0.7,7.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph$Estimate, graph$Factor_n, col = "cornflowerblue", pch = 20, cex = 2)
arrows(graph$Lower, graph$Factor_n,
       graph$Upper, graph$Factor_n,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line

#add axes
axis(side= 1, at = c(-0.04, 0, 0.04, 0.08,0.12), labels =  c(-0.04, 0, 0.04, 0.08, 0.12), 
     tick=T ,col = NA, col.ticks = 1, tck =-.015)

if("BLH_z" %in% graph$Factor){
axis(side = 2, at = c(1:7),
     labels = c("Wind support : Uplift : Migration year ", "Uplift : Migration year", "Wind support : Migration year", "Wind support : Uplift","Migration year", "Uplift",
                "Wind support"),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 
  mtext("Model A (boundary layer height as uplift)", side = 3, cex = 0.8, line = 1)
  
} else {
  axis(side = 2, at = c(1:7),
       labels = c("","", "", "", "", "", ""),
       tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
       tck=-.015 , #tick marks smaller than default by this proportion
       las=2) # text perpendicular to axis label 
  
  mtext("Model B (vertical velocity as uplift)", side = 3, cex = 0.8, line = 1)
}

}

dev.off()



#### STEP 4: individual variation plots ---------------------------------------------------------------------

#assign migratory year to each individual ID
load("/home/enourani/Desktop/Hester_HB/HB_annotated_data.RData") #annotated_data

age_id <- annotated_data %>% 
  group_by(id) %>% 
  summarize(age = as.character(head(category,1))) %>% #extract adult vs first year data
  arrange(age) %>% 
  as.data.frame()

#re-order names based on age
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

#plot
#X11(width = 9, height = 7)

png("/home/enourani/Desktop/Hester_HB/ind_var.png", units = "in", width = 9, height = 7, res = 300)

par(mfrow = c(1,2), 
    bty="n",
    cex = 0.7,
    oma = c(0,3.5,0,3),
    mar = c(3, 2, 3,1),
    font.axis = 3
)

for(var in c("i","v")){ #all files with i are related to the model with blh. all with v are related to vertical velocity
  
  variables <- str_sub(list.files("/home/enourani/Desktop/Hester_HB/", pattern = paste0("^tab.*", var), full.names = T), 39, -9) #extract variable names for this model
  
  files <- list.files("/home/enourani/Desktop/Hester_HB/", pattern = paste0("^tab.*", var), full.names = T) %>% #list files for this model
    lapply(readRDS)
  
  names(files) <- variables
  
  plot(0, bty = "l", labels = FALSE, tck = 0, xlim = c(-0.23,0.17), ylim = c(0,29.5), xlab = "", ylab = "")
  #add vertical line for zero
  abline(v = 0, col = "grey30",lty = 2)
  
  
  points(files[[1]]$mean, as.numeric(files[[1]]$ID) - 0.2, col = "darkgoldenrod2", pch = 19, cex = 1.3)
  arrows(files[[1]]$IClower, as.numeric(files[[1]]$ID) - 0.2,
         files[[1]]$ICupper, as.numeric(files[[1]]$ID) - 0.2,
         col = "darkgoldenrod2", code = 3, length = 0.03, angle = 90, lwd = 1) #angle of 90 to make the arrow head as straight as a line
  
  points(files[[2]]$mean, as.numeric(files[[2]]$ID) , col = "cornflowerblue", pch = 19, cex = 1.3)
  arrows(files[[2]]$IClower, as.numeric(files[[2]]$ID) ,
         files[[2]]$ICupper, as.numeric(files[[2]]$ID) ,
         col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 1) 
  
  axis(side= 1, at = c(-0.15,0,0.15), labels = c(-0.15,0,0.15), 
       tick=T ,col = NA, col.ticks = 1, tck=-.015)
  
  #add title
  if( var == "i"){
    axis(side= 2, at= c(1:28), 
         labels = files[[2]]$ID, 
         tick = T ,col = NA, col.ticks = 1, 
         tck = -.015 , 
         las = 2) 
    
    mtext("Model A (boundary layer height as uplift)", side = 3, cex = 0.8, line = 0.5)
  
    } else {
    mtext("Model B (vertical velocity as uplift)", side = 3, cex = 0.75, line = 0.5)
      axis(side= 2, at= c(1:28), 
           labels =  NA, 
           tick = T ,col = NA, col.ticks = 1, 
           tck = -.015 , 
           las = 2) 
      
      #add legend
      legend("topright", inset=c(-0.01,0), legend = c("Wind support", "Uplift", "Tagged as adult"), 
             col = c("cornflowerblue","darkgoldenrod1", "black"), 
             pch = 19, bg="white",bty="n", cex = 0.9)
      }
  
  points(x =rep(-.24,5), y = which(files[[2]]$ID %in% age_id[age_id$age == "Adult","id"]), pch = 20,cex = 1)
  
}

dev.off()

#### STEP 5: make maps for Mohammed ---------------------------------------------------------------------

m <- all_data %>% 
  filter(id == "Mohammed" & case_ == 1) %>% 
  arrange(timestamp)

ggTraj(m, xyi = c("lon2", "lat2"))

ggplot(m, aes(lon2,lat2)) +
  geom_path() +
  geom_tile(map.raster) +
  facet_wrap(.~year)


#https://r-spatial.org/r/2018/10/25/ggplot2-sf-3.html
library(rnaturalearth)
library(rnaturalearthdata)

base <- world <- ne_coastline(scale = 'medium', returnclass = 'sf')



ggplot(data = base) +
  geom_sf(col = "gray", fill = "gray") +
  coord_sf(xlim = c(-10, 38), ylim = c(0, 64), expand = FALSE) +
  geom_path(data = m, aes(x = lon2, y = lat2, col = tail), size = 2, lineend = "round") +
  #scale_colour_gradientn(colours = wes_palette("Zissou1", 3, type = "continuous"), limits = c(-17,22), 
  #                     na.value = "white", name = "Wind support\n (m/s)")+
  
  scale_colour_viridis(option = "mako", na.value = "white", name = "Wind support\n (m/s)") +
  theme_bw() +
  facet_wrap(.~year, nrow = 1)

ggplot(data = base) +
  geom_sf(col = "gray", fill = "gray") +
  coord_sf(xlim = c(-10, 38), ylim = c(0, 64), expand = FALSE) +
  geom_path(data = m, aes(x = lon2, y = lat2, col = BLH), size = 2, lineend = "round") +
  scale_colour_viridis(option = "mako", na.value = "white",  name = "Boundary layer\n height (m)")+
  theme_bw() +
  facet_wrap(.~year, nrow = 1)

ggplot(data = base) +
  geom_sf(col = "gray", fill = "gray") +
  coord_sf(xlim = c(-10, 38), ylim = c(0, 64), expand = FALSE) +
  geom_path(data = m, aes(x = lon2, y = lat2, col = vertical_pressure), size = 2.5, lineend = "round") +
  #scale_colour_gradientn(colours = wes_palette("Darjeeling1", 150, type = "continuous"), limits = c(-1,1.06), 
  #                       na.value = "white", name = "Vertical velocity\n (Pa/s)")+
  scale_colour_viridis(option = "mako", na.value = "white", name = "Vertical velocity\n (Pa/s)") +
  theme_bw() +
  facet_wrap(.~year, nrow = 1)



