###Age at maturity for Banded tulips
#
#ELEFAN methods fine tuning
#
#
#Use Alt+O to collapse all sections, Alt+Shift+O to expand all sections
#
#
.rs.restartR() #Restarts session (good if rerunning after working with other files)
graphics.off()  # turns off any plots from previous work session
rm(list=ls(all=TRUE)) # clears out environment 
#
#Load require packages (install as necessary)
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
if (!require("remotes")) install.packages("remotes")
remotes::install_github("rschwamborn/fishboot")
devtools::install_github("stefanedwards/lemon")
devtools::install_github("briandconnelly/growthcurve", build_vignettes = TRUE)
pacman::p_load(plyr, tidyverse, rstatix, #Df manipulation, basic summary
               zoo, FSA, nlstools, #length bins, growth tags
               ggmap, ggsn, ggpubr, ggfortify, scales, lemon, #Mapping and figures
               car, TropFishR, lubridate, MASS, fishboot, #ELEFAN
               flextable, janitor, officer)
#
#
#
#
####Load files####
#
#########Mark Recapture
Station_MR <- read.csv("../CSV/TB/TB_Pop_Target_final.csv", na.string = c("Z", "", "NA", " "))
glimpse(Station_MR)
#Limit to data required and reorder columns
Station_MR <- Station_MR %>% 
  dplyr::select(Date, Month, Year, SurveyMonth, Site, Station, Species, ID, SL, SW, Recaptured, Helper) %>%
  filter(Species == "BT") %>% drop_na(Station)
head(Station_MR)
#
#
Extra_MR <- read.csv("../CSV/TB/TB_Pop_Extra_final.csv", na.string = c("Z", "", "NA", " "))
glimpse(Extra_MR)
Extra_MR <- Extra_MR %>% dplyr::select(Date, Month, Year, Site:SW, Recaptured) %>%
  filter(Species == "BT") %>% drop_na(Station)
head(Extra_MR)
#
#
#########ELEFAN
Elefant <- read.csv("../CSV/TB/All_snails.csv", stringsAsFactors = FALSE, na.strings = c("NA", " ", "", "Z")) %>%
  filter(Species == "BT") %>% drop_na(Station, SL) %>% dplyr::select(Month, Year, SL) #limit to data needed
#
head(Elefant)
#
Elefant <- Elefant %>% mutate(Date = make_date(year = Year, month = Month),
                              SL_norm = scale(SL)) %>%
  subset(Date < "2021-06-01") #limit to similar timeframe as recapture data
#
#
#
####Length frequency - set up####
#
summary(Elefant)
(Ele_summ <- Elefant %>% get_summary_stats(SL, type = "full") %>% dplyr::select(n, min, max, mean, sd) %>% mutate(Type = "ELEFAN"))
#
Elefant %>% 
  ggplot(aes(SL))+
  geom_histogram(binwidth = 6, boundary = 0, closed = "left")+
  lemon::facet_rep_grid(vars(Month))+
  theme(panel.spacing = unit(0, "lines"), 
                          axis.text.y = element_text(size = 6, margin = margin(r = 5)), 
                          axis.text.x = element_text(size = 8, margin = margin(t = 5)))+
  #scale_y_continuous("Count", expand = c(0,0), limits = c(0, 48), breaks = seq(0, 48, 16))+
  scale_x_continuous("Shell length(mm)", expand = c(0,0), limits = c(0, 100), breaks = seq(0, 100, 20))
#
#Manually save:(file = "Output/Trials/LFQ_bin6.tiff", width = 1200) 
#
###Try: bin 2, MA = 9, bin 4, MA = 5; bin = 6, MA = 3 
#
#
####ELEFAN run - b2, MA9####
set_MA <- c(9)
###Create lfq structure
lfq2_9 <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3_9 <- lfqModify(lfq2_9, bin_size = 2) # modify bin size
lfq4_9 <- lfqRestructure(lfq3_9, MA = set_MA, addl.sqrt = FALSE) # restructuring process
#
#Plot bins 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_9, Fname = "catch", date.axis = "modern")
plot(lfq4_9, Fname = "rcounts", date.axis = "modern")
#Manually save:(file = "Output/Trials/Elefan_trials_bin2_MA9_LFrestructure.tiff", width = 1200) 
par(opar)
par(mfrow = c(1,1))
#
#
####Preliminary estimates of Linf
#Powell-Wetherall plot - 1:28
res_PW_9 <- powell_wetherall(param = lfq4_9,
                           catch_columns = 1:ncol(lfq2_9$catch))
#
paste("Linf =",round(res_PW_9$Linf_est), "+/-", round(res_PW_9$se_Linf)) # show results
Linf_pw_9 <- res_PW_9$Linf_est
#
#Max length from raw data
(Linf_ml <- max(Elefant$SL, na.rm = T))
# Froese and Binohlan formula
(Linf_pb <- exp(0.44 + 0.984*log(Linf_ml)))

#Linf range 1
Linf_1_9 <- c(Linf_pw_9, Linf_pb)
Linf_1_mid_9 <- Linf_1_9[1] + ((Linf_1_9[2] - Linf_1_9[1])/2)

#Linf range 2
Linf_2 <- c(Linf_ml, Linf_pb)
Linf_2_mid <- Linf_2[1] + ((Linf_2[2] - Linf_2[1])/2)

# tanchor may be from Mar-Jul
t_an <- c(3/12, 4/12, 5/12, 6/12, 7/12)
t_an_d <- c(0.162, 0.247, 0.329, 0.411, 0.496)
#
#
#### ELEFAN analysis with simulated annealing
# run ELEFAN_SA with range 1, NOT seasonalized
set.seed(1)
res_SA_1_9 <- ELEFAN_SA(lfq4_9, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                      MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                      init_par = list(Linf = Linf_1_mid_9, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                      low_par = list(Linf = Linf_1_9[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_1_9[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2_9 <- ELEFAN_SA(lfq4_9, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                      MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                      init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                      low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S_9 <- ELEFAN_SA(lfq4_9, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_1_mid_9, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_1_9[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_9[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S_9 <- ELEFAN_SA(lfq4_9, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1_9 <- ELEFAN_GA(lfq4_9, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                      addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                      low_par = list(Linf = Linf_1_9[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_1_9[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2_9 <- ELEFAN_GA(lfq4_9, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                      addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                      low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S_9 <- ELEFAN_GA(lfq4_9, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_1_9[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_9[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S_9 <- ELEFAN_GA(lfq4_9, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### plot LFQ and growth curves
opar <- par(mfrow = c(4,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_9, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq4_9, par = res_SA_1_9$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq4_9, par = res_GA_1_9$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_9, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_9, par = res_SA_2_9$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_9, par = res_GA_2_9$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_9, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_9, par = res_SA_1_S_9$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_9, par = res_GA_1_S_9$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_9, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_9, par = res_SA_2_S_9$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_9, par = res_GA_2_S_9$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
#Save file manually: tiff("Trials/Elefan_b2M9_model fits.tiff", width = 1200)
par(opar)
par(mfrow = c(1,1))
#
#Combine results for comparison
res_all_9 <- as.data.frame(rbind(c(paste("SA 1"), paste(round(Linf_1_9[1], 3), round(Linf_1_9[2],3)), paste("NS"), 
                                 paste("SA"), round(res_SA_1_9$par$Linf,3), round(res_SA_1_9$par$t_anchor,3), 
                                 round(res_SA_1_9$par$K,3), round(res_SA_1_9$Rn_max,3)),
                               c(paste("SA 1S"), paste(round(Linf_1_9[1], 3), round(Linf_1_9[2], 3)), paste("S"), 
                                 paste("SA"), round(res_SA_1_S_9$par$Linf,3), round(res_SA_1_S_9$par$t_anchor,3), 
                                 round(res_SA_1_S_9$par$K,3), round(res_SA_1_S_9$Rn_max,3)),
                               c(paste("SA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                 paste("SA"), round(res_SA_2_9$par$Linf,3), round(res_SA_2_9$par$t_anchor,3), 
                                 round(res_SA_2_9$par$K,3), round(res_SA_2_9$Rn_max,3)),
                               c(paste("SA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                 paste("SA"), round(res_SA_2_S_9$par$Linf,3), round(res_SA_2_S_9$par$t_anchor,3), 
                                 round(res_SA_2_S_9$par$K,3), round(res_SA_2_S_9$Rn_max,3)),
                               c(paste("GA 1"), paste(round(Linf_1_9[1], 3), round(Linf_1_9[2], 3)), paste("NS"), 
                                 paste("GA"), round(res_GA_1_9$par$Linf,3), round(res_GA_1_9$par$t_anchor,3), 
                                 round(res_GA_1_9$par$K,3), round(res_GA_1_9$Rn_max,3)),
                               c(paste("GA 1S"), paste(round(Linf_1_9[1], 3), round(Linf_1_9[2], 3)), paste("S"), 
                                 paste("GA"), round(res_GA_1_S_9$par$Linf,3), round(res_GA_1_S_9$par$t_ancho,3), 
                                 round(res_GA_1_S_9$par$K,3), round(res_GA_1_S_9$Rn_max,3)),
                               c(paste("GA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                 paste("GA"), round(res_GA_2_9$par$Linf,3), round(res_GA_2_9$par$t_anchor,3), 
                                 round(res_GA_2_9$par$K,3), round(res_GA_2_9$Rn_max,3)),
                               c(paste("GA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                 paste("GA"), round(res_GA_2_S_9$par$Linf,3), round(res_GA_2_S_9$par$t_anchor,3), 
                                 round(res_GA_2_S_9$par$K,3), round(res_GA_2_S_9$Rn_max,3))))
names(res_all_9) <- c("Model", "Linf_range", "Seasonalized", "SA_GA", "Linf", "t_anchor", "K", "Rn_max")
print(res_all_9) #GA1S
#
plot(lfq4_9, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_9, par = res_GA_1_S_9$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
#Save manually: (path = "Output/Trials", filename = "Elefan_b2m9_GA1S_fit.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK_9 <- vector("list", length(lfq4_9$dates))
for(i in 1:length(lfq4_9$dates)){
  loop_data <- list(dates = lfq4_9$dates[-i],
                    midLengths = lfq4_9$midLengths,
                    catch = lfq4_9$catch[,-i])
  tmp <- ELEFAN_GA(loop_data, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                   addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                   low_par = list(Linf = Linf_1_9[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_1_9[2], K = 1, t_anchor = 1, C = 1, ts = 1))
  
  
  JK_9[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres_9 <- do.call(cbind, JK_9)
# mean
JKmeans_9 <- apply(as.matrix(JKres_9), MARGIN = 1, FUN = mean)
#confidence intervals
JKconf_9 <- apply(as.matrix(JKres_9), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JK_table_9 <- t(rbind(JKmeans_9, JKconf_9))
colnames(JK_table_9) <- c("means", "lower","upper")
# show results
JK_table_9
#
#### estimation of Z and M
Ms_9 <- M_empirical(Linf = res_GA_1_S_9$par$Linf, K_l = res_GA_1_S_9$par$K, method = "Then_growth")
lfq4_9$M <- as.numeric(Ms_9)
paste("M =", as.numeric(Ms_9))
#> [1] "M = 0.422"
res_GA_1_S_9$agemax
#
#Figure of model output
print(res_GA_1_S_9$par) #print to copy numbers into figure
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=63.10317, K=0.2874028, C=0.3966056, ts = 0.4193841, t0 = 0.3277789)),
           "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 25))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_GA_1_S_9$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+ theme_classic()
#
ggsave(path = "Output/Trials/", filename = "b2m9_GA1S_VB_63.1_0.287_0.470.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4_9, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_9, par = res_GA_1_S_9$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
L_9 <- seq(0,100,0.1)
t_9 <- VBGF(L = seq(0,100,0.1), list(Linf=63.10317, K=0.2874028, C=0.3966056, ts = 0.4193841, t0 = 0.3277789)) 
plot(t_9, L_9, t="l", xlim = c(0, (res_GA_1_S_9$agemax +  5)), ylim = c(0, 100))
abline(v = res_GA_1_S_9$agemax , col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Trials", filename = "b2m9_GA1S_model_Age_11.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
#
#

####ELEFAN run - b4, MA5####
set_MA <- c(5)
###Create lfq structure
lfq2_5 <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3_5 <- lfqModify(lfq2_5, bin_size = 4) # modify bin size
lfq4_5 <- lfqRestructure(lfq3_5, MA = set_MA, addl.sqrt = FALSE) # restructuring process
#
#Plot bins 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_5, Fname = "catch", date.axis = "modern")
plot(lfq4_5, Fname = "rcounts", date.axis = "modern")
#Manually save:(file = "Output/Trials/Elefan_trials_bin4_MA5_LFrestructure.tiff", width = 1200) 
par(opar)
par(mfrow = c(1,1))
#
#
####Preliminary estimates of Linf
#Powell-Wetherall plot - 1:15
res_PW_5 <- powell_wetherall(param = lfq4_5,
                             catch_columns = 1:ncol(lfq2_5$catch))
#
paste("Linf =",round(res_PW_5$Linf_est), "+/-", round(res_PW_5$se_Linf)) # show results
Linf_pw_5 <- res_PW_5$Linf_est
#
#Max length from raw data
(Linf_ml <- max(Elefant$SL, na.rm = T))
# Froese and Binohlan formula
(Linf_pb <- exp(0.44 + 0.984*log(Linf_ml)))

#Linf range 1
Linf_1_5 <- c(Linf_pw_5, Linf_pb)
Linf_1_mid_5 <- Linf_1_5[1] + ((Linf_1_5[2] - Linf_1_5[1])/2)

#Linf range 2
Linf_2 <- c(Linf_ml, Linf_pb)
Linf_2_mid <- Linf_2[1] + ((Linf_2[2] - Linf_2[1])/2)

# tanchor may be from Mar-Jul
t_an <- c(3/12, 4/12, 5/12, 6/12, 7/12)
t_an_d <- c(0.162, 0.247, 0.329, 0.411, 0.496)
#
#
#### ELEFAN analysis with simulated annealing
# run ELEFAN_SA with range 1, NOT seasonalized
set.seed(1)
res_SA_1_5 <- ELEFAN_SA(lfq4_5, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_1_mid_5, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_1_5[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_5[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2_5 <- ELEFAN_SA(lfq4_5, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S_5 <- ELEFAN_SA(lfq4_5, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = Linf_1_mid_5, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                          low_par = list(Linf = Linf_1_5[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1_5[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S_5 <- ELEFAN_SA(lfq4_5, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                          low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1_5 <- ELEFAN_GA(lfq4_5, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_1_5[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_5[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2_5 <- ELEFAN_GA(lfq4_5, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S_5 <- ELEFAN_GA(lfq4_5, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                          addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                          low_par = list(Linf = Linf_1_5[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1_5[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S_5 <- ELEFAN_GA(lfq4_5, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                          addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                          low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### plot LFQ and growth curves
opar <- par(mfrow = c(4,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_5, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq4_5, par = res_SA_1_5$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq4_5, par = res_GA_1_5$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_5, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_5, par = res_SA_2_5$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_5, par = res_GA_2_5$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_5, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_5, par = res_SA_1_S_5$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_5, par = res_GA_1_S_5$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_5, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_5, par = res_SA_2_S_5$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_5, par = res_GA_2_S_5$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
#Save file manually: tiff("Trials/Elefan_b4M5_model fits.tiff", width = 1200)
par(opar)
par(mfrow = c(1,1))
#
#Combine results for comparison
res_all_5 <- as.data.frame(rbind(c(paste("SA 1"), paste(round(Linf_1_5[1], 3), round(Linf_1_5[2],3)), paste("NS"), 
                                   paste("SA"), round(res_SA_1_5$par$Linf,3), round(res_SA_1_5$par$t_anchor,3), 
                                   round(res_SA_1_5$par$K,3), round(res_SA_1_5$Rn_max,3)),
                                 c(paste("SA 1S"), paste(round(Linf_1_5[1], 3), round(Linf_1_5[2], 3)), paste("S"), 
                                   paste("SA"), round(res_SA_1_S_5$par$Linf,3), round(res_SA_1_S_5$par$t_anchor,3), 
                                   round(res_SA_1_S_5$par$K,3), round(res_SA_1_S_5$Rn_max,3)),
                                 c(paste("SA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                   paste("SA"), round(res_SA_2_5$par$Linf,3), round(res_SA_2_5$par$t_anchor,3), 
                                   round(res_SA_2_5$par$K,3), round(res_SA_2_5$Rn_max,3)),
                                 c(paste("SA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                   paste("SA"), round(res_SA_2_S_5$par$Linf,3), round(res_SA_2_S_5$par$t_anchor,3), 
                                   round(res_SA_2_S_5$par$K,3), round(res_SA_2_S_5$Rn_max,3)),
                                 c(paste("GA 1"), paste(round(Linf_1_5[1], 3), round(Linf_1_5[2], 3)), paste("NS"), 
                                   paste("GA"), round(res_GA_1_5$par$Linf,3), round(res_GA_1_5$par$t_anchor,3), 
                                   round(res_GA_1_5$par$K,3), round(res_GA_1_5$Rn_max,3)),
                                 c(paste("GA 1S"), paste(round(Linf_1_5[1], 3), round(Linf_1_5[2], 3)), paste("S"), 
                                   paste("GA"), round(res_GA_1_S_5$par$Linf,3), round(res_GA_1_S_5$par$t_ancho,3), 
                                   round(res_GA_1_S_5$par$K,3), round(res_GA_1_S_5$Rn_max,3)),
                                 c(paste("GA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                   paste("GA"), round(res_GA_2_5$par$Linf,3), round(res_GA_2_5$par$t_anchor,3), 
                                   round(res_GA_2_5$par$K,3), round(res_GA_2_5$Rn_max,3)),
                                 c(paste("GA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                   paste("GA"), round(res_GA_2_S_5$par$Linf,3), round(res_GA_2_S_5$par$t_anchor,3), 
                                   round(res_GA_2_S_5$par$K,3), round(res_GA_2_S_5$Rn_max,3))))
names(res_all_5) <- c("Model", "Linf_range", "Seasonalized", "SA_GA", "Linf", "t_anchor", "K", "Rn_max")
print(res_all_5) #GA1S
#
plot(lfq4_5, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_5, par = res_GA_2_S_5$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
#Save manually: (path = "Output/Trials", filename = "Elefan_b4m5_GA2S_fit.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK_5 <- vector("list", length(lfq4_5$dates))
for(i in 1:length(lfq4_5$dates)){
  loop_data <- list(dates = lfq4_5$dates[-i],
                    midLengths = lfq4_5$midLengths,
                    catch = lfq4_5$catch[,-i])
  tmp <- ELEFAN_GA(loop_data, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                   addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                   low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
  
  
  JK_5[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres_5 <- do.call(cbind, JK_5)
# mean
JKmeans_5 <- apply(as.matrix(JKres_5), MARGIN = 1, FUN = mean)
#confidence intervals
JKconf_5 <- apply(as.matrix(JKres_5), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JK_table_5 <- t(rbind(JKmeans_5, JKconf_5))
colnames(JK_table_5) <- c("means", "lower","upper")
# show results
JK_table_5
#
#### estimation of Z and M
Ms_5 <- M_empirical(Linf = res_GA_2_S_5$par$Linf, K_l = res_GA_2_S_5$par$K, method = "Then_growth")
lfq4_5$M <- as.numeric(Ms_5)
paste("M =", as.numeric(Ms_5))
#> [1] "M = 0.422"
res_GA_2_S_5$agemax
#
#Figure of model output
print(res_GA_2_S_5$par) #print to copy numbers into figure
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=100.7117, K=0.2089022, C=0.7227814, ts = 0.6174155, t0 = 0.7415839)),
           "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 25))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_GA_2_S_5$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+ theme_classic()
#
ggsave(path = "Output/Trials/", filename = "b4m5_GA2S_VB_100.7_0.209_0.741.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4_5, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_5, par = res_GA_2_S_5$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
L_5 <- seq(0,100,0.1)
t_5 <- VBGF(L = seq(0,100,0.1), list(Linf=100.7117, K=0.2089022, C=0.7227814, ts = 0.6174155, t0 = 0.7415839)) 
plot(t_5, L_5, t="l", xlim = c(0, (res_GA_2_S_5$agemax +  5)), ylim = c(0, 100))
abline(v = res_GA_2_S_5$agemax , col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Trials", filename = "b4m5_GA2S_model_Age_15.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
#
#
####ELEFAN run - b6, MA3####
set_MA <- c(3)
###Create lfq structure
lfq2_3 <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3_3 <- lfqModify(lfq2_3, bin_size = 6) # modify bin size
lfq4_3 <- lfqRestructure(lfq3_3, MA = set_MA, addl.sqrt = FALSE) # restructuring process
#
#Plot bins 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_3, Fname = "catch", date.axis = "modern")
plot(lfq4_3, Fname = "rcounts", date.axis = "modern")
#Manually save:(file = "Output/Trials/Elefan_trials_bin6_MA3_LFrestructure.tiff", width = 1200) 
par(opar)
par(mfrow = c(1,1))
#
#
####Preliminary estimates of Linf
#Powell-Wetherall plot - 1:12
res_PW_3 <- powell_wetherall(param = lfq4_3,
                             catch_columns = 1:ncol(lfq2_3$catch))
#
paste("Linf =",round(res_PW_3$Linf_est), "+/-", round(res_PW_3$se_Linf)) # show results
Linf_pw_3 <- res_PW_3$Linf_est
#
#Max length from raw data
(Linf_ml <- max(Elefant$SL, na.rm = T))
# Froese and Binohlan formula
(Linf_pb <- exp(0.44 + 0.984*log(Linf_ml)))

#Linf range 1
Linf_1_3 <- c(Linf_pw_3, Linf_pb)
Linf_1_mid_3 <- Linf_1_3[1] + ((Linf_1_3[2] - Linf_1_3[1])/2)

#Linf range 2
Linf_2 <- c(Linf_ml, Linf_pb)
Linf_2_mid <- Linf_2[1] + ((Linf_2[2] - Linf_2[1])/2)

# tanchor may be from Mar-Jul
t_an <- c(3/12, 4/12, 5/12, 6/12, 7/12)
t_an_d <- c(0.162, 0.247, 0.329, 0.411, 0.496)
#
#
#### ELEFAN analysis with simulated annealing
# run ELEFAN_SA with range 1, NOT seasonalized
set.seed(1)
res_SA_1_3 <- ELEFAN_SA(lfq4_3, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_1_mid_3, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_1_3[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_3[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2_3 <- ELEFAN_SA(lfq4_3, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S_3 <- ELEFAN_SA(lfq4_3, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = Linf_1_mid_3, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                          low_par = list(Linf = Linf_1_3[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1_3[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S_3 <- ELEFAN_SA(lfq4_3, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                          low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1_3 <- ELEFAN_GA(lfq4_3, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_1_3[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_3[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2_3 <- ELEFAN_GA(lfq4_3, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S_3 <- ELEFAN_GA(lfq4_3, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                          addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                          low_par = list(Linf = Linf_1_3[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1_3[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S_3 <- ELEFAN_GA(lfq4_3, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                          addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                          low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### plot LFQ and growth curves
opar <- par(mfrow = c(4,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_3, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq4_3, par = res_SA_1_3$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq4_3, par = res_GA_1_3$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_3, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3, par = res_SA_2_3$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_3, par = res_GA_2_3$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_3, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3, par = res_SA_1_S_3$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_3, par = res_GA_1_S_3$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_3, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3, par = res_SA_2_S_3$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_3, par = res_GA_2_S_3$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
#Save file manually: tiff("Trials/Elefan_b6M3_model fits.tiff", width = 1200)
par(opar)
par(mfrow = c(1,1))
#
#Combine results for comparison
res_all_3 <- as.data.frame(rbind(c(paste("SA 1"), paste(round(Linf_1_3[1], 3), round(Linf_1_3[2],3)), paste("NS"), 
                                   paste("SA"), round(res_SA_1_3$par$Linf,3), round(res_SA_1_3$par$t_anchor,3), 
                                   round(res_SA_1_3$par$K,3), round(res_SA_1_3$Rn_max,3)),
                                 c(paste("SA 1S"), paste(round(Linf_1_3[1], 3), round(Linf_1_3[2], 3)), paste("S"), 
                                   paste("SA"), round(res_SA_1_S_3$par$Linf,3), round(res_SA_1_S_3$par$t_anchor,3), 
                                   round(res_SA_1_S_3$par$K,3), round(res_SA_1_S_3$Rn_max,3)),
                                 c(paste("SA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                   paste("SA"), round(res_SA_2_3$par$Linf,3), round(res_SA_2_3$par$t_anchor,3), 
                                   round(res_SA_2_3$par$K,3), round(res_SA_2_3$Rn_max,3)),
                                 c(paste("SA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                   paste("SA"), round(res_SA_2_S_3$par$Linf,3), round(res_SA_2_S_3$par$t_anchor,3), 
                                   round(res_SA_2_S_3$par$K,3), round(res_SA_2_S_3$Rn_max,3)),
                                 c(paste("GA 1"), paste(round(Linf_1_3[1], 3), round(Linf_1_3[2], 3)), paste("NS"), 
                                   paste("GA"), round(res_GA_1_3$par$Linf,3), round(res_GA_1_3$par$t_anchor,3), 
                                   round(res_GA_1_3$par$K,3), round(res_GA_1_3$Rn_max,3)),
                                 c(paste("GA 1S"), paste(round(Linf_1_3[1], 3), round(Linf_1_3[2], 3)), paste("S"), 
                                   paste("GA"), round(res_GA_1_S_3$par$Linf,3), round(res_GA_1_S_3$par$t_ancho,3), 
                                   round(res_GA_1_S_3$par$K,3), round(res_GA_1_S_3$Rn_max,3)),
                                 c(paste("GA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                   paste("GA"), round(res_GA_2_3$par$Linf,3), round(res_GA_2_3$par$t_anchor,3), 
                                   round(res_GA_2_3$par$K,3), round(res_GA_2_3$Rn_max,3)),
                                 c(paste("GA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                   paste("GA"), round(res_GA_2_S_3$par$Linf,3), round(res_GA_2_S_3$par$t_anchor,3), 
                                   round(res_GA_2_S_3$par$K,3), round(res_GA_2_S_3$Rn_max,3))))
names(res_all_3) <- c("Model", "Linf_range", "Seasonalized", "SA_GA", "Linf", "t_anchor", "K", "Rn_max")
print(res_all_3) #SA1S
#
plot(lfq4_3, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3, par = res_SA_1_S_3$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
#Save manually: (path = "Output/Trials", filename = "Elefan_b6m3_SA1S_fit.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK_3 <- vector("list", length(lfq4_3$dates))
for(i in 1:length(lfq4_3$dates)){
  loop_data <- list(dates = lfq4_3$dates[-i],
                    midLengths = lfq4_3$midLengths,
                    catch = lfq4_3$catch[,-i])
  tmp <- ELEFAN_SA(loop_data, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                   MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                   init_par = list(Linf = Linf_1_mid_3, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                   low_par = list(Linf = Linf_1_3[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_1_3[2], K = 1, t_anchor = 1, C = 1, ts = 1))
  
  
  JK_3[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres_3 <- do.call(cbind, JK_3)
# mean
JKmeans_3 <- apply(as.matrix(JKres_3), MARGIN = 1, FUN = mean)
#confidence intervals
JKconf_3 <- apply(as.matrix(JKres_3), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JK_table_3 <- t(rbind(JKmeans_3, JKconf_3))
colnames(JK_table_3) <- c("means", "lower","upper")
# show results
JK_table_3
#
#### estimation of Z and M
Ms_3 <- M_empirical(Linf = res_SA_1_S_3$par$Linf, K_l = res_SA_1_S_3$par$K, method = "Then_growth")
lfq4_3$M <- as.numeric(Ms_3)
paste("M =", as.numeric(Ms_3))
#> [1] "M = 0.455
res_SA_1_S_3$agemax
#
#Figure of model output
print(res_SA_1_S_3$par) #print to copy numbers into figure
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=94.31666, K=0.3818931, C=0.6742018, ts = 0.6676365, t0 = 0.1958189)),
           "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 25))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_SA_1_S_3$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+ theme_classic()
#
ggsave(path = "Output/Trials/", filename = "b6m3_SA1S_VB_94.32_0.381_0.196.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4_3, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3, par = res_SA_1_S_3$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
L_3 <- seq(0,100,0.1)
t_3 <- VBGF(L = seq(0,100,0.1), list(Linf=94.31666, K=0.3818931, C=0.6742018, ts = 0.6676365, t0 = 0.1958189)) 
plot(t_3, L_3, t="l", xlim = c(0, (res_SA_1_S_3$agemax +  5)), ylim = c(0, 100))
abline(v = res_SA_1_S_3$agemax , col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Trials", filename = "b6m3_SA1S_model_Age_8.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
#
#
####ELEFAN run - b4, MA3####
set_MA <- c(3)
###Create lfq structure
lfq2_3b <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3_3b <- lfqModify(lfq2_3b, bin_size = 4) # modify bin size
lfq4_3b <- lfqRestructure(lfq3_3b, MA = set_MA, addl.sqrt = FALSE) # restructuring process
#
#Plot bins 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_3b, Fname = "catch", date.axis = "modern")
plot(lfq4_3b, Fname = "rcounts", date.axis = "modern")
#Manually save:(file = "Output/Trials/Elefan_trials_bin4_MA3_LFrestructure.tiff", width = 1200) 
par(opar)
par(mfrow = c(1,1))
#
#
####Preliminary estimates of Linf
#Powell-Wetherall plot - 1:16
res_PW_3b <- powell_wetherall(param = lfq4_3b,
                             catch_columns = 1:ncol(lfq2_3b$catch))
#
paste("Linf =",round(res_PW_3b$Linf_est), "+/-", round(res_PW_3b$se_Linf)) # show results
Linf_pw_3b <- res_PW_3b$Linf_est
#
#Max length from raw data
(Linf_ml <- max(Elefant$SL, na.rm = T))
# Froese and Binohlan formula
(Linf_pb <- exp(0.44 + 0.984*log(Linf_ml)))

#Linf range 1
Linf_1_3b <- c(Linf_pw_3b, Linf_pb)
Linf_1_mid_3b <- Linf_1_3b[1] + ((Linf_1_3b[2] - Linf_1_3b[1])/2)

#Linf range 2
Linf_2 <- c(Linf_ml, Linf_pb)
Linf_2_mid <- Linf_2[1] + ((Linf_2[2] - Linf_2[1])/2)

# tanchor may be from Mar-Jul
t_an <- c(3/12, 4/12, 5/12, 6/12, 7/12)
t_an_d <- c(0.162, 0.247, 0.329, 0.411, 0.496)
#
#
#### ELEFAN analysis with simulated annealing
# run ELEFAN_SA with range 1, NOT seasonalized
set.seed(1)
res_SA_1_3b <- ELEFAN_SA(lfq4_3b, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_1_mid_3b, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_1_3b[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_3b[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2_3b <- ELEFAN_SA(lfq4_3b, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S_3b <- ELEFAN_SA(lfq4_3b, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = Linf_1_mid_3b, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                          low_par = list(Linf = Linf_1_3b[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1_3b[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S_3b <- ELEFAN_SA(lfq4_3b, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                          low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1_3b <- ELEFAN_GA(lfq4_3b, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_1_3b[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1_3b[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2_3b <- ELEFAN_GA(lfq4_3b, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S_3b <- ELEFAN_GA(lfq4_3b, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                          addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                          low_par = list(Linf = Linf_1_3b[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1_3b[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S_3b <- ELEFAN_GA(lfq4_3b, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                          addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                          low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### plot LFQ and growth curves
opar <- par(mfrow = c(4,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_3b, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq4_3b, par = res_SA_1_3b$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq4_3b, par = res_GA_1_3b$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_3b, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3b, par = res_SA_2_3b$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_3b, par = res_GA_2_3b$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_3b, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3b, par = res_SA_1_S_3b$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_3b, par = res_GA_1_S_3b$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_3b, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3b, par = res_SA_2_S_3b$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_3b, par = res_GA_2_S_3b$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
#Save file manually: tiff("Trials/Elefan_b4M3_model fits.tiff", width = 1200)
par(opar)
par(mfrow = c(1,1))
#
#Combine results for comparison
res_all_3b <- as.data.frame(rbind(c(paste("SA 1"), paste(round(Linf_1_3b[1], 3), round(Linf_1_3b[2],3)), paste("NS"), 
                                   paste("SA"), round(res_SA_1_3b$par$Linf,3), round(res_SA_1_3b$par$t_anchor,3), 
                                   round(res_SA_1_3b$par$K,3), round(res_SA_1_3b$Rn_max,3)),
                                 c(paste("SA 1S"), paste(round(Linf_1_3b[1], 3), round(Linf_1_3b[2], 3)), paste("S"), 
                                   paste("SA"), round(res_SA_1_S_3b$par$Linf,3), round(res_SA_1_S_3b$par$t_anchor,3), 
                                   round(res_SA_1_S_3b$par$K,3), round(res_SA_1_S_3b$Rn_max,3)),
                                 c(paste("SA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                   paste("SA"), round(res_SA_2_3b$par$Linf,3), round(res_SA_2_3b$par$t_anchor,3), 
                                   round(res_SA_2_3b$par$K,3), round(res_SA_2_3b$Rn_max,3)),
                                 c(paste("SA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                   paste("SA"), round(res_SA_2_S_3b$par$Linf,3), round(res_SA_2_S_3b$par$t_anchor,3), 
                                   round(res_SA_2_S_3b$par$K,3), round(res_SA_2_S_3b$Rn_max,3)),
                                 c(paste("GA 1"), paste(round(Linf_1_3b[1], 3), round(Linf_1_3b[2], 3)), paste("NS"), 
                                   paste("GA"), round(res_GA_1_3b$par$Linf,3), round(res_GA_1_3b$par$t_anchor,3), 
                                   round(res_GA_1_3b$par$K,3), round(res_GA_1_3b$Rn_max,3)),
                                 c(paste("GA 1S"), paste(round(Linf_1_3b[1], 3), round(Linf_1_3b[2], 3)), paste("S"), 
                                   paste("GA"), round(res_GA_1_S_3b$par$Linf,3), round(res_GA_1_S_3b$par$t_ancho,3), 
                                   round(res_GA_1_S_3b$par$K,3), round(res_GA_1_S_3b$Rn_max,3)),
                                 c(paste("GA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                   paste("GA"), round(res_GA_2_3b$par$Linf,3), round(res_GA_2_3b$par$t_anchor,3), 
                                   round(res_GA_2_3b$par$K,3), round(res_GA_2_3b$Rn_max,3)),
                                 c(paste("GA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                   paste("GA"), round(res_GA_2_S_3b$par$Linf,3), round(res_GA_2_S_3b$par$t_anchor,3), 
                                   round(res_GA_2_S_3b$par$K,3), round(res_GA_2_S_3b$Rn_max,3))))
names(res_all_3b) <- c("Model", "Linf_range", "Seasonalized", "SA_GA", "Linf", "t_anchor", "K", "Rn_max")
print(res_all_3b) #SA1S
#
plot(lfq4_3b, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3b, par = res_GA_1_S_3b$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
#Save manually: (path = "Output/Trials", filename = "Elefan_b4m3_GA1S_fit.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK_3b <- vector("list", length(lfq4_3b$dates))
for(i in 1:length(lfq4_3b$dates)){
  loop_data <- list(dates = lfq4_3b$dates[-i],
                    midLengths = lfq4_3b$midLengths,
                    catch = lfq4_3b$catch[,-i])
  tmp <- ELEFAN_GA(loop_data, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                   addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                   low_par = list(Linf = Linf_1_3b[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_1_3b[2], K = 1, t_anchor = 1, C = 1, ts = 1))
  
  
  JK_3b[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres_3b <- do.call(cbind, JK_3b)
# mean
JKmeans_3b <- apply(as.matrix(JKres_3b), MARGIN = 1, FUN = mean)
#confidence intervals
JKconf_3b <- apply(as.matrix(JKres_3b), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JK_table_3b <- t(rbind(JKmeans_3b, JKconf_3b))
colnames(JK_table_3b) <- c("means", "lower","upper")
# show results
JK_table_3b
#
#### estimation of Z and M
Ms_3b <- M_empirical(Linf = res_GA_1_S_3b$par$Linf, K_l = res_GA_1_S_3b$par$K, method = "Then_growth")
lfq4_3b$M <- as.numeric(Ms_3b)
paste("M =", as.numeric(Ms_3b))
#> [1] "M = 0.455
res_GA_1_S_3b$agemax
#
#Figure of model output
print(res_GA_1_S_3b$par) #print to copy numbers into figure
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=104.6447, K=0.1052683, C=0.7116502, ts = 0.7863794, t0 = 0.523758)),
           "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 35))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_GA_1_S_3b$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+ theme_classic()
#
ggsave(path = "Output/Trials/", filename = "b4m3_GA1S_VB_104.64_0.105_0.524.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4_3b, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_3b, par = res_GA_1_S_3b$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
L_3b <- seq(0,100,0.1)
t_3b <- VBGF(L = seq(0,100,0.1), list(Linf=104.6447, K=0.1052683, C=0.7116502, ts = 0.7863794, t0 = 0.523758)) 
plot(t_3b, L_3b, t="l", xlim = c(0, (res_GA_1_S_3b$agemax +  5)), ylim = c(0, 100))
abline(v = res_GA_1_S_3b$agemax , col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Trials", filename = "b4m3_GA1S_model_Age_29.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
#
#
####ELEFAN run - b2, MA7####
set_MA <- c(7)
###Create lfq structure
lfq2_7 <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3_7 <- lfqModify(lfq2_7, bin_size = 2) # modify bin size
lfq4_7 <- lfqRestructure(lfq3_7, MA = set_MA, addl.sqrt = FALSE) # restructuring process
#
#Plot bins 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_7, Fname = "catch", date.axis = "modern")
plot(lfq4_7, Fname = "rcounts", date.axis = "modern")
#Manually save:(file = "Output/Trials/Elefan_trials_bin2_MA7_LFrestructure.tiff", width = 1200) 
par(opar)
par(mfrow = c(1,1))
#
#
####Preliminary estimates of Linf
#Powell-Wetherall plot - 1:31
res_PW_7 <- powell_wetherall(param = lfq4_7,
                              catch_columns = 1:ncol(lfq2_7$catch))
#
paste("Linf =",round(res_PW_7$Linf_est), "+/-", round(res_PW_7$se_Linf)) # show results
Linf_pw_7 <- res_PW_7$Linf_est
#
#Max length from raw data
(Linf_ml <- max(Elefant$SL, na.rm = T))
# Froese and Binohlan formula
(Linf_pb <- exp(0.44 + 0.984*log(Linf_ml)))

#Linf range 1
Linf_1_7 <- c(Linf_pw_7, Linf_pb)
Linf_1_mid_7 <- Linf_1_7[1] + ((Linf_1_7[2] - Linf_1_7[1])/2)

#Linf range 2
Linf_2 <- c(Linf_ml, Linf_pb)
Linf_2_mid <- Linf_2[1] + ((Linf_2[2] - Linf_2[1])/2)

# tanchor may be from Mar-Jul
t_an <- c(3/12, 4/12, 5/12, 6/12, 7/12)
t_an_d <- c(0.162, 0.247, 0.329, 0.411, 0.496)
#
#
#### ELEFAN analysis with simulated annealing
# run ELEFAN_SA with range 1, NOT seasonalized
set.seed(1)
res_SA_1_7 <- ELEFAN_SA(lfq4_7, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                         MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                         init_par = list(Linf = Linf_1_mid_7, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                         low_par = list(Linf = Linf_1_7[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                         up_par = list(Linf = Linf_1_7[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2_7 <- ELEFAN_SA(lfq4_7, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                         MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                         init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                         low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                         up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S_7 <- ELEFAN_SA(lfq4_7, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                           MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                           init_par = list(Linf = Linf_1_mid_7, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                           low_par = list(Linf = Linf_1_7[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                           up_par = list(Linf = Linf_1_7[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S_7 <- ELEFAN_SA(lfq4_7, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                           MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                           init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                           low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                           up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1_7 <- ELEFAN_GA(lfq4_7, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                         addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                         low_par = list(Linf = Linf_1_7[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                         up_par = list(Linf = Linf_1_7[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2_7 <- ELEFAN_GA(lfq4_7, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                         addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                         low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                         up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S_7 <- ELEFAN_GA(lfq4_7, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                           addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                           low_par = list(Linf = Linf_1_7[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                           up_par = list(Linf = Linf_1_7[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S_7 <- ELEFAN_GA(lfq4_7, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                           addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                           low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                           up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### plot LFQ and growth curves
opar <- par(mfrow = c(4,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4_7, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq4_7, par = res_SA_1_7$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq4_7, par = res_GA_1_7$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_7, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_7, par = res_SA_2_7$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_7, par = res_GA_2_7$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_7, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_7, par = res_SA_1_S_7$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_7, par = res_GA_1_S_7$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4_7, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_7, par = res_SA_2_S_7$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2_7, par = res_GA_2_S_7$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
#Save file manually: tiff("Trials/Elefan_b2M7_model fits.tiff", width = 1200)
par(opar)
par(mfrow = c(1,1))
#
#Combine results for comparison
res_all_7 <- as.data.frame(rbind(c(paste("SA 1"), paste(round(Linf_1_7[1], 3), round(Linf_1_7[2],3)), paste("NS"), 
                                    paste("SA"), round(res_SA_1_7$par$Linf,3), round(res_SA_1_7$par$t_anchor,3), 
                                    round(res_SA_1_7$par$K,3), round(res_SA_1_7$Rn_max,3)),
                                  c(paste("SA 1S"), paste(round(Linf_1_7[1], 3), round(Linf_1_7[2], 3)), paste("S"), 
                                    paste("SA"), round(res_SA_1_S_7$par$Linf,3), round(res_SA_1_S_7$par$t_anchor,3), 
                                    round(res_SA_1_S_7$par$K,3), round(res_SA_1_S_7$Rn_max,3)),
                                  c(paste("SA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                    paste("SA"), round(res_SA_2_7$par$Linf,3), round(res_SA_2_7$par$t_anchor,3), 
                                    round(res_SA_2_7$par$K,3), round(res_SA_2_7$Rn_max,3)),
                                  c(paste("SA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                    paste("SA"), round(res_SA_2_S_7$par$Linf,3), round(res_SA_2_S_7$par$t_anchor,3), 
                                    round(res_SA_2_S_7$par$K,3), round(res_SA_2_S_7$Rn_max,3)),
                                  c(paste("GA 1"), paste(round(Linf_1_7[1], 3), round(Linf_1_7[2], 3)), paste("NS"), 
                                    paste("GA"), round(res_GA_1_7$par$Linf,3), round(res_GA_1_7$par$t_anchor,3), 
                                    round(res_GA_1_7$par$K,3), round(res_GA_1_7$Rn_max,3)),
                                  c(paste("GA 1S"), paste(round(Linf_1_7[1], 3), round(Linf_1_7[2], 3)), paste("S"), 
                                    paste("GA"), round(res_GA_1_S_7$par$Linf,3), round(res_GA_1_S_7$par$t_ancho,3), 
                                    round(res_GA_1_S_7$par$K,3), round(res_GA_1_S_7$Rn_max,3)),
                                  c(paste("GA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                    paste("GA"), round(res_GA_2_7$par$Linf,3), round(res_GA_2_7$par$t_anchor,3), 
                                    round(res_GA_2_7$par$K,3), round(res_GA_2_7$Rn_max,3)),
                                  c(paste("GA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                    paste("GA"), round(res_GA_2_S_7$par$Linf,3), round(res_GA_2_S_7$par$t_anchor,3), 
                                    round(res_GA_2_S_7$par$K,3), round(res_GA_2_S_7$Rn_max,3))))
names(res_all_7) <- c("Model", "Linf_range", "Seasonalized", "SA_GA", "Linf", "t_anchor", "K", "Rn_max")
print(res_all_7) #SA1S
#
plot(lfq4_7, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_7, par = res_SA_2_S_7$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
#Save manually: (path = "Output/Trials", filename = "Elefan_b2m7_SA2S_fit.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK_7 <- vector("list", length(lfq4_7$dates))
for(i in 1:length(lfq4_7$dates)){
  loop_data <- list(dates = lfq4_7$dates[-i],
                    midLengths = lfq4_7$midLengths,
                    catch = lfq4_7$catch[,-i])
  tmp <- ELEFAN_SA(loop_data, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                  MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                  init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                  low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                  up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
  
  
  JK_7[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres_7 <- do.call(cbind, JK_7)
#
# mean
JKmeans_7 <- apply(as.matrix(JKres_7), MARGIN = 1, FUN = mean)
#confidence intervals
JKconf_7 <- apply(as.matrix(JKres_7), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JK_table_7 <- t(rbind(JKmeans_7, JKconf_7))
colnames(JK_table_7) <- c("means", "lower","upper")
# show results
 
#
#### estimation of Z and M
Ms_7 <- M_empirical(Linf = res_SA_2_S_7$par$Linf, K_l = res_SA_2_S_7$par$K, method = "Then_growth")
lfq4_7$M <- as.numeric(Ms_7)
paste("M =", as.numeric(Ms_7))
#> [1] "M = 0.455
res_SA_2_S_7$agemax
#
#Figure of model output
print(res_SA_2_S_7$par) #print to copy numbers into figure
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=101.7031, K=0.05894798, C=0.8672958, ts = 0.697694, t0 = 0.2866962)),
           "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 60))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_SA_2_S_7$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+ theme_classic()
#
ggsave(path = "Output/Trials/", filename = "b2m7_SA2S_VB_101.70_0.059_0.287.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4_7, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2_7, par = res_SA_2_S_7$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
L_7 <- seq(0,100,0.1)
t_7 <- VBGF(L = seq(0,100,0.1), list(Linf=101.7031, K=0.05894798, C=0.8672958, ts = 0.697694, t0 = 0.2866962)) 
plot(t_7, L_7, t="l", xlim = c(0, (res_SA_2_S_7$agemax +  5)), ylim = c(0, 100))
abline(v = res_SA_2_S_7$agemax , col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Trials", filename = "b2m7_SA2S_model_Age_51.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
#
#