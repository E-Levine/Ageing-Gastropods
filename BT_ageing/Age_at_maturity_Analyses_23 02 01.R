###Age at maturity for Banded tulips
#
#ELEFAN methods
#Mark-recapture 
#Wet lab growth
#Histological analyses
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
#########Wet Lab
WetLab <- read.csv("../CSV/WetLab/WetLab_BT_growth_final.csv", na.string = c("Z", "", "NA")) %>% drop_na(Date) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
         MonYr = as.yearmon(Date, format = "%m/%y")) %>% filter(Species == "BT")
#
head(WetLab)
summary(WetLab)
#
WetLab <- WetLab %>% subset(Date < "2021-02-28") %>% mutate(IDNumber = as.numeric(substr(ID, 3, nchar(ID)-1)))
#
#
#########Histology
histology <- read.csv("../CSV/Histology/Histo_combined_final.csv", na.string = c("Z", "", "NA")) %>%
  subset(Species == "BT")
head(histology)
#
histo <- histology %>% mutate(Date = as.Date(Collection_Date, format = "%m/%d/%Y"),
                              MonYr = as.yearmon(Date, format = "%m/%Y"),
                              SLclass = cut(histology$SL, breaks = seq(0, 5*ceiling(max(histology$SL, na.rm = T)/5), by = 5)),
                              OpL = as.numeric(Opercula_L),
                              OpW = as.numeric(Opercula_W),
                              Stage = as.factor(Stage)) %>%
  dplyr::select(Date, MonYr, Year, Month, Species:Penis, SLclass, MF_Final, Stage, Mature, OpL, OpW) %>% 
  drop_na(MF_Final) %>% #keep results to samples used for histology
  subset(Date < "2021-03-01")
#
#
#
###Figure formatting####
#
###Table properties
set_flextable_defaults(
  font.family = "Times New Roman",
  font.size = 12, theme_fun = theme_apa,
  padding = 0,
  background.color = "white")
#
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.5, height = 11
  ),
  type = "continuous",
  page_margins = page_mar(1)
)
#
#
Base <- theme(panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(size = 14, color = "black", family = "serif"),
              axis.text.x = element_text(size = 13, color = "black", 
                                         margin = unit(c(0.5, 0.5, 0, 0.5), "cm"), family = "serif"),
              axis.text.y = element_text(size = 14, color = "black", 
                                         margin = unit(c(0, 0.5, 0, 0), "cm"), family = "serif"),
              axis.ticks.length = unit(-0.15, "cm"))
Base2 <- theme(panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(size = 12, color = "black", family = "serif"),
              axis.text.x = element_text(size = 9, color = "black", 
                                         margin = unit(c(0.5, 0.5, 0, 0.5), "cm"), family = "serif"),
              axis.text.y = element_text(size = 8, color = "black", 
                                         margin = unit(c(0, 0.5, 0, 0), "cm"), family = "serif"),
              axis.ticks.length = unit(-0.15, "cm"))
XCate <- theme(axis.title.x = element_blank(),
               axis.text.x = element_text(color = "black", size = 14, family = "serif",
                                          margin = unit(c(0.5, 0.5, 0, 0.5), "cm")))
theme_f <- theme(strip.text.y = element_text(color = "black", size = 11, family = "serif", face = "bold"),
                 strip.background = element_rect(fill = "#CCCCCC"),
                 panel.spacing = unit(0.75, "lines"),
                 strip.text.x = element_text(size = 10, face = "bold", family = "serif"))
#
Sites <- c("WE" = "Weedon", "PP" = "Pinellas")
Stations <- c("OY" = "Oyster Reef", "SG" = "Seagrass", "SS" = "Soft Sediment")
Spp <- c("BT" = "Banded T", "CC" = "Crown C", "GM" = "Green M", "LW" = "Lightning W", "MC" = "Other 1", "MG" = "Other 2", "PW" = "Pear W", "TT" = "True T")
#
Stages <- c("0" = "Undeveloped", "1" = "Early Developing", "2" = "Late Developing", "3" = "Mature", "4" = "Recovering")
Sex <- c("M" = "Male", "F" = "Female", "U" = "Undetermined")
Months <- c("1" = "Jan", "2" = "Feb", "3" = "Mar", "4" = "Apr", "5" = "May", "6" = "Jun",
            "7" = "Jul", "8" = "Aug", "9" = "Sep", "10" = "Oct", "11" = "Nov", "12" = "Dec")
#
cbPalette <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")
#Map color to Stage
names(cbPalette) <- levels(histo$Stage)
StaFill <- scale_fill_manual(name = "Stage", labels = Stages, values = cbPalette, na.value = "#999999")
#
Model_color <- scale_color_manual(name = "Model", labels = c("ELEFAN", "Recaptures", "Laboratory"),
                                  values = c("#000000", "#0072B2", "#E69F00"))
#
#
####Length frequency####
#
summary(Elefant)
(Ele_summ <- Elefant %>% get_summary_stats(SL, type = "full") %>% dplyr::select(n, min, max, mean, sd) %>% mutate(Type = "ELEFAN"))
#
Elefant %>% 
  ggplot(aes(SL))+
  geom_histogram(binwidth = 4, boundary = 0, closed = "left")+
  lemon::facet_rep_wrap(vars(Month), ncol = 2, labeller = labeller(Month = Months))+
  Base2 + theme_f + theme(panel.spacing = unit(0, "lines"), 
                          axis.text.y = element_text(size = 6, margin = margin(r = 5)), 
                          axis.text.x = element_text(size = 8, margin = margin(t = 5)))+
  scale_y_continuous("Count", expand = c(0,0), limits = c(0, 48), breaks = seq(0, 48, 16))+
  scale_x_continuous("Shell length(mm)", expand = c(0,0), limits = c(0, 100), breaks = seq(0, 100, 20))
#
ggsave(path = "Output/Figures/2023 02/", filename = "A_LFD_all_2023 02.tiff", dpi = 1000, height = 5, width = 5, unit = "in")
#
set_MA <- c(5)
###Create lfq structure
lfq2 <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3 <- lfqModify(lfq2, bin_size = 4) # modify bin size
lfq4 <- lfqRestructure(lfq3, MA = set_MA, addl.sqrt = FALSE) # restructuring process
#
#Plot bins 
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4, Fname = "catch", date.axis = "modern")
plot(lfq4, Fname = "rcounts", date.axis = "modern")
#Manually save:(file = "Output/Figures/A_Elefan_LFcatches_LFrestructure.tiff", width = 1200) 
par(opar)
par(mfrow = c(1,1))
#
#
#
####Preliminary estimates of Linf
#Powell-Wetherall plot 
res_PW <- powell_wetherall(param = lfq4,
                           catch_columns = 1:ncol(lfq2$catch))
#
paste("Linf =",round(res_PW$Linf_est), "?", round(res_PW$se_Linf)) # show results
Linf_pw <- res_PW$Linf_est
#
#Max length from raw data
(Linf_ml <- max(Elefant$SL, na.rm = T))
# Froese and Binohlan formula
(Linf_pb <- exp(0.44 + 0.984*log(Linf_ml)))

#Linf range 1
Linf_1 <- c(Linf_pw, Linf_pb)
Linf_1_mid <- Linf_1[1] + ((Linf_1[2] - Linf_1[1])/2)

#Linf range 2
Linf_2 <- c(Linf_ml, Linf_pb)
Linf_2_mid <- Linf_2[1] + ((Linf_2[2] - Linf_2[1])/2)

# tanchor may be from Mar-Jul
t_an <- c(3/12, 4/12, 5/12, 6/12, 7/12)
#
#
#### ELEFAN analysis with simulated annealing
# run ELEFAN_SA with range 1, NOT seasonalized
set.seed(1)
res_SA_1 <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 500,
                      MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                      init_par = list(Linf = Linf_1_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                      low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2 <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 500,
                      MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                      init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                      low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 500,
                        MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_1_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 500,
                        MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1 <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = FALSE, maxiter = 50, pmutation = 0.2, popSize = 60, 
                      addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                      low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2 <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = FALSE, maxiter = 50, pmutation = 0.2, popSize = 60, 
                      addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                      low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = TRUE, maxiter = 50, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = TRUE, maxiter = 50, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### plot LFQ and growth curves
opar <- par(mfrow = c(4,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq4, par = res_SA_1$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq4, par = res_GA_1$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2, par = res_SA_2$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2, par = res_GA_2$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2, par = res_SA_1_S$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2, par = res_GA_1_S$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2, par = res_SA_2_S$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq2, par = res_GA_2_S$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
#
#Save file manually: tiff("A_Elefan_b4M5_model fits.tiff", width = 1200)
par(opar)
par(mfrow = c(1,1))
#
#Combine results for comparison
res_all <- as.data.frame(rbind(c(paste("SA 1"), paste(round(Linf_1[1], 3), round(Linf_1[2],3)), paste("NS"), 
                                 paste("SA"), round(res_SA_1$par$Linf,3), round(res_SA_1$par$t_anchor,3), 
                                 round(res_SA_1$par$K,3), round(res_SA_1$Rn_max,3)),
                               c(paste("SA 1S"), paste(round(Linf_1[1], 3), round(Linf_1[2], 3)), paste("S"), 
                                 paste("SA"), round(res_SA_1_S$par$Linf,3), round(res_SA_1_S$par$t_anchor,3), 
                                 round(res_SA_1_S$par$K,3), round(res_SA_1_S$Rn_max,3)),
                               c(paste("SA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                 paste("SA"), round(res_SA_2$par$Linf,3), round(res_SA_2$par$t_anchor,3), 
                                 round(res_SA_2$par$K,3), round(res_SA_2$Rn_max,3)),
                               c(paste("SA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                 paste("SA"), round(res_SA_2_S$par$Linf,3), round(res_SA_2_S$par$t_anchor,3), 
                                 round(res_SA_2_S$par$K,3), round(res_SA_2_S$Rn_max,3)),
                               c(paste("GA 1"), paste(round(Linf_1[1], 3), round(Linf_1[2], 3)), paste("NS"), 
                                 paste("GA"), round(res_GA_1$par$Linf,3), round(res_GA_1$par$t_anchor,3), 
                                 round(res_GA_1$par$K,3), round(res_GA_1$Rn_max,3)),
                               c(paste("GA 1S"), paste(round(Linf_1[1], 3), round(Linf_1[2], 3)), paste("S"), 
                                 paste("GA"), round(res_GA_1_S$par$Linf,3), round(res_GA_1_S$par$t_ancho,3), 
                                 round(res_GA_1_S$par$K,3), round(res_GA_1_S$Rn_max,3)),
                               c(paste("GA 2"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("NS"), 
                                 paste("GA"), round(res_GA_2$par$Linf,3), round(res_GA_2$par$t_anchor,3), 
                                 round(res_GA_2$par$K,3), round(res_GA_2$Rn_max,3)),
                               c(paste("GA 2S"), paste(round(Linf_2[1], 3), round(Linf_2[2], 3)), paste("S"), 
                                 paste("GA"), round(res_GA_2_S$par$Linf,3), round(res_GA_2_S$par$t_anchor,3), 
                                 round(res_GA_2_S$par$K,3), round(res_GA_2_S$Rn_max,3))))
names(res_all) <- c("Model", "Linf_range", "Seasonalized", "SA_GA", "Linf", "t_anchor", "K", "Rn_max")
print(res_all) #GA2S is best (m5,b4)
#
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2, par = res_GA_2_S$par,
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
#Save manually: (path = "Output/Figures", filename = "A_Elefan_final model fit_GA2S.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK <- vector("list", length(lfq4$dates))
for(i in 1:length(lfq4$dates)){
  loop_data <- list(dates = lfq4$dates[-i],
                    midLengths = lfq4$midLengths,
                    catch = lfq4$catch[,-i])
  tmp <- ELEFAN_GA(loop_data, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 100, 
                   addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                   low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
  
  
  JK[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres <- do.call(cbind, JK)
# mean
JKmeans <- apply(as.matrix(JKres), MARGIN = 1, FUN = mean)
#confidence intervals
JKconf <- apply(as.matrix(JKres), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JK_table <- t(rbind(JKmeans, JKconf))
colnames(JK_table) <- c("means", "lower","upper")
# show results
JK_table
#
#### estimation of Z and M
Ms <- M_empirical(Linf = res_GA_2_S$par$Linf, K_l = res_GA_2_S$par$K, method = "Then_growth")
lfq4$M <- as.numeric(Ms)
paste("M =", as.numeric(Ms))
#> [1] "M = 0.286"
res_GA_2_S$agemax
#
#Figure of model output
print(res_GA_2_S$par) #print to copy numbers into figure
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=100.6423, K=0.2087536, C=0.7337414, ts = 0.6103143, t0 = 0.7354528)),
              "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 25))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_GA_2_S$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+
  Base
#
ggsave(path = "Output/Figures/2023 02/", filename = "A_JKGA2S_VB_100.642_0.209_0.735_2023 02.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2, par = res_GA_1_S$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
plot(t, L, t="l", xlim = c(0, (res_GA_1_S$agemax +  5)), ylim = c(0, 100))
abline(v = 15, col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Figures", filename = "A_Elefan_GA2S_fit_model output.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
###Bootstrapped values
BT <- ELEFAN_GA_boot(lfq4, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 100, 
                   addl.sqrt = FALSE, parallel = FALSE, nresamp = 500, 
                   low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#Get means and CIs for variables
(BT_table <- t(rbind(cbind(apply(matrix(BT$bootRaw$Linf, 1, byrow = T) , MARGIN = 1, FUN = mean),
              apply(matrix(BT$bootRaw$K, 1, byrow = T) , MARGIN = 1, FUN = mean),
              apply(matrix(BT$bootRaw$t_anchor, 1, byrow = T) , MARGIN = 1, FUN = mean),
              apply(matrix(BT$bootRaw$C, 1, byrow = T) , MARGIN = 1, FUN = mean),
              apply(matrix(BT$bootRaw$ts, 1, byrow = T) , MARGIN = 1, FUN = mean)), 
        cbind(apply(matrix(BT$bootRaw$Linf, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
              apply(matrix(BT$bootRaw$K, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
              apply(matrix(BT$bootRaw$t_anchor, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
              apply(matrix(BT$bootRaw$C, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
              apply(matrix(BT$bootRaw$ts, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))))))
#
rownames(BT_table) <- c("Linf", "K", "t_anchor", "C", "ts")
colnames(BT_table) <- c("means", "lower","upper")
BT_table
#
#
#
##Age-length estimations
#Save age-length estimations using bootstrapped values
BT_table
t <- VBGF(L = seq(0,100,0.1), list(Linf=101.2241442, K=0.1840055, C=0.6099577, ts = 0.5401681, t0 = 0.5166941))
best_ageLength <- data.frame(Age = t,
                             SL = seq(0,100,0.1),
                             Age_n = round(t, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_ranges <- best_ageLength %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
#
#
#
####Mark Recapture####
#
#####################################SUMMARY
MarkRecap <- full_join(Station_MR, Extra_MR) %>% 
  mutate(ID = as.factor(ID), #Make sure ID is factor for summarizing
         Helper = ifelse(is.na(Helper), paste(Site, Station, sep = ""), Helper),
         Date = as.Date(Date, format = "%m/%d/%Y")) #Add missing helper values
head(MarkRecap)
#
##Data summary - sizes 
(BT_summ <- rbind(cbind(MarkRecap %>% distinct(ID, .keep_all = TRUE) %>% get_summary_stats(SL, type = "full") %>% 
                          dplyr::select(min, max, mean, sd) %>% mutate(Type = "Tagged"), #Summary of tagged BT
                        MarkRecap %>% distinct(ID, .keep_all = TRUE) %>% tally()) %>% #Number of tagged BT
                    dplyr::select(Type, n, everything()),
                  cbind(MarkRecap %>% filter(is.na(ID)) %>% get_summary_stats(SL, type = "full") %>% 
                          dplyr::select(min, max, mean, sd) %>% mutate(Type = "Not Tagged"),
                        MarkRecap %>% filter(is.na(ID)) %>% tally()) %>% #Number of not tagged BT
                    dplyr::select(Type, n, everything()),
                  cbind(MarkRecap %>% get_summary_stats(SL, type = "full") %>% 
                          dplyr::select(min, max, mean, sd) %>% mutate(Type = "Total"),
                        MarkRecap %>% tally()) %>% #Number of all BT
                    dplyr::select(Type, n, everything())))
#
##Data summary - monthly counts
(BT_months <- MarkRecap %>% group_by(Year, Month, Site, Station) %>% tally() %>% 
    mutate(Station_Site = paste(Station, "_", Site, sep = "")) %>%
    ungroup() %>% dplyr::select(-Site, -Station) %>% 
    spread(Station_Site, n) %>% arrange(Year, Month))
#
##Recapture summary
(IDs <- MarkRecap %>% group_by(ID) %>% drop_na(ID) %>% arrange(Date, ID) %>%
    mutate_at(.vars = vars(Date), funs(diff = .-lag(.)))) #Calculate days between capture events
#
(Recap_summ <- rbind(IDs %>% filter(Recaptured == "Y") %>% ungroup() %>% tally() %>% mutate(Type = "Total Recaps"),
                     IDs %>% filter(Recaptured == "Y") %>% tally() %>% filter(n == 1) %>% tally() %>% mutate(Type = "1 Recap"),
                     IDs %>% filter(Recaptured == "Y") %>% tally() %>% filter(n == 2) %>% tally() %>% mutate(Type = "2 Recaps"),
                     IDs %>% filter(Recaptured == "Y") %>% tally() %>% filter(n == 3) %>% tally() %>% mutate(Type = "3 Recaps")) %>%
    dplyr::select(Type, n))
#
IDs %>% filter(!is.na(diff)) %>% arrange(desc(diff)) #Arrange by descending days to recapture
IDs %>% filter(!is.na(diff)) %>% ungroup() %>% summarise(mean(diff)) #Mean number of days to recapture - 114.6
(ID_Days <- IDs %>% filter(!is.na(diff)) %>% mutate(binDays = lencat(as.numeric(diff), 0, w = 10)) %>% 
    group_by(binDays) %>% tally())

#
IDs %>% filter(!is.na(diff)) %>%
  ggplot(aes(x = diff))+
  geom_histogram(binwidth = 10, boundary = 0, closed = "right") +
  Base +
  scale_x_continuous("Number of days", expand = c(0,0), limits = c(0, 442), breaks = seq(0, 440, 40))+
  scale_y_continuous("Number of recaptures", expand = c(0,0), limits = c(0, 12))
#
ggsave(path = "Output/Figures/2023 02", filename = "B_Recaptures_daysSince.tiff", dpi = 1000)
#
#######################################ANALYSES
#
#Recaptured snails
(Recaps <- MarkRecap %>% drop_na(ID) %>% group_by(ID) %>% 
    filter(n() > 1) %>% arrange(ID, Date) %>% #Filter to all IDs with recapture and organize by date
    mutate_at(.vars = vars(Date, SL), list(diff = ~.-lag(.))) %>% #Calculate difference in days and SL
    mutate(Date_diff = ifelse(is.na(Date_diff), 0, Date_diff), #Add in 0s for starting date
           Rate = as.numeric(SL_diff)/as.numeric(Date_diff),
           binSL = lencat(SL, 0, w = 5)))#Calculate rate (mm/day) per indiviudal
#
ungroup(Recaps) %>% filter(Rate <= 0) %>% tally() #Number with 0 or negative growth - 21
ungroup(Recaps) %>% summarize(mean = mean(Rate, na.rm = T)) #All average - 0.0306 mm/day
ungroup(Recaps) %>% filter(Rate > 0) %>% summarize(mean = mean(Rate, na.rm = T)) #>0 (non-negative) - 0.0459 mm/day
#
#
#Plot observed growth
Recaps %>% 
  ggplot(aes(Date, SL, color = ID, group = ID))+
  geom_point()+
  geom_line()+
  scale_x_date("", labels = date_format("%b %Y"), breaks = date_breaks("12 weeks"),
               limits = c(as.Date("2018-10-01"), as.Date("2020-10-31")), expand = c(0,0))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 100)) +
  Base + XCate + theme(legend.position = "none")
#
ggsave(path = "Output/Figures/2023 02", filename = "B_Recaptures_Observed growth.tiff", dpi = 1000)
#
#
##Binned growth
Re_bin_growth <- Recaps %>% group_by(binSL) %>% #Group bins
  mutate(Rate = ifelse(Rate < 0, 0, Rate)) %>% drop_na(binSL) %>% #Limit to actual growth/negative as zero growth
  summarise(meanRate = mean(Rate, na.rm = T),
            sdRate = sd(Rate,na.rm = T))
#
Re_bin_growth %>%
  ggplot(aes(binSL, meanRate))+
  geom_point(position = position_nudge(x = -2.5))+
  geom_errorbar(aes(ymin = meanRate - sdRate, ymax = meanRate + sdRate),
                width = 1, position = position_nudge(x = -2.5))+
  scale_y_continuous(name = "Average growth rate (mm/day)", limits = c(-0.04, 0.12), expand = c(0,0))+
  scale_x_continuous(name = "Shell length (mm)", limits = c(0, 100), expand = c(0,0),
                     breaks = seq(0, 100, 5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  Base 
#
ggsave(path = "Output/Figures/2023 02", filename = "B_Recaptures_Ave growth binned.tiff", dpi = 1000)
#
#
#
##Get dataframe of starting SL, ending SL, and difference between dates
(FirstLast <- full_join(Recaps %>% slice(1) %>% dplyr::select(ID, SL, Date) %>% 
                          rename(Start_SL = SL, Start_t = Date),
                        Recaps %>% slice(n()) %>% dplyr::select(ID, SL, Date) %>% 
                          rename(End_SL = SL)) %>%
    arrange(ID, Date) %>% mutate(diff = as.numeric(Date - Start_t)) %>%
    ungroup())
#
#
#
##Round 1- All data
####Fit model to determine parameters estimates:
#
###Fit Wang model 
#Starting parameters
max(FirstLast$End_SL) #Linf = 76.5
with(FirstLast, mean((log(End_SL) - log(Start_SL))/diff)) #K = 0.0006030654 but 'minFactor'=0.000976562
Wang.start <- list(Linf = 76.5, K = 0.001, b = 0)
#
#Fit model
mvb <- vbFuns("Wang")
set.seed(54321)
MR_model_all <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start, data = FirstLast)
#
summary(MR_model_all)
1-((sum(residuals(MR_model_all)^2))/sum((FirstLast$End_SL-mean(FirstLast$End_SL, na.rm = T))^2)) #R squared = 0.8272819
AIC(MR_model_all) #AIC = 324.465
#
(MR_boot_all <- nlsBoot(MR_model_all, niter = 1000))
confint(MR_boot_all, plot = TRUE)
plot(residuals(MR_model_all)) #Variability assessment
hist(residuals(MR_model_all)) #Normality assessment
#
#Dataframe of estimates and CI
(Recap_table_all <- left_join(tidy(MR_model_all)[,1:2], 
                              data.frame(MR_boot_all$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(Recap_table_all) <- c("Param", "means", "lower", "upper")
(Recap_table_all <- Recap_table_all %>% mutate(Method = "Recaps_all"))
#
##Predict end SL of test set
predict2 <- function(x) predict(x, FirstLast)
(FirstLast <- cbind(FirstLast,
                    predict(MR_model_all, FirstLast), 
                    confint(Boot(MR_model_all, f = predict2))) %>%
    rename(predSL = 'predict(MR_model_all, FirstLast)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((FirstLast$End_SL - FirstLast$predSL)^2)) #2.982281
#Plot obs vs pred
FirstLast %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 02", filename = "B_Recaptures_All_EndSL_preds.tiff", dpi = 1000)
#
##Estimated SL ranges
t_all <- VBGF(L = seq(0,100,0.1), list(Linf=109.5512136, K=0.1761356, t0 = 0.5166941))
best_ageLength_all <- data.frame(Age = t_all,
                                SL = seq(0,100,0.1),
                                Age_n = round(t_all, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_MR_all <- best_ageLength_all %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
##Round 2 - Train and test dfs
#Divide data into training (70%) and test (30%) sets
set.seed(54321)
sample <- sample(nrow(FirstLast[1:6]), size = nrow(FirstLast[1:6]) * 0.7)
MR_train <- FirstLast[sample, 1:6]
MR_test <- FirstLast[-sample, 1:6]
#
####Fit model to determine parameters estimates:
#
###Fit Wang model 
mvb <- vbFuns("Wang")
set.seed(54321)
MR_model_tt <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start, data = MR_train)
#
summary(MR_model_tt)
1-((sum(residuals(MR_model_tt)^2))/sum((MR_train$End_SL-mean(MR_train$End_SL, na.rm = T))^2)) #R squared = 0.8201797
AIC(MR_model_tt) #AIC = 231.9138
#
(MR_boot_tt <- nlsBoot(MR_model_tt, niter = 1000))
confint(MR_boot_tt, plot = TRUE)
plot(residuals(MR_model_tt)) #Variability assessment
hist(residuals(MR_model_tt)) #Normality assessment
#
#Dataframe of estimates and CI
(Recap_table_tt <- left_join(tidy(MR_model_tt)[,1:2], 
                          data.frame(MR_boot_tt$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(Recap_table_tt) <- c("Param", "means", "lower", "upper")
(Recap_table_tt <- Recap_table_tt %>% mutate(Method = "Recaps_tt"))
#
##Predict end SL of test set
predict2 <- function(x) predict(x, MR_test)
(MR_test <- cbind(MR_test,
                 predict(MR_model_tt, MR_test), 
                 confint(Boot(MR_model_tt, f = predict2))) %>%
    rename(predSL = 'predict(MR_model_tt, MR_test)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((MR_test$End_SL - MR_test$predSL)^2)) #3.40198
#Plot obs vs pred
MR_test %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 02", filename = "B_Recaptures_tt_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges
t_tt <- VBGF(L = seq(0,100,0.1), list(Linf=109.0244963, K=0.1586193, t0 = 0.5166941))
best_ageLength_tt <- data.frame(Age = t_tt,
                                SL = seq(0,100,0.1),
                                Age_n = round(t_tt, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_MR_tt <- best_ageLength_tt %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
#
#
####Wet Lab####
#
(WL_summ <- rbind(WetLab %>% 
                    filter(Date == "2019-02-13") %>% summarise(n = n(), 
                                                               min = min(SL, na.rm = T),
                                                               max = max(SL, na.rm = T),
                                                               mean = mean(SL, na.rm = T)) %>%  
                    mutate(Type = "Start",
                           IDs = paste(WetLab %>% filter(Date == "2019-02-13") %>% distinct(IDNumber))),
                  WetLab %>% 
                    filter(Date == "2021-02-11") %>% summarise(n = n(), 
                                                               min = min(SL, na.rm = T),
                                                               max = max(SL, na.rm = T),
                                                               mean = mean(SL, na.rm = T)) %>%  
                    mutate(Type = "End",
                           IDs = paste(WetLab %>% filter(Date == "2021-02-11") %>% distinct(IDNumber)))) %>% dplyr::select(Type, everything()))
#
#
WL_IDs <- rbind(WetLab %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(1) %>% ungroup(),
                WetLab %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(n()) %>% ungroup()) %>%
  group_by(ID) %>% arrange(ID, Date) %>%
  mutate_at(.vars = vars(Date, SL), funs(diff = .-lag(.))) #Calculate days in wet lab
#
ungroup(WL_IDs) %>% distinct(IDNumber) %>% tally() #number of snails during course of study
WL_IDs %>% filter(!is.na(Date_diff)) %>% arrange(desc(Date_diff)) #Arrange by descending days in lab
WL_IDs %>% filter(!is.na(Date_diff)) %>% ungroup() %>% summarise(mean(Date_diff,na.rm = T)) #Mean number of days in wet lab - 373.67 days
(WL_IDs_Days <- WL_IDs %>% filter(!is.na(Date_diff)) %>% mutate(binDays = lencat(as.numeric(Date_diff), 0, w = 10)) %>% 
    group_by(binDays) %>% tally())
#
#
#
#######################################ANALYSES
#
#Individual growth
WL_IDs_growth <- WL_IDs %>% mutate(Date_diff = ifelse(is.na(Date_diff), 0, Date_diff), #Add in 0s for starting date
                                   Rate = as.numeric(SL_diff)/as.numeric(Date_diff),
                                   binSL = lencat(SL, 0, w = 5)) #Calculate rate (mm/day) per indiviudal
#
ungroup(WL_IDs_growth) %>% filter(Rate <= 0) %>% tally() #Number with 0 or negative growth - 9
ungroup(WL_IDs_growth) %>% summarize(mean = mean(Rate, na.rm = T)) #All average - 0.0121 mm/day
ungroup(WL_IDs_growth) %>% filter(Rate > 0) %>% summarize(mean = mean(Rate, na.rm = T)) #>0 (non-negative) - 0.0154 mm/day
#
#
##Compare treatments - SL, SL_diff, Rate
WL_treats <- WL_IDs_growth %>% mutate(Treatment = str_sub(ID, -1))
#
#Summary by treatment
WL_treats %>% group_by(Treatment) %>% summarise(count = n(),
                                                meanSL = mean(SL, na.rm = T),
                                                sdSL = sd(SL, na.rm = T))
ggboxplot(WL_treats, "Treatment", "SL") #Compare means
(SL_KW  <- kruskal.test(SL ~ Treatment, data = WL_treats)) #Kruskal-Wallis
#
#
#Summary by treatment
WL_treats %>% group_by(Treatment) %>% summarise(count = n(),
                                                meanDiff = mean(SL_diff, na.rm = T),
                                                sdDiff = sd(SL_diff, na.rm = T))
ggboxplot(WL_treats, "Treatment", "SL_diff") #Compare means
(Growth_KW <- kruskal.test(SL_diff ~ Treatment, data = WL_treats)) #Kruskal-Wallis
#
#
#Summary by treatment
WL_treats %>% group_by(Treatment) %>% summarise(count = n(),
                                                meanRate = mean(Rate, na.rm = T),
                                                sdRate = sd(Rate, na.rm = T))
ggboxplot(WL_treats, "Treatment", "Rate") #Compare means
(Rate_KW <- kruskal.test(Rate ~ Treatment, data = WL_treats)) #Kruskal-Wallis
(WL_test_summs <- rbind(data.frame(Test = "Shell Length", 
                                   tidy(SL_KW)),
                        data.frame(Test = "SL difference",
                                   tidy(Growth_KW)),
                        data.frame(Test = "Rate difference",
                                   tidy(Rate_KW))) %>%
    rename("Chi-squared" = statistic,
           "p-value" = p.value,
           "df" = parameter)) 
#
##Get dataframe of starting SL, ending SL, and difference between dates
(WL_FL <- full_join(WL_IDs %>% slice(1) %>% dplyr::select(ID, SL, Date) %>% 
                      rename(Start_SL = SL, Start_t = Date),
                    WL_IDs %>% slice(n()) %>% dplyr::select(ID, SL, Date) %>% 
                      rename(End_SL = SL)) %>%
    arrange(ID, Date) %>% mutate(diff = as.numeric(Date - Start_t)) %>%
    ungroup())
#
WL_IDs %>% 
  ggplot(aes(Date, SL, color = ID, group = ID))+
  geom_point()+
  geom_line()+
  scale_x_date("", labels = date_format("%b %Y"), breaks = date_breaks("12 weeks"),
               limits = c(as.Date("2018-12-5"), as.Date("2021-04-12")), expand = c(0,0.1))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 100)) +
  Base + XCate + theme(legend.position = "none") 
#
ggsave(path = "Output/Figures/2023 02", filename = "C_WetLab_Observed growth.tiff", dpi = 1000)
#
#
##Binned growth
WL_bin_growth <- WL_IDs_growth %>% mutate(binSL = lencat(SL, 0, w = 5)) %>% group_by(binSL) %>% #Group bins
  mutate(Rate = ifelse(Rate < 0, 0, Rate)) %>% drop_na(binSL) %>% #Limit to actual growth/negative as zero growth
  summarise(meanRate = mean(Rate, na.rm = T),
            sdRate = sd(Rate,na.rm = T))
#
WL_bin_growth %>%
  ggplot(aes(binSL, meanRate))+
  geom_point(position = position_nudge(x = -2.5))+
  geom_errorbar(aes(ymin = meanRate - sdRate, ymax = meanRate + sdRate),
                width = 1, position = position_nudge(x = -2.5))+
  scale_y_continuous(name = "Average growth rate (mm/day)", limits = c(-0.01, 0.06), expand = c(0,0))+
  scale_x_continuous(name = "Shell length (mm)", limits = c(0, 100), expand = c(0,0),
                     breaks = seq(0, 100, 5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  Base 
#
ggsave(path = "Output/Figures/2023 02", filename = "C_WetLab_Ave growth binned.tiff", dpi = 1000)
#
#
##Round 1- All data
####Fit model to determine parameters estimates:
#
###Fit Wang model 
#Starting parameters
max(WL_FL$End_SL) #Linf = 83.1
with(WL_FL, mean((log(End_SL) - log(Start_SL))/diff)) #K = 0.0006030654 but 'minFactor'=0.000976562
Wang.start_WL <- list(Linf = 83.1, K = 0.001, b = 0)
#
#Fit model
mvb <- vbFuns("Wang")
set.seed(54321)
WL_model_all <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start_WL, data = WL_FL)
#
summary(WL_model_all)
1-((sum(residuals(WL_model_all)^2))/sum((WL_FL$End_SL-mean(WL_FL$End_SL, na.rm = T))^2)) #R squared = 0.795302
AIC(WL_model_all) #AIC = 315.6769
#
(WL_boot_all <- nlsBoot(WL_model_all, niter = 1000))
confint(WL_boot_all, plot = TRUE)
plot(residuals(WL_model_all)) #Variability assessment
hist(residuals(WL_model_all)) #Normality assessment
#
#Dataframe of estimates and CI
(WL_table_all <- left_join(tidy(WL_model_all)[,1:2], 
                              data.frame(WL_boot_all$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(WL_table_all) <- c("Param", "means", "lower", "upper")
(WL_table_all <- WL_table_all %>% mutate(Method = "WetLab_all"))
#
##Predict end SL of test set
predict2_WL <- function(x) predict(x, WL_FL)
(WL_FL <- cbind(WL_FL,
                predict(WL_model_all, WL_FL), 
                confint(Boot(WL_model_all, f = predict2_WL))) %>%
    rename(predSL = 'predict(WL_model_all, WL_FL)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((WL_FL$End_SL - WL_FL$predSL)^2)) #7.386796
#Plot obs vs pred
WL_FL %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 02", filename = "C_WetLab_All_EndSL_preds.tiff", dpi = 1000)
#
##Estimated SL ranges
t_WL_all <- VBGF(L = seq(0,100,0.1), list(Linf=104.47654441, K=0.04120467, t0 = 0.5166941))
WL_ageLength_all <- data.frame(Age = t_WL_all,
                                 SL = seq(0,100,0.1),
                                 Age_n = round(t_WL_all, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_WL_all <- WL_ageLength_all %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
##Round 2 - Train and test dfs
#Divide data into training (70%) and test (30%) sets
set.seed(54321)
sample <- sample(nrow(WL_FL[1:6]), size = nrow(WL_FL[1:6]) * 0.7)
WL_train <- WL_FL[sample, 1:6]
WL_test <- WL_FL[-sample, 1:6]
#
####Fit model to determine parameters estimates:
#
###Fit Wang model 
mvb <- vbFuns("Wang")
set.seed(54321)
WL_model_tt <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start_WL, data = WL_train)
#
summary(WL_model_tt)
1-((sum(residuals(WL_model_tt)^2))/sum((WL_train$End_SL-mean(WL_train$End_SL, na.rm = T))^2)) #R squared = 0.6546132
AIC(WL_model_tt) #AIC = 225.6294
#
(WL_boot_tt <- nlsBoot(WL_model_tt, niter = 1000))
confint(WL_boot_tt, plot = TRUE)
plot(residuals(WL_model_tt)) #Variability assessment
hist(residuals(WL_model_tt)) #Normality assessment
#
#Dataframe of estimates and CI
(WL_table_tt <- left_join(tidy(WL_model_tt)[,1:2], 
                             data.frame(WL_boot_tt$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(WL_table_tt) <- c("Param", "means", "lower", "upper")
(WL_table_tt <- WL_table_tt %>% mutate(Method = "WetLab_tt"))
#
##Predict end SL of test set
predict2_WL <- function(x) predict(x, WL_test)
(WL_test <- cbind(WL_test,
                  predict(WL_model_tt, WL_test), 
                  confint(Boot(WL_model_tt, f = predict2_WL))) %>%
    rename(predSL = 'predict(WL_model_tt, WL_test)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((WL_test$End_SL - WL_test$predSL)^2)) #7.503713
#Plot obs vs pred
WL_test %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 02", filename = "C_WetLab_tt_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges
t_WL_tt <- VBGF(L = seq(0,100,0.1), list(Linf=106.22441897, K=0.03838359, t0 = 0.5166941))
WL_ageLength_tt <- data.frame(Age = t_WL_tt,
                                SL = seq(0,100,0.1),
                                Age_n = round(t_WL_tt, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_WL_tt <- WL_ageLength_tt %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
#
#
####Combined####
#
#
#Binned growth rates
Comb_bin_growth <- rbind(WL_IDs_growth %>% dplyr::select(Date, ID, SL, binSL, Rate),
      Recaps %>% dplyr::select(Date, ID, SL, binSL, Rate)) %>% group_by(binSL) %>% #Group bins
  mutate(Rate = ifelse(Rate < 0, 0, Rate)) %>% drop_na(binSL) %>% #Limit to actual growth/negative as zero growth
  summarise(meanRate = mean(Rate, na.rm = T),
            sdRate = sd(Rate,na.rm = T))
#
Comb_bin_growth %>%
  ggplot(aes(binSL, meanRate))+
  geom_point(position = position_nudge(x = -2.5))+
  geom_errorbar(aes(ymin = meanRate - sdRate, ymax = meanRate + sdRate),
                width = 1, position = position_nudge(x = -2.5))+
  scale_y_continuous(name = "Average growth rate (mm/day)", limits = c(-0.04, 0.12), expand = c(0,0))+
  scale_x_continuous(name = "Shell length (mm)", limits = c(0, 100), expand = c(0,0),
                     breaks = seq(0, 100, 5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  Base 
#
ggsave(path = "Output/Figures/2023 02", filename = "D_Combined_Ave growth binned.tiff", dpi = 1000)
#
#
##DF of all Start and End SLs
Comb_FL <- rbind(FirstLast[1:6], WL_FL[1:6])
#
#
##Round 1- All data
####Fit model to determine parameters estimates:
#
###Fit Wang model 
#Starting parameters
max(Comb_FL$End_SL) #Linf = 83.1
with(Comb_FL, mean((log(End_SL) - log(Start_SL))/diff)) #K = 0.0006030654 but 'minFactor'=0.000976562
Wang.start_C <- list(Linf = 83.1, K = 0.001, b = 0)
#
#Fit model
mvb <- vbFuns("Wang")
set.seed(54321)
C_model_all <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start_C, data = Comb_FL)
#
summary(C_model_all)
1-((sum(residuals(C_model_all)^2))/sum((Comb_FL$End_SL-mean(Comb_FL$End_SL, na.rm = T))^2)) #R squared = 0.7918797
AIC(C_model_all) #AIC = 315.6769
#
(C_boot_all <- nlsBoot(C_model_all, niter = 1000))
confint(C_boot_all, plot = TRUE)
plot(residuals(C_model_all)) #Variability assessment
hist(residuals(C_model_all)) #Normality assessment
#
#Dataframe of estimates and CI
(C_table_all <- left_join(tidy(C_model_all)[,1:2], 
                           data.frame(C_boot_all$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(C_table_all) <- c("Param", "means", "lower", "upper")
(C_table_all <- C_table_all %>% mutate(Method = "Combined_all"))
#
##Predict end SL of test set
predict2_C <- function(x) predict(x, Comb_FL)
(Comb_FL <- cbind(Comb_FL,
                predict(C_model_all, Comb_FL), 
                confint(Boot(C_model_all, f = predict2_C))) %>%
    rename(predSL = 'predict(C_model_all, Comb_FL)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((Comb_FL$End_SL - Comb_FL$predSL)^2)) #5.424732
#Plot obs vs pred
Comb_FL %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 02", filename = "D_Combined_All_EndSL_preds.tiff", dpi = 1000)
#
##Estimated SL ranges
t_C_all <- VBGF(L = seq(0,100,0.1), list(Linf=107.3617893, K=0.1434032, t0 = 0.5166941))
C_ageLength_all <- data.frame(Age = t_C_all,
                               SL = seq(0,100,0.1),
                               Age_n = round(t_C_all, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_C_all <- C_ageLength_all %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
##Round 2 - Train and test dfs
#Divide data into training (70%) and test (30%) sets
set.seed(54321)
sample <- sample(nrow(Comb_FL[1:6]), size = nrow(Comb_FL[1:6]) * 0.72)
C_train <- Comb_FL[sample, 1:6]
C_test <- Comb_FL[-sample, 1:6]
#
####Fit model to determine parameters estimates:
#
###Fit Wang model 
mvb <- vbFuns("Wang")
set.seed(54321)
C_model_tt <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start_C, data = C_train)
#
summary(C_model_tt)
1-((sum(residuals(C_model_tt)^2))/sum((C_train$End_SL-mean(C_train$End_SL, na.rm = T))^2)) #R squared = 0.8354909
AIC(C_model_tt) #AIC = 479.5319
#
(C_boot_tt <- nlsBoot(C_model_tt, niter = 1000))
confint(C_boot_tt, plot = TRUE)
plot(residuals(C_model_tt)) #Variability assessment
hist(residuals(C_model_tt)) #Normality assessment
#
#Dataframe of estimates and CI
(C_table_tt <- left_join(tidy(C_model_tt)[,1:2], 
                          data.frame(C_boot_tt$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(C_table_tt) <- c("Param", "means", "lower", "upper")
(C_table_tt <- C_table_tt %>% mutate(Method = "Combined_tt"))
#
##Predict end SL of test set
predict2_C <- function(x) predict(x, C_test)
(C_test <- cbind(C_test,
                  predict(C_model_tt, C_test), 
                  confint(Boot(C_model_tt, f = predict2_C))) %>%
    rename(predSL = 'predict(C_model_tt, C_test)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((C_test$End_SL - C_test$predSL)^2)) #6.094726
#Plot obs vs pred
C_test %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 02", filename = "D_Combined_tt_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges
t_C_tt <- VBGF(L = seq(0,100,0.1), list(Linf=107.284857, K=0.154844, t0 = 0.5166941))
C_ageLength_tt <- data.frame(Age = t_C_tt,
                              SL = seq(0,100,0.1),
                              Age_n = round(t_C_tt, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_C_tt <- C_ageLength_tt %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
#
#
#
####Histology####
#
histo <- histo %>% drop_na(Stage)
#
summary(histo)
(Histo_summ <- rbind(histo %>% drop_na(SL) %>% summarise(n = n(), 
                                                         min = min(SL),
                                                         max = max(SL),
                                                         mean = mean(SL)) %>% mutate(MF_Final = "Both"),
                     histo %>% drop_na(SL)%>% group_by(MF_Final)%>% summarise(n = n(), 
                                                                              min = min(SL),
                                                                              max = max(SL),
                                                                              mean = mean(SL))))
histo$bin10 <- lencat(histo$SL, 0, w = 10)
#
(H_samples <- histo %>% group_by(Year, Month, bin10) %>% drop_na(SL) %>% tally() %>% 
    mutate(Year_Month = as.factor(paste(Year, ".", Month, sep = ""))) %>% ungroup() %>% dplyr::select(-Year, -Month) %>%
    mutate(Year_Month = fct_relevel(Year_Month, c("2019.2", "2019.4", "2019.6", "2019.8", "2019.10", "2019.12", 
                                                  "2020.1", "2020.2", "2020.3", "2020.5", "2020.7", "2020.9",
                                                  "2020.10", "2020.11", "2021.1"))))
#
(Monthly <- full_join(histo %>% group_by(Month) %>% 
                        summarise(Sample = n(), aveSL = round(mean(SL, na.rm = T),1)) %>% mutate(MF_Final = "B"),
                      histo %>% group_by(Month, MF_Final) %>% 
                        summarise(Sample = n(), aveSL = round(mean(SL, na.rm = T),1))) %>%
    arrange(Month))
#
#Number of males and females per SL bin
histo %>%  
  ggplot(aes(x= SL, fill = MF_Final))+
  geom_histogram(aes(y = ..count..), breaks = seq(0, 90, by = 10), alpha = 0.8)+ 
  ylab("Number of Snails")+
  lemon::facet_rep_grid(MF_Final~., labeller = labeller(MF_Final = Sex)) +
  Base + scale_fill_grey(start = 0, end = 0.4)+
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_continuous(name = "Shell length (mm)", expand = c(0,0), 
                     limits = c(0, 90), breaks = seq(0, 90, by = 10))+
  scale_y_continuous(expand = c(0,0), limits = c(0,40))
#
ggsave(path = "Output/Figures/2023 02/", filename = "E_BT_MF_SLbins_BW.tiff", dpi = 1000)
#
#Stages per month
histo %>% mutate(Month = as.factor(Month)) %>% 
  ggplot(aes(Month, fill = Stage)) +
  geom_bar(position = "fill")+
  scale_y_continuous("Percent", labels = scales::percent_format(), expand = c(0,0))+
  Base + StaFill +
  scale_fill_grey(start = 0, end = 0.9)
  
#
ggsave(path = "Output/Figures/2023 02/", filename = "E_BT_Monthly_Stages_BW.tiff", dpi = 1000)
#
##Stages per sex
histo %>% mutate(Month = as.factor(Month)) %>% 
  ggplot(aes(Month, fill = Stage)) +
  geom_bar(position = "fill")+
  lemon::facet_rep_grid(MF_Final~.)+
  scale_y_continuous("Percent", labels = scales::percent_format(), expand = c(0,0))+
  Base + StaFill +
  scale_fill_grey(start = 0, end = 0.9)
#
ggsave(path = "Output/Figures/2023 02/", filename = "E_BT_Monthly_Stages_MaleFemale_BW.tiff", dpi = 1000)
#
##Stages per sex per SL bin
histo %>% drop_na(SLclass, Stage) %>%
  ggplot(aes(SLclass, fill = Stage)) +
  geom_bar(position = "fill")+
  lemon::facet_rep_grid(MF_Final~.)+
  scale_y_continuous("Percent", labels = scales::percent_format(), expand = c(0,0))+
  Base + StaFill + theme(axis.text.x = element_text(size = 10)) +
  scale_fill_grey(start = 0, end = 0.9)
#
ggsave(path = "Output/Figures/2023 02/", filename = "E_BT_SLClass_Stages_MaleFemale_BW.tiff", dpi = 1000)
#
#
####Size at maturity figure
#Need histo file
#Type 1 = all M and F in one, extra == NA; 
#Type 2 = either Male or Female, must specify extra == "M" or "F"; 
#Type 3 = Male and female with individual propMature, specify extra == 1 for 1 plot, extra == 2 for faceted
matureSL <- function(df, Spp, proportionMature, Type, extra){
  #Functions needed
  round10 <- function(x){10*ceiling(x/10)} #Rounds x up to nearest 10
  round5 <- function(x){5*ceiling(x/5)} #Rounds x up to nearest 5
  
  #df to work with - if FC, need to include "U" and use minSL M or F for maturity size
  mat_all <- df %>%  mutate(Mature = as.factor(Mature), MF_Final = as.factor(MF_Final)) %>% 
    subset(Species == Spp) %>% 
    dplyr::select(SL, SLclass, MF_Final, Mature) %>% drop_na()
  mat_M <- mat_all %>% filter(MF_Final != "F") %>% mutate(MF_Final = replace(MF_Final, MF_Final == "U", "M"))
  mat_F <- mat_all %>% filter(MF_Final != "M") %>% mutate(MF_Final = replace(MF_Final, MF_Final == "U", "F"))
  
  ##Color to sex
  Sex <- c("M" = "Male", "F" = "Female", "U" = "Undetermined")
  color_og <- c("#000000", "#666666", "#CCCCCC")
  #Map color to Sex
  names(color_og) <- c("F","M","U")
  MFCol <- scale_color_manual(name = "", labels = Sex, values = color_og) 
  
  ###ALL
  #Fit model
  lrSL <- glm(Mature ~ SL, family = binomial, data = mat_all)
  output <- data.frame(SL = seq(0, round10(max(mat_all$SL)), 5))
  output$Mature <- predict(lrSL, newdata = output, type = "response")
  #Get x where y = 0.5
  LD50 <- MASS::dose.p(lrSL, p = proportionMature)
  
  ###Males
  lrSLM <- glm(Mature ~ SL, family = binomial, data = mat_M)
  outputM <- data.frame(SL = seq(0, round10(max(mat_M$SL)), 5))
  outputM$Mature <- predict(lrSLM, newdata = outputM, type = "response")
  #Get x where y = p
  LD50M <- MASS::dose.p(lrSLM, p = proportionMature)
  
  ###Females
  lrSLF <- glm(Mature ~ SL, family = binomial, data = mat_F)
  outputF <- data.frame(SL = seq(0, round10(max(mat_F$SL)), 5))
  outputF$Mature <- predict(lrSLF, newdata = outputF, type = "response")
  #Get x where y = p
  LD50F <- MASS::dose.p(lrSLF, p = proportionMature)
  
  ##Plots
  All <- mat_all %>%
    ggplot(aes(SL, as.numeric(Mature)-1))+
    geom_point(aes(color = MF_Final), size = 3, alpha = 0.6)+
    stat_smooth(method = "glm", se = FALSE, fullrange = TRUE, 
                method.args = list(family = binomial), size = 1.25)+
    Base + XCate + MFCol +
    geom_vline(xintercept = LD50[[1]],linetype = "dashed", color = "black", size = 1)+
    scale_x_continuous(name = "Shell length (mm)", expand = c(0,0), limits = c(0, round10(max(mat_all$SL))), breaks = seq(0, round10(max(mat_all$SL)), by = 10))+
    scale_y_continuous(name = "Proportion mature", expand = c(0.025,0.025), limits = c(0,1))
  Single <- mat_all %>% subset(MF_Final == extra | MF_Final == "U") %>%
    ggplot(aes(SL, as.numeric(Mature)-1))+
    geom_point(aes(color = MF_Final), size = 3, alpha = 0.6)+
    stat_smooth(method = "glm", se = FALSE, fullrange = TRUE, 
                method.args = list(family = binomial), size = 1.25)+
    Base + XCate + MFCol +
    geom_vline(xintercept = ifelse(extra == "M", LD50M[[1]], LD50F[[1]]),linetype = "dashed", color = "black", size = 1)+
    scale_x_continuous(name = "Shell length (mm)", expand = c(0,0), limits = c(0, round10(max(mat_all$SL))), breaks = seq(0, round10(max(mat_all$SL)), by = 10))+
    scale_y_continuous(name = "Proportion mature", expand = c(0.025,0.025), limits = c(0,1))
  Both <- mat_all %>%
    ggplot(aes(SL, as.numeric(Mature)-1))+
    geom_point(aes(color = MF_Final), size = 3, alpha = 0.6)+
    stat_smooth(method = "glm", se = FALSE, fullrange = TRUE, 
                method.args = list(family = binomial), size = 1.25)+
    Base + XCate + MFCol +
    geom_vline(xintercept = LD50[[1]],linetype = "dashed", color = "black", size = 1)+
    geom_vline(xintercept = LD50M[[1]],linetype = "dashed", color = "#E69F00", size = 1)+
    geom_vline(xintercept = LD50F[[1]],linetype = "dashed", color = "#009E73", size = 1)+
    scale_x_continuous(name = "Shell length (mm)", expand = c(0,0), limits = c(0, round10(max(mat_all$SL))), breaks = seq(0, round10(max(mat_all$SL)), by = 10))+
    scale_y_continuous(name = "Proportion mature", expand = c(0.025,0.025), limits = c(0,1))
  Facet <- rbind(mat_M, mat_F) %>%
    ggplot(aes(SL, as.numeric(Mature)-1))+
    geom_point(aes(color = MF_Final), size = 3, alpha = 0.6)+
    stat_smooth(method = "glm", se = FALSE, fullrange = TRUE, 
                method.args = list(family = binomial), size = 1.25)+
    lemon::facet_rep_grid(MF_Final~., labeller = labeller(MF_Final = Sex))+
    Base + XCate + MFCol + theme_f+ theme(legend.position = "none")+
    geom_vline(data = filter(mat_all, MF_Final == "M"), aes(xintercept = LD50M[[1]]),linetype = "dashed", color = "black", size = 1)+
    geom_vline(data = filter(mat_all, MF_Final == "F"), aes(xintercept = LD50F[[1]]),linetype = "dashed", color = "black", size = 1)+
    scale_x_continuous(name = "Shell length (mm)", expand = c(0,0), limits = c(0, round10(max(mat_all$SL))), breaks = seq(0, round10(max(mat_all$SL)), by = 10))+
    scale_y_continuous(name = "Proportion mature", expand = c(0.025,0.025), limits = c(0,1))
  
  #OUTPUT
  ifelse(Type == 1,
         return(list(All, paste("All:", proportionMature,"=",LD50[1],",", "SE =",attr(LD50, "SE")))),
         ifelse(Type == 2,
                return(list(Single, ifelse(extra == "M", 
                                           paste("Males:", proportionMature,"=",LD50M[1],",", "SE =",attr(LD50M, "SE")), 
                                           paste("Females:", proportionMature,"=",LD50F[1],",", "SE =",attr(LD50F, "SE"))))),
                ifelse(Type == 3,
                       if(extra == 1){return(list(Both, 
                                                  paste("All:", proportionMature,"=",LD50[1],",", "SE =",attr(LD50, "SE")), 
                                                  paste("Males:", proportionMature,"=",LD50M[1],",", "SE =",attr(LD50M, "SE")), 
                                                  paste("Females", proportionMature,"=",LD50F[1],",", "SE =",attr(LD50F, "SE"))))}
                       else if(extra == 2){return(list(Facet, 
                                                       paste("Males:", proportionMature,"=",LD50M[1],",", "SE =",attr(LD50M, "SE")), 
                                                       paste("Females", proportionMature,"=",LD50F[1],",", "SE =",attr(LD50F, "SE"))))}
                       else{print("1 or 2 plots?")},
                       print("Need more info"))))
  
}
#
#Name SP_maturity_Prop#_PropSL_YYYY MM
matureSL(histo, "BT", 0.50, 3, 2)
#
ggsave(path = "Output/Figures/2023 02/", filename = "E_BT_sizeMaturity_M38.24_F42.16.tiff", dpi = 1000)
#
#
#
####Table outputs####
#
#(Ele_table <- Ele_summ %>% flextable())
#
(Res_tbl <- res_all %>% flextable() %>% autofit())
#
(JK_tble <- data.frame(JK_table) %>% rownames_to_column() %>% rename("Param" = rowname) %>% 
    mutate(Method = "ELEFAN-JK") %>% flextable() %>% autofit())
#
(BT_tble <- data.frame(BT_table) %>% rownames_to_column() %>% rename("Param" = rowname) %>% 
    mutate(Method = "ELEFAN-boot") %>% flextable() %>% autofit())
#
(Age_est_tble <- Age_SL_ranges %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Elefan_initial"))
#
#
#
#
##Mark-Recapture Snail summary
#(BT_summ_tbl <- BT_summ %>% 
#   flextable() %>% hline(i = 2) %>% autofit())
#
#Monthly summary
(BT_month_tbl <- BT_months %>% adorn_totals(c("col", "row"),,,,c(-Month, -Year)) %>%
    flextable() %>%
    separate_header(opts = c("center-hspan")) %>%
    colformat_num(j = c(1, 9), big.mark = "", digits = 0) %>%
    colformat_num(j = c(3:8), na_str = "0") %>%
    merge_v(j = 1) %>% vline(j = c(2, 4, 6, 8)) %>% hline(i = c(3, 15)) %>%
    autofit())
#
##Recap summary - Counts
(Recap_tbl <- Recap_summ %>% flextable() %>% autofit())
#
#Days
(MR_ID_Days_tbl <- ID_Days %>% flextable() %>% 
    colformat_num(j = 1, digits = 0) %>%
    autofit())
#
##Binned rates
(MR_bin_growth_tbl <- Re_bin_growth %>% flextable() %>%
    set_formatter("binSL" = function(x){formatC(x, format = "f", digits = 0)},
                  "meanRate" = function(x){formatC(x, format = "f", digits = 3)},
                  "sdRate" = function(x){formatC(x, format = "f", digits = 3)}) %>% autofit())
#
##Mark recapture SL model
(RecapModel_tbl <- rbind(left_join(tidy(MR_model_all) %>% mutate(Model = "Recaps-all", R2 = "0.827", AIC = "324.465", RMSE = "2.982281"),
                             data.frame(term = c("Linf", "K", "b"), LCI = confint(MR_boot_all)[,1], UCI = confint(MR_boot_all)[,2])) %>% 
                           dplyr::select(Model, everything()),
                         left_join(tidy(MR_model_tt) %>% mutate(Model = "Recaps-tt", R2 = "0.820", AIC = "231.914", RMSE = "3.40198"),
                                   data.frame(term = c("Linf", "K", "b"), LCI = confint(MR_boot_tt)[,1], UCI = confint(MR_boot_tt)[,2])) %>% 
                           dplyr::select(Model, everything())) %>%
    flextable() %>% 
    set_formatter("estimate" = function(x){formatC(x, format = "f", digits = 3)},
                  "std.error" = function(x){formatC(x, format = "f", digits = 3)},
                  "statistic" = function(x){formatC(x, format = "f", digits = 3)},
                  "p.value" = function(x){formatC(x, format = "e", digits = 3)}) %>%
    merge_v(j = c(1, 7, 8)) %>% hline(i = 3) %>%
    autofit())
#
(MR_Age_all_tble <- Age_SL_MR_all %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Recaps-All"))
#
#
(MR_Age_tt_tble <- Age_SL_MR_tt %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Recaps-tt"))
#
#
#
#
##Wet lab summaries
#(WL_tbl <- WL_summ %>% flextable() %>% autofit())
#
#Days
(WL_ID_Days_tbl <- WL_IDs_Days %>% flextable() %>% 
    colformat_num(j = 1, digits = 0) %>%
    autofit())
#
(Treatment_tble <- WL_test_summs %>% flextable() %>% autofit())
#
##Binned rates
(WL_bin_growth_tbl <- WL_bin_growth %>% flextable() %>%
    set_formatter("binSL" = function(x){formatC(x, format = "f", digits = 0)},
                  "meanRate" = function(x){formatC(x, format = "f", digits = 3)},
                  "sdRate" = function(x){formatC(x, format = "f", digits = 3)}) %>% autofit())
#
#
##Mark recapture SL model
(WLModel_tbl <- rbind(left_join(tidy(WL_model_all) %>% mutate(Model = "Wet Lab-all", R2 = "0.795", AIC = "315.677", RMSE = "7.386796"),
                                    data.frame(term = c("Linf", "K", "b"), LCI = confint(WL_boot_all)[,1], UCI = confint(WL_boot_all)[,2])) %>% 
                            dplyr::select(Model, everything()),
                          left_join(tidy(WL_model_tt) %>% mutate(Model = "Wet Lab-t", R2 = "0.655", AIC = "225.6294", RMSE = "7.503713"),
                                    data.frame(term = c("Linf", "K", "b"), LCI = confint(WL_boot_tt)[,1], UCI = confint(WL_boot_tt)[,2])) %>% 
                            dplyr::select(Model, everything())) %>%
    flextable() %>% 
    set_formatter("estimate" = function(x){formatC(x, format = "f", digits = 3)},
                  "std.error" = function(x){formatC(x, format = "f", digits = 3)},
                  "statistic" = function(x){formatC(x, format = "f", digits = 3)},
                  "p.value" = function(x){formatC(x, format = "e", digits = 3)}) %>%
    merge_v(j = c(1, 7, 8)) %>% hline(i = 3) %>%
    autofit())
#
(WL_Age_all_tble <- Age_SL_WL_all %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "WetLab-All"))
#
(WL_Age_tt_tble <- Age_SL_WL_tt %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "WetLab-tt"))
#
#
#
#
#
##Combined
##Binned rates
(Comb_bin_growth_tbl <- Comb_bin_growth %>% flextable() %>%
    set_formatter("binSL" = function(x){formatC(x, format = "f", digits = 0)},
                  "meanRate" = function(x){formatC(x, format = "f", digits = 3)},
                  "sdRate" = function(x){formatC(x, format = "f", digits = 3)}) %>% autofit())
#
##Mark recapture SL model
(CombModel_tbl <- rbind(left_join(tidy(C_model_all) %>% mutate(Model = "Combined-all", R2 = "0.792", AIC = "679.740", RMSE = "5.4247"),
                                   data.frame(term = c("Linf", "K", "b"), LCI = confint(C_boot_all)[,1], UCI = confint(C_boot_all)[,2])) %>% 
                           dplyr::select(Model, everything()),
                         left_join(tidy(C_model_tt) %>% mutate(Model = "Combined-tt", R2 = "0.835", AIC = "479.5319", RMSE = "6.094726"),
                                   data.frame(term = c("Linf", "K", "b"), LCI = confint(C_boot_tt)[,1], UCI = confint(C_boot_tt)[,2])) %>% 
                           dplyr::select(Model, everything())) %>%
    flextable() %>% 
    set_formatter("estimate" = function(x){formatC(x, format = "f", digits = 3)},
                  "std.error" = function(x){formatC(x, format = "f", digits = 3)},
                  "statistic" = function(x){formatC(x, format = "f", digits = 3)},
                  "p.value" = function(x){formatC(x, format = "e", digits = 3)}) %>%
    merge_v(j = c(1, 7, 8)) %>% hline(i = 3) %>%
    autofit())
#
(C_Age_all_tble <- Age_SL_C_all %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Combined-All"))
#
(C_Age_tt_tble <- Age_SL_C_tt %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Combined-tt"))
#
#
#
##Histology

#Ovearll size summary
(H_summ_tbl <- Histo_summ %>% dplyr::select(MF_Final, everything()) %>%
    flextable() %>% autofit())
#
#Monthly samples
(H_sam_tbl <- H_samples %>% spread(Year_Month, n) %>% 
    rename("Shell length" = bin10) %>%  adorn_totals(c('row', "col")) %>%
    flextable() %>% separate_header(opts = c("center-hspan")) %>%
    colformat_num(j = 17, digits = 0) %>%
    hline(i = c(7)) %>% vline(j = c(1, 7, 15, 16)) %>% autofit())
#
(His_mon_tble <- Monthly %>% flextable() %>% autofit() %>% merge_v(j = 1) %>%
    hline(i = c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33)))
#
#
#
#
##Final Combines
(Data_summs <- rbind(Ele_summ, BT_summ, 
                    WL_summ %>% dplyr::select(-IDs) %>% mutate(sd = NA), 
                    Histo_summ %>% rename("Type" = MF_Final) %>% mutate(sd = NA)) %>% 
  dplyr::select(Type, everything()) %>%
  flextable() %>% hline(i = c(1,4, 6)))

#
#
#
####Save output to Word####
#
save_as_docx("Data summaries" = Data_summs,
             "ELEFAN Model summaries" = Res_tbl,
             "Model means, 95% CI" = JK_tble,
             "Bootstrapped means, 95% CI" = BT_tble,
             "Elefan estimated SLs at age" = Age_est_tble,
             "Monthly Counts" = BT_month_tbl,
             "Number of days to recapture" = MR_ID_Days_tbl,
             "Recapture binned growth" = MR_bin_growth_tbl,
             "Recapture Models" = RecapModel_tbl,
             "Recapture All Data SL ranges" = MR_Age_all_tble,
             "Recapture T/T Data SL ranges" = MR_Age_tt_tble,
             "Number of days in Wet lab" = WL_ID_Days_tbl,
             "Kruskal-Wallis tests" = Treatment_tble,
             "Wet lab binned growth" = WL_bin_growth_tbl,             
             "Wet lab Models" = WLModel_tbl,
             "Wet lab All Data SL ranges" = WL_Age_all_tble,
             "Wet lab T/T Data SL ranges" = WL_Age_tt_tble,
             "Combined binned growth" = Comb_bin_growth_tbl,             
             "Combined Models" = CombModel_tbl,
             "Combined All Data SL ranges" = C_Age_all_tble,
             "Combined T/T Data SL ranges" = C_Age_tt_tble,
             "Histolgy monthly samples by SL" = H_sam_tbl,
             "Histology monthly size ranges" = His_mon_tble,
             path = "Output/Age Size at Maturity_Tables_2023 02.docx", 
             pr_section = sect_properties)
#

#
####Figures####
#
##Model estimated maxSL/age
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), aes(Age, maxSL, color = "ELEFAN", linetype = "All"))+
  #geom_line(data = Age_SL_MR_all %>% rename(Age = 'round(Age_n, 0)'), aes(Age, maxSL, color = "Recapture", linetype = "All"))+
  #geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)'), aes(Age, maxSL, color = "Recapture", linetype = "Parital"))+ 
  #geom_line(data = Age_SL_WL_all %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), aes(Age, maxSL, color = "Wet Lab", linetype = "All"))+
  #geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), aes(Age, maxSL, color = "Wet Lab", linetype = "Parital"))+
  #geom_line(data = Age_SL_C_all %>% rename(Age = 'round(Age_n, 0)'), aes(Age, maxSL, color = "Combined", linetype = "All"))+ 
  #geom_line(data = Age_SL_C_tt %>% rename(Age = 'round(Age_n, 0)'), aes(Age, maxSL, color = "Combined", linetype = "Parital"))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 100))+
  geom_hline(yintercept = 38.24, linetype = "dashed", color =  "red")+
  geom_hline(yintercept = 42.16, linetype = "dashed", color = "blue")+
  Base
#
#
#Run WL models with double growth rate. Decide All v Parital. Leaning Partial.
#Figure of raw output. Figure with each model added. Figure of all models. Figure all models w/ updated WL.
#
#
##Individual figures - progression - with or without size at maturity
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), aes(Age, maxSL, color = "ELEFAN"), linewidth = 1)+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)'), aes(Age, maxSL, color = "Recaptures"), linewidth = 1)+ 
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), aes(Age, maxSL, color = "Laboratory"), linewidth = 1)+
  geom_line(data = Age_SL_C_all %>% rename(Age = 'round(Age_n, 0)'), aes(Age, maxSL, color = "Combined"), linewidth = 1)+ 
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+ 
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base + theme(legend.key = element_blank(), legend.title = element_blank()) + 
  scale_color_manual(breaks = c("ELEFAN", "Recaptures", "Laboratory", "Combined"),
                     labels = c("ELEFAN", "Recaptures", "Laboratory", "Combined"),
                     values = c("#000000", "#0072B2", "#E69F00", "#009E73"))+
  geom_hline(yintercept = 38.24, linetype = "dashed", color =  "red")+
  geom_hline(yintercept = 42.16, linetype = "dashed", color = "blue")
  
#
ggsave(path = "Output/Figures/2023 02/", filename = "P_ELEFAN_Age_SLs_2023 02.tiff", dpi = 1000)
ggsave(path = "Output/Figures/2023 02/", filename = "P_ELE_Rec_Age_SLs_2023 02.tiff", dpi = 1000)
ggsave(path = "Output/Figures/2023 02/", filename = "P_ELE_Rec_WL_Age_SLs_2023 02.tiff", dpi = 1000)
ggsave(path = "Output/Figures/2023 02/", filename = "P_ELE_Rec_WL_Comb_Mat_2023 02.tiff", dpi = 1000)
#
##Difference between max SLs averaged over 365 days - growth rate per age
left_join(Age_SL_MR_tt %>% mutate(MR_Rate = (maxSL - lag(maxSL))/365) %>% dplyr::select(-minSL, -maxSL),
          Age_SL_WL_tt %>% mutate(WL_Rate = (maxSL - lag(maxSL))/365) %>% dplyr::select(-minSL, -maxSL))
#
WL_bin_growth %>%
  ggplot(aes(binSL, meanRate))+
  geom_point(position = position_nudge(x = -3), color = "#E69F00")+
  geom_errorbar(aes(ymin = meanRate - sdRate, ymax = meanRate + sdRate),
                width = 1, position = position_nudge(x = -3), color = "#E69F00")+
  geom_point(data = Re_bin_growth, aes(binSL, meanRate), position = position_nudge(x = -2), color = "#0072B2")+
  geom_errorbar(data = Re_bin_growth, aes(ymin = meanRate - sdRate, ymax = meanRate + sdRate),
                width = 1, position = position_nudge(x = -2), color = "#0072B2")+
  scale_y_continuous(name = "Average growth rate (mm/day)", limits = c(-0.03, 0.12), expand = c(0,0))+
  scale_x_continuous(name = "Shell length (mm)", limits = c(0, 100), expand = c(0,0),
                     breaks = seq(0, 100, 5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  Base 
#
ggsave(path = "Output/Figures/2023 02/", filename = "P_GrowthRates_Re_WL_2023 02.tiff", dpi = 1000)
