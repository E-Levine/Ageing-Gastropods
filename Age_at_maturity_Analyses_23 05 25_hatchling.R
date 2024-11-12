###Age at maturity for Banded tulips
#
#ELEFAN methods v LFD methods
#Mark-recapture (train/test) v MLE methods
#Wet lab growth - train/test
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
if (!require("pacman")) {install.packages("pacman")}
if (!require("devtools")) {install.packages("devtools")}
install.packages("https://cran.r-project.org/src/contrib/Archive/raster/raster_3.6-23.tar.gz", repos = NULL, type = "source")
devtools::install_github("stefanedwards/lemon")
devtools::install_github("briandconnelly/growthcurve", build_vignettes = TRUE)
#devtools::install_github("jonathansmart/BayesGrowth")
if (!require("remotes")) {install.packages("remotes")}
remotes::install_github("rschwamborn/fishboot")
pacman::p_load(plyr, tidyverse, rstatix, #Df manipulation, basic summary
               zoo, FSA, nlstools, #length bins, growth tags
               ggmap, ggsn, ggpubr, ggfortify, scales, lemon, raster, maps, mapdata, extrafont, #Mapping and figures
               lmPerm, car, TropFishR, lubridate, MASS, fishboot, #ELEFAN
               fishmethods, mclust, ClusterR, TMB, tmbstan, #Alternatives
               broom, flextable, janitor, officer)
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
WetLab <- read.csv("../CSV/WetLab/WetLab_BT_growth_final_all.csv", na.string = c("Z", "", "NA")) %>% drop_na(Date) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
         MonYr = as.yearmon(Date, format = "%m/%y")) %>% filter(Species == "BT")
#
head(WetLab)
summary(WetLab)
#
WetLab <- WetLab %>% subset(Date < "2021-02-28") #%>% mutate(IDNumber = as.numeric(substr(ID, 3, 5)))
#
Hatchlings <- read.csv("CSV/2017_Hatchling_growth.csv", na.string = c("Z", "", "NA")) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
         MonYr = as.yearmon(Date, format = "%m/%y"))
head(Hatchlings)
#
Juvs <- read.csv("CSV/2014_BTjuveniles.csv", na.string = c("Z", "", "NA")) %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
         MonYr = as.yearmon(Date, format = "%m/%y"))
head(Juvs)
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
##Get Arial font
font_import()
loadfonts(device = "win")
fonts()
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
##Regression equation
makeRegEqnLabel <- function(fit, x_char, y_char) {
  x_lab <- noquote(x_char)
  y_lab <- noquote(y_char)
  # Isolate coefficients (and control decimals)
  cfs <- coef(fit)
  m <- formatC(cfs[[2]], format="f", digits=3)
  # Handle b differently because of minus in the equation
  b <- cfs[[1]]
  b <- paste0(ifelse(b<0, " - ", " + "), formatC(abs(b), format="f", digits=3))
  # Put together and return
  paste0(y_lab," = ",m, x_lab, b)
}
#
Base <- theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(size = 14, color = "black", family = "Arial"),
              axis.text.x = element_text(size = 13, color = "black", 
                                         margin = unit(c(0.5, 0.5, 0, 0.5), "cm"), family = "Arial"),
              axis.text.y = element_text(size = 14, color = "black", 
                                         margin = unit(c(0, 0.5, 0, 0), "cm"), family = "Arial"),
              axis.ticks.length = unit(-0.15, "cm"), plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"))
Base2 <- theme(panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
               axis.line = element_line(color = "black"),
               axis.title = element_text(size = 12, color = "black", family = "Arial"),
               axis.text.x = element_text(size = 9, color = "black", 
                                          margin = unit(c(0.5, 0.5, 0, 0.5), "cm"), family = "Arial"),
               axis.text.y = element_text(size = 8, color = "black", 
                                          margin = unit(c(0, 0.5, 0, 0), "cm"), family = "Arial"),
               axis.ticks.length = unit(-0.15, "cm"), plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"))
XCate <- theme(axis.title.x = element_blank(),
               axis.text.x = element_text(color = "black", size = 14, family = "Arial",
                                          margin = unit(c(0.5, 0.5, 0, 0.5), "cm")),
               plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"))
theme_f <- theme(strip.text.y = element_text(color = "black", size = 11, family = "Arial", face = "bold"),
                 strip.background = element_rect(fill = "#CCCCCC"),
                 panel.spacing = unit(0.75, "lines"),
                 strip.text.x = element_text(size = 10, face = "bold", family = "Arial"))
#
Sites <- c("WE" = "Weedon Island", "PP" = "Pinellas Point")
Stations <- c("OY" = "Oyster reef", "SG" = "Seagrass", "SS" = "Soft sediment")
Spp <- c("BT" = "Banded T", "CC" = "Crown C", "GM" = "Green M", "LW" = "Lightning W", "MC" = "Other 1", "MG" = "Other 2", "PW" = "Pear W", "TT" = "True T")
Treatments <- c("A" = "Calcein", "T" = "Tag", "C" = "Control")
#
Stages <- c("0" = "Undeveloped", "1" = "Early Developing", "2" = "Late Developing", "3" = "Mature", "4" = "Recovering")
Sex <- c("M" = "Male", "F" = "Female", "U" = "Undetermined")
Months <- c("1" = "Jan", "2" = "Feb", "3" = "Mar", "4" = "Apr", "5" = "May", "6" = "Jun",
            "7" = "Jul", "8" = "Aug", "9" = "Sep", "10" = "Oct", "11" = "Nov", "12" = "Dec")
#
cbPalette <- c("#DD0000", "#E69F00", "#F0E442", "#009E73", "#56B4E9")
#Map color to Stage
names(cbPalette) <- levels(histo$Stage)
StaFill <- scale_fill_manual(name = "Stage", labels = Stages, values = cbPalette, na.value = "#999999")
#
Model_color <- scale_color_manual(name = "Model", labels = c("ELEFAN", "Recaptures", "Laboratory"),
                                  values = c("#000000", "#0072B2", "#E69F00"))
#
#
####Survey summary####
#
All_data <- read.csv("CSV/All_snails.csv", na.string = c("Z", "", "NA", " ")) %>% 
  filter(Species == "BT" & Site != "FtD") %>% drop_na(Station) %>%
  mutate(MonYr = as.factor(paste(Year, Month, sep = "/")),
         Site = as.factor(Site),
         Station = as.factor(Station))
All_summ_data <- All_data %>% group_by(MonYr, Site, Station) %>% summarise(n = n()) %>% 
  ungroup() %>% 
  complete(MonYr, nesting(Site, Station), fill = list(n = 0))
head(All_summ_data)
#Visualize site and station counts
ggboxplot(All_summ_data, x = "Site", y = "n", color = "Station")
#
##2-way Anova - does count change over sites or stations
All_summ_data %>% ggplot(aes(n)) + geom_histogram(aes(y = ..count..)) #skewed and 0 inflated
#
set.seed(54321)
All_anova <- aovp(n ~ Site + Station, data = All_summ_data, perm = "", nperm = 10000)
summary(All_anova)
All_A_tidy <- tidy(All_anova)
names(All_A_tidy) <- c("Factors", "df", "SS", "MS", "Iter", "Pr")
#
(All_TUK <- TukeyHSD(All_anova))
TUK_tab <- tidy(All_TUK) %>% dplyr::select(-c("term", "null.value"))
names(TUK_tab) <- c("Contrast", "Estimate", "CIlower", "CIupper", "adj.p-value")
#Get letters, load summary table, append letters and save with letters
(SiteLett <- merge(All_data %>% group_by(Site, MonYr) %>% summarise(n = n()) %>% ungroup() %>% 
                     group_by(Site) %>% summarise(mean = mean(n, na.rm = T), sd = sd(n, na.rm = T)),
  data.frame(Letters = multcompView::multcompLetters4(All_anova, comp = All_TUK)$Site$Letters) %>%
  tibble::rownames_to_column("Site")))
#
(StationLett <- merge(All_data %>% group_by(Station, MonYr) %>% summarise(n = n()) %>% ungroup() %>% 
                        group_by(Station) %>% summarise(mean = mean(n, na.rm = T), sd = sd(n, na.rm = T)),
                   data.frame(Letters = multcompView::multcompLetters4(All_anova, comp = All_TUK)$Station$Letters) %>%
                     tibble::rownames_to_column("Station")))
#
#
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
ggsave(path = "Output/Figures/2023 03/", filename = "A_LFD_all_2023 03.tiff", dpi = 1000, height = 5, width = 5, unit = "in")
#
set_MA <- c(3)
###Create lfq structure
lfq2 <- lfqCreate(Elefant, Lname = "SL", Dname = "Date")
#Plot length-frequency data
lfq3 <- lfqModify(lfq2, bin_size = 6) # modify bin size
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
res_SA_1 <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                      MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                      init_par = list(Linf = Linf_1_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                      low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, NOT seasonalized
set.seed(1)
res_SA_2 <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                      MA = set_MA, seasonalised = FALSE, addl.sqrt = FALSE,
                      init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                      low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 1, seasonalized
set.seed(1)
res_SA_1_S <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_1_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_SA with range 2, seasonalized
set.seed(1)
res_SA_2_S <- ELEFAN_SA(lfq4, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                        MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                        init_par = list(Linf = Linf_2_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                        low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
#
#### ELEFAN analysis with genetic algorithm 
# run ELEFAN_GA with range 1, not seasonalized
set.seed(1)
res_GA_1 <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                      addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                      low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, not seasonalized
set.seed(1)
res_GA_2 <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = FALSE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                      addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                      low_par = list(Linf = Linf_2[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = Linf_2[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 1, seasonalized
set.seed(1)
res_GA_1_S <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
                        addl.sqrt = FALSE, parallel = TRUE, monitor = FALSE, 
                        low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                        up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
# run ELEFAN_GA with range 2, seasonalized
set.seed(1)
res_GA_2_S <- ELEFAN_GA(lfq4, MA = set_MA, seasonalised = TRUE, maxiter = 100, pmutation = 0.2, popSize = 60, 
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
#Save manually: (path = "Output/Figures", filename = "A_Elefan_final model fit_SA1S.tiff", dpi = 1000)
#
#Create confidence intervals for best solution
JK <- vector("list", length(lfq4$dates))
for(i in 1:length(lfq4$dates)){
  loop_data <- list(dates = lfq4$dates[-i],
                    midLengths = lfq4$midLengths,
                    catch = lfq4$catch[,-i])
  tmp <- ELEFAN_SA(loop_data, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                   MA = set_MA, seasonalised = TRUE, addl.sqrt = FALSE,
                   init_par = list(Linf = Linf_1_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                   low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1))
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
Ms <- M_empirical(Linf = res_SA_1_S$par$Linf, K_l = res_SA_1_S$par$K, method = "Then_growth")
lfq4$M <- as.numeric(Ms)
paste("M =", as.numeric(Ms))
#> [1] "M = 0.455"
res_SA_1_S$agemax
#
#Figure of model output
print(res_SA_1_S$par) #print to copy numbers into figure
print(JK_table) #Use in figure if reporting these numbers
#
data.frame("t" = VBGF(L = seq(0,100,0.1), list(Linf=95.5825525, K=0.3847507, C=0.7161238, ts = 0.6769432, t0 = 0.2579785)),
           "L" = seq(0,100,0.1)) %>%
  ggplot(aes(t, L))+
  geom_line()+
  scale_x_continuous("Age", expand = c(0,0), limits = c(0, 15))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0),limits = c(0, 100))+
  geom_vline(xintercept = res_GA_2_S$agemax) +
  geom_hline(yintercept = mean(Elefant$SL, na.rm = T))+
  geom_hline(yintercept = max(Elefant$SL, na.rm = T), linetype = "dashed")+
  Base
#
ggsave(path = "Output/Figures/2023 03/", filename = "A_JKSA1S_VB_95.583_k0.385_t0.258_2023 03.tiff", dpi = 1000)
#
#
par(mfrow = c(2,1), mai = c(0.5, 0.5, 0.5, 0.5))
plot(lfq4, Fname = "rcounts",date.axis = "modern")
lt <- lfqFitCurves(lfq2, par = res_SA_1_S$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
#
L <- seq(0,100,0.1)
t <- VBGF(L = seq(0,100,0.1), list(Linf=95.5825525, K=0.3847507, C=0.7161238, ts = 0.6769432, t0 = 0.2579785))
plot(t, L, t="l", xlim = c(0, (res_SA_1_S$agemax +  5)), ylim = c(0, 100))
abline(v = 8, col = "red") #Add age max
abline(h = mean(Elefant$SL, na.rm = T), col = "blue")
abline(h = max(Elefant$SL, na.rm = T), col = "blue", lty = 2)
#
#Manually save:(path = "Output/Figures", filename = "A_Elefan_SA1S_fit_model output.tiff", width = 1500)
#
par(opar)
par(mfrow = c(1,1))
#
#
###Bootstrapped values 
BT_boot <- ELEFAN_SA_boot(lfq4, seasonalised = TRUE, SA_time = 60*5, SA_temp = 6e5, maxit = 1000,
                          MA = set_MA,  addl.sqrt = FALSE, parallel = FALSE, nresamp = 500, 
                          low_par = list(Linf = Linf_1[1], K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = Linf_1[2], K = 1, t_anchor = 1, C = 1, ts = 1),
                          init_par = list(Linf = Linf_1_mid, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5))

#
#
#
#Get means and CIs for variables
(BT_table <- t(rbind(cbind(apply(matrix(BT_boot$bootRaw$Linf, 1, byrow = T) , MARGIN = 1, FUN = mean),
                           apply(matrix(BT_boot$bootRaw$K, 1, byrow = T) , MARGIN = 1, FUN = mean),
                           apply(matrix(BT_boot$bootRaw$t_anchor, 1, byrow = T) , MARGIN = 1, FUN = mean),
                           apply(matrix(BT_boot$bootRaw$C, 1, byrow = T) , MARGIN = 1, FUN = mean),
                           apply(matrix(BT_boot$bootRaw$ts, 1, byrow = T) , MARGIN = 1, FUN = mean)), 
                     cbind(apply(matrix(BT_boot$bootRaw$Linf, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
                           apply(matrix(BT_boot$bootRaw$K, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
                           apply(matrix(BT_boot$bootRaw$t_anchor, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
                           apply(matrix(BT_boot$bootRaw$C, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975))),
                           apply(matrix(BT_boot$bootRaw$ts, 1, byrow = T) , MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))))))
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
t <- VBGF(L = seq(0,120,0.1), list(Linf=102.9655824, K=0.2630187, C=0.7422185, ts = 0.6119861, t0 = 0.4223677))
best_ageLength <- data.frame(Age = t,
                             SL = seq(0,120,0.1),
                             Age_n = round(t, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_ranges <- best_ageLength %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "A_EBoot_Age_maxSL.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
##Single plot with data
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1)+
  geom_jitter(data = Elefant %>% inner_join(best_ageLength, by = "SL") %>% mutate(Age_n = round(Age_n, 0)), 
              aes(Age_n, SL), width = 0.25)+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "A_EBoot_Age_maxSL_data.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
rm(lfq2, lfq3, lfq4, set_MA, opar, res_PW, res_SA_1, res_SA_2, res_SA_2_S, res_GA_1, res_GA_2, res_GA_1_S, res_GA_2_S, lt, L, t)
#
#
####LFD Other####
#
#####Repeat for each month
####January
#Estimate clusters
Month1 <- (Elefant %>% subset(Month == 1)) %>% rename(x = SL)
Month1 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit1 <- Mclust(Month1$x)
GMM_fit1$parameters$mean
#
Month1 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit1$parameters$mean)), color = "red")
#
#cluster parameters
(Age_means <- data.frame(GMM_fit1$parameters$mean) %>% 
    rowid_to_column() %>%
    rename("Age" = rowid, "SL" = 'GMM_fit1.parameters.mean') %>% 
    mutate(Month = 1))
#
#
#
####February
#Estimate clusters
Month2 <- (Elefant %>% subset(Month == 2)) %>% rename(x = SL)
Month2 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit2 <- Mclust(Month2$x)
GMM_fit2$parameters$mean
#
Month2 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit2$parameters$mean)), color = "red")
#
#cluster parameters
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit2$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit2.parameters.mean') %>% 
                      mutate(Month = 2)))
#
#
#
####March
#Estimate clusters
Month3 <- (Elefant %>% subset(Month == 3)) %>% rename(x = SL)
Month3 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit3 <- Mclust(Month3$x)
GMM_fit3$parameters$mean
#
Month3 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit3$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit3$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit3.parameters.mean') %>% 
                      mutate(Month = 3)))
#
#
#
####April
#Estimate clusters
Month4 <- (Elefant %>% subset(Month == 4)) %>% rename(x = SL)
Month4 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit4 <- Mclust(Month4$x)
GMM_fit4$parameters$mean
#
Month4 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit4$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit4$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit4.parameters.mean') %>% 
                      mutate(Month = 4)))
#
#
#
####May
#Estimate clusters
Month5 <- (Elefant %>% subset(Month == 5)) %>% rename(x = SL)
Month5 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit5 <- Mclust(Month5$x)
GMM_fit5$parameters$mean
#
Month5 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit5$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit5$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit5.parameters.mean') %>% 
                      mutate(Month = 5)))
#
#
#
#####June
#Estimate clusters
Month6 <- (Elefant %>% subset(Month == 6)) %>% rename(x = SL)
Month6 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit6 <- Mclust(Month6$x)
GMM_fit6$parameters$mean
#
Month6 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit6$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit6$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit6.parameters.mean') %>% 
                      mutate(Month = 6)))
#
#
#
####July
#Estimate clusters
Month7 <- (Elefant %>% subset(Month == 7)) %>% rename(x = SL)
Month7 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit7 <- Mclust(Month7$x)
GMM_fit7$parameters$mean
#
Month7 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit7$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit7$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit7.parameters.mean') %>% 
                      mutate(Month = 7)))
#
#
#
####August
#Estimate clusters
Month8 <- (Elefant %>% subset(Month == 8)) %>% rename(x = SL)
Month8 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit8 <- Mclust(Month8$x)
GMM_fit8$parameters$mean
#
Month8 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit8$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit8$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit8.parameters.mean') %>% 
                      mutate(Month = 8)))
#
#
#
#####September
#Estimate clusters
Month9 <- (Elefant %>% subset(Month == 9)) %>% rename(x = SL)
Month9 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit9 <- Mclust(Month9$x)
GMM_fit9$parameters$mean
#
Month9 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit9$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit9$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit9.parameters.mean') %>% 
                      mutate(Month = 9)))
#
#
#
####October
#Estimate clusters
Month10 <- (Elefant %>% subset(Month == 10)) %>% rename(x = SL)
Month10 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit10 <- Mclust(Month10$x)
GMM_fit10$parameters$mean
#
Month10 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit10$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit10$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit10.parameters.mean') %>% 
                      mutate(Month = 10)))
#
#
#
####November
#Estimate clusters
Month11 <- (Elefant %>% subset(Month == 11)) %>% rename(x = SL)
Month11 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit11 <-Mclust(Month11$x)
GMM_fit11$parameters$mean
#
Month11 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit11$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit11$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit11.parameters.mean') %>% 
                      mutate(Month = 11)))
#
#
#
####December
#Estimate clusters
Month12 <- (Elefant %>% subset(Month == 12)) %>% rename(x = SL)
Month12 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3)
#
set.seed(54321)
GMM_fit12 <- Mclust(Month12$x)
GMM_fit12$parameters$mean
#
Month12 %>% ggplot() + geom_histogram(aes(x = x, y = ..density..), binwidth = 3) +
  geom_vline(xintercept = (unlist(GMM_fit12$parameters$mean)), color = "red")
#
(Age_means <- rbind(Age_means, 
                    data.frame(GMM_fit12$parameters$mean) %>% 
                      rowid_to_column() %>%
                      rename("Age" = rowid, "SL" = 'GMM_fit12.parameters.mean') %>% 
                      mutate(Month = 12)))
#
#
#
ggplot() + 
  geom_point(data = Age_means, aes(Age, SL))
#
#
###Fit VB curve
#Starting parameters
max(Age_means$SL) #Linf = 66.7
vb <- vbFuns("Typical")
set.seed(54321)
vb.starts <- vbStarts(SL ~ Age, data = Age_means)
#
#Fit model
set.seed(54321)
LFD_vb <- nls(SL ~ vb(Age, Linf, K, t0), start = vb.starts, data = Age_means)
#
summary(LFD_vb)
1-((sum(residuals(LFD_vb)^2))/sum((Age_means$SL-mean(Age_means$SL, na.rm = T))^2)) #R squared = 0.482354
AIC(LFD_vb) #AIC = 170.6317
#
(LFD_boot <- nlsBoot(LFD_vb, niter = 1000))
confint(LFD_boot, plot = TRUE)
plot(residuals(LFD_vb)) #Variability assessment
hist(residuals(LFD_vb)) #Normality assessment
#
#
#Dataframe of estimates and CI
(LFD_table <- left_join(tidy(LFD_vb)[,1:2], 
                        data.frame(LFD_boot$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(LFD_table) <- c("Param", "means", "lower", "upper")
(LFD_table <- LFD_table %>% mutate(Method = "LFD"))
#
##Predict end SL
predict2 <- function(x) predict(x, Age_means)
(Age_means <- cbind(Age_means,
                    predict(LFD_vb, Age_means), 
                    confint(Boot(LFD_vb, f = predict2))) %>%
    rename(predSL = 'predict(LFD_vb, Age_means)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((Age_means$SL - Age_means$predSL)^2)) #9.748699
#Plot obs vs pred
Age_means %>% ggplot(aes(SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed mean shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  scale_y_continuous("Predicted mean shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "A_GMM_LFD_meanSL_preds_mclust.tiff", dpi = 1000)
#
#
###Compare to ELEFAN methods
#
##Age-length estimations
#Save age-length estimations using bootstrapped values
data.frame(LFD_table)
t_LFD <- VBGF(L = seq(0,120,0.1), list(Linf=98.9378240, K=0.3233681, t0 = -0.4691356))
best_ageLength_LFD <- data.frame(Age = t_LFD,
                                 SL = seq(0,120,0.1),
                                 Age_n = round(t_LFD, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_LFD_ranges <- best_ageLength_LFD %>% drop_na(Age) %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "dashed")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "A_GMM_Age_maxSL_mclust.tiff", dpi = 1000)
#
##Single with data
ggplot()+
  geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "dashed")+
  geom_jitter(data = Age_means, aes(Age, SL), width = 0.25)+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "A_GMM_Age_maxSL_mclust_data.tiff", dpi = 1000)
#
##Combined
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "solid", color = "#000000")+
  geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL),linewidth = 1, linetype = "dashed", color = "#000000")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "P_A_LFD_both_Age_maxSL.tiff", dpi = 1000)
#
#
rm(Month1, GMM_fit1, Month2, GMM_fit2, Month3, GMM_fit3, Month4, GMM_fit4, Month5, GMM_fit5, Month6, GMM_fit6, Month7, 
   GMM_fit7, Month8, GMM_fit8, Month9, GMM_fit9, Month10, GMM_fit10, Month11, GMM_fit11, Month12, GMM_fit12, vb, vb.starts)
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
ggsave(path = "Output/Figures/2023 03", filename = "B_Recaptures_daysSince.tiff", dpi = 1000)
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
#ggsave(path = "Output/Figures/2023 03", filename = "B_Recaptures_Observed growth.tiff", dpi = 1000)
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
ggsave(path = "Output/Figures/2023 03", filename = "B_Recaptures_Ave growth binned.tiff", dpi = 1000)
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
###vonBert
#
##Train and test dfs
#Divide data into training (70%) and test (30%) sets
set.seed(54321)
sample <- sample(nrow(FirstLast[1:6]), size = nrow(FirstLast[1:6]) * 0.7)
MR_train <- FirstLast[sample, 1:6]
MR_test <- FirstLast[-sample, 1:6]
#
####Fit model to determine parameters estimates:
#
max(FirstLast$End_SL)
###Fit Wang model 
Wang.start <- list(Linf = 76.5, K = 0.001, b = 0)
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
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "B_Recaptures_tt_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges (from Recap_table_tt, t0 from ELEFAN)
t_tt <- VBGF(L = seq(0,120,0.1), list(Linf=109.0244963, K=0.1586193, t0 = 0.4223677))
best_ageLength_tt <- data.frame(Age = t_tt,
                                SL = seq(0,120,0.1),
                                Age_n = round(t_tt, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_MR_tt <- best_ageLength_tt %>% drop_na(Age) %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1)+#, color = "#0072B2")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "B_MR_Age_maxSL_tt_P.tiff", dpi = 1000)
#
#
#
#
#
#
###Gaussian likelihood approach
#
source('GrowthEstimation_CapRecapSim.r')
source('GrowthEstimation_Methods.r')  
##
compile("FabensBayesian.cpp") #Only need to run the first time
dyn.load(dynlib("FabensBayesian"))
#
##Set base info
ITER <- 10000 
WARM <- 5000 
CHAINS <- 4
Exclude<-"N" 
MinDays<-14 ### Min number of days between a mark and a recapture event (minimum deltaT)
### Define prior distribution for the parameters Choices are: uniform, normal or lognormal
LinfPriorDist<-'lognormal' 
KPriorDist<-'uniform' 
SigmaPriorDist<-'uniform' 
#
Lmax <- 8.32 ### Define the maximum length (in cm) 
UpperLmax <- 12.0 ### Define an upper value (in cm) 
#
### Plot the prior distribution and generate the input mean and sd
LinfPrior <- hp.lognormal(Lmax, UpperLmax, plot=T) # for a normal distributed prior on Linf use LinfPrior<-hp.normal(Lmax, UpperLmax, plot=T)

### Specify the priors here, when prior distribution is uniform specify min and max, when otherwise specify mean and sd. 
LinfPr1<-as.numeric(LinfPrior[1]) #Linf prior, if uniform than this is the minimum bound, else it is the mean
LinfPr2<-as.numeric(LinfPrior[2]) #Linf prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
KPr1<-10^-10 #k prior, if uniform than this is the minimum bound, else it is the mean
KPr2<-100 #k prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
SPr1<-10^-10 #sigma prior, if uniform than this is the minimum bound, else it is the mean
SPr2<-10^3 #sigma prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
#
priorlist<-list(c(LinfPr1,LinfPr2),
                c(KPr1,KPr2),
                c(SPr1, SPr2))
#
### Convert train and test to cm:
MR_train_raw <- MR_train %>% dplyr::select(Start_SL, End_SL, diff) %>% 
  mutate(Start_SL = Start_SL/10, #Convert to cm
         End_SL = End_SL/10) %>% #Convert to cm
  rename(L1 = Start_SL, L2 = End_SL, deltaT = diff)
MR_test_raw <- MR_test %>% dplyr::select(Start_SL, End_SL, diff) %>% 
  mutate(Start_SL = Start_SL/10, #Convert to cm
         End_SL = End_SL/10) %>% #Convert to cm
  rename(L1 = Start_SL, L2 = End_SL, deltaT = diff)
rawdata<- MR_train_raw
#
#
###Highlight through line 1113 and run as one script
#
# Delete observations with short times at liberty
ndays <- 365.2425 # nb of days in 1 year, for conversion
rawdata2<-subset(rawdata, deltaT>(MinDays/ndays))
rawdata2<-droplevels(rawdata2)
#
### Calculate growth rate cm per year
rawdata2$Growth_rate_cm_Year<-(rawdata2$L2-rawdata2$L1)/rawdata2$deltaT
#
#Not removing outliers or excluding negative growth due to low sample size and neg growth included in Wang model 
boxplot(rawdata2$Growth_rate_cm_Year,horizontal = TRUE, xlab="Growth per year in cm")
outs<-boxplot.stats(rawdata2$Growth_rate_cm_Year)
RemovedN<-length(outs$out)
outlier<-outs$out
rawdata2$MATCH<-match(rawdata2$Growth_rate_cm_Year, outlier)
rawdata3<-subset(rawdata2, is.na(MATCH)=="TRUE")
rawdata3<-droplevels(rawdata3)
#Exclude negative growth rates
rawdata3$Length_change<-rawdata3$L2-rawdata3$L1
data<-subset(rawdata3, Length_change>=0)
data<-droplevels(data)

### Use original dataset (no data points excluded) or approach as described in the manuscript (Exclude = Y):
invisible(ifelse(Exclude=="Y", data<-data,
                 data<-rawdata2))
### Transform dataframe to list
datalist <- as.list(as.data.frame(data))
### Generate starting values
Linfstart<-mean(as.numeric(Lmax)/0.99)
Kstart <- as.numeric(median(-log((Linfstart-data$L2)/(Linfstart-data$L1))/data$deltaT, na.rm=T))
starts<-c(Linfstart, Kstart)
### Run BFa:
set.seed(54321)
Bfa_MR <- Bfa65(par=starts,L1=datalist$L1,L2=datalist$L2, #Using Fabens since it's recommended from lit and is Gaussian likelihood like LFD data
             deltaT=datalist$deltaT,
             priordist.Linf=LinfPriorDist,
             priordist.K=KPriorDist,
             priordist.sigma=SigmaPriorDist,
             hyperpar=priorlist,
             meth='nlminb',compute.se=T,
             onlyTMB=F,output.post.draws=F,
             mcmc.control=list('nchains'=CHAINS,'iter'=ITER,'warmup'=WARM,
                               'adapt_delta'=.8,'max_treedepth'=20))
# Initial sample size before data point removal (if selected)
InitialN<-length(rawdata$L1)
InitialN 
# Get final sample after data point removal
FinalN_MR<-length(data$L1)
FinalN_MR
# Get priorlist 
# [[1]]== Linf min max or mean and sd for uniform or normal/lognormal priors respectively
# [[2]]== K min max or mean and sd for uniform or normal/lognormal priors respectively
# [[3]]== Sigma min max or mean and sd for uniform or normal/lognormal priors respectively
priorlist
#Get parameter estimates and 95% credible intervals
Bfa_MR$par
Bfa_MR$cred.int
#
##Output table
Bfa_MR_table <- left_join(as.data.frame(Bfa_MR$par) %>% rename("means" = 'Bfa_MR$par') %>% rownames_to_column(),
                          as.data.frame(Bfa_MR$cred.int) %>% dplyr::select(Linf, K) %>% 
                            t() %>% data.frame() %>% rownames_to_column()) 
colnames(Bfa_MR_table) <- c("Param", "means", "lower", "upper")
(Bfa_MR_table <- Bfa_MR_table %>% mutate(Method = "Fabens"))
Bfa_MR_table[1, 1] <- "Linf (cm)"
#
#
#
fvb <- vbFuns("Fabens")
set.seed(54321)
MR_fa_mod <- nls(L2 ~ fvb(L1, deltaT, Linf, K),start = list(Linf = 76.5, K = 0.001), data = MR_train_raw)
#
summary(MR_fa_mod)
##Predict end SL of test set
predict2 <- function(x) predict(x, MR_test_raw)
(MR_test_raw <- cbind(MR_test_raw,
                  predict(MR_fa_mod, MR_test_raw), 
                  confint(Boot(MR_fa_mod, f = predict2))) %>%
    rename(predSL = 'predict(MR_fa_mod, MR_test_raw)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((MR_test_raw$L2 - MR_test_raw$predSL)^2)) #1.634393
#Plot obs vs pred
MR_test_raw %>% ggplot(aes((L2*10), (predSL*10)))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "B_Recaptures_f_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges (from Bfa_MR, t0 from ELEFAN)
MR_f_tt <- VBGF(L = seq(0,12.0,0.1), list(Linf=6.7887369780, K=0.001568173, t0 = 0.4223677))  
best_ageLength_MRf <- data.frame(Age = MR_f_tt,
                                SL = seq(0,12.0,0.1),
                                Age_n = round(MR_f_tt, 2)) %>%
  mutate(Age_n = as.numeric(Age_n/365.2425), #Convert from days back to years
         SL = SL * 10) #Convert back into mm
#
(Age_SL_MR_f <- best_ageLength_MRf %>% drop_na(Age) %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_MR_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "dashed")+#, color = "#0072B2")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "B_MR_Age_maxSL_F_P.tiff", dpi = 1000)
#
##Combined
ggplot()+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "solid", linewidth = 1)+#, color = "#0072B2")+
  geom_line(data = Age_SL_MR_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dashed", linewidth = 1)+#, color = "#0072B2")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "B_MR_both_Age_maxSL.tiff", dpi = 1000)
#
#
rm(ITER, WARM, CHAINS, Exclude, MinDays, LinfPriorDist, KPriorDist, SigmaPriorDist, Lmax, UpperLmax, LinfPrior, 
   LinfPr1, LinfPr2, KPr1, KPr2, SPr1, SPr2, priorlist, rawdata, ndays, rawdata2, outs, RemovedN, outlier, rawdata3, 
   data, datalist, Linfstart, Kstart, starts, InitialN, fvb)
#
#
#
####Wet Lab####
#
(WL_summ <- rbind(WetLab %>% 
                    filter(Date == "2019-02-13") %>% summarise(n = n(), 
                                                               min = min(SL, na.rm = T),
                                                               max = max(SL, na.rm = T),
                                                               mean = mean(SL, na.rm = T),
                                                               sd = sd(SL, na.rm = T)) %>%  
                    mutate(Type = "Start",
                           IDs = paste(WetLab %>% filter(Date == "2019-02-13") %>% distinct(ID))),
                  WetLab %>% 
                    filter(Date == "2021-02-11") %>% summarise(n = n(), 
                                                               min = min(SL, na.rm = T),
                                                               max = max(SL, na.rm = T),
                                                               mean = mean(SL, na.rm = T),
                                                               sd = sd(SL, na.rm = T)) %>%  
                    mutate(Type = "End",
                           IDs = paste(WetLab %>% filter(Date == "2021-02-11") %>% distinct(ID))),
                  WetLab %>% 
                    filter(Date >= "2019-02-13") %>% summarise(n = n(), 
                                                               min = min(SL, na.rm = T),
                                                               max = max(SL, na.rm = T),
                                                               mean = mean(SL, na.rm = T),
                                                               sd = sd(SL, na.rm = T)) %>%  
                    mutate(Type = "All",
                           IDs = paste(WetLab %>% filter(Date >= "2019-02-13") %>% distinct(ID) %>% summarise(n())))) %>% 
   dplyr::select(Type, everything()))
#
#
WL_IDs <- rbind(WetLab %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(1) %>% ungroup(),
                WetLab %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(n()) %>% ungroup()) %>%
  group_by(ID) %>% arrange(ID, Date) %>% filter(Project == "SWG2018") %>%
  mutate_at(.vars = vars(Date, SL), funs(diff = .-lag(.))) #Calculate days in wet lab
#
ungroup(WL_IDs) %>% distinct(ID) %>% tally() #number of snails during course of study
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
WL_treats <- WL_IDs_growth %>% mutate(Treatment = str_sub(ID, -1)) %>%
  mutate(Treatment = ifelse(Treatment == "A", "A", ifelse(Treatment == "T", "T", "C")))
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
                                                meanRate = mean(Rate, na.rm = T),
                                                sdRate = sd(Rate, na.rm = T))
ggboxplot(WL_treats, "Treatment", "Rate") #Compare means
(Rate_KW <- kruskal.test(Rate ~ Treatment, data = WL_treats)) #Kruskal-Wallis
(WL_test_summs <- rbind(data.frame(Test = "Shell Length", 
                                   tidy(SL_KW)),
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
    ungroup() %>% drop_na(End_SL))
#
WL_IDs %>% 
  ggplot(aes(Date, SL, color = ID, group = ID))+
  geom_point()+
  geom_line()+
  scale_x_date("", labels = date_format("%b %Y"), breaks = date_breaks("12 weeks"),
               limits = c(as.Date("2019-02-01"), as.Date("2021-04-12")), expand = c(0,0.1))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 100)) +
  Base + XCate + theme(legend.position = "none") 
#
ggsave(path = "Output/Figures/2023 04", filename = "C_WetLab_Observed growth_all.tiff", dpi = 1000)
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
ggsave(path = "Output/Figures/2023 04", filename = "C_WetLab_Ave growth binned_all.tiff", dpi = 1000)
#
#
##Train and test dfs
#Divide data into training (70%) and test (30%) sets
set.seed(54321)
sample <- sample(nrow(WL_FL[1:6]), size = nrow(WL_FL[1:6]) * 0.7)
WL_train <- WL_FL[sample, 1:6]
WL_test <- WL_FL[-sample, 1:6]
#
###Fit Wang model 
#Starting parameters
max(WL_FL$End_SL) #Linf = 83.1
with(WL_FL, mean((log(End_SL) - log(Start_SL))/diff)) #K = 0.0006030654 but 'minFactor'=0.000976562
Wang.start_WL <- list(Linf = 83.1, K = 0.00098, b = 0)
#
####Fit Wang model to determine parameters estimates:
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
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 04", filename = "C_WetLab_tt_EndSL_preds_tt.tiff", dpi = 1000)
#
#
##Estimated SL ranges
t_WL_tt <- VBGF(L = seq(0,120,0.1), list(Linf=106.22441670, K=0.03838364, t0 = 0.4223677))
WL_ageLength_tt <- data.frame(Age = t_WL_tt,
                              SL = seq(0,120,0.1),
                              Age_n = round(t_WL_tt, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_WL_tt <- WL_ageLength_tt %>% drop_na(Age) %>%  group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, color = "#E69F00")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 03", filename = "C_WL_Age_maxSL_tt_P_all.tiff", dpi = 1000)
#
#
#
#
#
#
#
###Gaussian likelihood approach
#
source('GrowthEstimation_CapRecapSim.r')
source('GrowthEstimation_Methods.r')  
##
compile("FabensBayesian.cpp") #Only need to run the first time
dyn.load(dynlib("FabensBayesian"))
#
##Set base info
ITER <- 10000 
WARM <- 5000 
CHAINS <- 4 
Exclude<-"N" 
MinDays<-14 ### Min number of days between a mark and a recapture event (minimum deltaT)
### Define prior distribution for the parameters Choices are: uniform, normal or lognormal
LinfPriorDist<-'lognormal' 
KPriorDist<-'uniform' 
SigmaPriorDist<-'uniform' 
#
Lmax <- 8.3 ### Define the maximum length (cm)
UpperLmax <- 12.0 ### Define an upper value for your best guess on Linf (cm)
#
### Plot the prior distribution and generate the input mean and sd
LinfPrior <- hp.lognormal(Lmax, UpperLmax, plot=T) # for a normal distributed prior on Linf use LinfPrior<-hp.normal(Lmax, UpperLmax, plot=T)

### Specify the priors here, when prior distribution is uniform specify min and max, when otherwise specify mean and sd. 
LinfPr1<-as.numeric(LinfPrior[1]) #Linf prior, if uniform than this is the minimum bound, else it is the mean
LinfPr2<-as.numeric(LinfPrior[2]) #Linf prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
KPr1<-10^-10 #k prior, if uniform than this is the minimum bound, else it is the mean
KPr2<-100 #k prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
SPr1<-10^-10 #sigma prior, if uniform than this is the minimum bound, else it is the mean
SPr2<-10^3 #sigma prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
#
priorlist<-list(c(LinfPr1,LinfPr2),
                c(KPr1,KPr2),
                c(SPr1, SPr2))
#
### Convert train and test to cm:
WL_train_raw <- WL_train %>% dplyr::select(Start_SL, End_SL, diff) %>% 
  mutate(Start_SL = Start_SL/10, #Convert to cm
         End_SL = End_SL/10) %>% #Convert to cm
  rename(L1 = Start_SL, L2 = End_SL, deltaT = diff)
WL_test_raw <- WL_test %>% dplyr::select(Start_SL, End_SL, diff) %>% 
  mutate(Start_SL = Start_SL/10, #Convert to cm
         End_SL = End_SL/10) %>% #Convert to cm
  rename(L1 = Start_SL, L2 = End_SL, deltaT = diff)
rawdata<- WL_train_raw
#
#
###Highlight through line 1483 and run as one script
#
# Delete observations with short times at liberty
ndays <- 365.2425 # nb of days in 1 year, for conversion
rawdata2<-subset(rawdata, deltaT>(MinDays/ndays))
rawdata2<-droplevels(rawdata2)
#
### Calculate growth rate per year
rawdata2$Growth_rate_cm_Year<-(rawdata2$L2-rawdata2$L1)/rawdata2$deltaT
#
#remove outliers
boxplot(rawdata2$Growth_rate_cm_Year,horizontal = TRUE, xlab="Growth per year in cm")
outs<-boxplot.stats(rawdata2$Growth_rate_cm_Year)
RemovedN<-length(outs$out)
outlier<-outs$out
rawdata2$MATCH<-match(rawdata2$Growth_rate_cm_Year, outlier)
rawdata3<-subset(rawdata2, is.na(MATCH)=="TRUE")
rawdata3<-droplevels(rawdata3)
#Exclude negative growth rates
rawdata3$Length_change<-rawdata3$L2-rawdata3$L1
data<-subset(rawdata3, Length_change>=0)
data<-droplevels(data)

### Use original dataset (no data points excluded) or approach as described in the manuscript (Exclude = Y):
invisible(ifelse(Exclude=="Y", data<-data,
                 data<-rawdata2))
### Transform dataframe to list
datalist <- as.list(as.data.frame(data))
### Generate starting values
Linfstart<-mean(as.numeric(Lmax)/0.99)
Kstart <- as.numeric(median(-log((Linfstart-data$L2)/(Linfstart-data$L1))/data$deltaT, na.rm=T))
starts<-c(Linfstart, Kstart)
### Run BFa:
set.seed(54321)
Bfa_WL <- Bfa65(par=starts,L1=datalist$L1,L2=datalist$L2, #Using Fabens since it's recommended from lit and is Gaussian likelihood like LFD data
                deltaT=datalist$deltaT,
                priordist.Linf=LinfPriorDist,
                priordist.K=KPriorDist,
                priordist.sigma=SigmaPriorDist,
                hyperpar=priorlist,
                meth='nlminb',compute.se=T,
                onlyTMB=F,output.post.draws=F,
                mcmc.control=list('nchains'=CHAINS,'iter'=ITER,'warmup'=WARM,
                                  'adapt_delta'=.8,'max_treedepth'=20))
# Initial sample size before data point removal (if selected)
InitialN<-length(rawdata$L1)
InitialN 
# Get final sample after data point removal
FinalN_MR<-length(data$L1)
FinalN_MR
# Get priorlist 
# [[1]]== Linf min max or mean and sd for uniform or normal/lognormal priors respectively
# [[2]]== K min max or mean and sd for uniform or normal/lognormal priors respectively
# [[3]]== Sigma min max or mean and sd for uniform or normal/lognormal priors respectively
priorlist
#Get parameter estimates and 95% credible intervals
Bfa_WL$par
Bfa_WL$cred.int
#
##Output table
Bfa_WL_table <- left_join(as.data.frame(Bfa_WL$par) %>% rename("means" = 'Bfa_WL$par') %>% rownames_to_column(),
                          as.data.frame(Bfa_WL$cred.int) %>% dplyr::select(Linf, K) %>% 
                            t() %>% data.frame() %>% rownames_to_column())
colnames(Bfa_WL_table) <- c("Param", "means", "lower", "upper")
(Bfa_WL_table <- Bfa_WL_table %>% mutate(Method = "Fabens"))
Bfa_WL_table[1, 1] <- "Linf (cm)"
#
#
#
fvb <- vbFuns("Fabens")
set.seed(54321)
WL_fa_mod <- nls(L2 ~ fvb(L1, deltaT, Linf, K),start = list(Linf = 83.1, K = 0.001), data = WL_train_raw)
#
summary(WL_fa_mod)
##Predict end SL of test set
predict2 <- function(x) predict(x, WL_test_raw)
(WL_test_raw <- cbind(WL_test_raw,
                      predict(WL_fa_mod, WL_test_raw), 
                      confint(Boot(WL_fa_mod, f = predict2))) %>%
    rename(predSL = 'predict(WL_fa_mod, WL_test_raw)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((WL_test_raw$L2 - WL_test_raw$predSL)^2)) #1.634393
#Plot obs vs pred
WL_test_raw %>% ggplot(aes((L2*10), (predSL*10)))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 04", filename = "C_WetLab_f_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges (from Bfa_WL, t0 from ELEFAN)
WL_f_tt <- VBGF(L = seq(0,12.0,0.1), list(Linf=7.3921768741, K=0.0009640975, t0 = 0.4223677))
best_ageLength_WLf <- data.frame(Age = WL_f_tt,
                                 SL = seq(0,12.0,0.1),
                                 Age_n = round(WL_f_tt, 2)) %>%
  mutate(Age_n = as.numeric(Age_n/365.2425), #Convert from days back to years
         SL = SL * 10) #Convert back into mm
#
(Age_SL_WL_f <- best_ageLength_WLf %>% drop_na(Age) %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "dashed", color = "#E69F00")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 04", filename = "C_WL_Age_maxSL_F_P.tiff", dpi = 1000)
#
##Combined
ggplot()+
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "solid", linewidth = 1, color = "#E69F00")+
  geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dashed", linewidth = 1, color = "#E69F00")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 04", filename = "P_C_WL_both_Age_maxSL_all.tiff", dpi = 1000)
#
#
rm(ITER, WARM, CHAINS, Exclude, MinDays, LinfPriorDist, KPriorDist, SigmaPriorDist, Lmax, UpperLmax, LinfPrior, 
   LinfPr1, LinfPr2, KPr1, KPr2, SPr1, SPr2, priorlist, rawdata, ndays, rawdata2, outs, RemovedN, outlier, rawdata3, 
   data, datalist, Linfstart, Kstart, starts, InitialN, fvb)
#
#
#
#
####Wet Lab + Hatchlings####
#
(WL_summ) 
(Hatch_summ <- rbind(Hatchlings %>% group_by(ID) %>% filter(Date == min(Date)) %>%
                       summarise(n = n(), 
                                 min = min(SL, na.rm = T),
                                 max = max(SL, na.rm = T),
                                 mean = mean(SL, na.rm = T),
                                 sd = sd(SL, na.rm = T)),
                     Hatchlings %>% group_by(ID) %>% filter(Date == max(Date)) %>%
                       summarise(n = n(), 
                                 min = min(SL, na.rm = T),
                                 max = max(SL, na.rm = T),
                                 mean = mean(SL, na.rm = T),
                                 sd = sd(SL, na.rm = T)))) %>%
  arrange(ID)
#
#
#Juvenile
(Juv_summ <- rbind(Juvs %>% group_by(ID) %>% filter(Date == min(Date)) %>%
                    summarise(n = n(), 
                              min = min(SL, na.rm = T),
                              max = max(SL, na.rm = T),
                              mean = mean(SL, na.rm = T)),
                  Juvs %>% group_by(ID) %>% filter(Date == max(Date)) %>%
                    summarise(n = n(), 
                              min = min(SL, na.rm = T),
                              max = max(SL, na.rm = T),
                              mean = mean(SL, na.rm = T))) %>%
    arrange(ID))
#
#
#######################################ANALYSES
#
#Hatchling group growth
Hatch_growth <- rbind(Hatchlings %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(1) %>% ungroup(),
                      Hatchlings %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(n()) %>% ungroup()) %>%
  ungroup() %>% arrange(ID, Date) %>% group_by(ID) %>%
  mutate_at(.vars = vars(Date, SL), list(diff = ~.-lag(.))) %>% #Calculate days in wet lab
  mutate(Date_diff = ifelse(is.na(Date_diff), 0, Date_diff), #Add in 0s for starting date
                                   Rate = as.numeric(SL_diff)/as.numeric(Date_diff), #Calculate rate (mm/day) per indiviudal
                                   binSL = lencat(SL, 0, w = 5))
#
#
ungroup(Hatch_growth) %>% filter(Rate < 0) %>% tally() #Number with 0 or negative growth - 6
ungroup(Hatch_growth) %>% filter(Date_diff != 0) %>% summarize(mean = mean(Rate, na.rm = T)) #All average - 0.0156 mm/day
ungroup(Hatch_growth) %>% filter(Date_diff != 0 & Rate > 0) %>% summarize(mean = mean(Rate, na.rm = T)) #>0 (non-negative) - 0.0156 mm/day
#
#
#
#Juvenile group growth
Juv_growth <- rbind(Juvs %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(1) %>% ungroup(),
                    Juvs %>% group_by(ID) %>% drop_na(ID) %>% arrange(ID, Date) %>% slice(n()) %>% ungroup()) %>%
  group_by(ID) %>% arrange(ID, Date) %>% 
  mutate_at(.vars = vars(Date, SL), funs(diff = .-lag(.))) %>%
  mutate(Date_diff = ifelse(is.na(Date_diff), 0, Date_diff), #Add in 0s for starting date
                                Rate = as.numeric(SL_diff)/as.numeric(Date_diff),
                                binSL = lencat(SL, 0, w = 5)) #Calculate rate (mm/day) per individual
#
ungroup(Juv_growth) %>% filter(Rate <= 0) %>% tally() #Number with 0 or negative growth - 9
ungroup(Juv_growth) %>% summarize(mean = mean(Rate, na.rm = T)) #All average - 0.0124 mm/day
ungroup(Juv_growth) %>% filter(Rate > 0) %>% summarize(mean = mean(Rate, na.rm = T)) #>0 (non-negative) - 0.0160 mm/day
#
#
#
##Get dataframe of starting SL, ending SL, and difference between dates
head(WL_FL)
(WL_FL_h <- rbind(WL_FL,
                    full_join(Juv_growth %>% slice(1) %>% dplyr::select(ID, SL, Date) %>% 
                                rename(Start_SL = SL, Start_t = Date),
                              Juv_growth %>% slice(n()) %>% dplyr::select(ID, SL, Date) %>% 
                                rename(End_SL = SL)) %>%
                    arrange(ID, Date) %>% mutate(diff = as.numeric(Date - Start_t)) %>%
                    ungroup(),
                  full_join(Hatchlings %>% group_by(ID) %>% slice(1) %>% dplyr::select(ID, SL, Date) %>% 
                              rename(Start_SL = SL, Start_t = Date),
                            Hatchlings %>% group_by(ID) %>% slice(n()) %>% dplyr::select(ID, SL, Date) %>% 
                              rename(End_SL = SL)) %>%
                    arrange(ID, Date) %>% mutate(diff = as.numeric(Date - Start_t)) %>%
                    ungroup()))
#
##Binned growth
WL_bin_growth_h <- rbind(WL_IDs_growth, Hatch_growth) %>% rbind(Juv_growth) %>%
  mutate(binSL = lencat(SL, 0, w = 5)) %>% group_by(binSL) %>% #Group bins
  mutate(Rate = ifelse(Rate < 0, 0, Rate)) %>% drop_na(binSL) %>% #Limit to actual growth/negative as zero growth
  summarise(meanRate = mean(Rate, na.rm = T),
            sdRate = sd(Rate,na.rm = T))
#
WL_bin_growth_h %>%
  ggplot(aes(binSL, meanRate))+
  geom_point(position = position_nudge(x = 2.5))+
  geom_errorbar(aes(ymin = meanRate - sdRate, ymax = meanRate + sdRate),
                width = 1, position = position_nudge(x = 2.5))+
  scale_y_continuous(name = "Average growth rate (mm/day)", limits = c(-0.02, 0.06), seq(-0.02, 0.06, 0.02), expand = c(0,0))+
  scale_x_continuous(name = "Shell length (mm)", limits = c(0, 100), expand = c(0,0),
                     breaks = seq(0, 100, 5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  Base 
#
ggsave(path = "Output/Figures/2023 08", filename = "C2_WetLab_Ave growth binned_all_hatchlings.tiff", dpi = 1000)
#
#
##Train and test dfs
#Divide data into training (70%) and test (30%) sets
set.seed(54321)
sample_h <- sample(nrow(WL_FL_h[1:6]), size = nrow(WL_FL_h[1:6]) * 0.7)
WL_h_train <- WL_FL_h[sample_h, 1:6]
WL_h_test <- WL_FL_h[-sample_h, 1:6]
#
###Fit Wang model 
#Starting parameters
max(WL_FL_h$End_SL) #Linf = 83.1
with(WL_FL_h %>% filter(ID != "BT8"), mean((log(End_SL) - log(Start_SL))/diff)) #NaN - with BT8 removed - 0.0006386899
Wang.start_WL_h <- list(Linf = 83.1, K = 0.00097, b = 0)
#
####Fit Wang model to determine parameters estimates:
mvb <- vbFuns("Wang")
set.seed(54321)
WL_h_model_tt <- nls(End_SL ~ mvb(Start_SL, diff, Linf, K, b), start = Wang.start_WL_h, data = WL_h_train)
#
summary(WL_h_model_tt)
1-((sum(residuals(WL_h_model_tt)^2))/sum((WL_h_train$End_SL-mean(WL_h_train$End_SL, na.rm = T))^2)) #Rsquared = 0.6586337
AIC(WL_h_model_tt) #AIC = 322.9564
#
(WL_h_boot_tt <- nlsBoot(WL_h_model_tt, niter = 1000))
confint(WL_h_boot_tt, plot = TRUE)
plot(residuals(WL_h_model_tt)) #Variability assessment
hist(residuals(WL_h_model_tt)) #Normality assessment
#
#Dataframe of estimates and CI
(WL_h_table_tt <- left_join(tidy(WL_h_model_tt)[,1:2], 
                          data.frame(WL_h_boot_tt$bootCI[,2:3]) %>% rownames_to_column() %>% rename(term = rowname))) 
colnames(WL_h_table_tt) <- c("Param", "means", "lower", "upper")
(WL_h_table_tt <- WL_h_table_tt %>% mutate(Method = "WetLab_h_tt"))
#
##Predict end SL of test set
predict2_WL_h <- function(x) predict(x, WL_h_test)
(WL_h_test <- cbind(WL_h_test,
                  predict(WL_h_model_tt, WL_h_test), 
                  confint(Boot(WL_h_model_tt, f = predict2_WL_h), type = "perc")) %>%
    rename(predSL = 'predict(WL_h_model_tt, WL_h_test)', lower = '2.5 %', upper = '97.5 %'))
#
#RMSE
sqrt(mean((WL_h_test$End_SL - WL_h_test$predSL)^2)) #7.675782
#Plot obs vs pred
WL_h_test %>% ggplot(aes(End_SL, predSL))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 08", filename = "C2_WetLab_tt_EndSL_preds_h_tt.tiff", dpi = 1000)
#
#
##Estimated SL ranges
t_WL_h_tt <- VBGF(L = seq(0,120,0.1), list(Linf=94.85314, K=0.15396, t0 = 0.4223677))
WL_h_ageLength_tt <- data.frame(Age = t_WL_h_tt,
                              SL = seq(0,120,0.1),
                              Age_n = round(t_WL_h_tt, 2)) %>%
  mutate(Age_n = as.numeric(ifelse(Age_n > 19.99, 20, Age_n))) 
#
(Age_SL_WL_h_tt <- WL_h_ageLength_tt %>% drop_na(Age) %>%  group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_WL_h_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, color = "#999999")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 08", filename = "C2_WL_Age_maxSL_h_tt_P_all.tiff", dpi = 1000)
#
#
#
#
#
#
#
###Gaussian likelihood approach
#
source('GrowthEstimation_CapRecapSim.r')
source('GrowthEstimation_Methods.r')  
##
compile("FabensBayesian.cpp") #Only need to run the first time
dyn.load(dynlib("FabensBayesian"))
#
##Set base info
ITER <- 10000 
WARM <- 5000 
CHAINS <- 4 
Exclude<-"N" 
MinDays<-14 ### Min number of days between a mark and a recapture event (minimum deltaT)
### Define prior distribution for the parameters Choices are: uniform, normal or lognormal
LinfPriorDist<-'lognormal' 
KPriorDist<-'uniform' 
SigmaPriorDist<-'uniform' 
#
Lmax <- 8.3 ### Define the maximum length (cm)
UpperLmax <- 12.0 ### Define an upper value for your best guess on Linf (cm)
#
### Plot the prior distribution and generate the input mean and sd
LinfPrior <- hp.lognormal(Lmax, UpperLmax, plot=T) # for a normal distributed prior on Linf use LinfPrior<-hp.normal(Lmax, UpperLmax, plot=T)

### Specify the priors here, when prior distribution is uniform specify min and max, when otherwise specify mean and sd. 
LinfPr1<-as.numeric(LinfPrior[1]) #Linf prior, if uniform than this is the minimum bound, else it is the mean
LinfPr2<-as.numeric(LinfPrior[2]) #Linf prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
KPr1<-10^-10 #k prior, if uniform than this is the minimum bound, else it is the mean
KPr2<-100 #k prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
SPr1<-10^-10 #sigma prior, if uniform than this is the minimum bound, else it is the mean
SPr2<-10^3 #sigma prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
#
priorlist<-list(c(LinfPr1,LinfPr2),
                c(KPr1,KPr2),
                c(SPr1, SPr2))
#
### Convert train and test to cm:
WL_h_train_raw <- WL_h_train %>% filter(ID != "BT8") %>% dplyr::select(Start_SL, End_SL, diff) %>% 
  mutate(Start_SL = Start_SL/10, #Convert to cm
         End_SL = End_SL/10) %>% #Convert to cm
  rename(L1 = Start_SL, L2 = End_SL, deltaT = diff)
WL_h_test_raw <- WL_h_test %>% dplyr::select(Start_SL, End_SL, diff) %>% 
  mutate(Start_SL = Start_SL/10, #Convert to cm
         End_SL = End_SL/10) %>% #Convert to cm
  rename(L1 = Start_SL, L2 = End_SL, deltaT = diff)
rawdata_h<- WL_h_train_raw
#
#
###Highlight through line 1910 and run as one script
#
# Delete observations with short times at liberty
ndays <- 365.2425 # nb of days in 1 year, for conversion
rawdata2_h<-subset(rawdata_h, deltaT>(MinDays/ndays))
rawdata2_h<-droplevels(rawdata2_h)
#
### Calculate growth rate per year
rawdata2_h$Growth_rate_cm_Year<-(rawdata2_h$L2-rawdata2_h$L1)/rawdata2_h$deltaT
#
#remove outliers
boxplot(rawdata2_h$Growth_rate_cm_Year,horizontal = TRUE, xlab="Growth per year in cm")
outs_h<-boxplot.stats(rawdata2_h$Growth_rate_cm_Year)
RemovedN_h<-length(outs_h$out)
outlier_h<-outs_h$out
rawdata2_h$MATCH<-match(rawdata2_h$Growth_rate_cm_Year, outlier_h)
rawdata3_h<-subset(rawdata2_h, is.na(MATCH)=="TRUE")
rawdata3_h<-droplevels(rawdata3_h)
#Exclude negative growth rates
rawdata3_h$Length_change<-rawdata3_h$L2-rawdata3_h$L1
data_h<-subset(rawdata3_h, Length_change>=0)
data_h<-droplevels(data_h)

### Use original dataset (no data points excluded) or approach as described in the manuscript (Exclude = Y):
invisible(ifelse(Exclude=="Y", data_h<-data_h,
                 data_h<-rawdata2_h))
### Transform dataframe to list
datalist_h <- as.list(as.data.frame(data_h))
### Generate starting values
Linfstart<-mean(as.numeric(Lmax)/0.99)
Kstart_h <- as.numeric(median(-log((Linfstart-data_h$L2)/(Linfstart-data_h$L1))/data_h$deltaT, na.rm=T))
starts_h<-c(Linfstart, Kstart_h)
### Run BFa:
set.seed(54321)
Bfa_WL_h <- Bfa65(par=starts_h,L1=datalist_h$L1,L2=datalist_h$L2, #Using Fabens since it's recommended from lit and is Gaussian likelihood like LFD data
                deltaT=datalist_h$deltaT,
                priordist.Linf=LinfPriorDist,
                priordist.K=KPriorDist,
                priordist.sigma=SigmaPriorDist,
                hyperpar=priorlist,
                meth='nlminb',compute.se=T,
                onlyTMB=F,output.post.draws=F,
                mcmc.control=list('nchains'=CHAINS,'iter'=ITER,'warmup'=WARM,
                                  'adapt_delta'=.8,'max_treedepth'=20))
# Initial sample size before data point removal (if selected)
InitialN_h<-length(rawdata_h$L1)
InitialN_h 
# Get final sample after data point removal
FinalN_MR_h<-length(data_h$L1)
FinalN_MR_h
# Get priorlist 
# [[1]]== Linf min max or mean and sd for uniform or normal/lognormal priors respectively
# [[2]]== K min max or mean and sd for uniform or normal/lognormal priors respectively
# [[3]]== Sigma min max or mean and sd for uniform or normal/lognormal priors respectively
priorlist
#Get parameter estimates and 95% credible intervals
Bfa_WL_h$par
Bfa_WL_h$cred.int
#
##Output table
Bfa_WL_h_table <- left_join(as.data.frame(Bfa_WL_h$par) %>% rename("means" = 'Bfa_WL_h$par') %>% rownames_to_column(),
                          as.data.frame(Bfa_WL_h$cred.int) %>% dplyr::select(Linf, K) %>% 
                            t() %>% data.frame() %>% rownames_to_column())
colnames(Bfa_WL_h_table) <- c("Param", "means", "lower", "upper")
(Bfa_WL_h_table <- Bfa_WL_h_table %>% mutate(Method = "Fabens-hatch"))
Bfa_WL_h_table[1, 1] <- "Linf (cm)"
#
#
#
fvb2 <- vbFuns("Fabens2")
set.seed(54321)
WL_h_fa_mod <- nls(L2 ~ fvb2(L1, deltaT, Linf, K),start = list(Linf = 83.1, K = 0.001), data = WL_h_train_raw, 
                   #control = nls.control(maxiter = 300), 
                  trace = TRUE)
#
summary(WL_h_fa_mod)
##Predict end SL of test set
predict2_h <- function(x) predict(x, WL_h_test_raw)
(WL_h_test_raw <- cbind(WL_h_test_raw,
                      predict(WL_h_fa_mod, WL_h_test_raw)) %>%#, Dropping CI since not used 
                      #confint(Boot(WL_h_fa_mod, f = predict2_h), type = "perc")) %>%
    rename(predSL = 'predict(WL_h_fa_mod, WL_h_test_raw)'))#, lower = '95% LCI', upper = '95% UCI'))
#
#RMSE
sqrt(mean((WL_h_test_raw$L2 - WL_h_test_raw$predSL)^2)) #0.4914922
#Plot obs vs pred
WL_h_test_raw %>% ggplot(aes((L2*10), (predSL*10)))+
  geom_point()+
  geom_abline(intercept = 0, linetype = "dashed")+
  scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
  Base
#
ggsave(path = "Output/Figures/2023 08", filename = "C2_WetLab_h_f_EndSL_preds.tiff", dpi = 1000)
#
#
##Estimated SL ranges (from Bfa_WL_h, t0 from ELEFAN)
WL_h_f_tt <- VBGF(L = seq(0,12.0,0.1), list(Linf=7.8128319356, K=0.0007219363, t0 = 0.4223677))
best_ageLength_WLf_h <- data.frame(Age = WL_h_f_tt,
                                 SL = seq(0,12.0,0.1),
                                 Age_n = round(WL_h_f_tt, 2)) %>%
  mutate(Age_n = as.numeric(Age_n/365.2425), #Convert from days back to years
         SL = SL * 10) #Convert back into mm
#
(Age_SL_WL_h_f <- best_ageLength_WLf_h %>% drop_na(Age) %>% group_by(round(Age_n, 0)) %>% 
    summarise(minSL = min(SL), maxSL= max(SL)))
#
##Single plot
ggplot()+
  geom_line(data = Age_SL_WL_h_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "solid", color = "#999999")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 08", filename = "C2_WL_h_Age_maxSL_F_P.tiff", dpi = 1000)
#
##Combined
ggplot()+
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "solid", linewidth = 1, color = "#E69F00")+
  geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dashed", linewidth = 1, color = "#E69F00")+
  geom_line(data = Age_SL_WL_h_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "solid", linewidth = 1, color = "#999999")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base
#
ggsave(path = "Output/Figures/2023 08", filename = "P_C2_WL_both_Age_maxSL_all.tiff", dpi = 1000)
#
#
rm(ITER, WARM, CHAINS, Exclude, MinDays, LinfPriorDist, KPriorDist, SigmaPriorDist, Lmax, UpperLmax, LinfPrior, 
   LinfPr1, LinfPr2, KPr1, KPr2, SPr1, SPr2, priorlist, rawdata, ndays, rawdata2, outs, RemovedN, outlier, rawdata3, 
   data, datalist, Linfstart, Kstart, starts, InitialN, fvb)
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
                                                         mean = mean(SL),
                                                         sd = sd(SL)) %>% mutate(MF_Final = "Both"),
                     histo %>% drop_na(SL)%>% group_by(MF_Final)%>% summarise(n = n(), 
                                                                              min = min(SL),
                                                                              max = max(SL),
                                                                              mean = mean(SL),
                                                                              sd = sd(SL))))
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
  ylab("Number of snails")+
  lemon::facet_rep_grid(MF_Final~., labeller = labeller(MF_Final = Sex)) +
  Base + scale_fill_grey(start = 0, end = 0.4)+
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_continuous(name = "Shell length (mm)", expand = c(0,0), 
                     limits = c(0, 100), breaks = seq(0, 90, by = 10))+
  scale_y_continuous(expand = c(0,0), limits = c(0,40))
#
ggsave(path = "Output/Figures/2023 03/", filename = "D_BT_MF_SLbins_BW.tiff", dpi = 1000)
#
#Stages per month
histo %>% mutate(Month = as.factor(Month)) %>% 
  ggplot(aes(Month, fill = Stage)) +
  geom_bar(position = "fill")+
  lemon::facet_rep_grid(MF_Final~.)+
  scale_y_continuous("Percent", labels = scales::percent_format(), expand = c(0,0))+
  Base + StaFill #+
  scale_fill_grey(start = 0, end = 0.9)

#
ggsave(path = "Output/Figures/2023 03/", filename = "D_BT_Monthly_Stages_BW.tiff", dpi = 1000)
#
##Stages per sex
histo %>% mutate(Month = as.factor(Month)) %>% 
  ggplot(aes(Month, fill = Stage)) +
  geom_bar(position = "fill")+
  lemon::facet_rep_grid(MF_Final~.)+
  scale_y_continuous("Percent", labels = scales::percent_format(), expand = c(0,0))+
  Base + theme_f + StaFill #+scale_fill_grey(start = 0, end = 0.9)
#
ggsave(path = "Output/Figures/2023 03/", filename = "D_BT_Monthly_Stages_MaleFemale.tiff", dpi = 1000)
#
##Stages per sex per SL bin
histo %>% drop_na(SLclass, Stage) %>%
  ggplot(aes(SLclass, fill = Stage)) +
  geom_bar(position = "fill")+
  lemon::facet_rep_grid(MF_Final~.)+
  scale_y_continuous("Percent", labels = scales::percent_format(), expand = c(0,0))+
  xlab("Shell length (mm)") +
  Base + StaFill + theme_f + theme(axis.text.x = element_text(size = 9)) #+ scale_fill_grey(start = 0, end = 0.9)
#
ggsave(path = "Output/Figures/2023 03/", filename = "D_BT_SLClass_Stages_MaleFemale.tiff", dpi = 1000)
#
#
####Size at maturity figure
#Need histo file, "Species code", Proportion maturity, Type, Extra
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
  Sex <- c("M" = "Male", "F" = "Female")#, "U" = "Undetermined")
  color_og <- c("#009E73", "#E69F00")#, "#CCCCCC")
  #Map color to Sex
  names(color_og) <- c("F","M")#,"U")
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
(Out <- matureSL(histo, "BT", 0.50, 3, 2))
#Use following for axis title, BW
Out[[1]] + theme(axis.title.x = element_text(size = 14, color = "black", family = "sans")) #+ scale_color_manual(name = "", labels = c("M" = "Male", "F" = "Female"), values = c("#000000", "#666666")) 
#
ggsave(path = "Output/Figures/2023 03/", filename = "D_BT_sizeMaturity_M38.24_F42.16_BW.tiff", dpi = 1000)
#
#
#
####Histology regressions####
#
#
head(histology)
#
###SL v Wet Weight
#
WetWeights <- histology %>% dplyr::select(SL,WW) %>% drop_na()
#
#Check data
boxplot(WW ~ SL, data = WetWeights)
gghistogram(WetWeights$WW)
#Unequal variances, ~normal distribution
#
WetWeights <- WetWeights %>% mutate(logWW = log(WW))
#
gghistogram(WetWeights$logWW)
boxplot(logWW ~ SL, data = WetWeights)
#
set.seed(54321)
SL_WW <- lm(logWW ~ SL, data = WetWeights)
summary(SL_WW)
#Residuals
plot(fitted(SL_WW), resid(SL_WW))
#
#Get weights and run weighted model
WW_wt <- 1/lm(abs(SL_WW$residuals) ~ SL_WW$fitted.values)$fitted.values^2
set.seed(54321)
wt_WW <- lm(logWW ~ SL, data = WetWeights, weights = WW_wt)
#
#
(S_WW_sum <- summary(wt_WW))
(S_WW_tab <- tidy(wt_WW) %>% mutate(across(2:4, round, 3)) %>% 
    mutate(p.value = formatC(p.value, format = "e", digits = 3)))
names(S_WW_tab) <- c("term", "Est.", "SE", "t", "p-value")
(S_WW_sum_tab <- glance(wt_WW) %>% dplyr::select(r.squared:df, deviance:df.residual) %>%
    mutate(across(1:5, round, 3)) %>% 
    mutate(p.value = formatC(p.value, format = "e", digits = 3)))
names(S_WW_sum_tab) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
#
confint(wt_WW)
#
S_WW_fill <- WetWeights[complete.cases(WetWeights), ] 
S_WWmeans <- WetWeights %>% mutate(binSL = lencat(SL, 0, w = 5)) %>%
  group_by(binSL) %>% summarise(aveWW = mean(WW, na.rm = T), sdWW = sd(WW, na.rm = T)) 
#
predicted_S_WW <- data.frame(predWW = exp(predict(wt_WW, S_WW_fill)), SL = S_WW_fill$SL)
S_WW_CI <- exp(predict(wt_WW, interval = "predict"))
model_S_WW <- cbind(predicted_S_WW, S_WW_CI) %>% 
  dplyr::select(-fit) %>% dplyr::select(SL, predWW, everything()) %>%
  mutate(binSL = lencat(SL, 0, w = 5)) %>% group_by(binSL) %>% 
  dplyr::summarise(predWW = mean(predWW, na.rm = T), lwr = mean(lwr), upr = mean(upr))
#
head(model_S_WW)
#
(S_WWwCI  <- ggplot()+
    geom_jitter(data = histology, aes(SL, WW, color = "Observed"), width = 0.15, alpha = 0.4)+
    geom_point(data = S_WWmeans, aes(binSL, aveWW, color = "Mean"), shape = 15, size = 4.5)+
    geom_line(data = model_S_WW, aes(binSL, predWW, color = "Predict"), size = 1)+
    geom_line(data = model_S_WW, aes(binSL, lwr, color = "95% CI"), linetype = "dashed")+
    geom_line(data = model_S_WW, aes(binSL, upr, color = "95% CI"), linetype = "dashed")+
    ylab("Wet weight (g)")+ 
    xlab("Shell length (mm)")+
    scale_x_continuous(expand = c(0,0), limits = c(0, 90), breaks = seq(0, 90, by = 30)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 90), breaks = seq(0, 90, by = 30)) +
    scale_color_manual(name = "",
                       breaks = c("Observed", "Mean", "Predict", "95% CI"),
                       values = c("#000000", "#000000", "#DD0000", "#999999"),
                       labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                       guide = guide_legend(override.aes = list(
                         linetype = c("blank", "blank", "solid", "dashed"),
                         shape = c(19, 15, NA, NA))))+
    Base + theme(legend.position = "none"))
#
#Check equation and add to figure
makeRegEqnLabel(wt_WW, "SL", "WW")
exp(coef(wt_WW))
#
(S_WWwCI_eq <- S_WWwCI +
    annotate("text", x = 12, y = 80, label = "log(WW) = 0.06 * SL - 0.37",
             size = 3.5, color = "black", family = "sans")+
    annotate("text", x = 12, y = 72, label = paste("R^2 ==", round(S_WW_sum$r.squared, 2)), parse = TRUE,
             size = 3.5, color = "black", family = "sans"))
#
rm(WetWeights, WW_wt, S_WW_sum, S_WW_fill, predicted_S_WW, S_WW_CI, model_S_WW, S_WWwCI)
#
#
#
###SL v Shell Width
#
ShellWidth <- histology %>% dplyr::select(SL,SW) %>% drop_na()
#
#Check data
boxplot(SW ~ SL, data = ShellWidth)
gghistogram(ShellWidth$SW)
#Unequal variances, ~normal distribution
#
#
set.seed(54321)
SL_SW <- lm(SW ~ SL, data = ShellWidth)
summary(SL_SW)
#Residuals
plot(fitted(SL_SW), resid(SL_SW))
#
#Get weights and run weighted model
SW_wt <- 1/lm(abs(SL_SW$residuals) ~ SL_SW$fitted.values)$fitted.values^2
set.seed(54321)
wt_SW <- lm(SW ~ SL, data = ShellWidth, weights = SW_wt)
#
#
(S_SW_sum <- summary(wt_SW))
(S_SW_tab <- tidy(wt_SW) %>% mutate(across(2:4, round, 3)) %>% 
    mutate(p.value = formatC(p.value, format = "e", digits = 3)))
names(S_SW_tab) <- c("term", "Est.", "SE", "t", "p-value")
(S_SW_sum_tab <- glance(wt_SW) %>% dplyr::select(r.squared:df, deviance:df.residual) %>%
    mutate(across(1:5, round, 3)) %>% 
    mutate(p.value = formatC(p.value, format = "e", digits = 3)))
names(S_SW_sum_tab) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
#
confint(wt_SW)
#
S_SW_fill <- ShellWidth[complete.cases(ShellWidth), ] 
S_SWmeans <- ShellWidth %>% mutate(binSL = lencat(SL, 0, w = 5)) %>%
  group_by(binSL) %>% summarise(aveSW = mean(SW, na.rm = T), sdSW = sd(SW, na.rm = T)) 
#
predicted_S_SW <- data.frame(predSW = predict(wt_SW, S_SW_fill), SL = S_SW_fill$SL)
S_SW_CI <- predict(wt_SW, interval = "predict")
model_S_SW <- cbind(predicted_S_SW, S_SW_CI) %>% 
  dplyr::select(-fit) %>% dplyr::select(SL, predSW, everything()) %>%
  mutate(binSL = lencat(SL, 0, w = 5)) %>% group_by(binSL) %>% 
  dplyr::summarise(predSW = mean(predSW, na.rm = T), lwr = mean(lwr), upr = mean(upr))
#
head(model_S_SW)
#
(S_SWwCI  <- ggplot()+
    geom_jitter(data = histology, aes(SL, SW, color = "Observed"), width = 0.15, alpha = 0.4)+
    geom_point(data = S_SWmeans, aes(binSL, aveSW, color = "Mean"), shape = 15, size = 4.5)+
    geom_line(data = model_S_SW, aes(binSL, predSW, color = "Predict"), size = 1)+
    geom_line(data = model_S_SW, aes(binSL, lwr, color = "95% CI"), linetype = "dashed")+
    geom_line(data = model_S_SW, aes(binSL, upr, color = "95% CI"), linetype = "dashed")+
    ylab("Shell width (mm)")+ 
    xlab("Shell length (mm)")+
    scale_x_continuous(expand = c(0,0), limits = c(0, 90), breaks = seq(0, 90, by = 30)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 90), breaks = seq(0, 90, by = 30)) +
    scale_color_manual(name = "",
                       breaks = c("Observed", "Mean", "Predict", "95% CI"),
                       values = c("#000000", "#000000", "#DD0000", "#999999"),
                       labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                       guide = guide_legend(override.aes = list(
                         linetype = c("blank", "blank", "solid", "dashed"),
                         shape = c(19, 15, NA, NA))))+
    Base + theme(legend.position = "none"))
#
#Check equation and add to figure
makeRegEqnLabel(wt_SW, "SL", "SW")
coef(wt_SW)
#
(S_SWwCI_eq <- S_SWwCI +
    annotate("text", x = 12, y = 80, label = "SW = 0.45 * SL + 2.44",
             size = 3.5, color = "black", family = "sans")+
    annotate("text", x = 12, y = 72, label = paste("R^2 ==", round(S_SW_sum$r.squared, 2)), parse = TRUE,
             size = 3.5, color = "black", family = "sans"))
#
rm(ShellWidth, SW_wt, S_SW_sum, S_SW_fill, predicted_S_SW, S_SW_CI, model_S_SW, S_SWwCI)
#
#
#
#
#
####Table outputs####
#
(Anova_tbl <- All_A_tidy %>% mutate(SS = round(SS, 0),
                                   MS = round(MS, 0),
                                   Iter = round(Iter, 3)) %>%
  flextable() %>% 
  set_formatter("Pr" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
#
(Site_tbl <- SiteLett %>% mutate(mean = round(mean, 2),
                                 sd = round(sd, 2)) %>% flextable() %>% autofit())
(Station_tbl <- StationLett %>% mutate(mean = round(mean, 2),
                                       sd = round(sd, 2)) %>% flextable() %>% autofit())
#
#
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
(LFD_tble <- data.frame(LFD_table) %>% flextable() %>% autofit())
#
(Age_LFD_est_tble <- Age_SL_LFD_ranges %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "LFD"))
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
(RecapModel_tbl <- left_join(tidy(MR_model_tt) %>% mutate(Model = "Recaps-tt", R2 = "0.820", AIC = "231.914", RMSE = "3.40198"),
                             data.frame(term = c("Linf", "K", "b"), LCI = confint(MR_boot_tt)[,1], UCI = confint(MR_boot_tt)[,2])) %>% 
    dplyr::select(Model, everything()) %>%
    flextable() %>% 
    set_formatter("estimate" = function(x){formatC(x, format = "f", digits = 3)},
                  "std.error" = function(x){formatC(x, format = "f", digits = 3)},
                  "statistic" = function(x){formatC(x, format = "f", digits = 3)},
                  "p.value" = function(x){formatC(x, format = "e", digits = 3)}) %>%
    merge_v(j = c(1, 7, 8)) %>% 
    autofit())
#
#
(MR_Age_tt_tble <- Age_SL_MR_tt %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Recaps-tt"))
#
(Bfa_MR_tble <- Bfa_MR_table %>% dplyr::select(Method, everything()) %>% flextable() %>% autofit() %>%
  set_formatter("means" = function(x){formatC(x, format = "f", digits = 4)},
                "lower" = function(x){formatC(x, format = "f", digits = 4)},
                "upper" = function(x){formatC(x, format = "f", digits = 4)}))
#
(MR_Age_f_tble <- Age_SL_MR_f %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "Recaps-F"))
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
(WL_h_bin_growth_tbl <- WL_bin_growth_h %>% flextable() %>%
    set_formatter("binSL" = function(x){formatC(x, format = "f", digits = 0)},
                  "meanRate" = function(x){formatC(x, format = "f", digits = 3)},
                  "sdRate" = function(x){formatC(x, format = "f", digits = 3)}) %>% autofit())
#
##Wet lab SL model
(WLModel_tbl <- left_join(tidy(WL_model_tt) %>% mutate(Model = "Wet Lab-t", R2 = "0.655", AIC = "225.6294", RMSE = "7.503713"),
                          data.frame(term = c("Linf", "K", "b"), LCI = confint(WL_boot_tt)[,1], UCI = confint(WL_boot_tt)[,2])) %>% 
    dplyr::select(Model, everything()) %>%
    flextable() %>% 
    set_formatter("estimate" = function(x){formatC(x, format = "f", digits = 3)},
                  "std.error" = function(x){formatC(x, format = "f", digits = 3)},
                  "statistic" = function(x){formatC(x, format = "f", digits = 3)},
                  "p.value" = function(x){formatC(x, format = "e", digits = 3)}) %>%
    merge_v(j = c(1, 7, 8, 9)) %>%  autofit())
#
#
(WL_Age_tt_tble <- Age_SL_WL_tt %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "WetLab-tt"))
#
#
(Bfa_WL_tble <- Bfa_WL_table %>% dplyr::select(Method, everything()) %>% flextable() %>% autofit() %>%
    set_formatter("means" = function(x){formatC(x, format = "f", digits = 4)},
                  "lower" = function(x){formatC(x, format = "f", digits = 4)},
                  "upper" = function(x){formatC(x, format = "f", digits = 4)}))
#
(WL_Age_f_tble <- Age_SL_WL_f %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "WetLab-F"))
#
#
#Wet lab + hatchlings SL model
(WLModel_h_tbl <- left_join(tidy(WL_h_model_tt) %>% mutate(Model = "Wet Lab-t", R2 = "0.659", AIC = "322.956", RMSE = "7.676"),
                          data.frame(term = c("Linf", "K", "b"), LCI = confint(WL_h_boot_tt)[,1], UCI = confint(WL_h_boot_tt)[,2])) %>% 
    dplyr::select(Model, everything()) %>%
    flextable() %>% 
    set_formatter("estimate" = function(x){formatC(x, format = "f", digits = 3)},
                  "std.error" = function(x){formatC(x, format = "f", digits = 3)},
                  "statistic" = function(x){formatC(x, format = "f", digits = 3)},
                  "p.value" = function(x){formatC(x, format = "e", digits = 3)}) %>%
    merge_v(j = c(1, 7, 8, 9)) %>%  autofit())
#
#
(WL_h_Age_tt_tble <- Age_SL_WL_h_tt %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "WetLab_h-tt"))
#
#
(Bfa_WL_h_tble <- Bfa_WL_h_table %>% dplyr::select(Method, everything()) %>% flextable() %>% autofit() %>%
    set_formatter("means" = function(x){formatC(x, format = "f", digits = 4)},
                  "lower" = function(x){formatC(x, format = "f", digits = 4)},
                  "upper" = function(x){formatC(x, format = "f", digits = 4)}))
#
(WL_h_Age_f_tble <- Age_SL_WL_h_f %>% flextable() %>% autofit() %>% 
    colformat_num(j = 1, digits = 0) %>% set_header_labels('round(Age_n, 0)' = "WetLab_h-F"))
#
#
#
#
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
##Regressions
(Reg_means_tbl <- S_WWmeans %>% left_join(S_SWmeans) %>% mutate(across(2:5, round, 2)) %>%
  flextable() %>% vline(j = c(1, 3)) %>% autofit())
#
(Reg_tbl <- rbind(S_WW_tab, S_SW_tab) %>% flextable() %>% autofit() %>% hline(i = 2))
(Reg_sum_tbl <- rbind(S_WW_sum_tab, S_SW_sum_tab) %>% flextable() %>% autofit() %>% hline(i = 1))
#
#
#
##Final Combines
(Data_summs <- rbind(Ele_summ, BT_summ, 
                     WL_summ %>% dplyr::select(-IDs) %>% mutate(sd = NA), 
                     rbind(Hatchlings %>% filter(MonYr == min(MonYr)) %>%
                             summarise(n = n(), min = min(SL, na.rm = T), max = max(SL, na.rm = T), 
                                       mean = mean(SL, na.rm = T), sd = sd(SL, na.rm = T)) %>% mutate(Type = "H_Start"), 
                           Hatchlings %>% filter(MonYr == max(MonYr)) %>%
                             summarise(n = n(), min = min(SL, na.rm = T), max = max(SL, na.rm = T),
                                       mean = mean(SL, na.rm = T), sd = sd(SL, na.rm = T)) %>% mutate(Type = "H_End")),
                     rbind(Juvs %>% filter(MonYr == min(MonYr)) %>%
                             summarise(n = n(), min = min(SL, na.rm = T), max = max(SL, na.rm = T), 
                                       mean = mean(SL, na.rm = T), sd = sd(SL, na.rm = T)) %>% mutate(Type = "J_Start"), 
                           Juvs %>% filter(MonYr == max(MonYr)) %>%
                             summarise(n = n(), min = min(SL, na.rm = T), max = max(SL, na.rm = T),
                                       mean = mean(SL, na.rm = T), sd = sd(SL, na.rm = T)) %>% mutate(Type = "J_End")),
                     Histo_summ %>% rename("Type" = MF_Final) %>% mutate(sd = NA)) %>% 
    dplyr::select(Type, everything()) %>%
    flextable() %>% hline(i = c(1,4, 7, 9, 11)) %>% autofit() %>%
    set_formatter("min" = function(x){formatC(x, format = "f", digits = 3)},
                  "max" = function(x){formatC(x, format = "f", digits = 3)},
                  "mean" = function(x){formatC(x, format = "f", digits = 3)},
                  "sd" = function(x){formatC(x, format = "f", digits = 3)}))

#
#
#
####Save output to Word####
#
save_as_docx("Data summaries" = Data_summs,
             "Site, Station PermANOVA" = Anova_tbl,
             "Site Comps" = Site_tbl,
             "Station Comps" = Station_tbl,
             "ELEFAN Model summaries" = Res_tbl,
             "Bootstrapped means, 95% CI" = BT_tble,
             "Elefan estimated SLs at age" = Age_est_tble,
             "LFD means, 95% CI" = LFD_tble,
             "LFD estimated SLs at age" = Age_LFD_est_tble,
             "Monthly Counts" = BT_month_tbl,
             "Number of days to recapture" = MR_ID_Days_tbl,
             "Recapture binned growth" = MR_bin_growth_tbl,
             "Recapture Models" = RecapModel_tbl,
             "Recapture T/T Data SL ranges" = MR_Age_tt_tble,
             "Recapture Fabens" = Bfa_MR_tble,
             "Recapture Fabens SL ranges" = MR_Age_f_tble,
             "Number of days in Wet lab" = WL_ID_Days_tbl,
             "Kruskal-Wallis tests" = Treatment_tble,
             "Wet lab binned growth" = WL_bin_growth_tbl,             
             "Wet lab Models" = WLModel_tbl,
             "Wet lab T/T Data SL ranges" = WL_Age_tt_tble,
             "Wet lab Fabens" = Bfa_WL_tble,
             "Wet lab Fabens SL ranges" = WL_Age_f_tble,
             "Wet lab binned growth - with SWG 2016" = WL_h_bin_growth_tbl,             
             "Wet lab Models - with SWG 2016" = WLModel_h_tbl,
             "Wet lab T/T Data SL ranges - with SWG 2016" = WL_h_Age_tt_tble,
             "Wet lab Fabens - with SWG 2016" = Bfa_WL_h_tble,
             "Wet lab Fabens SL ranges - with SWG 2016" = WL_h_Age_f_tble,
             "Histolgy monthly samples by SL" = H_sam_tbl,
             "Histology monthly size ranges" = His_mon_tble,
             "Mean WW and SW per binSL" = Reg_means_tbl,
             "Weighted least squares" = Reg_tbl,
             "WLS models" = Reg_sum_tbl,
             path = "Output/Age Size at Maturity_Tables_2023 08 15.docx", 
             pr_section = sect_properties)
#
#
#
####Figures - <04/2023####
#
##Model estimated maxSL/age - 1 ELEFAN/Wang, 1 Bayesian
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "LFD"), linewidth = 1)+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "Recaptures"), linewidth = 1)+ 
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), 
            aes(Age, maxSL, color = "Laboratory"), linewidth = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+ 
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base + theme(legend.key = element_blank(), legend.title = element_blank()) + 
  scale_color_manual(breaks = c("LFD", "Recaptures", "Laboratory"),
                     labels = c("LFD", "Recaptures", "Laboratory"),
                     values = c("#000000", "#0072B2", "#E69F00"))+
  geom_hline(yintercept = 38.24, linetype = "dotted", color =  "blue", linewidth = 0.75)+
  geom_hline(yintercept = 42.16, linetype = "dotted", color = "red", linewidth = 0.75)
#
ggsave(path = "Output/Figures/2023 03/", filename = "Z_Standard_Age_SLs_2023 03.tiff", dpi = 1000)
#
ggplot()+
  geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "LFD"), linewidth = 1, linetype = "dashed")+
  geom_line(data = Age_SL_MR_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "Recaptures"), linewidth = 1, linetype = "dashed")+ 
  geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), 
            aes(Age, maxSL, color = "Laboratory"), linewidth = 1, linetype = "dashed")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+ 
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base + theme(legend.key = element_blank(), legend.title = element_blank()) + 
  scale_color_manual(breaks = c("LFD", "Recaptures", "Laboratory"),
                     labels = c("LFD", "Recaptures", "Laboratory"),
                     values = c("#000000", "#0072B2", "#E69F00"))+
  geom_hline(yintercept = 38.24, linetype = "dotted", color =  "blue", linewidth = 0.75)+
  geom_hline(yintercept = 42.16, linetype = "dotted", color = "red", linewidth = 0.75)
#
ggsave(path = "Output/Figures/2023 03/", filename = "Z_Bayesian_Age_SLs_2023 03.tiff", dpi = 1000)
#
#
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "ELEFAN"), linewidth = 1)+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "Recaptures"), linewidth = 1)+ 
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), 
            aes(Age, maxSL, color = "Laboratory"), linewidth = 1)+
  geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "ELEFAN"), linewidth = 1, linetype = "dashed")+
  geom_line(data = Age_SL_MR_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "Recaptures"), linewidth = 1, linetype = "dashed")+ 
  geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), 
            aes(Age, maxSL, color = "Laboratory"), linewidth = 1, linetype = "dashed")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+ 
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base + theme(legend.key = element_blank(), legend.title = element_blank()) + 
  scale_color_manual(breaks = c("ELEFAN", "Recaptures", "Laboratory"),
                     labels = c("ELEFAN", "Recaptures", "Laboratory"),
                     values = c("#000000", "#0072B2", "#E69F00"))
#
ggsave(path = "Output/Figures/2023 03/", filename = "Z_All_Age_SLs_2023 03.tiff", dpi = 1000)
##
#
#
#
#
#Figure of raw output. Figure with each model added. Figure of all models. Figure all models w/ updated WL.
#
##Individual figures - progression - with or without size at maturity
ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "ELEFAN"), linewidth = 1)+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL, color = "Recaptures"), linewidth = 1)+ 
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age <= 19), 
            aes(Age, maxSL, color = "Laboratory"), linewidth = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+ 
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base + theme(legend.key = element_blank(), legend.title = element_blank()) + 
  scale_color_manual(breaks = c("ELEFAN", "Recaptures", "Laboratory"),
                     labels = c("ELEFAN", "Recaptures", "Laboratory"),
                     values = c("#000000", "#0072B2", "#E69F00"))+
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

####Final Figures output####
#
#
#
###Figure 1-Station Map####
#
Station_list <- read.csv("CSV/TB_Pop_Stations.csv", na.string = c("Z", "", "NA", " "))
US <- readRDS("CSV/gadm36_USA_1_sp.rds")
FL <- US %>% subset(NAME_1 == "Florida")
FL_map <- fortify(FL)
#
FL_inset <- ggplot() +
  geom_blank(data = FL_map, aes(long, lat))+ 
  geom_map(data = FL_map, map = FL_map, aes(group = group, map_id = id), fill = "#999999", color = "black")+
  #geom_polygon(data = FL_map, aes(long, lat, group = group), fill = "#999999", color = "black")+
  geom_text(aes(-85.7, 28.9, label = "Gulf of Mexico", fontface = "italic", family = "Arial"), color = "black", size = 3.5)+
  geom_text(aes(-80.4, 30, label = "Atlantic \n Ocean", fontface = "italic", family = "Arial"), color = "black", size = 3.5)+
  coord_fixed()+ theme_void() + theme(panel.border = element_rect(color = "black", fill = NA),
                                   panel.background = element_rect(color = "white"))+
  geom_point(aes(-82.6, 27.75), size = 6.5, color = "black", shape = 15)
#
#
(TB_Gastro <- ggplot()+
    geom_blank(data = FL_map, aes(long, lat))+ 
    geom_map(data = FL_map, map = FL_map, aes(group = group, map_id = id), fill = "#CCCCCC", color = "black")+
    scale_x_continuous(limits = c(-82.78, -82.53), expand = c(0,0))+ #(limits = c(-82.78, -82.45)
    scale_y_continuous(limits = c(27.60, 27.884), expand = c(0,0))+ #limits = c(27.48, 27.884)
    coord_fixed()+
    geom_text(aes(x = -82.57, y = 27.7, label = "Tampa \n Bay", fontface = "italic", family = "Arial"), color = "black", size = 10)+
    geom_text(aes(x = -82.64, y = 27.684, label = "Pinellas \n Point  ", fontface = "bold", family = "Arial"), color = "black", size = 7)+
    geom_text(aes(x = -82.56, y = 27.78, label = "Weedon \n Island", fontface = "bold", family = "Arial"), color = "black", size = 7)+
    geom_segment(aes(x = -82.575, y = 27.795, xend = -82.597, yend = 27.84), color = "black", size = 1.25)+
    geom_point(data = Station_list, aes(Longitude, Latitude), color = "black", size = 4)+
    theme(panel.border = element_rect(color = "black", fill = NA))+
    theme_void() + 
    theme(panel.border = element_rect(color = "black", fill = NA)))
#
#
#Saving at 1400
TB_Gastro +  
  north(symbol = 12, x.min = -82.60, x.max = -82.57, y.min = 27.605, y.max = 27.62, scale = 1)+ 
  ggsn::scalebar(dist = 2, dist_unit = "km", transform = TRUE, model = "WGS84",
                 x.min = -82.64, x.max = -82.60, y.min = 27.605, y.max = 27.655,
                 height = 0.1, st.dist = 0.10, st.bottom = FALSE, st.size = 3.5)+
  patchwork::inset_element(FL_inset, 0, 0.56, 0.5, 1.08, align_to = "panel") #left, bottom, right, top
#
#
#
###Figure 2-Site and Habitat Comps####
#
SiteLett
StationLett
#
###Sites
#Use letters to differentiate groups
(Site_comps <- ggplot(SiteLett, aes(Site, mean, fill = "black"))+
  geom_bar(stat = "identity", size = 2)+
  geom_errorbar(aes(Site, ymin = mean, ymax = mean + sd), width = 0.2)+
  scale_fill_grey()+
  geom_text(aes(label = Letters, y = mean+sd), size = 5, vjust = -0.35)+
  Base + theme(legend.position = "none") +
  theme(axis.title.y = element_text(size = 14))+
  ylab("Average number of snails")+
  scale_x_discrete(labels = Sites)+ 
  scale_y_continuous(expand = c(0,0), limits = c(0,80))) 
#
###Stations
(Station_comps <- ggplot(StationLett, aes(Station, mean, fill = "black"))+
  geom_bar(stat = "identity", size = 2)+
  geom_errorbar(aes(Station, ymin = mean, ymax = mean+sd), width = 0.2)+
  scale_fill_grey()+
  geom_text(aes(label = Letters, y = mean+sd), size = 5, vjust = -0.35)+
  Base + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_text(size = 11)) +
  xlab("Habitat")+
  scale_x_discrete(labels = Stations)+ 
  scale_y_continuous(expand = c(0,0), limits = c(0,80)))
#
#
#
(Comps <- ggpubr::ggarrange(Site_comps + rremove("xlab"), Station_comps + rremove("xlab"), 
                            labels = c("A", "B"), 
                            nrow = 1, ncol = 2, hjust = -5.9, vjust = 1.55, align = "v"))

#
#
ggsave(path = "Output/Figures", filename = "Fig2_Site_Hab_Comps.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
#
###Figure 3-LFD histograms####
#

(Month_labels <- Elefant %>% mutate(SL = 8,
                                   count = 40) %>% distinct(Month, SL, count) %>%
  right_join(as.data.frame(Months) %>% rownames_to_column("Month") %>% mutate(Month = as.integer(Month))))
#
Elefant %>% 
  ggplot(aes(SL))+
  geom_histogram(binwidth = 4, boundary = 0, closed = "left")+
  geom_text(data = Month_labels, aes(SL, count, label = Months), size = 3.5, family = "sans")+
  lemon::facet_rep_wrap(vars(Month), ncol = 2)+
  theme(strip.text = element_blank())+
  Base2 +
  theme(panel.spacing = unit(0, "lines"), 
        axis.text.y = element_text(size = 6, margin = margin(r = 5)), 
        axis.text.x = element_text(size = 8, margin = margin(t = 5)))+
  scale_y_continuous("Count", expand = c(0,0), limits = c(0, 48), breaks = seq(0, 48, 16))+
  scale_x_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 100), breaks = seq(0, 100, 20))
#
ggsave(path = "Output/Figures/Final/", filename = "C_LFD_monthly.tiff", dpi = 1000, height = 5, width = 5, unit = "in")
#
#
#
###Figure 4-Recaps####
#
(Recaps_by_days <- IDs %>% group_by(ID) %>% #Group by ID
   filter(Recaptured == "Y") %>% mutate(N = as.factor(row_number()))  %>% #Filter to Recaps and count number of recaps/ID
   ggplot(aes(x = diff, fill = N))+ #Color by recap N
   geom_histogram(position = "stack", binwidth = 10, boundary = 0, closed = "right") +
   Base + #theme(legend.position = "none")+
   scale_fill_manual(values = c("#99CCFF", "#003399", "#000000"), labels = c(1, 2, 3))+
   #scale_fill_grey(start = 0.1, end = 0.8)+
   scale_x_continuous("Number of days", expand = c(0,0), limits = c(0, 442), breaks = seq(0, 440, 40))+
   scale_y_continuous("Number of recaptures", expand = c(0,0), limits = c(0, 10), breaks = seq(0, 10, 2)))
#
#
ggsave(path = "Output/Figures", filename = "4_Recaptures.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
#
###Figure6-
###Figure 5-Tag Treatments####
#
#
(Treat_SL <- ggboxplot(WL_treats, "Treatment", "SL", fill = "darkgray")+
   scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0,100))+
   scale_x_discrete("",labels = Treatments)+
   Base)
#
(Treat_Rate <- ggboxplot(WL_treats, "Treatment", "Rate", fill = "darkgray")+
    scale_y_continuous("Growth rate (mm/day)", expand = c(0,0), limits = c(0,0.06))+
    scale_x_discrete("", labels = Treatments)+
    Base)
#
#
(Treat_plots <- ggpubr::ggarrange(Treat_SL, Treat_Rate, labels = c("A", "B"), 
                                  nrow = 1, ncol = 2, hjust = -6.9, vjust = 1.55, align = "v"))
#
ggsave(path = "Output/Figures/Final", filename = "G_Tagging_Treatments.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
#
###Figure 6-Regressions####
#
S_WWwCI_eq
S_SWwCI_eq
#
#
#
(Reg_plots <- ggpubr::ggarrange(S_WWwCI_eq, S_SWwCI_eq, labels = c("A", "B"), 
                                nrow = 2, ncol = 1, hjust = -5.9, vjust = 1.25, align = "v", 
                                common.legend = FALSE))
#
ggsave(path = "Output/Figures", filename = "F6_Regressions_SL_vWW_SW.tiff", dpi = 1000, width = 8, height = 6, units = "in")

#
#
###Figure 7-Histology####
#
##Stages per sex
(Stages_MF <- histo %>% drop_na(SLclass, Stage) %>%
   ggplot(aes(SLclass, fill = Stage)) +
   geom_bar(position = "fill")+
   lemon::facet_rep_grid(MF_Final~.)+
   scale_y_continuous("Percent", labels = label_percent(suffix = ""), expand = c(0,0))+
   xlab("Shell length (mm)") + 
   Base + theme_f + theme(axis.text.x = element_text(size = 9.5),
                          axis.title.y = element_blank(),
                          legend.title = element_blank(),
                          panel.spacing.y = unit(0.05, "lines"), plot.margin = margin(t = 10, b = 5)) + 
   StaFill)#scale_fill_grey(start = 0.9, end = 0, labels = Stages))
#
#Stages per month
(Stages_Month <- histo %>% mutate(Month = as.factor(Month)) %>% 
    ggplot(aes(Month, fill = Stage)) +
    geom_bar(position = "fill")+
    lemon::facet_rep_grid(MF_Final~.)+
    scale_y_continuous("Percent", labels = label_percent(suffix = ""), expand = c(0,0))+
    scale_x_discrete(labels = Months)+
    Base + theme_f + theme(axis.text.x = element_text(size = 11),
                           axis.title.y = element_blank(),
                           legend.title = element_blank(),
                           panel.spacing.y = unit(0.05, "lines")) + 
    StaFill) #scale_fill_grey(start = 0.9, end = 0, labels = Stages))
#
(Stage_plots <- (ggpubr::ggarrange(Stages_MF, Stages_Month, labels = c("A", "B"), 
                                  nrow = 2, ncol = 1, hjust = -4.99, vjust = 0.8, align = "v", 
                                  common.legend = TRUE, legend = "bottom") + theme(legend.text = element_text(family = "Arial"), 
                                                                                   plot.margin = margin(t=2), legend.margin = margin(t = 25, unit = "lines"))) %>%
    annotate_figure(left = text_grob("Percent (%)", rot = 90, family = "Arial", size = 14, color = "black")))
#
ggsave(path = "Output/Figures", filename = "Fig7_Histo_Stages_SL_Month.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
#
#
###Figure 8-Size at maturity####
#
#
(Out <- matureSL(histo, "BT", 0.50, 3, 2))
#Use following for axis title, BW
Out[[1]] + Base + theme_f + 
  theme(axis.title.x = element_text(size = 14, color = "black", family = "sans"),
        legend.position = "none") + 
  scale_color_manual(name = "", labels = c("M" = "Male", "F" = "Female"), values = c("#000000", "#666666")) 
#
ggsave(path = "Output/Figures/Final", filename = "J_Size_Maturity.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#

###Figure 9-Models####
#
##Pop LFD
#ELEFAN
(Pop_LF <- ggplot()+
  geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1)+
  geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linewidth = 1, linetype = "dashed")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("", expand = c(0,0), limits = c(0, 120))+
  Base + theme(axis.title.x = element_blank()))
#
##MR
(MR_plot <- ggplot()+
  geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dotted", linewidth = 1)+
  geom_line(data = Age_SL_MR_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dotdash", linewidth = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
  Base + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16)))
#
##WL
(WL_plot <- ggplot()+
  geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dotted", linewidth = 1)+
  geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
            aes(Age, maxSL), linetype = "dotdash", linewidth = 1) +
  geom_line(data = Age_SL_WL_h_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL), linetype = "solid", linewidth = 1, color = "#999999")+
  geom_line(data = Age_SL_WL_h_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL), linetype = "dotdash", linewidth = 1, color = "#999999") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
  scale_y_continuous("", expand = c(0,0), limits = c(0, 120))+
  Base)
#
#
#
(Growth_models <- ggpubr::ggarrange(Pop_LF, MR_plot, WL_plot, labels = c("A", "B", "C"), 
                                    nrow = 3, ncol = 1, hjust = -6.3, vjust = 1.55, align = "v"))
#
#
ggsave(path = "Output/Figures/Final", filename = "D2_Growth_Models_SWG2016.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
# 
###Figure 9-Models (legend)####
#
##Pop LFD
#ELEFAN
(Pop_LF <- ggplot()+
   geom_line(data = Age_SL_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
             aes(Age, maxSL, linetype = "ELEFAN"), linewidth = 1)+
   geom_line(data = Age_SL_LFD_ranges %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
             aes(Age, maxSL, linetype = "GMM"), linewidth = 1)+
   scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
   scale_y_continuous("", expand = c(0,0), limits = c(0, 120))+
   Base + theme(axis.title.x = element_blank())+
   scale_linetype_manual(name = "",
                         breaks = c("ELEFAN", "GMM"),
                         values = c("solid", "dashed"),
                         labels = c("ELEFAN (W)", "GMM (W)"))+
   theme(legend.position = c(0.095, 0.948), legend.title = element_blank(),
         legend.text = element_text(family = "Arial"), 
         legend.background = element_rect(linetype = 1, color = NA)))
#
##MR
(MR_plot <- ggplot()+
    geom_line(data = Age_SL_MR_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL, linetype = "BFa"), linewidth = 1)+
    geom_line(data = Age_SL_MR_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL, linetype = "Wang"), linewidth = 1)+
    scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
    scale_y_continuous("Shell length (mm)", expand = c(0,0), limits = c(0, 120))+
    Base + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16))+
    scale_linetype_manual(name = "",
                          breaks = c("BFa", "Wang"),
                          values = c("dotdash", "dotted"),
                           labels = c("BFa (W)", "Wang (W)"))+
    theme(legend.position = c(0.090, 0.948), legend.title = element_blank(),
          legend.text = element_text(family = "sans"), 
          legend.background = element_rect(linetype = 1, color = NA)))
#
##WL
(WL_plot <- ggplot()+
    geom_line(data = Age_SL_WL_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL, color = "BFa"), linewidth = 1, linetype = "dotdash") +
    geom_line(data = Age_SL_WL_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL, color = "Wang"), linewidth = 1, linetype = "dotted")+
    geom_line(data = Age_SL_WL_h_f %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL, color = "BFa+"), linewidth = 1, linetype = "dotdash") +
    geom_line(data = Age_SL_WL_h_tt %>% rename(Age = 'round(Age_n, 0)') %>% filter(Age < 20), 
              aes(Age, maxSL, color = "Wang+"), linewidth = 1, linetype = "solid")+
    scale_x_continuous(expand = c(0,0), limits = c(0, 20))+
    scale_y_continuous("", expand = c(0,0), limits = c(0, 120))+
    Base+
    scale_color_manual(name = "",
                       breaks = c("BFa", "Wang", "BFa+", "Wang+"),
                       values = c("#0072B2", "#0072B2", "#E69F00", "#E69F00"),
                       labels = c("BFa (L)", "Wang (L)", "BFa+ (L+)", "Wang+ (L+)"),
                       guide = guide_legend(ncol = 2,
                         override.aes = list(
                         linetype = c("dotdash", "dotted", "dotdash", "solid"))))+
    theme(legend.position = c(0.145, 0.888), legend.title = element_blank(),
          legend.text = element_text(family = "sans"), 
          legend.background = element_rect(linetype = 1, color = NA)))
#
#
#
(Growth_models <- ggpubr::ggarrange(Pop_LF, MR_plot, WL_plot, labels = c("A", "B", "C"), 
                                    nrow = 3, ncol = 1, hjust = -6.6, vjust = 1.55, align = "v"))
#
#
ggsave(path = "Output/Figures/Final", filename = "D2_Growth_Models_SWG2016_legend.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#
# 
###Figure 10-Predicted RMSEs####
#
##MR
(BFa_MR <- MR_test %>% ggplot(aes(End_SL, predSL))+
   geom_point()+
   geom_abline(intercept = 0, linetype = "dashed")+
   annotate("text", x = 10, y = 94, label = "RMSE = 1.63", family = "sans", size = 3)+
   scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
   scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
   Base + theme(axis.title = element_blank(), axis.text.y = element_text(size = 14)))
#
(W_MR <- MR_test_raw %>% ggplot(aes((L2*10), (predSL*10)))+
    geom_point()+
    geom_abline(intercept = 0, linetype = "dashed")+
    annotate("text", x = 10, y = 94, label = "RMSE = 3.40", family = "sans", size = 3)+
    scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
    scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
    Base+ theme(axis.title = element_blank(), axis.text.y = element_text(size = 14)))
#
##WL
(BFa_WL <- WL_test %>% ggplot(aes(End_SL, predSL))+
    geom_point()+
    geom_point(data = WL_h_test, aes(End_SL, predSL), shape = 17, color = "#CC0000", size = 2)+
    geom_abline(intercept = 0, linetype = "dashed")+
    annotate("text", x = 10, y = 79, label = "RMSE = 4.28 \n RMSE2 = 0.49", family = "sans", size = 3)+
    scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
    scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
    Base+ theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 14), axis.title.y = element_text(hjust = -0.2)))
#
(W_WL <- WL_test_raw %>% ggplot(aes((L2*10), (predSL*10)))+
    geom_point()+
    geom_point(data = WL_h_test_raw, aes((L2*10), (predSL*10)), shape = 17, color = "#CC0000", size = 2)+
    geom_abline(intercept = 0, linetype = "dashed")+
    annotate("text", x = 10, y = 77, label = "RMSE = 7.50 \n RMSE2 = 7.68", family = "sans", size = 3)+
    scale_x_continuous("Observed shell length (mm)", expand = c(0,0), limits = c(0, 100))+
    scale_y_continuous("Predicted shell length (mm)", expand = c(0,0), limits = c(0, 100))+
    Base + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14)))
#
#
(RMSE_plots <- ggpubr::ggarrange(BFa_MR, W_MR, BFa_WL, W_WL, labels = c("A", "B", "C", "D"), 
                                    nrow = 4, ncol = 1, hjust = -6.5, vjust = 1.55, align = "v"))
#
#
ggsave(path = "Output/Figures", filename = "F10_Model_RMSE_SWG2016.tiff", dpi = 1000, width = 8, height = 6, units = "in")
#
#
#