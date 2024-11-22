##Operculum striae as an aging technique in Banded Tulips
#
##Data analyses
#
#Add set_formatter(p.adjust = function(x){formatC(x, format = "e", digits = 3)} to p values
#
#
#Use Alt+O to collapse all sections, Alt+Shift+O to expand all sections#
#Use Alt+L to collapse on section, Alt+Shift_L to expand section#
#
graphics.off()  # turns off any plots from previous work session
rm(list=ls(all=TRUE)) # clears out environment 
#
#
#Load require packages (install as necessary)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr, tidyverse, rstatix,
               vcd, #Cohen's Kappa - reader agreement
               FSA, AICcmodavg, propagate, lmPerm, #Length bins, AIC comps, nls predictions
               scales, broom, psych, #Re-scaling MG, cleaning summary tables, pairwise letters
               ggpubr, lemon, #Figure arrangement
               flextable, janitor, officer) #Tables, output to Word
if (!require("remotes")) install.packages("remotes")
remotes::install_github("GegznaV/biostat")
#
#
#
####Load files####
#
##Reader counts
Counts <- read.csv("Output/BT_final_counts_df.csv", na.string = c("Z", "", "NA"))
glimpse(Counts)
#
##Operculum: distance between bands, cumulative growth of bands
Growth <- read.csv("Output/BT_operculum_ave_growth_final_df.csv", na.string = c("Z", "", "NA"))
glimpse(Growth)
#
##Morphology of snails, MF, maturity, binned morphology, coded MF*maturity
Morphology <- read.csv("Output/BT_operculum_morphology_final_df.csv", na.string = c("Z", "", "NA"))
glimpse(Morphology)
Morphology <- Morphology %>% mutate_at(c("Quarter", "MF_Final", "Mature", "MF_Mat", 
                                         "Striae_fac", "Striae_VB"), as.factor)
#
##Final striae counts with MF, Maturity, SL, marginal growth (MG), and MG rate (MGR)
Margins <- read.csv("Output/BT_operculum_marginal_final_df.csv", na.string = c("Z", "", "NA"))
glimpse(Margins)
Margins <- Margins %>% mutate_at(c("Quarter", "MF_Final", "Mature", "MF_Mat"), as.factor)
#
#
#
####Tools####
#
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
#
#
####Figure Formatting####
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
##Titles and text both axis
basetheme <- theme_bw()+ 
  theme(panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(), #Plot background
        axis.line = element_line(color = "black"), 
        axis.title = element_text(size = 12, face = "bold", color = "black", family = "serif"), 
        axis.text = element_text(size = 11, color = "black", family = "serif", margin = unit(c(0.5, 0.5, 0, 0.5), "cm")), 
        axis.ticks.length = unit(-0.15, "cm"))
#
##No X Title (Categories (i.e. Month))
XCate <- theme(axis.title.x = element_blank(),
               axis.text.x = element_text(color = "black", size = 12, family = "serif",
                                          margin = unit(c(0.5, 0.5, 0, 0.5), "cm")))
#
##Legend
legends <- theme(legend.position = "right", 
                 legend.key = element_rect(color = "transparent", fill = "transparent"))
#
##Faceted format
theme_f <- theme(strip.text = element_blank(),
                 strip.background = element_blank(),
                 panel.spacing = unit(0.75, "lines"))
#
#
#
####Morphology summary and Reader agreement####
#
###Samples collected - all, Male, Female
Samples <- full_join(Morphology %>% dplyr::select(SL:Opercula_L) %>% 
                       get_summary_stats(type = "full", show = c("n", "min", "max", "mean", "sd", "se")) %>%
                       mutate(MF_Final = "All"),
                     Morphology %>% dplyr::select(SL:Opercula_L, MF_Final) %>% group_by(MF_Final) %>% 
                       get_summary_stats(type = "full", show = c("n", "min", "max", "mean", "sd", "se"))) %>%
  dplyr::select(MF_Final, everything())
#
(Sample_tab <- Samples %>% 
  flextable() %>% merge_v(j = ~MF_Final) %>% #Merge cells
  colformat_num(j = 3, digits = 0) %>% #Remove decimal places on n
  hline(i = c(4, 8)) %>% autofit()) #Add dividing lines

#
###Samples by Striae count
Sam_striae <- Morphology %>% dplyr::select(SL:Opercula_L, Striae) %>% 
  group_by(Striae) %>%
  get_summary_stats(type = "full", show = c("n", "min", "max", "mean", "sd", "se"))
#
(Sam_striae_tab <- Sam_striae %>% drop_na(Striae) %>%
    flextable() %>% merge_v(j = ~Striae) %>% #Merge cells
    colformat_num(j = 3, digits = 0) %>% #Remove decimal places on n
    hline(i = c(4, 8, 12, 16, 20)) %>% autofit()) #Add dividing lines
#
##Overall SL distribution
Morphology %>% dplyr::select(SL) %>% ggplot(aes(SL))+
  geom_histogram(aes(y=..count..)) 
#
#
#
###Reader agreement 
Counts <- Counts %>% filter(Damaged == "N") %>% #Determine agreement - initial and final
  mutate(Agree = ifelse(Count1 == Count2, "Y",
                        ifelse(is.na(Final), "N", "A")))  %>%
  mutate(Agree = ifelse(is.na(Agree), "N", Agree))
#
#Percent initial agreement (Y), agree after review (A), and disagreement
(Agree <- Counts %>%  group_by(Agree) %>% 
  summarise(n = n()) %>% mutate(Percent = n/sum(n)*100) %>%
  flextable() %>% autofit())
#
#Contingency table
(Tab <- xtabs(data = (Counts %>% filter(Damaged == "N") %>% dplyr::select(Count1, Count2))))
set.seed(54321)
print(Kappa(Tab), CI = TRUE)
#Copy over numbers from printed table for output
(Cohens_Kappa <- data.frame("Type" = c("Unweighted", "Weighted"),
                           "value" = c(0.7212, 0.8182), "ASE" = c(0.03523, 0.02548),
                           "z" = c(20.47, 32.11), "p" = c(4.009e-93, 3.614e-226),
                           "lower" = c(0.6522, 0.7682), "upper" = c(0.7903, 0.8681)) %>%
    flextable() %>% autofit())
#
#
##Figures
#Agreement direct counts
Agree_count <- Counts %>% ggplot()+
  geom_jitter(aes(Count1, Count2, shape = Agree), size =2, width = 0.25)+
  geom_abline(intercept = 0, slope = 1, size = 0.75)+
  geom_abline(intercept = -1, slope = 1, linetype = "dashed", size = 0.75) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", size = 0.75) +
  scale_shape_manual(values = c(10, 4, 16))+
  scale_x_continuous(name = "Reader 1 Counts", expand = c(0,0), limits = c(0, 8)) + 
  scale_y_continuous(name = "Reader 2 Counts", expand = c(0,0), limits = c(0, 8)) +
  basetheme + theme(legend.position = "none")
#
#Agreement frequencies
Agree_freq <- Counts %>% mutate(diff = as.numeric(Count1-Count2)) %>%
  ggplot(aes(x = diff))+
  geom_bar(aes(y = (..count..)/sum(..count..)), width = 0.85)+
  scale_x_continuous(name = "Difference between readers (R1 - R2)", expand = c(0,0), limits = c(-3.25, 3.25))+
  scale_y_continuous(name = "Frequency", expand = c(0,0), limits = c(0, 1))+
  basetheme
#
#Arrange and save figrue
ggarrange(Agree_count, Agree_freq, nrow = 1)
ggsave(path = "Output/Figures", 
       filename = paste("F1_Reader_agreement_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)

#
#
rm(Tab, Counts, Sam_striae, Samples) #Remove items to clean up workspace. Keep figures and tables.
#
#
#
####A. Regressions - continuous age v SL####
#
###Scale striae into continuous aging
#Within each striae group, scale marginal growth 0-0.95 (since 1 would be the next striae)
Morphology <- left_join(Morphology, 
                        rbind(Morphology %>% filter(Striae == 0) %>% dplyr::select(ID, MG) %>% drop_na() %>%
                                mutate(scaledMG = scales::rescale(MG, to = c(0,0.95))),
                              Morphology %>% filter(Striae == 1) %>% dplyr::select(ID, MG) %>% drop_na() %>%
                                mutate(scaledMG = scales::rescale(MG, to = c(0,0.95))),
                              Morphology %>% filter(Striae == 2) %>% dplyr::select(ID, MG) %>% drop_na() %>%
                                mutate(scaledMG = scales::rescale(MG, to = c(0,0.95))),
                              Morphology %>% filter(Striae == 3) %>% dplyr::select(ID, MG) %>% drop_na() %>%
                                mutate(scaledMG = scales::rescale(MG, to = c(0,0.95))),
                              Morphology %>% filter(Striae == 4) %>% dplyr::select(ID, MG) %>% drop_na() %>%
                                mutate(scaledMG = scales::rescale(MG, to = c(0,0.95))),
                              Morphology %>% filter(Striae == 5) %>% dplyr::select(ID, MG) %>% drop_na() %>%
                                mutate(scaledMG = scales::rescale(MG, to = c(0,0.95))))) %>%
  mutate(scaledAge = Striae+scaledMG)
#
Morphology %>% ggplot(aes(scaledAge, SL))+
  geom_point()
#
#
###Shell length
#
ShellLengths <- Morphology %>% dplyr::select(scaledAge, SL) %>% drop_na()
#
#Check data
boxplot(SL ~ scaledAge, data = Morphology)
gghistogram(Morphology$SL)
#Unequal variances, normal distribution
#
set.seed(54321)
S_SL <- lm(SL ~ scaledAge, data = ShellLengths)
summary(S_SL)
#Residuals
plot(fitted(S_SL), resid(S_SL))
#
#Get weights and run weighted lm
SL_wt <- 1/lm(abs(S_SL$residuals) ~ S_SL$fitted.values)$fitted.values^2
set.seed(54321)
wt_SL <- lm(SL ~ scaledAge, data = ShellLengths, weights = SL_wt)
#
#Summary and AIC comparisons
summary(wt_SL)
AIC(S_SL) #1555.896
AIC(wt_SL)#1516.496 - better
#
S_SL_sum <- summary(wt_SL)
#MLR table with test values
S_SL_tab <- tidy(wt_SL)
names(S_SL_tab) <- c("term", "Est.", "SE", "t", "p-value")
S_SL_sum_tab <- glance(wt_SL) %>% dplyr::select(r.squared:df, deviance:df.residual)
names(S_SL_sum_tab) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
confint(wt_SL)
#
S_SL_fill <- ShellLengths[complete.cases(ShellLengths), ] 
S_SLmeans <- Morphology %>% mutate(binAge = lencat(scaledAge, 0, w = 0.5)) %>%drop_na(binAge) %>%
  group_by(binAge) %>% summarise(aveSL = mean(SL, na.rm = T), sdSL = sd(SL, na.rm = T)) 
#
predicted_S_SL <- data.frame(predSL = predict(wt_SL, S_SL_fill), scaledAge = S_SL_fill$scaledAge)
S_SL_CI <- predict(wt_SL, interval = "predict")
model_S_SL <- cbind(predicted_S_SL, S_SL_CI) %>% 
  dplyr::select(-fit) %>% dplyr::select(scaledAge, predSL, everything()) %>%
  group_by(scaledAge) %>% dplyr::summarise(predSL = mean(predSL, na.rm = T), lwr = mean(lwr), upr = mean(upr))
#
#
(S_SLwCI <- ggplot()+
  #geom_point(aes(Striae_grp, ave, color = "Mean"))+
  geom_jitter(data = Morphology, aes(scaledAge, SL, color = "Observed"), width = 0.15, alpha = 0.4)+
  geom_point(data = S_SLmeans, aes(binAge, aveSL, color = "Mean"), shape = 15, size = 4.5)+
  geom_line(data = model_S_SL, aes(scaledAge, predSL, color = "Predict"))+
  geom_line(data = model_S_SL, aes(scaledAge, lwr, color = "95% CI"), linetype = "dashed")+
  geom_line(data = model_S_SL, aes(scaledAge, upr, color = "95% CI"), linetype = "dashed")+
  ylab("Shell length (mm)")+ 
  xlab("Number of striae")+
  scale_x_continuous(expand = c(0,0.1), limits = c(0, 6)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,125)) +
  scale_color_manual(name = "",
                     breaks = c("Observed", "Mean", "Predict", "95% CI"),
                     values = c("#000000", "#000000", "#FF0000", "#999999"),
                     labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "blank", "solid", "dashed"),
                       shape = c(19, 15, NA, NA))))+
  basetheme)
#
#Check equation and add to figure
makeRegEqnLabel(wt_SL, "Striae_grp", "SL")
(S_SLwCI_eq <- S_SLwCI +
  annotate("text", x = 1.2, y = 110, label = "SL = 10.39 * Striae + 25.71",
           size = 3.5, color = "black", family = "serif")+
  annotate("text", x = 0.9, y = 103, label = paste("R^2 ==", round(S_SL_sum$r.squared, 2)), parse = TRUE,
           size = 3.5, color = "black", family = "serif"))
#
##- Save figure -- Organize S_SL_tab, S_SL_sum_tab, and S_SLmeans for output
ggsave(path = "Output/Figures", 
       filename = paste("F2A_ShellLength_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
(S_SL_table <- S_SL_tab %>% flextable() %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(S_SL_sum_table <- S_SL_sum_tab %>% flextable() %>% 
    colformat_num(j = c(6, 8), digits = 0) %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>%autofit())
#
rm(S_SL, S_SL_CI, S_SL_sum, S_SL_fill, predicted_S_SL, ShellLengths, model_S_SL, wt_SL, SL_wt, 
   S_SLwCI, S_SL_tab, S_SL_sum_tab)
#
#
#
####A. Regressions - continuous age v SW####
#
#scaledAge vs SW
ShellWidths <- Morphology %>% dplyr::select(scaledAge, SW) %>% drop_na()
#
#Check data
boxplot(SW ~ scaledAge, data = Morphology)
gghistogram(Morphology$SW)
#Unequal variances, normal distribution
#
#Regression
set.seed(54321)
S_SW <- lm(SW ~ scaledAge, data = ShellWidths)
summary(S_SW)
#
#Residuals
plot(fitted(S_SW), resid(S_SW))
#
#Get weights and run weighted model
SW_wt <- 1/lm(abs(S_SW$residuals) ~ S_SW$fitted.values)$fitted.values^2
set.seed(54321)
wt_SW <- lm(SW ~ scaledAge, data = ShellWidths, weights = SW_wt)
#
(S_SW_sum <- summary(wt_SW))
#MLR table with test values
S_SW_tab <- tidy(wt_SW)
names(S_SW_tab) <- c("term", "Est.", "SE", "t", "p-value")
S_SW_sum_tab <- glance(wt_SW) %>% dplyr::select(r.squared:df, deviance:df.residual)
names(S_SW_sum_tab) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
confint(wt_SW)
#
S_SW_fill <- ShellWidths[complete.cases(ShellWidths), ] 
S_SWmeans <- ShellWidths %>% mutate(binAge = lencat(scaledAge, 0, w = 0.5)) %>%
  group_by(binAge) %>% summarise(aveSW = mean(SW, na.rm = T), sdSW = sd(SW, na.rm = T)) 
#
predicted_S_SW <- data.frame(predSW = predict(wt_SW, S_SW_fill), scaledAge = S_SW_fill$scaledAge)
S_SW_CI <- predict(wt_SW, interval = "predict")
model_S_SW <- cbind(predicted_S_SW, S_SW_CI) %>% 
  dplyr::select(-fit) %>% dplyr::select(scaledAge, predSW, everything()) %>%
  group_by(scaledAge) %>% dplyr::summarise(predSW = mean(predSW, na.rm = T), lwr = mean(lwr), upr = mean(upr))

#
(S_SWwCI  <- ggplot()+
  geom_jitter(data= Morphology, aes(scaledAge, SW, color = "Observed"), width = 0.15, alpha = 0.4)+
  geom_point(data = S_SWmeans, aes(binAge, aveSW, color = "Mean"), shape = 15, size = 4.5)+
  geom_line(data = model_S_SW, aes(scaledAge, predSW, color = "Predict"))+
  geom_line(data = model_S_SW, aes(scaledAge, lwr, color = "95% CI"), linetype = "dashed")+
  geom_line(data = model_S_SW, aes(scaledAge, upr, color = "95% CI"), linetype = "dashed")+
  ylab("Shell width (mm)")+ 
  xlab("Number of striae")+
  scale_x_continuous(expand = c(0,0.1), limits = c(0, 6)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
  scale_color_manual(name = "",
                     breaks = c("Observed" , "Mean", "Predict", "95% CI"),
                     values = c("#000000", "#000000", "#FF0000", "#999999"),
                     labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "blank", "solid", "dashed"),
                       shape = c(19, 15, NA, NA))))+
  basetheme)
#
makeRegEqnLabel(wt_SW, "Striae_grp", "SW")
#
(S_SWwCI_eq <- S_SWwCI + 
  annotate("text", x = 1.15, y = 52, label = "SW = 4.67 * Striae + 14.01",
           size = 3.5, color = "black", family = "serif")+
  annotate("text", x = 0.9, y = 47.5, label = paste("R^2 ==", round(S_SW_sum$r.squared, 2)), parse = TRUE,
           size = 3.5, color = "black", family = "serif"))
#
##- Save figure -- Organize S_SL_tab, S_SL_sum_tab, and S_SLmeans for output
ggsave(path = "Output/Figures", 
       filename = paste("F2B_ShellWidth_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
(S_SW_table <- S_SW_tab %>% flextable() %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(S_SW_sum_table <- S_SW_sum_tab %>% flextable() %>% colformat_num(j = c(6, 8), digits = 0) %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
#
rm(ShellWidths,S_SW, S_SW_CI, S_SW_tab, S_SW_sum, S_SW_sum_tab, S_SW_fill, predicted_S_SW, model_S_SW,
   wt_SW, SW_wt, S_SWwCI)
#
#
#
#
####A. Regressions - continuous age v WW####
#
#scaledAge vs WW
#
ShellWeight <- Morphology %>% dplyr::select(scaledAge, WW) %>% drop_na()
#Check data
boxplot(WW ~ scaledAge, data = ShellWeight)
gghistogram(ShellWeight$WW)
#Unequal variances, non-normal distribution
ShellWeight <- ShellWeight %>% mutate(logWW = log(WW+1),
                                      sqrtWW = sqrt(WW))
gghistogram(ShellWeight$sqrtWW)
gghistogram(ShellWeight$logWW)
boxplot(sqrtWW ~ scaledAge, data = ShellWeight)
#
#Regression
set.seed(54321)
S_WW <- lm(sqrtWW ~ scaledAge, data = ShellWeight)
summary(S_WW)
#
#Residuals
plot(fitted(S_WW), resid(S_WW))
#
#Get weights and run weighted model
WW_wt <- 1/lm(abs(S_WW$residuals) ~ S_WW$fitted.values)$fitted.values^2
set.seed(54321)
wt_WW <- lm(sqrtWW ~ scaledAge, data = ShellWeight, weights = WW_wt)
#
#
(S_WW_sum <- summary(wt_WW))
S_WW_tab <- tidy(wt_WW)
names(S_WW_tab) <- c("term", "Est.", "SE", "t", "p-value")
S_WW_sum_tab <- glance(wt_WW) %>% dplyr::select(r.squared:df, deviance:df.residual)
names(S_WW_sum_tab) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
#
confint(wt_WW)
#
S_WW_fill <- ShellWeight[complete.cases(ShellWeight), ] 
S_WWmeans <- ShellWeight %>% mutate(binAge = lencat(scaledAge, 0, w = 0.5)) %>%
  group_by(binAge) %>% summarise(aveWW = mean(WW, na.rm = T), sdWW = sd(WW, na.rm = T)) 
#
predicted_S_WW <- data.frame(predWW = (predict(wt_WW, S_WW_fill)^2), scaledAge = S_WW_fill$scaledAge)
S_WW_CI <- predict(wt_WW, interval = "predict")^2
model_S_WW <- cbind(predicted_S_WW, S_WW_CI) %>% 
  dplyr::select(-fit) %>% dplyr::select(scaledAge, predWW, everything()) %>%
  group_by(scaledAge) %>% dplyr::summarise(predWW = mean(predWW, na.rm = T), lwr = mean(lwr), upr = mean(upr))
#
#
(S_WWwCI  <- ggplot()+
  geom_jitter(data = Morphology, aes(scaledAge, WW, color = "Observed"), width = 0.15, alpha = 0.4)+
  geom_point(data = S_WWmeans, aes(binAge, aveWW, color = "Mean"), shape = 15, size = 4.5)+
  geom_line(data = model_S_WW, aes(scaledAge, predWW, color = "Predict"))+
  geom_line(data = model_S_WW, aes(scaledAge, lwr, color = "95% CI"), linetype = "dashed")+
  geom_line(data = model_S_WW, aes(scaledAge, upr, color = "95% CI"), linetype = "dashed")+
  ylab("Wet weight (g)")+ 
  xlab("Number of striae")+
  scale_x_continuous(expand = c(0,0.1), limits = c(0, 6)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,125), breaks = seq(0, 125, by = 25)) +
  scale_color_manual(name = "",
                     breaks = c("Observed", "Mean", "Predict", "95% CI"),
                     values = c("#000000", "#000000", "#FF0000", "#999999"),
                     labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "blank", "solid", "dashed"),
                       shape = c(19, 15, NA, NA))))+
  basetheme)
#
makeRegEqnLabel(wt_WW, "Striae_grp", "sqrtWW")
#
(S_WWwCI_eq <- S_WWwCI +
  annotate("text", x = 1.25, y = 110, label = "sqrtWW = 0.99 * Striae + 1.25",
           size = 3.5, color = "black", family = "serif")+
  annotate("text", x = 1.1, y = 103, label = paste("R^2 ==", round(S_WW_sum$r.squared, 2)), parse = TRUE,
           size = 3.5, color = "black", family = "serif"))
#
##- Save figure -- Organize S_SL_tab, S_SL_sum_tab, and S_SLmeans for output
ggsave(path = "Output/Figures", 
       filename = paste("F2C_ShellWeight_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
(S_WW_table <- S_WW_tab %>% flextable() %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(S_WW_sum_table <- S_WW_sum_tab %>% flextable() %>% colformat_num(j = c(6, 8), digits = 0) %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
#
rm(S_WW, S_WW_CI, S_WW_sum, S_WW_fill, predicted_S_WW, ShellWeight, model_S_WW, wt_WW, WW_wt, 
   S_WWwCI, S_WW_tab, S_WW_sum_tab)
#
#
#
####A. Regressions - continuous age v OL####
#
#scaledAge vs OL
#
OpLength <- Morphology %>% dplyr::select(scaledAge, Opercula_L) %>% drop_na()
#Check data
boxplot(Opercula_L ~ scaledAge, data = OpLength)
gghistogram(OpLength$Opercula_L)
#Unequal variances, normal distribution
#
#Regression
set.seed(54321)
S_OL <- lm(Opercula_L ~ scaledAge, data = OpLength)
summary(S_OL)
#
#Residuals
plot(fitted(S_OL), resid(S_OL))
#
#Get weights adn run weighted model
OL_wt <- 1/lm(abs(S_OL$residuals) ~ S_OL$fitted.values)$fitted.values^2
set.seed(54321)
wt_OL <- lm(Opercula_L ~ scaledAge, data = OpLength, weights = OL_wt)
#
#
(S_OL_sum <- summary(wt_OL))
#MLR table with test values
S_OL_tab <- tidy(wt_OL)
names(S_OL_tab) <- c("term", "Est.", "SE", "t", "p-value")
S_OL_sum_tab <- glance(wt_OL) %>% dplyr::select(r.squared:df, deviance:df.residual)
names(S_OL_sum_tab) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
#
confint(wt_OL, level = 0.95)
#
S_OL_fill <- OpLength[complete.cases(OpLength), ] 
S_OLmeans <- OpLength %>% mutate(binAge = lencat(scaledAge, 0, w = 0.5)) %>%
  group_by(binAge) %>% summarise(aveOL = mean(Opercula_L, na.rm = T), sdOL = sd(Opercula_L, na.rm = T)) 
#
predicted_S_OL <- data.frame(predOL = predict(wt_OL, S_OL_fill), scaledAge = S_OL_fill$scaledAge)
S_OL_CI <- predict(wt_OL, interval = "predict")
model_S_OL <- cbind(predicted_S_OL, S_OL_CI) %>% 
  dplyr::select(-fit) %>% dplyr::select(scaledAge, predOL, everything()) %>%
  group_by(scaledAge) %>% dplyr::summarise(predOL = mean(predOL, na.rm = T), lwr = mean(lwr), upr = mean(upr))
#
#
(S_OLwCI  <- ggplot()+
  geom_jitter(data = Morphology, aes(scaledAge, Opercula_L, color = "Observed"), width = 0.15, alpha = 0.4)+
  geom_point(data = S_OLmeans, aes(binAge, aveOL, color = "Mean"), shape = 15, size = 4.5)+
  geom_line(data = model_S_OL, aes(scaledAge, predOL, color = "Predict"))+
  geom_line(data = model_S_OL, aes(scaledAge, lwr, color = "95% CI"), linetype = "dashed")+
  geom_line(data = model_S_OL, aes(scaledAge, upr, color = "95% CI"), linetype = "dashed")+
  ylab("Operculum length (mm)")+ 
  xlab("Number of striae")+
  scale_x_continuous(expand = c(0,0.1), limits = c(0, 6)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
  scale_color_manual(name = "",
                     breaks = c("Observed", "Mean", "Predict", "95% CI"),
                     values = c("#000000", "#000000", "#FF0000", "#999999"),
                     labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "blank", "solid", "dashed"),
                       shape = c(19, 15, NA, NA))))+
  basetheme)
#
makeRegEqnLabel(wt_OL, "Striae_grp", "OL")
#
(S_OLwCI_eq <- S_OLwCI +
  annotate("text", x = 1.15, y = 52, label = "OL = 4.62 * Striae + 10.56",
           size = 3.5, color = "black", family = "serif")+
  annotate("text", x = 1.1, y = 47.5, label = paste("R^2 ==", round(S_OL_sum$r.squared, 2)), parse = TRUE,
           size = 3.5, color = "black", family = "serif"))
#
##- Save figure -- Organize S_SL_tab, S_SL_sum_tab, and S_SLmeans for output
ggsave(path = "Output/Figures", 
       filename = paste("F2D_OperculumLength_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
(S_OL_table <- S_OL_tab %>% flextable() %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(S_OL_sum_table <- S_OL_sum_tab %>% flextable() %>% colformat_num(j = c(6, 8), digits = 0) %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
#
rm(S_OL, S_OL_CI, S_OL_sum, S_OL_fill, predicted_S_OL, OpLength, model_S_OL, wt_OL, OL_wt, 
   S_OLwCI, S_OL_tab, S_OL_sum_tab)
#
#
#
####A. Regression summary####
#
#Combine regression table means for output
(Regression_sum_tab <- left_join(S_SLmeans, S_SWmeans) %>% 
   left_join(S_WWmeans) %>% left_join(S_OLmeans) %>%
   flextable() %>% vline(j = c(1, 3, 5, 7)) %>% autofit())
#
rm(S_SLmeans, S_SWmeans, S_WWmeans, S_OLmeans)
#
ggarrange(S_SLwCI_eq + rremove("x.title"), S_SWwCI_eq + rremove("x.title"), 
          S_WWwCI_eq, S_OLwCI_eq, 
          nrow = 2, ncol = 2, common.legend = TRUE)
##- Save figure 
ggsave(path = "Output/Figures", 
       filename = paste("F2_Regressions_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
rm(S_SLwCI_eq, S_SWwCI_eq, S_WWwCI_eq, S_OLwCI_eq)
#
#
####WORKING#####
#
summary(Morphology$SL)
quantile(Morphology$SL,probs=c(0.025,0.975)) #80.5625
Morphology %>% mutate(LW_r = (SL/WW)*-1,
                      scaledLW_r = scales::rescale(LW_r, to=c(0, 80.5625))) %>% #scale to SL with maximum being 97.5% quaantile as est for Linf (ave max SL)
  ggplot(aes(scaledAge, scaledLW_r))+ geom_point()+
  geom_point(aes(scaledAge, SL), color = "blue")
#
#
#
#
#
Morphology %>% filter(Striae_grp > 0) %>%
  ggplot(aes(as.factor(Month), MG, color = as.factor(Striae_grp)))+
  geom_point()+
  geom_point(data = Morphology %>% filter(Striae_grp > 0) %>%
               group_by(Striae_grp, Month) %>% summarise(Range_diff = max(MG)-min(MG),
                                                         Mean = mean(MG)),
             aes(as.factor(Month), Mean), color = "black")+
  lemon::facet_rep_grid(Striae_grp~., scales = "free")
#
temp <- Morphology %>% filter(Striae_grp > 0) %>%
  group_by(Striae_grp, Month) %>% summarise(Mean = mean(MG)) %>% drop_na(Mean)
#
Morphology %>% mutate(MonYr = zoo::as.yearmon(MonYr, format = "%b %Y")) %>%
  filter(Striae_grp > 0) %>%
  ggplot(aes(MonYr, MGR))+
  geom_point()+
  lemon::facet_rep_grid(Striae_grp~., scales = "free")

#
####B. Best model - continuous(age)####
#
ShellLengths <- Morphology %>% dplyr::select(scaledAge, SL) %>% drop_na()
#Check for major outliers
ShellLengths %>% ggplot(aes(scaledAge, SL))+
  geom_point()
boxplot(Morphology$SL ~ Morphology$scaledAge)
gghistogram(Morphology$scaledAge)
#Unequeal var but normal
#
Ages <- seq(0, 10, by = 0.25)
S_SLmeans <- Morphology %>% group_by(Striae_grp) %>% summarise(aveSL = mean(SL, na.rm = TRUE),
                                                               sdSL = sd(SL, na.rm = TRUE))
#
#
###WLS
set.seed(54321)
S_SL <- lm(SL ~ scaledAge, data = ShellLengths)
SL_wt <- 1/lm(abs(S_SL$residuals) ~ S_SL$fitted.values)$fitted.values^2
#
##Model and summary
set.seed(54321)
wt_SL <- lm(SL ~ scaledAge, data = ShellLengths, weights = SL_wt)
summary(wt_SL) #R2 = 0.8074
AIC(wt_SL) #AIC = 1516.496
WLS_sum <- data.frame("Model" = "Weighted", "R2" = 0.8074, "AIC" = 1516.496)
#
WLSpred <- data.frame(scaledAge = Ages,
                      predSL = predict(wt_SL, list(scaledAge = Ages)))
(WLS_est <- WLSpred %>% filter(scaledAge %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%
    mutate(Model = "Weighted"))
#
#
#
###Linear
S_SLi <- lm(SL ~ scaledAge, data = ShellLengths)
par(mfrow = c(2,2))
plot(S_SLi)
par(mfrow = c(1,1))
##Slight hetero
S_SL2 <- lm(log(SL) ~ scaledAge, data = ShellLengths)
par(mfrow = c(2,2))
plot(S_SL2)
par(mfrow = c(1,1))
summary(S_SL2)
#Residuals
plot(fitted(S_SL2), resid(S_SL2))
#
Lin_sum <- glance(S_SL2) %>%  dplyr::select(r.squared, AIC) %>% 
  rename("R2" = r.squared) %>% mutate(Model = "Linear")
#
set.seed(54321)
LinPred <- data.frame(scaledAge = Ages, 
                      predSL = exp(predict(S_SL2, list(scaledAge = Ages)))) 
(Lin_est <- LinPred %>% filter(scaledAge %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%
  mutate(Model = "Linear"))
#
#
#
#
#
####vonBert
vbTyp <- vbFuns(param = "Typical") #SL~Linf*(1-exp(-K*(Striae_grp-t0)))
#
(svTyp <- vbStarts(SL~scaledAge, data = ShellLengths))
svStart <- list(Linf = 85, K = 0.2444347, t0 = -2.019411)
#
VBfit <- nls(SL ~ vbTyp(scaledAge, Linf, K, t0), data = ShellLengths, start = svStart)
#
summary(VBfit)
glance(VBfit)
1-((sum(residuals(VBfit)^2))/sum((ShellLengths$SL-mean(ShellLengths$SL, na.rm = T))^2)) 
AIC(VBfit)
#
vonB_sum <- data.frame("Model" = "vonBert", "R2" = 0.7172, "AIC" = 1554.576)
#
set.seed(54321)
VBPred <- data.frame(scaledAge = Ages, predSL = predict(VBfit, list(scaledAge = Ages), type = "probs"))
(VB_est <- VBPred %>% filter(scaledAge %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%
    mutate(Model = "VonBert"))
#
#
#
#
#
####Gompertz
#
plot(SL ~ scaledAge, data = ShellLengths)
getInitial(SL ~ SSgompertz(scaledAge, Asym, b2, b3), data = ShellLengths)
#
gomfit <- nls(SL ~ SSgompertz(scaledAge, Asym, b2, b3), data = ShellLengths)
1-((sum(residuals(gomfit)^2))/sum((ShellLengths$SL-mean(ShellLengths$SL, na.rm = T))^2)) 
AIC(gomfit)
#
gomp_sum <- data.frame("Model" = "Gompertz", "R2" = 0.7182, "AIC" = 1553.832)
#
set.seed(54321)
GomPred <- data.frame(scaledAge = Ages, predSL = predict(gomfit, list(scaledAge = Ages)))
(Gomp_est <- GomPred %>% filter(scaledAge %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%
    mutate(Model = "Gompertz"))
#
#
#
#
####Logistic
#
Logfit <- nls(SL ~ SSlogis(scaledAge, Asym, xmid, scal), data = ShellLengths)
#
summary(Logfit)
1-((sum(residuals(Logfit)^2))/sum((ShellLengths$SL-mean(ShellLengths$SL, na.rm = T))^2)) 
AIC(Logfit)
#
log_sum <- data.frame("Model" = "Logistic", "R2" = 0.7186, "AIC" = 1553.508)
#
set.seed(54321)
LogPred <- data.frame(scaledAge = Ages, predSL = predict(Logfit, list(scaledAge = Ages)))
(Log_est <- LogPred %>% filter(scaledAge %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%
    mutate(Model = "Logistic"))
#
#
#
#Get model terms for output
(Model_terms <- rbind(tidy(wt_SL) %>% mutate(Model = "Weighted"), 
                      tidy(S_SL2) %>% mutate(Model = "Linear"),
                      tidy(VBfit) %>% mutate(Model = "von Bertalanffy"),
                      tidy(gomfit) %>% mutate(Model = "Gompertz"),
                      tidy(Logfit) %>% mutate(Model = "Logistic")) %>% 
    dplyr::select(Model, term, estimate))
#
#
#Comparison
(All_models <- ggplot(ShellLengths)+
  geom_jitter(aes(scaledAge, SL, color = "Mean"), width = 0.15, alpha = 0.4)+
  geom_line(data = WLSpred, aes(scaledAge, predSL, color = "WLS"), linetype = "solid")+
  geom_line(data = LinPred, aes(scaledAge, predSL, color = "Linear"), linetype = "solid")+
  geom_line(data = VBPred, aes(scaledAge, predSL, color = "vonBert"), linetype = "dashed")+
  geom_line(data = GomPred, aes(scaledAge, predSL, color = "Gompertz"), linetype = "dotted")+  
  geom_line(data = LogPred, aes(scaledAge, predSL, color = "Logistic"), linetype = "longdash")+
  ylab("Shell length (mm)")+ 
  xlab("Number of striae")+
  scale_x_continuous(expand = c(0,0.1), limits = c(0, 8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,125)) +
  geom_hline(yintercept = 85, linetype = "dotdash")+
  scale_color_manual(name = "",
                     breaks = c("Mean", "WLS", "Linear", "vonBert", "Gompertz", "Logistic"),
                     values = c("#000000", "#000000", "#999999", "#333333", "#333333", "#333333"),
                     labels = c("Observed", "WLS", "Linear", "vonBert", "Gompertz", "Logistic"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "solid", "solid", "dashed", "dotted", "longdash"))))+
  basetheme) #+ theme(legend.position = "none")
#
ggsave(path = "Output/Figures", 
       filename = paste("F3A_Model_comparisons_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
#
###Best model = Logistic
#
summary(Logfit)
#Predictions
Log_pred <- predictNLS(Logfit, data.frame(scaledAge = Ages), interval = "prediction", alpha = 0.05, nsim = 10000)
#
#Combine ages with predicted means and CI
model_S_SL <- cbind(LogPred$scaledAge, 
                    Log_pred$summary[,c(7, 11:12)]) %>% rename(Mean = "Sim.Mean", 
                                                               Pred_lwr = "Sim.2.5%", 
                                                               Pred_upr = "Sim.97.5%") %>%
  rename(scaledAge = "LogPred$scaledAge")
#
(Best_model <- ggplot(ShellLengths)+
  geom_jitter(aes(scaledAge, SL, color = "Observed"), width = 0.1, alpha = 0.4)+
  geom_point(data = S_SLmeans, aes(Striae_grp, aveSL, color = "Mean"), shape = 15, size = 4.5)+
  geom_line(data = model_S_SL, aes(scaledAge, mean, color = "Predict"))+
  geom_line(data = model_S_SL, aes(scaledAge, lwr, color = "95% CI"), linetype = "dashed")+
  geom_line(data = model_S_SL, aes(scaledAge, upr, color = "95% CI"), linetype = "dashed")+
  ylab("Shell length (mm)")+ 
  xlab("Number of striae")+
  scale_x_continuous(expand = c(0,0.1), limits = c(0, 8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,125)) +
  geom_hline(yintercept = 85, linetype = "dotdash")+
  scale_color_manual(name = "",
                     breaks = c("Observed", "Mean", "Predict", "95% CI"),
                     values = c("#000000", "#000000", "#FF0000", "#999999"),
                     labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "blank", "solid", "dashed"),
                       shape = c(19, 15, NA, NA))))+
  basetheme + theme(legend.position = "none"))
#
ggsave(path = "Output/Figures", 
       filename = paste("F3B_BestModel_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
rm(ShellLengths, Ages, S_SLmeans, S_SL, SL_wt, wt_SL, S_SLi, S_SL2, vbTyp, svTyp, svStart, VBfit, 
   gomfit, Logfit, Log_pred, model_S_SL)
#
#
#
#
####B. Best model summary####
#
##Model predictions
(Model_preds <- left_join(WLSpred %>% rename("Weighted" = predSL), 
                         LinPred %>% rename("Linear" = predSL))  %>% 
  left_join(VBPred %>% rename("vonBert" = predSL)) %>%
  left_join(GomPred %>% rename("Gompertz" = predSL)) %>%
   flextable() %>% autofit())
#
#Model evals
(Model_eval <- left_join(Model_terms, rbind(WLS_sum, Lin_sum, 
                                            vonB_sum %>% mutate(Model = "von Bertalanffy"), gomp_sum, log_sum)) %>% 
    flextable() %>% set_formatter(R2 = function(x){formatC(x, format = "e", digits = 2)}) %>%
    merge_v(j = c(1, 4, 5)) %>% hline(i = c(2, 4, 7, 10)) %>% autofit())
#
#Model figures
(Models <- ggarrange(All_models + rremove("legend"), 
                    Best_model + rremove("y.title"), nrow = 1))
#
ggsave(path = "Output/Figures", 
       filename = paste("F3_Models_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
rm(WLSpred, WLS_est, LinPred, Lin_est, VBPred, VB_est, GomPred, Gomp_est, WLS_sum, Lin_sum, 
   vonB_sum, gomp_sum, log_sum, All_models, Best_model)
#
#
#
#
#
####C. Length at age####
#
#
###Using SL vs OL relationship for back calculations - need model of relationship
###SL vs OL regression 
set.seed(54321)
SLOL_mod1 <- lm(SL ~ Opercula_L, data = Morphology)
OL_SL_sum1 <- summary(SLOL_mod1)
(OL_SL_tab1 <- tidy(SLOL_mod1))
names(OL_SL_tab1) <- c("term", "Est.", "SE", "t", "p-value")
(OL_SL_sum_tab1 <- glance(SLOL_mod1) %>% dplyr::select(r.squared:df, deviance:df.residual))
names(OL_SL_sum_tab1) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
#
OL_SL_fill1 <- (Morphology %>% dplyr::select(SL, Opercula_L))[complete.cases(Morphology %>% dplyr::select(SL, Opercula_L)),]
(SL_OLmeans <- Morphology %>% group_by(BinOL) %>% summarise(aveSL = mean(SL, na.rm = T)) %>% drop_na(BinOL))
#
predicted_OL_SL <- data.frame(predSL = predict(SLOL_mod1, OL_SL_fill1), BinOL = OL_SL_fill1$Opercula_L)
OL_SL_CI1 <- predict(SLOL_mod1, interval = "prediction")
model_OL_SL1 <- cbind(predicted_OL_SL, OL_SL_CI1) %>% 
  dplyr::select(-fit) %>% dplyr::select(BinOL, predSL, everything()) %>%
  group_by(BinOL) %>% dplyr::summarise(predSL = mean(predSL, na.rm = T), lwr = mean(lwr), upr = mean(upr))
#
(OL_SLwCI <- ggplot()+
  geom_jitter(data = Morphology, aes(Opercula_L, SL, color = "Observed"), width = 0.15, alpha = 0.4)+
  geom_point(data = SL_OLmeans, aes(BinOL, aveSL, color = "Mean"), shape = 15, size = 4.5)+
  geom_line(data = model_OL_SL1, aes(BinOL, predSL, color = "Predict"))+
  geom_line(data = model_OL_SL1, aes(BinOL, lwr, color = "95% CI"), linetype = "dashed")+
  geom_line(data = model_OL_SL1, aes(BinOL, upr, color = "95% CI"), linetype = "dashed")+
  ylab("Shell length (mm)")+ 
  xlab("Operculum length (mm)")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 50)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,125)) +
  scale_color_manual(name = "",
                     breaks = c("Observed", "Mean", "Predict", "95% CI"),
                     values = c("#000000", "#000000", "#FF0000", "#999999"),
                     labels = c("Observed", "Observed mean", "Predicted mean", "95% confidence limit"),
                     guide = guide_legend(override.aes = list(
                       linetype = c("blank", "blank", "solid", "dashed"),
                       shape = c(19, 15, NA, NA))))+
  basetheme + theme(legend.position = "none"))
#
makeRegEqnLabel(SLOL_mod1, "Opercula_L", "SL")
(OL_SLwCI +
  annotate("text", x = 12, y = 110, label = "OL = 2.22 * SL + 2.76",
           size = 3.5, color = "black", family = "serif")+
  annotate("text", x = 12, y = 103, label = paste("R^2 ==", round(OL_SL_sum1$r.squared, 2)), parse = TRUE,
           size = 3.5, color = "black", family = "serif"))
#
ggsave(path = "Output/Figures", 
       filename = paste("F4_Operculum vs SL_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
#
#
#
##Age length key construction
age_length <- Morphology %>% dplyr::select(ID, SL, Striae_grp, BinSL) %>% 
  drop_na(Striae_grp) %>% filter(ID != "BT20-217")
#
##Observed age-length key
#Frequency per size group
(al_freq <- xtabs(~BinSL+Striae_grp, data = age_length))
rowSums(al_freq)
#
alk <- prop.table(al_freq, margin = 1)
round(alk, 3)
#
#
##Age-Length Key Figure
alkPlot(alk, type = "area", showLegend = T, leg.cex = 0.7, xlab = "Shell length (mm)", col = hcl.colors(5, "Grays"))
#Replicate in ggplot to modify text
as.data.frame(alk) %>% 
  mutate(BinSL = as.numeric(paste(BinSL)), #Change to numeric for sequential plotting
         Striae_grp = fct_reorder(Striae_grp, desc(Striae_grp))) %>%  #Reorder levels
  ggplot(aes(BinSL, Freq, fill = Striae_grp))+
  geom_area(position = "fill")+
  scale_x_continuous("Shell length (mm)", expand = c(0,0), limits = c(20, 85), breaks = seq(20, 85, 5))+
  scale_y_continuous("Frequency", expand = c(0,0))+
  scale_fill_grey(limits = c("0", "1", "2", "3", "4"))+
  basetheme + theme(panel.border = element_rect(color = "black", fill = NA),
                    legend.position = "top", legend.box.spacing = unit(-0.1, "cm"),
                    legend.key.height = unit(1, 'cm'),
                    legend.key.width = unit (3.75, 'cm'), 
                    legend.title = element_blank(), legend.background = element_blank(),
                    legend.text = element_text(color= "black", size = 10, family = "serif"))
#
ggsave(path = "Output/Figures", 
       filename = paste("F5_Age Length Proportions_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
#
#
##Applying key
(len_n <- xtabs(~BinSL, data = age_length))
(temp <- sweep(alk, MARGIN = 1, FUN = "*", STATS = len_n))
(ad1 <- colSums(temp)) #Number at each age
(prop.table(ad1)) #Proportion at each age
##Creating ALK
#Age distribution with SEs
alkAgeDist(alk, lenA.n = rowSums(al_freq), len.n = len_n)
#Mean length-at-age
alkMeanVar(alk, SL~BinSL + Striae_grp, data = age_length, len.n = len_n)
#Assign ages to un-aged snails
unaged <- alkIndivAge(alk, Striae_grp ~ SL, data = Morphology %>% 
                        dplyr::select(ID, SL, Striae_grp, BinSL) %>% filter(is.na(Striae_grp)))
final <- rbind(age_length, unaged)
#Final Length at Age means
(AgeLength_sum <- final %>% group_by(Striae_grp) %>% summarise(n = n(), 
                                                              meanSL = mean(SL, na.rm = T),
                                                              sdSL = sd(SL, na.rm  = T),
                                                              seSL = se(SL, na.rm = T)))
#
plot(SL ~ Striae_grp, data = final, pch = 19, xlab = "Striae", ylab = "Shell length (mm)", ylim = c(0, 90))
lines(meanSL ~ Striae_grp, data = AgeLength_sum, lwd = 2, lty = 2)
#
ggplot(AgeLength_sum, aes(Striae_grp, meanSL))+
  geom_jitter(data = final, aes(Striae_grp, SL), color = "black", width = 0.1, alpha = 0.3)+
  geom_point(size = 3.5)+
  #geom_line()+
  geom_smooth(method = "lm", formula = y~x, se = F, color = "black")+
  geom_errorbar(aes(Striae_grp, ymin = meanSL-sdSL, ymax = (meanSL+sdSL)), width = 0.2, size = 1)+
  basetheme+
  scale_x_continuous(name = "Number of striae", expand = c(0,0), limits = c(-0.15,4.15))+
  scale_y_continuous(name = "Shell length (mm)", expand = c(0,0), 
                     limits = c(0,125), breaks = seq(0, 125, 25))
#
ggsave(path = "Output/Figures", 
       filename = paste("F6_Mean Length at Age_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
#
#
####Back calculated lengths
#Comparison of measured Operculum length vs ImageJ measured length
ggplot()+
  geom_smooth(data = Morphology, aes(SL, Opercula_L), se = F, color= "black")+
  geom_smooth(data = Morphology, aes(SL, OL2), se = F, color = "red")#+
#geom_smooth(data = OL_Striae, aes(predSL, Opercula_L), se = F)
#
#Convert data to wide format
OL_Wide <- full_join(
  ##Data frame of all BT with 1 or more striae
  Growth %>% group_by(ID) %>% 
     filter(Number > 1) %>% #Remove OLs to keep just measurements
     mutate(Striae_grp = Number - 1) %>% #Get striae number rather than mark number 
     slice(-n()), #remove last measurement = MG. Keep all >4 for predicting SL
  ##Data frame of all age 0s
  right_join(Growth, Morphology %>% subset(Striae == 0)) %>% 
     dplyr::select(ID, Measure, Number, Sum_mm) %>% filter(Number > 1) %>% 
     mutate(Striae_grp = Number - 2)) %>%
  rename(Opercula_L = Sum_mm) %>% dplyr::select(-Measure, -Ave_mm, -Number) %>% ##Modify output - only keep total OL and Striae counts
  spread(Striae_grp, Opercula_L) %>%  rename_at(3:8, ~ paste("Striae", ., sep = "_")) %>% #Divide into columns of length at each striae count
  #Add columns for Age at capture, SL at capture, and OL at capture
  merge(Morphology %>% dplyr::select(ID, Striae, SL, OL2) %>% rename(AgeCap = Striae, SLCap = SL, OLCap = OL2)) %>%
  dplyr::select(Species, ID, AgeCap, SLCap, OLCap, everything()) %>%
  filter(ID != "BT20-217")
#
##Long format
#Last measurement removed prior so only have striae measurements.
OL_Long <- OL_Wide %>% gather(key = Striae_grp, value = Opercula_L, -Species, -ID, - AgeCap, - SLCap, -OLCap)
#
##Visualize lengths at each age (striae) for all individuals
OL_Long %>% 
  ggplot()+
  geom_point(aes(Striae_grp, Opercula_L))+
  geom_line(aes(Striae_grp, Opercula_L, group = ID), alpha = 0.6)+
  #scale_x_discrete(expand = c(0,0), limits = c(0,6))+
  scale_y_continuous(expand = c(0,0), limits = c(0,45))+
  theme_classic()
#
#
###Regressions - get coefficients
#length-scale
set.seed(54321)
SLOL_mod2 <- lm(SLCap ~ OLCap, data = OL_Wide)
#Out put summary and tables
SLOL_tab2 <- tidy(SLOL_mod2) 
names(SLOL_tab2) <- c("term", "Est.", "SE", "t", "p-value")
SLOL_sum_tab2 <- glance(SLOL_mod2) %>% dplyr::select(r.squared:df, deviance:df.residual)
names(SLOL_sum_tab2) <- c("R2", "adjR2", "RSE", "F", "p-value", "df", "RSS", "Resid.df")
#
LO_coef <- coef(SLOL_mod2)[1]
LS_coef <- coef(SLOL_mod2)[2]
#
#Applying models
OL_Long$FL_lengths <- with(OL_Long, (Opercula_L/OLCap)*(SLCap - LO_coef) + LO_coef)
##Mean lengths - back-calculated
(OL_bc_means <- OL_Long %>% group_by(Striae_grp) %>% drop_na(FL_lengths) %>%
  summarise(n = n(), 
            meanSL = mean(FL_lengths, na.rm = T), 
            sdSL = sd(FL_lengths, na.rm = T)))
#Visually check back calculated means for odd numbers
OL_Long %>%
  ggplot()+
  geom_jitter(aes(Striae_grp, FL_lengths), width = 0.2, alpha = 0.4)+
  scale_x_discrete(name = "Number of Striae", 
                   breaks = c("Striae_0", "Striae_1", "Striae_2", "Striae_3", "Striae_4", "Striae_5"), 
                   labels = c("0", "1", "2", "3", "4", "5"))+
  scale_y_continuous(name = "Back calculated shell length (mm)", limits = c(0, 100), expand = c(0,0))+
  theme_classic()
#
#
#
#
####Compare observed vs back-calc (Beaty and Chen)
#
#Combine observed and back calculated values
Comb_ori_BC <- rbind(Morphology %>% dplyr::select(SL, Striae) %>% drop_na(Striae) %>%
                       mutate(Striae = paste("Obs_", Striae, sep = "")),
                     OL_Long %>% dplyr::select(Striae_grp, FL_lengths) %>% 
                       rename(Striae = Striae_grp, SL = FL_lengths) %>%
                       mutate(Striae = paste("BC_", substr(Striae, 8, 8), sep = "")) %>% drop_na(SL))
#Visual check
ggplot(Comb_ori_BC, aes(Striae, SL))+
  geom_boxplot()
#
###Welch's t-test comparing observed and back calculated SLs for each striae group
(Welchs <- rbind(
  #Age 0 - CohensD
  cbind(tidy(t.test((Comb_ori_BC %>% filter(Striae == "BC_0"))$SL, (Comb_ori_BC %>% filter(Striae == "Obs_0"))$SL)),
        (Comb_ori_BC %>% filter(grepl("0", Striae))) %>% cohens_d(SL ~ Striae, var.equal = FALSE) %>% dplyr::select(effsize:n2)),
  #Age 1
  cbind(tidy(t.test((Comb_ori_BC %>% filter(Striae == "BC_1"))$SL, (Comb_ori_BC %>% filter(Striae == "Obs_1"))$SL)),
        (Comb_ori_BC %>% filter(grepl("1", Striae))) %>% cohens_d(SL ~ Striae, var.equal = FALSE) %>% dplyr::select(effsize:n2)),
  #Age 2  
  cbind(tidy(t.test((Comb_ori_BC %>% filter(Striae == "BC_2"))$SL, (Comb_ori_BC %>% filter(Striae == "Obs_2"))$SL)),
        (Comb_ori_BC %>% filter(grepl("2", Striae))) %>% cohens_d(SL ~ Striae, var.equal = FALSE) %>% dplyr::select(effsize:n2)),
  #Age 3
  cbind(tidy(t.test((Comb_ori_BC %>% filter(Striae == "BC_3"))$SL, (Comb_ori_BC %>% filter(Striae == "Obs_3"))$SL)),
        (Comb_ori_BC %>% filter(grepl("3", Striae))) %>% cohens_d(SL ~ Striae, var.equal = FALSE) %>% dplyr::select(effsize:n2)),
  #Age 4
  cbind(tidy(t.test((Comb_ori_BC %>% filter(Striae == "BC_4"))$SL, (Comb_ori_BC %>% filter(Striae == "Obs_4"))$SL)),
        (Comb_ori_BC %>% filter(grepl("4", Striae))) %>% cohens_d(SL ~ Striae, var.equal = FALSE) %>% dplyr::select(effsize:n2)),
  #Age 5
  cbind(tidy(t.test((Comb_ori_BC %>% filter(Striae == "BC_5"))$SL, (Comb_ori_BC %>% filter(Striae == "Obs_5"))$SL)),
        (Comb_ori_BC %>% filter(grepl("5", Striae))) %>% cohens_d(SL ~ Striae, var.equal = FALSE) %>% dplyr::select(effsize:n2))) %>%
  dplyr::select(-method, -alternative) %>%
  rename(mean_BC = estimate1, mean_Obs = estimate2, t = statistic, df = parameter, Cohen_d = effsize) %>%
    mutate(across(1:9, round, 3)))
#
#
#Figure comparing observed means and back calculated means
Final_comps <- Comb_ori_BC %>% group_by(Striae) %>% summarise(n = n(), 
                                                              meanSL = mean(SL, na.rm = T),
                                                              sdSL = sd(SL, na.rm = T)) %>% drop_na(Striae) %>%
  mutate(Group = gsub("_.$", "", Striae), Age = substr(Striae, nchar(Striae)-1+1, nchar(Striae)))
#
Final_comps %>% 
  ggplot(aes(Age, meanSL, color = Group))+
  geom_point(position = position_dodge(0.25))+
  geom_errorbar(aes(ymin = meanSL - sdSL, ymax = meanSL + sdSL), width = 0.25, position = position_dodge(0.25))+
  xlab("Number of striae")+
  scale_y_continuous(name = "Shell length (mm)",  limits = c(0, 125), expand = c(0,0))+
  basetheme 
#
ggsave(path = "Output/Figures", 
       filename = paste("F7_Obs vs Back Calc Means_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
#
#
Ages <- seq(0, 10, by = 0.05)
##Fit logistic to original and BC data
set.seed(54321)
Logfit_or <- nls(SLCap ~ SSlogis(AgeCap, Asym, xmid, scal), data = (OL_Long %>% dplyr::select(AgeCap, SLCap) %>% drop_na()))
#
summary(Logfit_or)
glance(Logfit_or)
AIC(Logfit_or)
1-((sum(residuals(Logfit_or)^2))/sum((OL_Long$SLCap-mean(OL_Long$SLCap, na.rm = T))^2)) 
#
LogPred_or <- data.frame(AgeCap = Ages, predSL = predict(Logfit_or, list(AgeCap = Ages)))
#
#
#
set.seed(54321)
Logfit_bc <- nls(FL_lengths ~ SSlogis(AgeCap, Asym, xmid, scal), data = (OL_Long %>% dplyr::select(AgeCap, FL_lengths) %>% drop_na()))
#
summary(Logfit_bc)
glance(Logfit_bc)
AIC(Logfit_bc)
1-((sum(residuals(Logfit_bc)^2))/sum((OL_Long$FL_lengths-mean(OL_Long$FL_lengths, na.rm = T))^2)) 
#
LogPred_bc <- data.frame(AgeCap = Ages, predSL = predict(Logfit_bc, list(AgeCap = Ages)))
#
#
ggarrange(ggplot()+
            geom_jitter(data = OL_Long, aes(AgeCap, SLCap), width = 0.25, alpha = 0.4)+
            geom_line(data = LogPred_or, aes(AgeCap, predSL), size = 1)+
            scale_y_continuous(name = "Shell length (mm)", expand = c(0,0), limits = c(0, 125), 
                               breaks = seq(0, 125, 25))+ 
            scale_x_continuous(name = "Age", expand = c(0.01,0), limits = c(0, 10),
                               breaks = seq(0, 10, 2))+
            basetheme,
          ggplot()+
            geom_jitter(data = OL_Long, aes(AgeCap, FL_lengths), width = 0.25, alpha = 0.4, shape = 17)+
            geom_line(data = LogPred_bc, aes(AgeCap, predSL), size = 1)+
            scale_y_continuous(name = "Shell length (mm)", expand = c(0,0), limits = c(0, 125), 
                               breaks = seq(0, 125, 25))+ 
            scale_x_continuous(name = "Age", expand = c(0.01,0), limits = c(0, 10),
                               breaks = seq(0, 10, 2))+
            basetheme + theme(axis.title.y = element_blank()))
#
ggsave(path = "Output/Figures", 
       filename = paste("F8_Obs vand Back Calc Models_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
(Model_comps <- rbind(LogPred_or %>% filter(AgeCap %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>% 
  round(2) %>% mutate(Model = "Original"),
  LogPred_bc %>% filter(AgeCap %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))  %>% 
    round(2) %>% mutate(Model = "Back-calculated")))
#
#
#
rm(age_length, Comb_ori_BC, final, Logfit_bc, Logfit_or, LogPred_bc, LogPred_or, model_OL_SL1, OL_bc_means,
   OL_Long, OL_SL_CI1, OL_SL_fill1, OL_SL_sum1, OL_Wide, predicted_OL_SL, SL_OLmeans, SLOL_mod1, SLOL_mod2,
   unaged)
#
#
#
####C. Length at age summary####
#
#OL vs SL relationship
(OL_SL_tab1_out <- OL_SL_tab1 %>% flextable() %>% 
   set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(OL_SL_sum_table <- OL_SL_sum_tab1 %>% flextable() %>% colformat_num(j = 6, digits = 0) %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
rm(OL_SL_tab1, OL_SL_sum_tab1)
#
#Proportion of ages
(Age_prop_table <- as.data.frame(alkAgeDist(alk, lenA.n = rowSums(al_freq), len.n = len_n)) %>% 
    flextable() %>% colformat_num(j = 1, digits = 0) %>% autofit())
rm(alk, al_freq, len_n)
#
#Mean lengths at age
(AL_sum_table <- AgeLength_sum %>% flextable() %>% colformat_num(j = 1, digits = 0) %>% autofit())
rm(AgeLength_sum)
#
#Length at capture regression
(SLOL_tab <- SLOL_tab2 %>% flextable() %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(SLOL_sum_table <- SLOL_sum_tab2 %>% flextable() %>% colformat_num(j = c(6, 8), digits = 0) %>% 
    set_formatter("p-value" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
rm(SLOL_tab2, SLOL_sum_tab2)
#
#Original and Back calculated means
(OL_bc_table <- Final_comps %>% dplyr::select(Group, Age, n, meanSL, sdSL, - Striae) %>% 
    flextable() %>% merge_v(j = 1) %>% hline(i = 6) %>% autofit())
rm(Final_comps)
#
#Welchs test comps
(Welchs_table <- Welchs %>% flextable() %>% colformat_num(j = c(1:9), digits = 2) %>% autofit())
rm(Welchs)
#
#Model mean comps
(Model_comps_table <- Model_comps %>% 
  spread(Model, predSL) %>%
  flextable() %>% colformat_num(j = 1, digits = 0) %>% autofit())
rm(Model_comps)
#
#
#
####D. Marginal growth####
#
Marginal <- (Margins %>% mutate_at(c("Month"), as.factor) %>% 
  dplyr::select(Month, Total_Striae, Mature, MF_Final, MGR))[complete.cases(Margins %>% dplyr::select(Month, Total_Striae, Mature, MF_Final, MGR)),]
#
#Check normality of marginal growth rates
Marginal %>% ggplot(aes(MGR))+
  geom_histogram(aes(y = ..density..))+ #Not normal
  stat_function(fun = dnorm, args = list(mean = mean(Marginal$MGR), sd = sd(Marginal$MGR)), color = "red")
Marginal %>% ggplot(aes(sqrt(MGR)))+
  geom_histogram(aes(y = ..density..))+ #More normal
  stat_function(fun = dnorm, args = list(mean = mean(sqrt(Marginal$MGR)), sd = sd(sqrt(Marginal$MGR))), color = "red")
#
#Add column of sqrt transformed values
Marginal$sqMGR <- sqrt(Marginal$MGR)
#
#Visual assessment Sex and Maturity
Marginal %>% group_by(Month, Total_Striae, Mature, MF_Final) %>% 
  get_summary_stats(sqMGR, type = "mean_sd") %>%
  ggboxplot(x = "Month", y = "mean", color = "Mature")
#
#
###Permuation based anova
#
#Full model
set.seed(54321)
perm1 <- aovp(sqMGR ~ Month + Total_Striae + MF_Final + Mature, 
              data = Marginal, perm = "", nperm = 50000)
summary(perm1) #Maturity and striae significant
#
#Best model - keep Month for comparisons - use for both comparisons
set.seed(54321)
perm2 <- aovp(sqMGR ~ Month + Total_Striae + Mature, 
              data = Marginal, perm = "", nperm = 50000)
summary(perm2) #Maturity and striae, month
#
#
qqnorm(resid(perm2)) 
qqline(resid(perm2))
#
#Best model summary
(MGR_model <- summary(perm2))
(MGR_tidy <- tidy(perm2))
names(MGR_tidy) <- c("Factors", "df", "SS", "MS", "F", "Pr")
(MGR_summ <- glance(perm2)) 
names(MGR_summ) <- c("logLike", "AIC", "BIC", "deviance", "n", "R2")
#
###Pairwise comparisons
#
##Months
set.seed(54321)
(Month_pair <- Marginal %>% pairwise_t_test(sqMGR ~ Month, p.adjust.method = "bonferroni"))
(Month_tab <- dplyr::select(Month_pair, c("group1", "group2", "p", "p.adj")) %>%
  mutate(Comparison = paste(group1, group2, sep = "-")) %>%   #Add new column of grp v grp
  dplyr::select("Comparison", everything()) %>% dplyr::select(-c("group1", "group2")) %>% 
    rename(p.value = p, p.adjust = p.adj) %>% arrange(p.value))    #Move 'Comparison' to front and drop grp1 & grp2
Month_letters <- biostat::make_cld(Month_tab) %>% dplyr::select(-c("spaced_cld")) %>% 
  rename(Month = group, Letters = cld)
#All seasons are sig diff except 3 & 4 so seasonality should be removed from data
(Months <- merge(Marginal %>% group_by(Month) %>% rstatix::get_summary_stats(MGR, type = "mean_sd") %>% 
                   dplyr::select(-c("variable")) %>% transform(lower = mean-sd, upper = mean+sd), 
                Month_letters, by = "Month") %>% arrange(as.numeric(Month)))
#
#
##Maturity
set.seed(54321)
(Mat_pair <- Marginal %>% pairwise_t_test(sqMGR ~ Mature, p.adjust.method = "bonferroni"))
(Mat_tab <- dplyr::select(Mat_pair, c("group1", "group2", "p", "p.adj")) %>%
    mutate(Comparison = paste(group1, group2, sep = "-")) %>%   #Add new column of grp v grp
    dplyr::select("Comparison", everything()) %>% dplyr::select(-c("group1", "group2")) %>% 
    rename(p.value = p, p.adjust = p.adj) %>% arrange(p.value))    #Move 'Comparison' to front and drop grp1 & grp2
Mat_letters <- biostat::make_cld(Mat_tab) %>% dplyr::select(-c("spaced_cld")) %>% 
  rename(Mature = group, Letters = cld)
#All seasons are sig diff except 3 & 4 so seasonality should be removed from data
(Maturity <- merge(Marginal %>% group_by(Mature) %>% rstatix::get_summary_stats(MGR, type = "mean_sd") %>% 
                  dplyr::select(-c("variable")) %>% transform(lower = mean-sd, upper = mean+sd), 
                Mat_letters, by = "Mature"))
#
#
##Striae
set.seed(54321)
(Striae_pair <- Marginal %>% pairwise_t_test(sqMGR ~ Total_Striae, p.adjust.method = "bonferroni"))
(Striae_tab <- dplyr::select(Striae_pair, c("group1", "group2", "p", "p.adj")) %>%
    mutate(Comparison = paste(group1, group2, sep = "-")) %>%   #Add new column of grp v grp
    dplyr::select("Comparison", everything()) %>% dplyr::select(-c("group1", "group2")) %>% 
    rename(p.value = p, p.adjust = p.adj) %>% arrange(p.value))    #Move 'Comparison' to front and drop grp1 & grp2
Striae_letters <- biostat::make_cld(Striae_tab) %>% dplyr::select(-c("spaced_cld")) %>% 
  rename(Total_Striae = group, Letters = cld)
#All seasons are sig diff except 3 & 4 so seasonality should be removed from data
(Striae <- merge(Marginal %>% group_by(Total_Striae) %>% rstatix::get_summary_stats(MGR, type = "mean_sd") %>% 
                  dplyr::select(-c("variable")) %>% transform(lower = mean-sd, upper = mean+sd), 
                Striae_letters, by = "Total_Striae"))
#
Margins %>% group_by(Month, Total_Striae) %>% 
  get_summary_stats(MGR, type = "mean_sd") %>%
  ggplot(aes(Month, mean, color = as.factor(Total_Striae), group = as.factor(Total_Striae)))+
  geom_point(aes(shape = as.factor(Total_Striae)))+
  geom_line(aes(linetype = as.factor(Total_Striae)))+
  scale_x_continuous(expand = c(0,0.1), limits = c(1,12), breaks = seq(1, 12, 1))+
  scale_y_continuous("Average marginal growth rate (mm)", expand = c(0,0), limits = c(0, 7))+
  scale_shape_manual(name = "", values = c(16, 15, 17, 2))+
  scale_linetype_manual(name = "", values = c("solid", "solid", "dashed", "dotted"))+
  scale_color_manual(name = "", values = c("#000000", "#999999", "#000000", "#000000"))+
  basetheme + theme(legend.position = "none")
#
ggsave(path = "Output/Figures", 
       filename = paste("F9_Monthly Age MGR_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
Margins %>% group_by(Month, Mature) %>% 
  get_summary_stats(MGR, type = "mean_sd") %>%
  ggplot(aes(Month, mean, color = Mature, group = Mature))+
  geom_point(size = 2)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.15)+
  geom_line()+
  scale_x_continuous(expand = c(0,0.2), limits = c(0.9,12.1), breaks = seq(1, 12, 1))+
  scale_y_continuous("Average marginal growth rate (mm)", expand = c(0,0), limits = c(-0.5, 8.5))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  scale_color_manual(name = "", values = c("#000000", "#999999"))+
  basetheme + theme(legend.position = "none")
#
ggsave(path = "Output/Figures", 
       filename = paste("F10_Monthly Mature MGR_", format(Sys.Date(), "%Y_%m_%d"),".tiff", sep = ""), dpi = 1000)
#
#
rm(perm1, MGR_model, perm2, Month_pair, Month_letters, Mat_pair, Mat_letters, Striae_pair, Striae_letters)
#
#
#
####D. MGR summary####
#
#PermANOVA model
(MGR_table <- MGR_tidy %>% flextable() %>% colformat_num(j = 2, digits = 0) %>% 
   set_formatter("Pr" = function(x){formatC(x, format = "e", digits = 3)}) %>% autofit())
(MGR_summ_table <- MGR_summ %>% flextable() %>% autofit())
#
rm(MGR_tidy, MGR_summ)
#
#Pairwise T-tests
(Pairwise_tests <- full_join(Month_tab %>% mutate(Factor = "Month") %>% 
                               dplyr::select(Factor, everything()) %>% filter(p.adjust <= 0.1), 
                            Mat_tab %>% mutate(Factor = "Maturity")) %>% 
  full_join(Striae_tab %>% mutate(Factor = "Striae"))  %>% 
  flextable() %>% set_formatter(p.adjust = function(x){formatC(x, format = "e", digits = 3)},
                                p.value = function(x){formatC(x, format = "e", digits = 3)}) %>%
  merge_v(j = 1) %>% hline(i = c(1, 2)) %>% autofit())
#
rm(Month_tab, Mat_tab, Striae_tab)
#
#Means and letters
(Months_tab <- Months %>% flextable() %>% colformat_num(j = 2, digits = 0) %>% autofit())
(Mature_tab <- Maturity %>% flextable() %>% colformat_num(j = 2, digits = 0) %>% autofit())
(Striae_tab <- Striae %>% flextable() %>% colformat_num(j = 2, digits = 0) %>% autofit())
#
rm(Months, Maturity, Striae)
#
#
#
####Write tables to Word####
#
save_as_docx("Ave morphometrics by sex" = Sample_tab,
             "Ave morphometrics by Striae" = Sam_striae_tab,
             "Percent reader agreement" = Agree,
             "Cohens Kappa for Reader Agreement" = Cohens_Kappa,
             "Shell length model" = S_SL_table, 
             "Shell length R2" = S_SL_sum_table,
             "Shell width model" = S_SW_table,
             "Shell width R2" = S_SW_sum_table,
             "Shell weight model" = S_WW_table,
             "Shell weight R2" = S_WW_sum_table,
             "Regression average age values (raw data)" = Regression_sum_tab,
             "Model predictions" = Model_preds,
             "Model evaluations" = Model_eval,
             "Operculum Length v. Shell Length model" = OL_SL_tab1_out,
             "OL v SL R2" = OL_SL_sum_table,
             "Proportion of ages" = Age_prop_table,
             "Mean length at age" = AL_sum_table,
             "OL Capture v SL Capture model" = SLOL_tab,
             "OL Cap v SL Cap R2" = SLOL_sum_table,
             "Observed & Back calculated means" = OL_bc_table,
             "Welch's results" = Welchs_table,
             "Mean length at age, obs and bc" = Model_comps_table,
             "PermANOVA model" = MGR_table,
             "PermANOVA R2" = MGR_summ_table,
             "Pairwise T-tests" = Pairwise_tests,
             "Monthly MGR" = Months_tab,
             "Mature MGR" = Mature_tab,
             "Age MGR" = Striae_tab,
             path = "Output/BT_operculum_tables.docx", 
             pr_section = sect_properties)
#
#
#