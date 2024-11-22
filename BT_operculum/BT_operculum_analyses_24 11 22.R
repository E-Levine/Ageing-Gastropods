##Operculum striae as an aging technique in Banded Tulips
#
##Data analyses
#
#
#
#Use Alt+O to collapse all sections, Alt+Shift+O to expand all sections#
#Use Alt+L to collapse on section, Alt+Shift_L to expand section#
#
#
#Load require packages (install as necessary)
if (!require("pacman")) {install.packages("pacman")}
pacman::p_load(plyr, tidyverse, rstatix,
               vcd, #Cohen's Kappa - reader agreement
               FSA, AICcmodavg, propagate, lmPerm, #Length bins, AIC comps, nls predictions
               scales, broom, psych, #Re-scaling MG, cleaning summary tables, pairwise letters
               ggpubr, lemon, extrafont, #Figure arrangement
               flextable, janitor, officer) #Tables, output to Word
if (!require("remotes")) {install.packages("remotes")}
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
####Figures and Tools####
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
        axis.title = element_text(size = 12, face = "bold", color = "black", family = "Arial"), 
        axis.text = element_text(size = 11, color = "black", family = "serif", margin = unit(c(0.5, 0.5, 0, 0.5), "cm")), 
        axis.ticks.length = unit(-0.15, "cm"))
#
##No X Title (Categories (i.e. Month))
XCate <- theme(axis.title.x = element_blank(),
               axis.text.x = element_text(color = "black", size = 12, family = "extrafont",
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
#Contingency table - agreement
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
             size = 2.5, color = "black", family = "Arial")+
    annotate("text", x = 0.9, y = 103, label = paste("R^2 ==", round(S_SL_sum$r.squared, 2)), parse = TRUE,
             size = 2.5, color = "black", family = "Arial"))
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
             size = 2.5, color = "black", family = "Arial")+
    annotate("text", x = 0.9, y = 47.5, label = paste("R^2 ==", round(S_SW_sum$r.squared, 2)), parse = TRUE,
             size = 2.5, color = "black", family = "Arial"))
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
             size = 2.5, color = "black", family = "Arial")+
    annotate("text", x = 1.1, y = 103, label = paste("R^2 ==", round(S_WW_sum$r.squared, 2)), parse = TRUE,
             size = 2.5, color = "black", family = "Arial"))
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
             size = 2.5, color = "black", family = "Arial")+
    annotate("text", x = 1.1, y = 47.5, label = paste("R^2 ==", round(S_OL_sum$r.squared, 2)), parse = TRUE,
             size = 2.5, color = "black", family = "Arial"))
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
#
#
#