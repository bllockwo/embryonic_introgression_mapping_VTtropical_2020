library(ggplot2)
library(readr)
library(drc)
library(Rmisc)
library(cowplot)

# Load survival data
survival_data <- read_delim("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/20.WEATHER_PLOT_plus_Lockwood_2018_PHENOS/Raw_survival_data_VT_trop_only.txt")

# subset for just SK and VTECK_8
sk_vt8 <- survival_data[which(survival_data$genotype == "SK" | survival_data$genotype == "VTECK_8"),]

# set up for plot
newdata <- expand.grid(temperature=seq(20, 45, length=250), genotype=c("SK","VTECK_8"))

# run logistic model
drc <- drm(hatched/eggs~temperature, genotype, weights = eggs, data = sk_vt8, 
                fct = LL.3(names = c("slope", "upper limit", "LT50")), type="binomial")

# predict line fits (Don't worry about warmings)
pm.VT <- predict(drc, newdata=newdata, interval="confidence", type="response")

# new data with predictions
newdata$p <- pm.VT[,1]
newdata$pmin <- pm.VT[,2]
newdata$pmax <- pm.VT[,3]

# survival plot

ggplot(data = sk_vt8, aes(x = temperature, y = survival)) +
  geom_ribbon(data=newdata, aes(x=temperature, y=p, ymin=pmin, ymax=pmax, fill=factor(genotype, levels = c("SK","VTECK_8"))), alpha=0.5) +
  geom_line(data=newdata, size = 0.7, aes(x=temperature, y=p, color = factor(genotype, levels = c("SK","VTECK_8")), linetype = factor(genotype, levels = c("SK","VTECK_8")))) +
  geom_point(aes(color = factor(genotype, levels = c("SK","VTECK_8")), shape = factor(genotype, levels = c("SK","VTECK_8"))), alpha=0.2, size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
  xlab("Temperature (°C)") + 
  ylab("Proportion hatched") +
  scale_color_manual(values = c("firebrick", "steelblue")) +
  scale_fill_manual(values = c("firebrick", "steelblue")) +
  xlim(25,45) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.title = element_blank(), legend.position = c(0.85,0.88)) + 
  guides(colour = guide_legend(override.aes = list(alpha = .6)))

# Try a plot of survival (hatching success) at 36°C only with a boxplot
temp_36 <- sk_vt8[which(sk_vt8$temperature == 36),]

ggplot(temp_36, aes(x = genotype, y = survival)) +
  geom_boxplot()

# Create a summary data frame 
sk_vt8_summary <- Rmisc::summarySE(sk_vt8, measurevar="survival", groupvars=c("temperature","genotype"))

# Plot a line interaction plot at 36°
ggplot(sk_vt8_summary[which(sk_vt8_summary$temperature == 36),], aes(x=as.factor(genotype), y=survival)) +
  geom_point(size=3, alpha=0.75) +
  geom_errorbar(aes(ymin=survival-se, ymax=survival+se), width=.15) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  ylim(0,1) +
  theme_half_open() +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6), axis.text.x = element_blank(), axis.title.x = element_blank())

#Plot a line interaction plot at 34°
ggplot(sk_vt8_summary[which(sk_vt8_summary$temperature == 34),], aes(x=as.factor(genotype), y=survival)) +
  geom_point(size=3, alpha=0.75) +
  geom_errorbar(aes(ymin=survival-se, ymax=survival+se), width=.15) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  ylim(0,1) +
  theme_half_open() +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6), axis.text.x = element_blank(), axis.title.x = element_blank())

# Plot LT50 embryos
lt50 <- c(34.77, 35.628)
lt50_se <- c(0.195, 0.131)
genotype <- c("VT8", "SK")

sk_vt8_lt50 <- data.frame(genotype, lt50, lt50_se)

ggplot(sk_vt8_lt50, aes(x=as.factor(genotype), y=lt50)) +
  geom_point(size=4, alpha=0.75) +
  geom_errorbar(aes(ymin=lt50-lt50_se, ymax=lt50+lt50_se), width=.15) +
  ylim(34, 36) +
  ylab(expression("LT"["50"]~"("*degree*"C)")) +
  theme_half_open() +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6), axis.title.x = element_blank())


# Plot LT50 adults
lt50_adult <- c(39.535, 39.397)
lt50_se_adult <- c(0.106, 0.058)
genotype <- c("VT8", "SK")

sk_vt8_lt50_adult <- data.frame(genotype, lt50_adult, lt50_se_adult)

ggplot(sk_vt8_lt50_adult, aes(x=as.factor(genotype), y=lt50_adult)) +
  geom_point(size=4, alpha=0.75) +
  geom_errorbar(aes(ymin=lt50_adult-lt50_se_adult, ymax=lt50_adult+lt50_se_adult), width=.15) +
  ylim(37, 40) +
  ylab(expression("LT"["50"]~"("*degree*"C)")) +
  theme_half_open() +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6), axis.title.x = element_blank())

# Plot embryo and adult LT50 together
lt50 <- c(34.77, 35.628, 39.535, 39.397)
lt50_se <- c(0.195, 0.131, 0.106, 0.058)
genotype <- c("VT8", "SK", "VT8", "SK")
stage <- c("embryo", "embryo", "adult", "adult")
sk_vt8_lt50 <- data.frame(genotype, stage, lt50, lt50_se)

ggplot(sk_vt8_lt50, aes(x=as.factor(genotype), y=lt50, color = stage, shape = stage)) +
  geom_point(size=3, alpha=0.75) +
  geom_errorbar(aes(ymin=lt50-lt50_se, ymax=lt50+lt50_se), width=.15) +
  scale_shape_manual(values = c(17,19))+
  scale_color_manual(values = c("#009E73", "#56B4E9")) +
  scale_y_continuous(breaks = seq(34, 42, by = 1), limits = c(34, 40)) +
  ylab(expression("LT"["50"]~"("*degree*"C)")) +
  theme_minimal_hgrid() +
  theme(panel.grid.minor.y = element_blank(), plot.margin = margin(0,6,0,6), axis.title.x = element_blank(), legend.title = element_blank())

