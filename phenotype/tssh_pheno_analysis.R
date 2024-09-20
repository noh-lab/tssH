# Analysis of phenotypic data for "Symbiotic T6SS affects horizontal transmission of Paraburkholderia bonniea among Dictyostelium discoideum amoeba hosts" by Chen et al.
# Please contact Suegene Noh (suegene.noh@colby.edu) with any issues

library(car)
library(lme4)
library(parameters)
library(effectsize)
library(ggplot2)
library(emmeans)

manuscript_theme = ggplot2::theme_bw() + ggplot2::theme(text=element_text(size=11))


############ ############# HOST FITNESS ANALYSIS ############# #############
tssh_spore <- read.table("tssH.pheno_fitness.tsv", h=T, row.names=NULL)
head(tssh_spore)
str(tssh_spore)
colnames(tssh_spore) <- c("date","host","background","MOI","infection_prevalence","host_fitness","treatment")
tssh_spore$MOI <- factor(tssh_spore$MOI, levels=c("0%_0","0.05%_0.3","0.25%_1.5","1.25%_7.5"))
tssh_spore$date <- as.factor(tssh_spore$date)
tssh_spore$treatment <- as.factor(gsub(tssh_spore$treatment, pattern="wt", replacement="wildtype"))
tssh_spore$treatment <- relevel(tssh_spore$treatment, "wildtype")
tssh_spore$host <- as.factor(tssh_spore$host)
tssh_spore$background <- as.factor(tssh_spore$background)

str(tssh_spore)

# overall patterns in the two variables measured during the experiment
ggplot(subset(tssh_spore, MOI!="control"), aes(x=MOI, y=infection_prevalence, fill=treatment)) + geom_boxplot()
ggplot(subset(tssh_spore, MOI!="control"), aes(x=MOI, y=host_fitness, fill=treatment)) + geom_boxplot()


# what affects host fitness, including symbiont treatment, symbiont background
fitness.1 <- lmer(data=tssh_spore, host_fitness ~ infection_prevalence*background*treatment + (1|date) + (1|host))
Anova(fitness.1) 

fitness.2 <- lmer(data=tssh_spore, host_fitness ~ infection_prevalence*background*treatment - infection_prevalence:background:treatment + (1|date) + (1|host))
Anova(fitness.2) 

anova(fitness.1, fitness.2) # not significant; infection_prevalence:treatment is not significant (mutant does not have a different slope, but bb433 and bb859 do)

hist(residuals(fitness.2))
plot(fitted(fitness.2), residuals(fitness.2))    

#model_parameters(fitness.2, effects="fixed")
#model_parameters(anova(fitness.2))

epsilon_squared(fitness.2)

emtrends(fitness.2, pairwise ~ treatment|background, var = "infection_prevalence", adjust = "tukey")

fitness.3 <- lmer(data=tssh_spore, infection_prevalence ~ background*treatment + (1|date) + (1|host))
Anova(fitness.3) 

emmeans(fitness.3, pairwise ~ treatment|background, var = "infection_prevalence", adjust = "tukey")


# visualize effect of symbiont, and lack of effect of host 
# full factorial design
ggplot(tssh_spore, aes(x=infection_prevalence, y=host_fitness, fill=symbiont, shape=treatment)) +geom_point(size=2) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE) + facet_grid(host~background) + ylim(60,110)

# grouped by background 
ggplot(tssh_spore, aes(x=infection_prevalence, y=host_fitness, fill=symbiont, shape=treatment)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(.~background) + ylab("Host fitness") + xlab("Infection prevalence") + manuscript_theme


############ ############# SYMBIONT TRANSMISSION ANALYSIS ############# #############
tssh_horiz <- read.table("tssH.pheno_transmission.tsv", h=T)
head(tssh_horiz)
str(tssh_horiz)
colnames(tssh_horiz) <- c("date","host","background","MOI","infection_prevalence","symbiont_transmission","treatment")
tssh_horiz$MOI <- factor(tssh_horiz$MOI, levels=c("0%_0","0.05%_0.3","0.25%_1.5","1.25%_7.5"))
tssh_horiz$date <- as.factor(tssh_horiz$date)
tssh_horiz$treatment <- as.factor(gsub(tssh_horiz$treatment, pattern="wt", replacement="wildtype"))
tssh_horiz$treatment <- relevel(tssh_horiz$treatment, "wildtype")
tssh_horiz$host <- as.factor(tssh_horiz$host)
tssh_horiz$background <- as.factor(tssh_horiz$background)

str(tssh_horiz)


# overall patterns in the two variables measured during the experiment
ggplot(tssh_horiz, aes(x=background, y=infection_prevalence, fill=treatment)) +geom_boxplot() 
ggplot(tssh_horiz, aes(x=background, y=symbiont_transmission, fill=treatment)) +geom_boxplot() 


# what affects transmission
transmit.1 <- lmer(data=tssh_horiz, symbiont_transmission ~ infection_prevalence*background*treatment + (1|date) + (1|host)) 
Anova(transmit.1) 

transmit.2 <- lmer(data=tssh_horiz, symbiont_transmission ~ infection_prevalence*background*treatment - background:treatment + (1|date) + (1|host))
Anova(transmit.2, type=2) 

anova(transmit.1, transmit.2) # significant

hist(residuals(transmit.1))
plot(fitted(transmit.1), residuals(transmit.1))    

#model_parameters(transmit.1, effects="fixed")
#model_parameters(anova(transmit.1))

epsilon_squared(transmit.1)

emtrends(transmit.1, pairwise ~ treatment|background, var = "infection_prevalence", adjust = "tukey")

transmit.3 <- lmer(data=tssh_horiz, infection_prevalence ~ background*treatment + (1|date) + (1|host))
Anova(fitness.3) 

emmeans(transmit.3, pairwise ~ treatment|background, var = "infection_prevalence", adjust = "tukey")


# visualize effect of symbiont, and lack of effect of host type
# full factorial design
ggplot(tssh_horiz, aes(x=infection_prevalence, y=symbiont_transmission, fill=background, shape=treatment)) +geom_point(size=2) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=treatment, color=background)) + facet_grid(host~background) 

# grouped by background 
ggplot(tssh_horiz, aes(x=infection_prevalence, y=symbiont_transmission, fill=background, shape=treatment)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=treatment, color=background)) + facet_grid(.~background) + ylab("Horizontal transmission") + xlab("Infection prevalence") + manuscript_theme


############ ############# MUTANT RESCUE ANALYSIS ############# #############
tssh_rescue <- read.table("tssH.pheno_rescue.tsv", h=T)
head(tssh_rescue)
str(tssh_rescue)
colnames(tssh_rescue) <- c("date","host","background","MOI","infection_prevalence","symbiont_transmission","treatment")
tssh_rescue$MOI <- factor(tssh_rescue$MOI, levels=c("0%_0","0.05%_0.3","0.25%_1.5","1.25%_7.5"))
tssh_rescue$date <- as.factor(tssh_rescue$date)
tssh_rescue$treatment <- as.factor(gsub(tssh_rescue$treatment, pattern="wt", replacement="wildtype"))
tssh_rescue$treatment <- relevel(tssh_rescue$treatment, "wildtype")
tssh_rescue$host <- as.factor(tssh_rescue$host)
tssh_rescue$background <- as.factor(tssh_rescue$background)

ggplot(tssh_rescue, aes(x=infection_prevalence, y=symbiont_transmission, shape=background, fill=treatment)) + geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23, 24)) + geom_smooth(method="lm", se=FALSE, aes(group=treatment, color=treatment)) + facet_grid(.~background, scales="free") + theme_bw() + xlab("Infection prevalence") + ylab("Horizontal transmission") + theme(plot.title = element_text(size=11)) 


# what affects transmission
temp <- subset(tssh_rescue, treatment!="kanR")
tssh_rescue <- temp
rescue.1 <- lmer(data=tssh_rescue, symbiont_transmission ~ infection_prevalence*background*treatment + (1|date) ) 
Anova(rescue.1) 

rescue.2 <- lmer(data=tssh_rescue, symbiont_transmission ~ infection_prevalence*background*treatment - infection_prevalence:background:treatment + (1|date) ) 
Anova(rescue.2) 

anova(rescue.1, rescue.2) # not significant

hist(residuals(rescue.2))
plot(fitted(rescue.2), residuals(rescue.2))    

#model_parameters(rescue.2, effects="fixed")
#model_parameters(anova(rescue.2))

epsilon_squared(rescue.2)

emtrends(rescue.2, pairwise ~ treatment|background, var = "infection_prevalence", adjust = "tukey")

rescue.3 <- lmer(data=tssh_rescue, infection_prevalence ~ background*treatment + (1|date))
Anova(rescue.3) 

emmeans(rescue.3, pairwise ~ treatment|background, var = "infection_prevalence", adjust = "tukey")

postscript(file=paste("tssh_fitness",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=3)
ggplot(tssh_spore, aes(x=infection_prevalence, y=host_fitness, fill=treatment, shape=treatment)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + scale_fill_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + scale_color_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + geom_smooth(method = "lm", se=FALSE, aes(group=treatment, color=treatment)) + facet_grid(.~background) + ylab("Host fitness") + xlab("Infection prevalence") + manuscript_theme
dev.off()
postscript(file=paste("tssh_transmission",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=3)
ggplot(tssh_horiz, aes(x=infection_prevalence, y=symbiont_transmission, fill=treatment, shape=treatment)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + scale_fill_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + scale_color_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + geom_smooth(method = "lm", se=FALSE, aes(group=treatment, color=treatment)) + facet_grid(.~background) + ylab("Symbiont transmission") + xlab("Infection prevalence") + manuscript_theme
dev.off()
postscript(file=paste("tssh_rescue",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=3)
ggplot(tssh_rescue, aes(x=infection_prevalence, y=symbiont_transmission, shape=treatment, fill=treatment)) + geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23, 24)) + scale_fill_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + scale_color_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + geom_smooth(method="lm", se=FALSE, aes(group=treatment, color=treatment)) + facet_grid(.~background) + theme_bw() + xlab("Infection prevalence") + ylab("Symbiont transmission") + manuscript_theme
dev.off()


############ ############# difference between two experimental rounds ############# #############
library(gridExtra)
library(ggbeeswarm)

r1 <- ggplot(subset(tssh_spore, MOI!="control"), aes(x=background, y=infection_prevalence, fill=treatment, shape=treatment)) + geom_beeswarm(priority="density", cex=2) + scale_shape_manual(values=c(21, 22, 23)) + scale_fill_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + theme(legend.position="none") + ylab("Infection prevalence") + xlab("Symbiont") + ggtitle("(a) Host fitness exp.") + theme(plot.title = element_text(size=11)) + ylim(0,100) + facet_grid(.~background) + manuscript_theme

r2 <- ggplot(tssh_horiz, aes(x=background, y=infection_prevalence, fill=treatment, shape=treatment)) + geom_beeswarm(priority = "density", cex=2) + scale_shape_manual(values=c(21, 22, 23)) + scale_fill_manual(values=c("#31688EFF", "#35B779FF", "#F1605DFF")) + theme(legend.position="none") + ylab("") + xlab("Symbiont") + ggtitle("(b) Symbiont transmission exp.") + theme(plot.title = element_text(size=11)) + ylim(0,100) + facet_grid(.~background) + manuscript_theme

grid.arrange(r1, r2, ncol=2)

postscript(file=paste("tssh_compare",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=3)
grid.arrange(r1, r2, ncol=2)
dev.off()


r1.data <- subset(spore, MOI!="control")[,c(1,2,3,4,6,8)]
r2.data <- horiz[c(1,2,3,4,5,7)]
r1.data$stage <- c("one")
r2.data$stage <- c("two")
social <- rbind(r1.data, r2.data)
rm(r1.data, r2.data)

round <- lmer(flow_infected ~ stage*symbiont + (1|host), data=social)
Anova(round) # all significant
hist(residuals(round))
plot(fitted(round), residuals(round))    

epsilon_squared(round)

emmeans(round, pairwise ~ symbiont | stage, adjust = "tukey")

