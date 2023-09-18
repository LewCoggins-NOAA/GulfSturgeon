#####
###   Plot PVA results - 2023-08-12
#####

library(tidyverse)
library(gridExtra)
library(grid)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(rprojroot)

# set wd to project location
setwd(find_rstudio_root_file())

res.w <- read.csv("PVA_results_20230812.csv")

res <- gather(res.w, year, prob, X50.year.Probability:X200.year.Probability)

res$year[res$year == "X50.year.Probability"]  <- 50
res$year[res$year == "X100.year.Probability"] <- 100
res$year[res$year == "X200.year.Probability"] <- 200

res$year <- as.numeric(res$year)


c.100_plot <- ggplot(filter(res, Threat.Definition=='chronic' & Vulnerable.Abundance==100),
                              aes(x = year, y = prob, group=Adult.Mortality, colour=factor(Adult.Mortality))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "100")) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 1),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(c.100_plot)

c.500_plot <- ggplot(filter(res, Threat.Definition=='chronic' & Vulnerable.Abundance==500),
                              aes(x = year, y = prob, group=Adult.Mortality, colour=factor(Adult.Mortality))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "500")) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 0),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(c.500_plot)

c.1k_plot <- ggplot(filter(res, Threat.Definition=='chronic' & Vulnerable.Abundance==1000),
                              aes(x = year, y = prob, group=Adult.Mortality, colour=factor(Adult.Mortality))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size = 5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "1,000")) +
  theme(axis.text.x = element_text(face = "bold", size = 16), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 1),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(c.1k_plot)

c.10k_plot <- ggplot(filter(res, Threat.Definition=='chronic' & Vulnerable.Abundance==10000),
                              aes(x = year, y = prob, group=Adult.Mortality, colour=factor(Adult.Mortality))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5, show.legend  = F) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(name="Adult\nMortality", values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "10,000")) +
  theme(axis.text.x = element_text(face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 0),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = c(.915,0.583)); print(c.10k_plot)


pva.c.fig <- ggarrange(c.100_plot, NULL, c.500_plot, 
                       c.1k_plot,  NULL, c.10k_plot, 
                         ncol = 3,
                         nrow = 2,
                         widths = c(1,0.05,1,
                                    1,0.05,1)) +
  theme(plot.margin = margin(0.05,0.1,0.1,0.5, "cm"))

setwd(find_rstudio_root_file())

tiff("results/chronic.tiff", width = 18, height = 10, units = 'in', res = 300)
annotate_figure(pva.c.fig, left = textGrob("Extirpation Probability", rot = 90, vjust = .7, gp = gpar(cex = 2.1)),
                top = textGrob("PVA Scenarios: Chronic Mortality", vjust = .35, gp = gpar(cex = 2.4)))
dev.off()

#-----------------------------------------------------------------------------------------------------------------

e.100_plot <- ggplot(filter(res, Threat.Definition=='episodic' & Vulnerable.Abundance==100),
                     aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(50,25,10)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "100")) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 1),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(e.100_plot)

e.500_plot <- ggplot(filter(res, Threat.Definition=='episodic' & Vulnerable.Abundance==500),
                     aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(50,25,10)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "500")) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 0),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(e.500_plot)

e.1k_plot <- ggplot(filter(res, Threat.Definition=='episodic' & Vulnerable.Abundance==1000),
                    aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(50,25,10)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size = 5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "1,000")) +
  theme(axis.text.x = element_text(face = "bold", size = 16), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 1),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(e.1k_plot)

e.10k_plot <- ggplot(filter(res, Threat.Definition=='episodic' & Vulnerable.Abundance==10000),
                     aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(50,25,10)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5, show.legend  = F) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(name="Mean\nFrequency", labels=c('1/50 years', '1/25 years', '1/10 years'), values=c("#1b7837", "#e08214", "#b2182b"))+
  ggtitle(bquote(N[0]~ "=" ~ "10,000")) +
  theme(axis.text.x = element_text(face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 0),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = c(.88,0.44)); print(e.10k_plot)


pva.e.fig <- ggarrange(e.100_plot, NULL, e.500_plot, 
                       e.1k_plot,  NULL, e.10k_plot, 
                       ncol = 3,
                       nrow = 2,
                       widths = c(1,0.05,1,
                                  1,0.05,1)) +
  theme(plot.margin = margin(0.05,0.1,0.1,0.5, "cm"))

setwd(find_rstudio_root_file())

tiff("results/episodic.tiff", width = 18, height = 10, units = 'in', res = 300)
annotate_figure(pva.e.fig, left = textGrob("Extirpation Probability", rot = 90, vjust = .7, gp = gpar(cex = 2.1)),
                top = textGrob("PVA Scenarios: Episodic Mortality \u2013 M = 0.11", vjust = .35, gp = gpar(cex = 2.4)))
dev.off()


#-----------------------------------------------------------------------------------------------------------------

r.100_plot <- ggplot(filter(res, Threat.Definition=='rec_fail' & Vulnerable.Abundance==100),
                     aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(10,5)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214"))+
  ggtitle(bquote(N[0]~ "=" ~ "100")) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 1),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(r.100_plot)

r.500_plot <- ggplot(filter(res, Threat.Definition=='rec_fail' & Vulnerable.Abundance==500),
                     aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(10,5)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214"))+
  ggtitle(bquote(N[0]~ "=" ~ "500")) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 0),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(r.500_plot)

r.1k_plot <- ggplot(filter(res, Threat.Definition=='rec_fail' & Vulnerable.Abundance==1000),
                    aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(10,5)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size = 5) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(values=c("#1b7837", "#e08214"))+
  ggtitle(bquote(N[0]~ "=" ~ "1,000")) +
  theme(axis.text.x = element_text(face = "bold", size = 16), 
        axis.text.y = element_text(face = "bold", size = 16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 1),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = "none"); print(r.1k_plot)

r.10k_plot <- ggplot(filter(res, Threat.Definition=='rec_fail' & Vulnerable.Abundance==10000),
                     aes(x = year, y = prob, group=Average.Frequency, colour=factor(Average.Frequency, levels = c(10,5)))) +
  geom_point(size = 4) +
  geom_line(size = 1) + 
  geom_text(aes(label = round(prob, 2), y = ifelse(prob==0.113,prob+.035,prob+.008), fontface = "bold"), vjust = -0.9, size=5, show.legend  = F) +
  coord_cartesian(ylim = c(0,1),
                  clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.002)) +
  scale_color_manual(name="Mean\nFrequency", labels=c('1/10 years', '1/5 years'), values=c("#1b7837", "#e08214"))+
  ggtitle(bquote(N[0]~ "=" ~ "10,000")) +
  theme(axis.text.x = element_text(face = "bold", size = 16), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black",size = 1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3.5, size = 26),
        plot.margin = margin(t = 20,
                             r = 0,
                             b = 3.5,
                             l = 0),
        panel.background = element_rect(fill = "white", colour = "white", 
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(colour = "black", size = 0.5),
        legend.position = c(.88,0.44)); print(r.10k_plot)


pva.r.fig <- ggarrange(r.100_plot, NULL, r.500_plot, 
                       r.1k_plot,  NULL, r.10k_plot, 
                       ncol = 3,
                       nrow = 2,
                       widths = c(1,0.05,1,
                                  1,0.05,1)) +
  theme(plot.margin = margin(0.05,0.1,0.1,0.5, "cm"))

setwd(find_rstudio_root_file())

tiff("results/rec_fail.tiff", width = 18, height = 10, units = 'in', res = 300)
annotate_figure(pva.r.fig, left = textGrob("Extirpation Probability", rot = 90, vjust = .7, gp = gpar(cex = 2.1)),
                top = textGrob("PVA Scenarios: Recruitment Failure \u2013 M = 0.11", vjust = .35, gp = gpar(cex = 2.4)))
dev.off()
