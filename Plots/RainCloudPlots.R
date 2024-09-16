install.packages(c("agridat", "ggplot2", "ghibli", "ggdist")

#Line below selects the dataset "Birth weight and weaning weight of Dorper x Red Maasi lambs"
data <- agridat::ilri.sheep

data <- data |> mutate(birthwt=as.numeric(birthwt),
                       weanwt=as.numeric(weanwt),
                       weanage=as.numeric(weanage),
#Line below creates the variable "weight gain from birth to weaning" displayed in grams per day
weight_gain_gram=as.numeric(round((((weanwt-birthwt)/weanage)*1000),2),na.rm=T))

data <- subset(data, select=c(lamb,gen,weight_gain_gram)) #selecting the variables of interest for this exercise

head(data,10) #shows the first 10 rows only

#lamb gen weight_gain_gram
#1   627  DD           108.80
#2   629  DD           138.39
#3   635  DD           111.93

ggplot(data, aes(x = gen, y = weight_gain_gram, fill=gen)) +
# Line below sets the Studio Ghibli color pallete, for the sake of nostalgia =)
  scale_fill_ghibli_d("SpiritedMedium", direction = -1) +
  geom_boxplot(width = 0.1) +
  xlab('Lamb genotype') +
  ylab('Weight gain, in g/d') +
  ggtitle("Weight gain from birth to weaning in 4 lamb genotypes") +
  theme_classic(base_size=18, base_family="serif")+
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=0, hjust=.5, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position="none")+
  scale_y_continuous(breaks = seq(0, 180, by=20), limits=c(0,180), expand = c(0, 0)) +
# Line below adds dot plots from {ggdist} package 
  stat_dots(side = "left", justification = 1.12, binwidth = 1.9) +
# Line below adds half-violin from {ggdist} package
  stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA)
