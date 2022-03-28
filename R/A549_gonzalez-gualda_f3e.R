#create figure 3e for González-Gualda et al.

library(tidyverse)

l3 <- read_csv("../data/a549_low_serum_level_3.csv")

df <- l3 %>%
  filter(ECMp == "COL1_NA_NA")

p_full <-ggplot(df, aes(x=reorder(Ligand,Nuclei_PA_Gated_EdUPositiveProportion, FUN=median),
                   y=Nuclei_PA_Gated_EdUPositiveProportion,
                   color = Ligand))+
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,.5))+
  guides(color = "none") +
  labs(x = "Ligand",
       y="EdU High Proportion",
       title=paste("EdU High Proportion for",unique(df$CellLine), "cells by ECM Protein"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(1.2)),
        axis.title.x = element_text(size=rel(1.5)),
        plot.title = element_text(size = rel(1)),
        legend.text=element_text(size = rel(1)),
        legend.title=element_text(size = rel(1)),
        strip.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

p_full

p_paper <-ggplot(df, aes(x=reorder(Ligand,Nuclei_PA_Gated_EdUPositiveProportion, FUN=median),
                   y=Nuclei_PA_Gated_EdUPositiveProportion,
                   color = Ligand))+
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0,.5))+
  guides(color = "none") +
  labs(x = "Ligand",
       y="EdU High Proportion")+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size=rel(1.5)),
        plot.title = element_text(size = rel(1)),
        legend.text=element_text(size = rel(1)),
        legend.title=element_text(size = rel(1)),
        strip.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

p_paper


pdf("../plots/Gonzalez-gualda_fig3e_annotated.pdf", useDingbats = FALSE, width = 8, height = 6)
print(p_full)
res <- dev.off()

pdf("../plots/Gonzalez-gualda_fig3e.pdf", useDingbats = FALSE, width = 8, height = 2)
print(p_paper)
res <- dev.off()
