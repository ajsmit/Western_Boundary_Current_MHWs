# create plot arrangements for figures
# load plots from these files first:
# /Users/ajsmit/Dropbox/R/WBCs/setup/Fig00_Aviso+_MKE_EKE.R
# /Users/ajsmit/Dropbox/R/WBCs/eddies/Fig018_Eddy_trajectories.R
# /Users/ajsmit/Dropbox/R/WBCs/setup/Fig06_OISST_MeanInt.R

# Figure 1
Figure_1 <- ggarrange(AC.Fig00a + labs(title = "Mean kinetic energy"),
                       AC.Fig00b + labs(title = "Eddy kinetic energy"),
                       AC.Fig018 + labs(title = "Eddy trajectories"),
                       AC.Fig06 + labs(title = "Heatwave intensity"),
                       ncol = 2, nrow = 2, labels = list("a", "b", "c", "d"))
ggplot2::ggsave("publ_plots/Figure_1_AC.jpg",
                width = 7.0 * (1/3), height = 5.2 * (1/3), scale = 3.7)


# Figure 3
Figure_1 <- ggarrange(AC.fig05,
                      AC.fig09,
                      BC.fig05,
                      BC.fig09
                      EAC.fig05,
                      EAC.fig09
                      GS.fig05,
                      GS.fig09
                      KC.fig05,
                      KC.fig09,
                      ncol = 2, nrow = 5, labels = "auto")
ggplot2::ggsave("publ_plots/Figure_3.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)

