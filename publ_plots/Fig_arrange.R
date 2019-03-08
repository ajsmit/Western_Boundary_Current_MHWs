# create plot arrangements for figures
# load plots from these files first:
# /Users/ajsmit/Dropbox/R/WBCs/setup/Fig00_Aviso+_MKE_EKE.R
# /Users/ajsmit/Dropbox/R/WBCs/eddies/Fig018_Eddy_trajectories.R
# /Users/ajsmit/Dropbox/R/WBCs/setup/Fig06_OISST_MeanInt.R

# Figure 1 (a-d). MKE, EKE, eddy trajectories, and MHW intensity of the Agulhas Current.
Fig.1 <- ggarrange(AC.Fig00a + labs(title = "Mean kinetic energy"),
                   AC.Fig00b + labs(title = "Eddy kinetic energy"),
                   AC.Fig018 + labs(title = "Eddy trajectories"),
                   AC.Fig06 + labs(title = "Heatwave intensity"),
                   ncol = 2, nrow = 2, labels = list("a", "b", "c", "d"))
ggplot2::ggsave("publ_plots/Figure_1.jpg",
                width = 7.0 * (1/3), height = 5.2 * (1/3), scale = 3.7)


# •	Suppl. Fig. 1 (a-t). Full set of panels corresponding to Fig. 1.
Supp.1 <- ggarrange(AC.fig05,
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
ggplot2::ggsave("publ_plots/Suppl.Fig.1.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)

# •	Suppl. Fig. 2 (a-d). Full set of panels of SST mean trend.
Supp.2 <- ggarrange(AC.Fig04,
                    BC.Fig04,
                    EAC.Fig04,
                    GS.Fig04,
                    KC.Fig04,
                    ncol = 1, nrow = 5, labels = "auto")
ggplot2::ggsave("publ_plots/Suppl.Fig.2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

# Figure 2 (a-e). Masked regions of the Agulhas Current that correspond to MKE, EKE, eddy tracks, and MHW intensity ≥ 90th percentile.
Fig.2 <- ggarrange(AC.poly.plt,
                   BC.poly.plt,
                   EAC.poly.plt,
                   GS.poly.plt,
                   KC.poly.plt,
                   ncol = 1, nrow = 5, labels = "auto")
ggplot2::ggsave("publ_plots/Figure_2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

# Figure 3. Trends of MHW metrics mean intensity and frequency of selected WBCs… [duration trends not strongly associated with KE features or with regions where the MHWs are most intense.]
Fig.3 <- ggarrange(AC.Fig010 + labs(title = "Mean MHW intensity trend"), AC.Fig016_v2 + labs(title = "MHW count trend"),
                   BC.Fig010 + labs(title = ""), BC.Fig016_v2 + labs(title = ""),
                   EAC.Fig010 + labs(title = ""), EAC.Fig016_v2 + labs(title = ""),
                   GS.Fig010 + labs(title = ""), GS.Fig016_v2 + labs(title = ""),
                   KC.Fig010 + labs(title = ""), KC.Fig016_v2 + labs(title = ""),
                   ncol = 2, nrow = 5, labels = "auto")
ggplot2::ggsave("publ_plots/Figure_3.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)

# •	Suppl. Fig. 3. Full set of panels matching Fig. 3, including also trend in duration.
Supp.3 <- ggarrange(AC.Fig010 + labs(title = "MHW mean int. trend"), AC.Fig016_v2 + labs(title = "MHW count trend"), AC.Fig013_v2 + labs(title = "MHW duration trend"),
                    BC.Fig010 + labs(title = ""), BC.Fig016_v2 + labs(title = ""), BC.Fig013_v2 + labs(title = ""),
                    EAC.Fig010 + labs(title = ""), EAC.Fig016_v2 + labs(title = ""), EAC.Fig013_v2 + labs(title = ""),
                    GS.Fig010 + labs(title = ""), GS.Fig016_v2 + labs(title = ""), GS.Fig013_v2 + labs(title = ""),
                    KC.Fig010 + labs(title = ""), KC.Fig016_v2 + labs(title = ""), KC.Fig013_v2 + labs(title = ""),
                    ncol = 3, nrow = 5, labels = "auto")
ggplot2::ggsave("publ_plots/Suppl.Fig.3.jpg",
                width = 10.5 * (1/3), height = 13 * (1/3), scale = 3.7)

# Figure 4 | Pixel-by-pixel time series correlations between (a, c, e, g, i) MKE vs. mean MHW intensity, and (b, d, f, h, j) EKE vs. mean MHW intensity.
Fig.4 <- ggarrange(AC.fig0x[[1]] + labs(title = "Mean kinetic energy"), AC.fig0x[[2]] + labs(title = "Eddy kinetic energy"),
                   BC.fig0x[[1]] + labs(title = ""), BC.fig0x[[2]] + labs(title = ""),
                   EAC.fig0x[[1]] + labs(title = ""), EAC.fig0x[[2]] + labs(title = ""),
                   GS.fig0x[[1]] + labs(title = ""), GS.fig0x[[2]] + labs(title = ""),
                   KC.fig0x[[1]] + labs(title = ""), KC.fig0x[[2]] + labs(title = ""),
                   ncol = 2, nrow = 5, labels = "auto")
ggplot2::ggsave("publ_plots/Figure_4.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)




