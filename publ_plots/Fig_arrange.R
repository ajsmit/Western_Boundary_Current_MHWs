# create plot arrangements for figures
# load plots from these files first:
# /Users/ajsmit/Dropbox/R/WBCs/setup/Fig00_Aviso+_MKE_EKE.R
# /Users/ajsmit/Dropbox/R/WBCs/eddies/Fig018_Eddy_trajectories.R
# /Users/ajsmit/Dropbox/R/WBCs/setup/Fig06_OISST_MeanInt.R


# Figure 1 ----------------------------------------------------------------

# •	Figure 1 (a-d). MKE, EKE, eddy trajectories, and MHW intensity of the Agulhas Current.
# Use figures produced by Fig00_Aviso_MKE_EKE_v2.R
Fig.1 <- ggarrange(
  AC.Fig00a + labs(title = "Mean kinetic energy"),
  AC.Fig00b + labs(title = "Eddy kinetic energy"),
  AC.Fig018 + labs(title = "Eddy trajectories"),
  AC.Fig06 + labs(title = "Heatwave intensity"),
  ncol = 2,
  nrow = 2,
  labels = "AUTO"
)
ggplot2::ggsave(
  "publ_plots/Figure_1.jpg",
  width = 7.0 * (1 / 3),
  height = 5.2 * (1 / 3),
  scale = 3.7
)

# Supplement Figure 1 -----------------------------------------------------

# •	Suppl. Fig. 1 (a-t). Full set of panels corresponding to Fig. 1.
# Use figures produced by Fig00_Aviso_MKE_EKE_v2.R
Supp.1 <- ggarrange(
  AC.Fig00a,
  AC.Fig00b,
  AC.Fig018,
  AC.Fig06,
  BC.Fig00a,
  BC.Fig00b,
  BC.Fig018,
  BC.Fig06,
  EAC.Fig00a,
  EAC.Fig00b,
  EAC.Fig018,
  EAC.Fig06,
  GS.Fig00a,
  GS.Fig00b,
  GS.Fig018,
  GS.Fig06,
  KC.Fig00a,
  KC.Fig00b,
  KC.Fig018,
  KC.Fig06,
  ncol = 4,
  nrow = 5,
  labels = "AUTO"
)
ggplot2::ggsave(
  "publ_plots/Suppl.Fig.1.jpg",
  width = 14.0 * (1 / 3),
  height = 13 * (1 / 3),
  scale = 3.7
)

# Supplement Figure 2 -----------------------------------------------------

# # •	Suppl. Fig. 2 (a-d). Full set of panels of SST mean trend.
# Supp.2 <- ggarrange(AC.Fig04,
#                     BC.Fig04,
#                     EAC.Fig04,
#                     GS.Fig04,
#                     KC.Fig04,
#                     ncol = 1, nrow = 5, labels = "AUTO")
# ggplot2::ggsave("publ_plots/Suppl.Fig.2.jpg",
#                 width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

# •	Figure 2 (a-e). Masked regions of the Agulhas Current that correspond to MKE, EKE, eddy tracks, and MHW intensity ≥ 90th percentile.
# In /Volumes/GoogleDrive/My Drive/R/WBCs/masks/1.create_MKE_EKE_event_mask_poly.R
Fig.2 <- ggarrange(
  AC.poly.plt,
  BC.poly.plt,
  EAC.poly.plt,
  GS.poly.plt,
  KC.poly.plt,
  ncol = 1,
  nrow = 5,
  labels = "AUTO"
)
ggplot2::ggsave(
  "publ_plots/Figure_2.jpg",
  width = 3.5 * (1 / 3),
  height = 13 * (1 / 3),
  scale = 3.7
)

# Figure 3 old ------------------------------------------------------------

# •	Figure 3 | Pixel-by-pixel time series correlations between (a, c, e, g, i) MKE vs. mean MHW intensity, and (b, d, f, h, j) EKE vs. mean MHW intensity.
# no longer used
Fig.3_old <-
  ggarrange(
    AC.fig0x[[1]] + labs(title = "Mean kinetic energy"),
    AC.fig0x[[2]] + labs(title = "Eddy kinetic energy"),
    BC.fig0x[[1]] + labs(title = ""),
    BC.fig0x[[2]] + labs(title = ""),
    EAC.fig0x[[1]] + labs(title = ""),
    EAC.fig0x[[2]] + labs(title = ""),
    GS.fig0x[[1]] + labs(title = ""),
    GS.fig0x[[2]] + labs(title = ""),
    KC.fig0x[[1]] + labs(title = ""),
    KC.fig0x[[2]] + labs(title = ""),
    ncol = 2,
    nrow = 5,
    labels = "AUTO"
  )
ggplot2::ggsave(
  "publ_plots/Figure_3.jpg",
  width = 7.0 * (1 / 3),
  height = 13 * (1 / 3),
  scale = 3.7
)

# Figure 3 correlations ---------------------------------------------------

# these data were prepared in
# /Volumes/GoogleDrive/My Drive/R/WBCs/correlate/correlate_rolling_means_xx_v3.R

load(file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/AC_corr_plt_v1.RData")
load(file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/BC_corr_plt_v1.RData")
load(file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/EAC_corr_plt_v1.RData")
load(file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/GS_corr_plt_v1.RData")
load(file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/KC_corr_plt_v1.RData")

l <- get_legend(AC.corr.plt)
Fig.3 <- ggarrange(AC.corr.plt + theme(legend.position = 'none'),
                              BC.corr.plt + theme(legend.position = 'none'),
                              EAC.corr.plt + theme(legend.position = 'none'),
                              KC.corr.plt + theme(legend.position = 'none'),
                              GS.corr.plt + theme(legend.position = 'none'),
                              l,
                              ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", NULL))

ggplot2::ggsave("publ_plots/Figure_3_correlations.jpg",
                width = 5.0, height = 4.5, scale = 3.5, units = "cm")

# Figure 5 ----------------------------------------------------------------

# •	Figure 5. Trends of MHW metrics mean intensity and frequency of selected WBCs… [duration trends not strongly associated with KE features or with regions where the MHWs are most intense.]
Fig.5 <-
  ggarrange(
    AC.Fig010 + labs(title = "Mean MHW intensity trend"),
    AC.Fig016_v2 + labs(title = "MHW count trend"),
    BC.Fig010 + labs(title = ""),
    BC.Fig016_v2 + labs(title = ""),
    EAC.Fig010 + labs(title = ""),
    EAC.Fig016_v2 + labs(title = ""),
    GS.Fig010 + labs(title = ""),
    GS.Fig016_v2 + labs(title = ""),
    KC.Fig010 + labs(title = ""),
    KC.Fig016_v2 + labs(title = ""),
    ncol = 2,
    nrow = 5,
    labels = "AUTO"
  )
ggplot2::ggsave(
  "publ_plots/Figure_5.jpg",
  width = 7.0 * (1 / 3),
  height = 13 * (1 / 3),
  scale = 3.7
)

# Supplement Figure 3 -----------------------------------------------------

# •	Suppl. Fig. 3. Full set of panels matching Fig. 4, including also trend in duration.
Supp.3 <-
  ggarrange(
    AC.Fig04 + labs(title = "Mean SST trend"),
    AC.Fig010 + labs(title = "MHW mean int. trend"),
    AC.Fig016_v2 + labs(title = "MHW count trend"),
    AC.Fig013_v2 + labs(title = "MHW duration trend"),
    BC.Fig04 + labs(title = ""),
    BC.Fig010 + labs(title = ""),
    BC.Fig016_v2 + labs(title = ""),
    BC.Fig013_v2 + labs(title = ""),
    EAC.Fig04 + labs(title = ""),
    EAC.Fig010 + labs(title = ""),
    EAC.Fig016_v2 + labs(title = ""),
    EAC.Fig013_v2 + labs(title = ""),
    GS.Fig04 + labs(title = ""),
    GS.Fig010 + labs(title = ""),
    GS.Fig016_v2 + labs(title = ""),
    GS.Fig013_v2 + labs(title = ""),
    KC.Fig04 + labs(title = ""),
    KC.Fig010 + labs(title = ""),
    KC.Fig016_v2 + labs(title = ""),
    KC.Fig013_v2 + labs(title = ""),
    ncol = 4,
    nrow = 5,
    labels = "AUTO"
  )
ggplot2::ggsave(
  "publ_plots/Suppl.Fig.3.jpg",
  width = 14.0 * (1 / 3),
  height = 13 * (1 / 3),
  scale = 3.7
)




