source("Settings/utilities.R")
source("Settings/Colorpalette.R")
output_folder <- "Image/Fused_plot/"
createFolder(output_folder)


stackbar.updating <- function(pl, paletteFixed) {
  library(cowplot)
  library(ggplotify)

  # ---- COSTRUISCO IL PLOT BASE (con legenda) ----
  p1 <- pl +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, size = 7, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      strip.text = element_text(face = "bold", size = 20),
      legend.key.size = unit(0.3, "cm"),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 20, face = "bold")
    ) +
    guides(fill = guide_legend(ncol = 3, title = "Species"))

  # ---- ORDINO I LIVELLI ----
  all_levels <- names(paletteFixed)
  p1$data$OTU_id <- as.character(p1$data$OTU_id)
  p1$data$OTU_id[p1$data$OTU_id == "NA"] <- NA_character_
  if (!("Other" %in% all_levels)) {
    p1$data$OTU_id[p1$data$OTU_id == "Other"] <- NA_character_
  }
  p1$data$OTU_id <- factor(p1$data$OTU_id, levels = all_levels)

  # ---- PALETTE FISSA ----
  p1 <- p1 + scale_fill_manual(values = paletteFixed, drop = FALSE, na.translate = FALSE)

  # ---- ESTRAI LEGENDA (QUI è il punto giusto!) ----
  legend <- get_legend(p1)

  # ---- RIMUOVO LEGENDA PER IL PLOT FINALE ----
  p1 <- p1 + theme(legend.position = "none")

  # ---- MODIFICA STRIP (come fai tu) ----
  g <- ggplot_gtable(ggplot_build(p1))
  stripr <- which(grepl("strip-t", g$layout$name))

  fills <- paletteStatus
  k <- 1
  for (i in stripr) {
    j <- which(grepl(
      "strip.background",
      g$grobs[[i]]$grobs[[1]]$childrenOrder
    ))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }

  final_plot <- ggplotify::as.ggplot(g)

  # ---- OUTPUT ----
  return(
    patchwork::wrap_plots(
      final_plot,
      patchwork::wrap_elements(legend, ignore_tag = TRUE),
      ncol = 1,
      heights = c(0.9, 0.1)
    )
  )
}

#### patch1

p1 = readRDS(file = "./Output/NEW_PLOT_RDS/Bacteria_StackedBar_Species.rds") + 
  theme(plot.margin = margin(t = 10, r = 10, b = 80, l = 10))  # Aumentare il margine inferiore
p1 = stackbar.updating(p1, paletteFixed = paletteBact)
p1

p2 = readRDS(file = "./Output/NEW_PLOT_RDS/Bacteria_alpha_category.rds") + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = paletteStatus)+labs(fill = "Disease Status")

p3 = readRDS(file = "./Output/NEW_PLOT_RDS/Bacteria_beta_cat.rds") + 
  theme(axis.text.x = element_text(size=24,face = "bold"))+
  scale_color_manual(values = unname(paletteStatus), breaks = c("HD","MS"),labels = c("HEALTHY","MS") )

p4 = readRDS(file = "./Output/NEW_PLOT_RDS/Archaea_StackedBar_Species.rds") + theme(plot.margin = margin(t = 10, r = 10, b = 80, l = 10),axis.text.x = element_text(angle = 45, hjust = 1))

p4 = stackbar.updating(p4, paletteFixed = paletteArchaea)
p4
p5 = readRDS(file = "./Output/NEW_PLOT_RDS/Archaea_alpha_category.rds") +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = paletteStatus)+labs(fill = "Disease Status")

p6 = readRDS(file = "./Output/NEW_PLOT_RDS/Archaea_beta_cat.rds") + theme(axis.text.x = element_text(size=24, face = "bold"))+
  scale_color_manual(values = unname(paletteStatus), breaks = c("HD","MS"),labels = c("HEALTHY","MS") )
p7 = readRDS(file = "./Output/NEW_PLOT_RDS/Eukaryote_StackedBar_Species.rds") + theme(plot.margin = margin(t = 10, r = 10, b = 80, l = 10),axis.text.x = element_text(angle = 45, hjust = 1)) 

p7 = stackbar.updating(p7, paletteFixed = paletteEuk)
p7

p8 = readRDS(file = "Output/NEW_PLOT_RDS/Eukaryote_alpha_category.rds") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = paletteStatus)+labs(fill = "Disease Status")

p9 = readRDS(file = "./Output/NEW_PLOT_RDS/Eukaryote_beta_cat.rds") +
  theme(axis.text.x = element_text(size=24, face = "bold"))+
  scale_color_manual(values = unname(paletteStatus), breaks = c("HD","MS"),labels = c("HEALTHY","MS") )

patchwork = (
  patchwork::wrap_elements(full = p1) +
  patchwork::wrap_elements(full = p4) +
  patchwork::wrap_elements(full = p7) +
  p2 + p5 + p8 + p3 + p6 + p9
) + plot_layout(
  ncol = 3,
  widths = c(1, 1, 1),
  heights = c(1, 0.5, 1),
  guides = "collect"
)

patchwork = patchwork & plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 30, face = "bold"),
        text = element_text(size=24, face = "bold"),
        plot.title = element_text(size=24, face = "bold"),
        legend.title = element_text(size=24, face = "bold"),
        legend.text = element_text(size=24, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(colour = "black", size=24, face = "bold"),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.ticks.x = element_blank())

ggsave(gsub(" ","",paste(output_folder,"patchworkplot_last.pdf")), plot = patchwork, width = 30, height = 20, dpi = 300)

p2 = readRDS(file = "Output/NEW_PLOT_RDS/Bacteria_alpha_gc_treatment.rds") + theme(axis.text.x = element_blank(),
                                                                          axis.title.x = element_blank(),
                                                                          axis.ticks.x = element_blank())
p2
p3 = readRDS(file = "Output/NEW_PLOT_RDS/Bacteria_beta_gc.rds") + theme(axis.text.x = element_text(size=24,face = "bold"))

p5 = readRDS(file = "Output/NEW_PLOT_RDS/Archaea_alpha_gc_treatment.rds") + theme(axis.text.x = element_blank(), 
                                                                          axis.title.x = element_blank(), 
                                                                          axis.ticks.x = element_blank())
p6 = readRDS(file = "Output/NEW_PLOT_RDS/Archaea_beta_gc.rds") + theme(axis.text.x = element_text(size=24, face = "bold"))

p8 = readRDS(file = "Output/NEW_PLOT_RDS/Eukaryote_alpha_gc_treatment.rds") + theme(axis.text.x = element_blank(), 
                                                                         axis.title.x = element_blank(),
                                                                         axis.ticks.x = element_blank())
p9 = readRDS(file = "Output/NEW_PLOT_RDS/Eukaryote_beta_gc.rds") + theme(axis.text.x = element_text(size=24, face = "bold"))

patchwork = (p2 / p3)| (p5 / p6) | (p8 / p9)
patchwork
patchwork & plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 42, face = "bold"),
        text = element_text(size=24, face = "bold"),
        axis.title.y = element_text(size=24, face = "bold"),
        axis.title.x = element_text(size=24, face = "bold"),
        axis.text.y = element_text(size=24, face = "bold"),
        plot.title = element_text(size=24, face = "bold"),
        legend.title = element_text(size=24, face = "bold"),
        legend.text = element_text(size=24, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(colour = "black", size=24, face = "bold"),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.ticks.x = element_blank())
patchwork

ggsave(gsub(" ","",paste(output_folder,"patchworkplot_GCtreatment.pdf")), plot = patchwork, width = 30, height = 20, dpi = 300)
