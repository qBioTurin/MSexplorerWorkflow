source("Settings/utilities.R")
source("Settings/Colorpalette.R")
output_folder = "Image/Fused_plot/"
createFolder(output_folder)


stackbar.updating = function(pl, paletteFixed){
  p1 = pl  +
    theme_bw()+
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 90, size = 6, face = "bold"),
          axis.text.y = element_text(size=12, face = "bold"),
          strip.text = element_text(face = "bold",size=20),
          legend.key.size = unit(0.3, "cm"),
          axis.title.x = element_text(size=18, face = "bold"),
          axis.title.y = element_text(size=18, face = "bold"),
          legend.text = element_text(size=10)
    )+
    guides(fill=guide_legend(ncol=3,title = "Species"))
   s1<-(data.frame(
   OTU_id = levels(p1$data$OTU_id),
   fill   = paletteFixed[levels(p1$data$OTU_id)]
 ))
#  p1<-s1$OTU_id
#   # forziamo i livelli dei dati a quelli della palette fissa
#   p1$data$OTU_id <- factor(p1$data$OTU_id, levels = names(paletteFixed))

  # usiamo la palette fissa, anche se non tutti i livelli sono presenti
  p1 <- p1 + scale_fill_manual(values = paletteFixed, drop = FALSE)

  g <- ggplot_gtable(ggplot_build(p1))
  stripr <- which(grepl('strip-t', g$layout$name))
  
  fills <- paletteStatus
  k <- 1
  for (i in stripr) {
    j <- which(grepl('strip.background', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  return(ggplotify::as.ggplot(g))
}

####patch1

p1 = readRDS(file = "./Output/NEW_PLOT_RDS/Bacteria_StackedBar_Species.rds") + 
  theme(plot.margin = margin(t = 10, r = 10, b = 80, l = 10))  # Aumentare il margine inferiore
p1
p1 = stackbar.updating(pl=p1, paletteFixed = paletteBact)
p1

p2 = readRDS(file = "./Output/NEW_PLOT_RDS/Bacteria_alpha_category.rds") + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = paletteStatus)+labs(fill = "Disease Status")
p2
p3 = readRDS(file = "./Output/NEW_PLOT_RDS/Bacteria_beta_cat.rds") + 
  theme(axis.text.x = element_text(size=15))+
  scale_color_manual(values = unname(paletteStatus), breaks = c("HD","MS"),labels = c("HEALTHY","MS") )

p4 = readRDS(file = "./Output/NEW_PLOT_RDS/Archaea_StackedBar_Species.rds") + theme(plot.margin = margin(t = 10, r = 10, b = 80, l = 10),axis.text.x = element_text(angle = 45, hjust = 1))
p4
p4 = stackbar.updating(pl=p4, paletteFixed = paletteArchaea)
p4
p5 = readRDS(file = "./Output/NEW_PLOT_RDS/Archaea_alpha_category.rds") +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = paletteStatus)+labs(fill = "Disease Status")

p6 = readRDS(file = "./Output/NEW_PLOT_RDS/Archaea_beta_cat.rds") + theme(axis.text.x = element_text(size=15))+
  scale_color_manual(values = unname(paletteStatus), breaks = c("HD","MS"),labels = c("HEALTHY","MS") )
p7 = readRDS(file = "./Output/NEW_PLOT_RDS/Eukaryote_StackedBar_Species.rds") + theme(plot.margin = margin(t = 10, r = 10, b = 80, l = 10),axis.text.x = element_text(angle = 45, hjust = 1)) 
p7
p7 = stackbar.updating(pl=p7, paletteFixed = paletteEuk)
p7
p8 = readRDS(file = "Output/NEW_PLOT_RDS/Eukaryote_alpha_category.rds") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = paletteStatus)+labs(fill = "Disease Status")

p9 = readRDS(file = "./Output/NEW_PLOT_RDS/Eukaryote_beta_cat.rds") +
  theme(axis.text.x = element_text(size=15))+
  scale_color_manual(values = unname(paletteStatus), breaks = c("HD","MS"),labels = c("HEALTHY","MS") )
lapply(list(p1,p2,p3,p4,p5,p6,p7,p8,p9), class)

patchwork = (p1 + p4 + p7 + p2 + p5 + p8 + p3 + p6 + p9) + plot_layout(ncol = 3, 
                                                                         widths = c(1, 1, 1), 
                                                                         heights = c(1, 0.5, 1),
                                                                       guide = "collect") 

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

ggsave(gsub(" ","",paste(output_folder,"patchworkplot.pdf")), plot = patchwork, width = 30, height = 25, dpi = 300)

p2 = readRDS(file = "Output/NEW_PLOT_RDS/Bacteria_alpha_gc_treatment.rds") + theme(axis.text.x = element_blank(),
                                                                          axis.title.x = element_blank(),
                                                                          axis.ticks.x = element_blank())
p2
p3 = readRDS(file = "Output/NEW_PLOT_RDS/Bacteria_beta_gc.rds") + theme(axis.text.x = element_text(size=15))
p3
p5 = readRDS(file = "Output/NEW_PLOT_RDS/Archaea_alpha_gc_treatment.rds") + theme(axis.text.x = element_blank(), 
                                                                          axis.title.x = element_blank(), 
                                                                        axis.ticks.x = element_blank())
 p5                                                                 
p6 = readRDS(file = "Output/NEW_PLOT_RDS/Archaea_beta_gc.rds") + theme(axis.text.x = element_text(size=15))
p6
p8 = readRDS(file = "Output/NEW_PLOT_RDS/Eukaryote_alpha_gc_treatment.rds") + theme(axis.text.x = element_blank(), 
                                                                         axis.title.x = element_blank(),
                                                                         axis.ticks.x = element_blank())
p8
p9 = readRDS(file = "Output/NEW_PLOT_RDS/Eukaryote_beta_gc.rds") + theme(axis.text.x = element_text(size=15))
p9
 patchwork = (p2 / p3)| (p5 / p6) | (p8 / p9)
# patchwork
# patchwork = (p2 | p5 + p8 + p3 + p6 + p9) + plot_layout(ncol = 3, widths = c(1, 1, 1), heights = c(1, 1))
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

ggsave(gsub(" ","",paste(output_folder,"patchworkplot_GCtreatment.pdf")), plot = patchwork, width = 30, height = 20, dpi = 300)
