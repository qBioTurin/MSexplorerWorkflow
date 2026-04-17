library(tidyverse)
library(phyloseq)
library(MicrobiomeStat)
library(topicmodels)
library(tidytext)
library(ldatuning)
library(cowplot)
library(curatedMetagenomicData)

ps = readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
#ps= subset_samples(ps, category == "MS")
(ps_g <- tax_glom(ps, taxrank = "Genus"))

taxa_names(ps_g) <- tax_table(ps_g)[, 6]
head(taxa_names(ps_g))

count_matrix <- data.frame(t(data.frame(otu_table(ps_g)))) ###conte, giusto??

result = readRDS("Output/Topic_modeling/HDvsMS_findtopicnumber.RDS")
#result = readRDS("Output/Topic_modeling/MS_findtopicnumber.RDS")
result <- FindTopicsNumber(
  count_matrix,
  topics = seq(from = 3, to = 50, by = 3),
  metrics = c("CaoJuan2009", "Arun2010"),
  method = "VEM",
  control = list(seed = 243),
  mc.cores = 4,
  verbose = TRUE
)

#saveRDS(result, "Output/Topic_modeling/HDvsMS_findtopicnumber.RDS")
#saveRDS(result, "Output/Topic_modeling/MS_findtopicnumber.RDS")

my_plot <- function (values){  #updating to allow plotting via cowplot::plot_grid
  if ("LDA_model" %in% names(values)) {
    values <- values[!names(values) %in% c("LDA_model")]
  }
  columns <- base::subset(values, select = 2:ncol(values))
  values <- base::data.frame(values["topics"], base::apply(columns, 
                                                           2, function(column) {
                                                             scales::rescale(column, to = c(0, 1), from = range(column))
                                                           }))
  values <- reshape2::melt(values, id.vars = "topics", na.rm = TRUE)
  values$group <- values$variable %in% c("Griffiths2004", "Deveaud2014")
  values$group <- base::factor(values$group, levels = c(FALSE, 
                                                        TRUE), labels = c("Minimize", "Maximize"))
  p <- ggplot(values, aes_string(x = "topics", y = "value", 
                                 group = "variable"))
  p <- p + geom_line()
  p <- p + geom_point(aes_string(shape = "variable"), size = 3)
  p <- p + guides(size = FALSE, shape = guide_legend(title = "Metrics:"))
  p <- p + scale_x_continuous(breaks = values$topics)
  p <- p + labs(x = "\nNumber of Topics", y = NULL)
  p <- p + facet_grid(group ~ .)
  p <- p + theme_bw() %+replace% theme(panel.grid.major.y = element_blank(), 
                                       panel.grid.minor.y = element_blank(), panel.grid.major.x = element_line(colour = "grey70"), 
                                       panel.grid.minor.x = element_blank(), legend.key = element_blank(), 
                                       strip.text.y = element_text(angle = 90),
                                       legend.position = "bottom")
  return(p)
}

(p1 <- my_plot(result)) ## choose 12, for genus HD vs MS; choose 7, for MS
ggsave("Output/Topic_modeling/MS_choosetopic.png", plot = p1, width = 6, height = 12, dpi = 300)

lda_k36 <- LDA(count_matrix, k = 12, method = "VEM", control = list(seed = 243)) # anche sani - genus
#lda_k36 <- LDA(count_matrix, k = 7, method = "VEM", control = list(seed = 243)) # solo MS - genus
#lda_k36 <- LDA(count_matrix, k = 13, method = "VEM", control = list(seed = 243)) # solo MS - species
#lda_k36 <- LDA(count_matrix, k = 42, method = "VEM", control = list(seed = 243)) # HD vs MS - species
#lda_k36 <- LDA(count_matrix, k = 6, method = "VEM", control = list(seed = 243)) # MS pos 
#lda_k36 <- LDA(count_matrix, k = 8, method = "VEM", control = list(seed = 243)) # MS neg

b_df <- data.frame(tidy(lda_k36, matrix = "beta"))

g_df <- data.frame(tidy(lda_k36, matrix = "gamma")) %>%
  arrange(document, topic)

head(b_df)
head(g_df)

lib_size_df <- sample_sums(ps_g) %>%
  # sample_sums returns a named vector; make a tibble with explicit names
  tibble::enframe(name = "document", value = "read_count")

# ensure g_df has the needed columns
stopifnot(all(c("document", "topic", "gamma") %in% colnames(g_df)))

tm_df <- lib_size_df %>%
  left_join(g_df, by = "document") %>%
  mutate(topic_count = round(read_count * gamma, 0)) %>%
  select(document, topic, topic_count) %>%
  pivot_wider(names_from = topic,
              values_from = topic_count,
              values_fill = 0) %>%
  # prefix topic columns; keep 'document' untouched
  rename_with(.cols = -document, .fn = ~ paste0("Topic_", .)) %>%
  # explicitly use tibble version to avoid masking issues
  tibble::column_to_rownames(var = "document") %>%
  t() %>%
  as.data.frame()

(ps_topic_g <- phyloseq(
  sample_data(ps_g),
  otu_table(tm_df, taxa_are_rows = TRUE)))

meta = data.frame(sample_data(ps_topic_g))
# meta <- meta %>%
#   mutate(gc_treatment = as.factor(gc_treatment)
#         gadolinium_contrast = as.factor(gadolinium_contrast), 
#          lesion_burden = as.factor(lesion_burden), 
#          spinal_cord_lesion = as.factor(spinal_cord_lesion), 
#          subtentorial_lesions = as.factor(subtentorial_lesions), 
#          sex = as.factor(sex), 
#          age = as.numeric(age), 
#          bmi = as.numeric(bmi)) # solo MS
meta <- meta %>%
   mutate(category = as.factor(category),
          sex = as.factor(sex), 
          age = as.numeric(age), 
          bmi = as.numeric(bmi))# anche healthy


###### 
#save.image(file='Output/Topic_modeling/topic_HDMS.RData')
#load('/Users/doratortarolo/Documents/Microbiome/MS/MS_HD/Rebuttal/MSexplorerWorkflow-main 2/Output/TopicModeling/topic_HDMS.RData')
######

linda <- LinDA::linda(otu.tab = data.frame(otu_table(ps_topic_g)), meta = meta, 
               formula = '~ category + age + sex + bmi', alpha = 0.05, winsor.quan = 0.97,
               prev.cut = 0.3, lib.cut = 0, n.cores = 1)

linda_res_df <- data.frame(linda$output$category) %>%
  dplyr::arrange(padj) %>%
  rownames_to_column(var = "Topic")

fdr_linda <- linda_res_df %>%
  mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
  separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
  mutate(Topic = gsub("_", " ", Topic)) %>%
  arrange(desc(t_n)) %>%
  mutate(Topic = factor(Topic, levels = Topic)) # MS: no significative padj (not even close); signif in healthy!!

(p2 <- ggplot(data = fdr_linda, aes(x = Topic, y = log2FoldChange, fill = reject)) +
    geom_col(alpha = .8) +
    labs(y = "\nLog2 Fold-Change", x = "") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    coord_flip() +
    scale_fill_manual(values = c("grey", "dodgerblue")) +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme(legend.position = "none")) # HD vs MS: 2 signif topics!!
ggsave("Output/Topic_modeling/MS_topics_subtentorial_lesions.png", plot = p2, width = 6, height = 12, dpi = 300)
# subtentorial_lesions: 1 signif
# check topic 4
b_plot <- b_df %>%
  dplyr::filter(topic == 11) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 20) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))

(p3 <- ggplot(data = b_plot, aes(x = beta, y = term, color = term)) +
    geom_point(aes(size = beta)) +
    labs(x = "\nTopic 11: HD vs MS: Genus-Topic Probabilities (beta)", y = "") +
    theme(legend.position = "none"))
ggsave("Output/Topic_modeling/MS_subtentorial_lesions_signiftopic4_probabilities.png", plot = p3, width = 6, height = 12, dpi = 300)

# Selecting topic five as an example, we can use the information in the per-topic-per-word probabilities matrix to examine which genera have the highest probabilities of assignment to this Topic_modeling/community type.
# Below we can see that Segatella has the highest probability of assignment followed by several others with lower values.

linda_g_all <- LinDA::linda(otu.tab = data.frame(otu_table(ps_g)), meta = data.frame(sample_data(ps_g)), 
                     formula = '~ category + age + sex + bmi', alpha = 0.05, winsor.quan = 0.97,
                     prev.cut = 0, lib.cut = 0, n.cores = 1)
linda_res_g_df <- data.frame(linda_g_all$output$category) %>%
  dplyr::arrange(padj) %>%
  rownames_to_column(var = "Genus")

keep_g <- b_plot$term

log2_df <- linda_res_g_df %>%
  filter(Genus %in% keep_g) %>%
  arrange(desc(Genus)) %>%
  mutate(Genus = factor(Genus, levels = Genus))

(p4 <- ggplot(data = log2_df, aes(x = Genus, y = log2FoldChange, fill = Genus)) +
    geom_col(alpha = .8) +
    labs(y = "HD vs MS, Topic 11: \nLog2 Fold-Change", x = "") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme(legend.position = "none"))
ggsave("Output/Topic_modeling/HDMS_signiftopic11_FC.png", plot = p4, width = 6, height = 12, dpi = 300)
