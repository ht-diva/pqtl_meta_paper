
library(readxl)
library(dplyr)
library(ggplot2)

files <- read_excel("/Users/giulia.pontali/Desktop/drug_rep/supplementary_table_7.xlsx", sheet=2)
files <- as.data.frame(files)

files$drug_category <- ifelse(
  files$drug_rep == "Drug repurposing",
  "Drug repurposing",
  ifelse(files$drug_rep %in% c(
    "Drug exists (Exact EFO_ID)",
    "Drug exists (Genetic_phenotype_phecode)",
    "Drug exists (Matched_genetic_phenotype)",
    "Drug exists (Parent EFO)"
    ),
    "Drug rediscovery", NA)
)

df <- subset(files, !is.na(files$drug_category))
df <- df[which(df$PHENOTYPE_class=="diseases"),]

efo_counts <- df %>%
  group_by(EFO_Term, drug_category) %>%
  summarise(Count = n(), .groups = "drop") %>%
    rename(Source = drug_category)


efo_counts %>%
  group_by(EFO_Term) %>%
  dplyr::mutate(cumulative_count = sum(Count)) %>%
  ungroup() %>%
  dplyr::mutate(
    efo_sorted = fct_reorder(EFO_Term, cumulative_count),
    source_reord = factor(Source, levels = c("Drug repurposing", "Drug rediscovery"))
    ) %>%
  ggplot(aes(x = efo_sorted, y = Count, fill = source_reord)) +
  geom_bar(
    stat = "identity",
    width = .8,
    position = position_dodge(0.9, preserve = "single")
    ) +
  scale_fill_manual(
    name = "Source",
    values = c(
    "Drug rediscovery" = "#0072B2",
    "Drug repurposing" = "#E69F00"
  )) +
  labs(
    y = "Count of entries",
    x = "EFO terms"
  ) +
  coord_flip() +
  theme_bw(base_size = 10) +
  theme(
    axis.ticks.length = unit(2.2, "mm"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = 1),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_blank()
  )

ggsave(filename = "26-Nov-25_efo_count_terms.png", 
       last_plot(), height = 12, width = 9, dpi = 500)

 
