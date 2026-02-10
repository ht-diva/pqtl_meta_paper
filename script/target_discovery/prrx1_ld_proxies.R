
# LD proxies from LD Link (GRCh37)
# reference variant is  , the GTEx eQTL with largest PIP
path_prrx1_proxies <- "/scratch/dariush.ghasemi/projects/pqtl_meta_paper/data/prrx1_proxies_1000genom_europeans_build_37.txt"

# lead pQTL in meta-analyis
index_pqtl <- 170630259 #rs601938
indep_pqtl <- c(
  170136068, 170565756, 170580655, 170583760, 170600133, 170600737,
  170617420, 170618504, 170628097, 170634487, 170641733, 170641803, 170652683
  )

prrx1_proxies <- fread(path_prrx1_proxies)

# Check LD between GTEx eQTL and Interval pQTLs
prrx1_proxies %>%
  mutate(
    pos = str_remove(Coord, "chr1:") %>% as.integer(),
    pqtl = ifelse(pos %in% c(index_pqtl, indep_pqtl), RS_Number, "")
  ) %>%
  arrange(pos) %>%
  #select(RS_Number:R2) %>% gt::gt()
  ggplot(aes(x = pos, y = R2, fill = R2)) +
  geom_point(shape = 24, size = 2.5, color = "grey30") +
  scale_x_continuous(
    breaks = seq(170428381, 170848260, 50000),
    labels = c(seq(170.4, 170.8, .05)),
    #expand = c(0,0)
    ) +
  ggrepel::geom_text_repel(aes(label = pqtl), box.padding = 0.55, max.overlaps = Inf) +
  scale_fill_gradient2(
    low = "grey", 
    mid = "orange", 
    high = "brown", 
    midpoint = .6
  ) +
  labs(x = "\nGenomic position (Mb, GRCh37)", y = "LD (r2)") +
  theme_light() +
  theme(
    legend.position = c(0.91, 0.75),
    legend.background = element_blank(),
    axis.title = element_text(size = 12, face = 2),
    axis.text = element_text(size = 11)
  )
  


ggsave(
  filename = "02-Feb-26_prrx1_eqtls_gtex_ld_with_pqtls.png", 
  height = 6.5, width = 11, dpi = 300
)

