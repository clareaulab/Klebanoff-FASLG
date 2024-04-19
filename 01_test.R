library(data.table)
library(dplyr)
library(BuenColors)
library(ggbeeswarm)

df <- fread("Fig1_CART_global_obs.csv")
umap <- fread("umap_coords.csv")
u <- fread("umis.csv")

df$FASLG <- u[,2]
df$UMAP1 <- umap[,2]
df$UMAP2 <- umap[,3]
df$FASLGscaled <- log1p(df$FASLG)

if(FALSE){
ggplot(df, aes(x = UMAP1, y = UMAP2, color = FASLGscaled)) + 
  geom_point() +
  facet_wrap(product~CAR)
}

ggplot(df, aes(x = cell_type, y = FASLGscaled)) + 
  geom_violin() + theme_bw() + 
  facet_grid(generic ~ timepoint)

df %>%
  group_by(generic, timepoint, barcode, response, CAR) %>%
  summarize(pct_exp = mean(FASLG > 0)*100, n = n()) -> g_df

# baseline analysis
df %>% filter(timepoint == "D7") %>%
  group_by(cell_type, generic, barcode) %>%
  summarize(pct_exp = mean(FASLG > 0)*100, n = n()) %>%
  filter(n > 50) %>%
  ggplot(aes(x = cell_type, y = pct_exp)) +
  geom_quasirandom() + theme_bw() +
  labs(y = "% expressing FASLG", x = "D7 celltype") +
  facet_wrap(~generic) + pretty_plot(fontsize = 8) + 
  geom_boxplot(color = "dodgerblue", fill = NA, outlier.shape = NA)

g_df %>%
  ggplot(aes(x = interaction(CAR, timepoint), y = pct_exp, color = response, group = interaction(barcode, timepoint))) +
  geom_quasirandom() + theme_bw() +
  geom_line(color = "black") +
  labs(y = "% expressing FASLG", x = "Timepoint") +
  facet_wrap(~generic)

collapse_7 <- g_df %>% 
  mutate(what2 = case_when(
    timepoint == "D7" & CAR ~ "D7CAR",
    timepoint == "D7-CART" ~ "D7CAR",
    timepoint == "D7" & !CAR ~ "D7noCAR",
    timepoint == "Infusion" ~ "cInfusion",
    TRUE ~ timepoint
  )) %>%
  group_by(generic, barcode, what2) %>%
  summarize(totalP = sum(pct_exp*n)/100, totalN = sum(n)) %>%
  mutate(new_pct = totalP/totalN*100)

write.table(collapse_7, file = "data_for_boxplots_ck.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

library(lme4)
collapse_7lm_axi <- collapse_7 %>% filter(what2 %in% c("D7CAR", "D7noCAR")) %>% filter(generic == "Axi-cel")
collapse_7lm_tisa <- collapse_7 %>% filter(what2 %in% c("D7CAR", "D7noCAR")) %>% filter(generic == "Tisa-cel")

glm_axi <- glmer(new_pct/100 ~ what2 + (1|barcode), weights=totalN, family=binomial, data = collapse_7lm_axi)
summary(glm_axi)

glm_tisa <- glmer(new_pct/100 ~ what2 + (1|barcode), weights=totalN, family=binomial, data = collapse_7lm_axi)
summary(glm_tisa)

pnorm(q=-26.55, lower.tail=TRUE)

t.test(collapse_7lm_axi$new_pct[collapse_7lm_axi$what2 == "D7CAR"],
       collapse_7lm_axi$new_pct[collapse_7lm_axi$what2 == "D7noCAR"], paired = TRUE)

t.test(collapse_7lm_tisa$new_pct[collapse_7lm_tisa$what2 == "D7CAR"],
       collapse_7lm_tisa$new_pct[collapse_7lm_tisa$what2 == "D7noCAR"], paired = TRUE)



pX <- ggplot(collapse_7,aes(x = what2, y = new_pct)) +
  geom_quasirandom() + theme_bw() +
  labs(y = "% expressing FASLG", x = "Timepoint") +
  facet_wrap(~generic) + pretty_plot(fontsize = 8) + 
  geom_boxplot(color = "dodgerblue", fill = NA, outlier.shape = NA)
cowplot::ggsave2(pX, file = "cleaner_boxplot.pdf", width = 3.5, height = 1.8)

axi_df <- collapse_7 %>% filter(grepl("^Axi", barcode)) 


t.test(axi_df[[2]], axi_df[[3]], paired = TRUE)
