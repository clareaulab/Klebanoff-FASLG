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

ggplot(df, aes(x = UMAP1, y = UMAP2, color = FASLGscaled)) + 
  geom_point() +
  facet_wrap(product~CAR)

ggplot(df, aes(x = cell_type, y = FASLGscaled)) + 
  geom_violin() + theme_bw() + 
  facet_grid(generic ~ timepoint)

df %>%
  group_by(generic, timepoint, barcode, response) %>%
  summarize(pct_exp = mean(FASLG > 0)*100, n = n()) %>%
  ggplot(aes(x = timepoint, y = pct_exp, color = response)) +
  geom_quasirandom() + theme_bw() +
  labs(y = "% expressing FASLG", x = "Timepoint") +
  facet_wrap(~generic)
