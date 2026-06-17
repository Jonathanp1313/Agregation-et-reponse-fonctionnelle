library(ggplot2)

library(ggplot2)

df_solitaires <- read.table("solitaires.txt", comment.char = "#", header = FALSE,
                            col.names = c("N", "capture"))
df_proie_banc <- read.table("proie_banc.txt", comment.char = "#", header = FALSE,
                            col.names = c("N", "capture"))
df_pred_banc  <- read.table("pred_banc.txt",  comment.char = "#", header = FALSE,
                            col.names = c("N", "capture"))
df_banc_banc  <- read.table("banc_banc.txt",  comment.char = "#", header = FALSE,
                            col.names = c("N", "capture"))

df_solitaires$groupe  <- "Proies indépendantes" ; df_solitaires$facette <- "Prédateurs indépendants"
df_proie_banc$groupe  <- "Proies en banc"        ; df_proie_banc$facette <- "Prédateurs indépendants"
df_pred_banc$groupe   <- "Proies indépendantes" ; df_pred_banc$facette  <- "Prédateurs en banc"
df_banc_banc$groupe   <- "Proies en banc"        ; df_banc_banc$facette  <- "Prédateurs en banc"

df <- rbind(df_solitaires, df_proie_banc, df_pred_banc, df_banc_banc)
 
x11()
ggplot(df, aes(x = N, y = capture, color = groupe)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Proies indépendantes" = "tomato",
                                "Proies en banc"       = "steelblue")) +
  facet_wrap(~ facette) +
  theme_minimal() +
  labs(x     = "N (nombre de proies)",
       y     = "Taux de capture",
       title = "Réponse fonctionnelle",
       color = NULL)