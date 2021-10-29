# RUV-III - PRPS((Pseudo replicate of pseudo sample)) map

sample.info$biology <- sample(letters[1:4], 431, replace = TRUE)
sample.info$new.batch <- paste0(
  sample.info$year_mda, #sample.info$Year,
  '_',
  sample.info$PlateId_mda #sample.info$Plates
)

library(tidyverse)
df_count <- sample.info %>%
  dplyr::count(new.batch, biology)

df_count$use <- 'Un-selected'
df_count$use[df_count$n > 2] <- 'Selected'

ggplot(df_count, aes(x = new.batch, y = biology)) +
  geom_count(aes(color = use)) +
  geom_text(aes(
    label = n,
    hjust = 0.5,
    vjust = 0.5
  )) +
  xlab('Years-plates') +
  ylab('Biological groups') +
  theme_bw()+
  theme(
    axis.line = element_line(colour = 'black', size = .85),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 10),
    legend.position = 'none'
  )
