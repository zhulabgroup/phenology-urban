library(ggplot2)
library(ggpmisc)
library(future)

scatters_lst_phe <- function(
    rs_lst, rs_phe,
    xlabel,
    ylabel,
    output_file
) {
  # data
  value.phe <- raster::getValues(rs_phe)
  value.lst <- raster::getValues(rs_lst)
  d <- data.table::data.table(x = value.lst, y = value.phe)
  d <- d[complete.cases(d)]
  d <- d[y > 0 & y <= 366]
  #d <- d[x > 0]  #for tcc
  
  # make bins
  # for (k in 0:100) d[(x >= min(d$x) + k) & (x < min(d$x) + k + 1), bin := (k + 1)]
  d[, bin := ceiling(x)]
  db <- d[, .(bin_x = median(x), bin_y = median(y)), by = bin]
  colnames(db) <- c('bin', 'x', 'y')
  

  # Create a ggplot
  # Create a ggplot with hexagonal bins
  plot <- ggplot(d, aes(x = x, y = y)) +
    geom_hex(bins = 60, color = NA) +  
    scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
    theme_classic() +
    labs(x = xlabel, y = ylabel, fill = "Point Count") +
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    scale_x_continuous(breaks = pretty(d$x)) +
    scale_y_continuous(breaks = pretty(d$y)) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 18)
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", formula = y ~ x) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                 formula = y ~ x,
                 parse = TRUE,
                 label.x = "right", label.y = 0.9,
                 size = 6, color = "blue")
  ggsave(output_file, plot, width = 8, height = 6, units = "in")
  
  
    print(plot)}
  
  
  