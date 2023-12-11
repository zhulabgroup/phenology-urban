library(ggplot2)
library(ggpmisc)

scatters_tcc_lst <- function(
    tree, lst, xlabel, ylabel,output_file
) {
  # data
  value.tree <- raster::getValues(tree)
  value.lst <- raster::getValues(lst)
  d <- data.table::data.table(x = value.tree, y = value.lst)
  d <- d[complete.cases(d)]
  d <- d[x > 0 & x <= 100]
  d <- d[y > 0]

    # make bins
  # for (k in 0:100) d[(x >= min(d$x) + k) & (x < min(d$x) + k + 1), bin := (k + 1)]
  #d[, bin := ceiling(x)]
  #db <- d[, .(bin_x = median(x), bin_y = median(y)), by = bin]
  #colnames(db) <- c('bin', 'x', 'y')
  
  # Non-linear regression
  #fit <- nls(y ~ a * bin ^ b, data = db, start = list(a = 50, b = -0.01))
  
  # Create a ggplot
  # Apply log transformation to x and y
  d$log_x <- log(d$x)
  d$log_y <- log(d$y)
  
  ggplot(d, aes(x = log_x, y = log_y)) +
    geom_point(color = rgb(0.25, 0.5, 0.75, 1), size = 0.5) +
    #geom_hex(bins = 60, color = NA) +  
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_classic() +
    labs(x = xlabel, y = ylabel, fill = "Point Count") +
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
    scale_x_continuous(breaks = pretty(d$log_x)) +
    scale_y_continuous(breaks = pretty(d$log_y)) +
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
                 size = 6, color = "blue")+
    geom_smooth(
      method = "nls",
      formula = y ~ a * x ^ b,
      se = FALSE,
      method.args = list(start = c(a = 50, b = -0.01)),
      col = 'red',
      lwd = 1.7
    )
  
  
  # Save the plot to a file
  # ggsave(output_file, plot, width = 8, height = 6, units = "in")
  

}

