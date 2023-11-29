library(ggplot2)

addeq = function(fit, decimal = 2) {
  eq = substitute(
    italic('y') ~ '=' ~ b * italic('x') ^ a,
    list(
      a = sprintf(paste0('%.', decimal, 'f'), as.numeric(stats::coef(fit))[2]),
      b = sprintf(paste0('%.', decimal, 'f'), as.numeric(stats::coef(fit))[1])
    )
  )
  return(eq)
}

addr2 <- function(actual, predict) {
  differences <- predict - actual
  squared_differences <- differences^2
  sum_squared_differences <- sum(squared_differences)
  return(sum_squared_differences)
}

scatters <- function(
    tree, lst,
    xlabel = expression(bold(paste('Tree cover (%)'))),
    ylabel = expression(bold(paste('LST (', degree, 'C)', sep = ''))),
    output_file
) {
  # data
  value.tree <- raster::getValues(tree)
  value.lst <- raster::getValues(lst)
  d <- data.table::data.table(x = value.tree, y = value.lst)
  d <- d[complete.cases(d)]
  d <- d[x > 0]

    # make bins
  # for (k in 0:100) d[(x >= min(d$x) + k) & (x < min(d$x) + k + 1), bin := (k + 1)]
  d[, bin := ceiling(x)]
  db <- d[, .(bin_x = median(x), bin_y = median(y)), by = bin]
  colnames(db) <- c('bin', 'x', 'y')
  
  # Non-linear regression
  fit <- nls(y ~ a * bin ^ b, data = db, start = list(a = 50, b = -0.01))
  
  # Create a ggplot
  plot <- ggplot(d, aes(x = x, y = y)) +
    geom_point(color = rgb(0.25, 0.5, 0.75, 1), size = 0.5) +
    geom_smooth(
      method = "nls",
      formula = y ~ a * x ^ b,
      se = FALSE,
      method.args = list(start = c(a = 50, b = -0.01)),
      col = 'red',
      lwd = 1.7
    ) +
    theme_classic() +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
    scale_x_continuous(breaks = pretty(d$x)) +
    scale_y_continuous(breaks = pretty(d$y)) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    ggtitle(addeq(fit)) +
    labs(subtitle = addr2(db$y, predict(fit))) +
    annotate(
      "text",x = 80,y = 20,label = xlabel,hjust = 1,vjust = 1,size = 5
    ) +
    annotate(
      "text",x = 0,y = 60,label = ylabel,hjust = 0,vjust = 0,size = 5
    )
  # Save the plot to a file
  ggsave(output_file, plot, width = 8, height = 6, units = "in")
  
  print(plot)
}

