#' Plot a matrix in Hi-C style
#' 
#' Style can be Juicebox or higlass.
#' Huge matrices can take a long time to plot!!!
#' 
#' @param mat the matrix to plot; will use row and colnames for labeling
#' @param max the maximum value for plotting
#' @param style can be 'juicebox' or 'higlass' (default)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
#' @export

plot_mat <- function(mat, max = NULL, style = 'higlass', 
                     main = '', xlab = '', ylab = '') {
  if (!class(mat)[1]=='matrix') {
    stop('input should be a matrix')
  }
  gg <- melt(mat)
  if (!is.null(max)) {
    gg$value <- ifelse(gg$value > max, max, gg$value)
  }
  gg <- gg |>
    ggplot(aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_y_reverse() +
    theme_void() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(main)
  if (style == 'higlass') {
    gg <- gg +
      scale_fill_distiller(palette='YlOrRd', direction = 1)
  } else if (style == 'juicebox') {
    gg <- gg +
      scale_fill_gradient(low = 'white', high = 'red')
  }
  gg  
}