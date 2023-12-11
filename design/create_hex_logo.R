pkg_color <- "red"
fill_color <- "lightgrey"
out_path <- file.path(getwd(), "man/figures/logo-large")

idm_plot <- function() {
  prodlim::plotIllnessDeathModel(
    style = 1,
    box1.label = "0",
    box2.label = "1",
    box3.label = "2",
    box1.col = "green",
    box2.col = "yellow",
    box3.col = "red",
    cex = 2
  )
}
idm_plot()

tmp <- tempdir()
idm_fig <- file.path(tmp, "idm_plot.png")
png(filename = idm_fig, res = 100, width = 500, height = 500, bg = "transparent")
idm_plot()
dev.off()

plot(magick::image_read(idm_fig))

make_hexplot <- function() {
  require(hexSticker)
  require(showtext)

  for (out_file in paste(out_path, c("svg", "png"), sep = ".")) {
    hexSticker::sticker(
      subplot = idm_fig,
      package = "simIDM",
      p_size = 140,
      p_color = pkg_color,
      p_y = 1.5,
      s_y = 0.8,
      s_x = 0.9,
      s_width = 0.48,
      s_height = 0.48,
      h_fill = fill_color,
      h_color = pkg_color,
      h_size = 2,
      url = "github.com/insightsengineering/simIDM",
      u_size = 24,
      u_color = pkg_color,
      filename = out_file,
      p_family = "mono",
      white_around_sticker = TRUE,
      dpi = 2000
    )
  }
}

make_hexplot()
usethis::use_logo(paste0(out_path, ".png"))
