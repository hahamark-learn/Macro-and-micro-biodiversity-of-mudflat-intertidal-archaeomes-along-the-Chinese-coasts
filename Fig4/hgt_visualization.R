#hgtector 
library(statnet)
library(circlize)
#测试数据
category <- paste0("sample", "_", 1:9)
percent <- sort(sample(40:80, 9))
color <- rev(rainbow(length(percent)))

circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 100))

circos.track(
  ylim = c(0.5, length(percent)+0.5), track.height = 0.8, 
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    # 添加圆形中心线
    circos.segments(rep(xlim[1], 9), 1:9,
                    rep(xlim[2], 9), 1:9,
                    col = "#CCCCCC")
    # 添加 9 个圆形矩形
    circos.rect(rep(0, 9),
                1:9 - 0.45,
                percent,
                1:9 + 0.45,
                col = color,
                border = "white")
    # 添加文本信息
    circos.text(
      rep(xlim[1], 9),
      1:9,
      paste(category, " - ", percent, "%"),
      facing = "downward",
      adj = c(1.05, 0.5),
      cex = 0.8
    )
    # 添加轴信息
    breaks = seq(0, 85, by = 5)
    circos.axis(
      h = "top",
      major.at = breaks,
      labels = paste0(breaks, "%"),
      labels.cex = 0.6
    )
  })
