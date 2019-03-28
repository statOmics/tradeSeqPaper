library(eulerr)
library(RColorBrewer)
library(magick)
counts <- c(139 + 263, 263, 69 + 263)
names(counts) <- c("edgeR", "edgeR&tradeR", "tradeR")
png("Venn.png", width = 4, height = 2.5, units = "in", res = 1000)
plot(euler(counts, shape = "ellipse", input = "union"),
  quantities = T,
  fills = c("#ff7f00",  "#4292C6", "#4DAF4A")
)
dev.off()

Venn <- image_read("Venn.png")
Venn <- image_border(Venn, color = "black", geometry = "20x20")
Venn <- image_annotate(Venn, "1533", size = 200, location = "+100+100",
                       color = "black", gravity = "southeast")
Venn
image_write(Venn, path = "Venn.png", format = "png")
