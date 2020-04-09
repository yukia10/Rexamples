# facet coloring example

library(ggplot2)

data(mpg)

ggplot(mpg, aes(displ, hwy)) +
  facet_wrap(~ drv) +
  geom_rect(aes(fill=drv), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# scale_fill_manual("legend", values=c("4"="#FFCCCC", "f"="#FFCCFF", "r"="#CCCCFF")) +
# theme(legend.position="none") +
  geom_point(shape=1)
