# smoker.R - Smoking rate in Germany

# Data source:
#   30 years after the fall of the Berlin Wall: Regional health differences in Germany,
#   Journal of Health Monitoring S2/2019, Robert Koch Institute

smokers <- data.frame(rbind(
  Women=c(
    17.4, 16.6, 21.3, 20.2, 24.2, 18.9, 17.5, 22.1,
    19.2, 19.4, 18.6, 17.7, 16.6, 20.0, 20.0, 21.7),
  Men=c(
    25.1, 24.6, 29.9, 29.0, 30.9, 27.8, 24.8, 33.4,
    26.7, 26.0, 24.9, 23.5, 26.5, 29.8, 27.1, 30.8),
  Total=c(
    21.2, 20.5, 25.5, 24.5, 27.4, 23.2, 21.1, 27.7,
    22.9, 22.6, 21.7, 20.6, 21.4, 24.8, 23.5, 26.2)))

# Federal States of Germany
state <- c(
  "BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV",
  "NI", "NW", "RP", "SL", "SN", "ST", "SH", "TH")

ew <- c(
  "w", "w", "b", "e", "w", "w", "w", "e",
  "w", "w", "w", "w", "e", "e", "w", "e")
names(ew) <- state

names(smokers) <- state

col <- c(e="red", w="blue", b="purple")
fg <- col[ew]
names(fg) <- state

plot(
  as.numeric(smokers[2,]), as.numeric(smokers[1,]), type="n",
  main="Smoking rate in Germany (2017)", xlab="Men [%]", ylab="Women [%]")
text(as.numeric(smokers[2,]), as.numeric(smokers[1,]), state, col=fg)
legend(
  "topleft", c("Former East Germany", "Former West Germany", "Berlin"),
  text.col=col)

library(tidyverse)
library(geofacet)

# library(cowplot)
# theme_set(theme_cowplot())

tb <- as_tibble(smokers, rownames="sex") %>%
  pivot_longer(-sex, names_to="state", values_to="rate") %>%
  mutate(sex=factor(sex, levels=c("Women", "Men", "Total")))

fg <- c(Women="#FF9999", Men="#9999FF", Total="lightgray")
bg <- c(e="#FFCCCC", w="#CCCCFF", b="#FFCCFF")[ew]
names(bg) <- state

ggplot(tb, aes(sex, rate)) +
  facet_geo(~ state, grid="de_states_grid1") +
  geom_rect(aes(fill=state), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_col(aes(fill=sex), color="black") +
  geom_text(aes(label=sprintf("%.1f", rate)), vjust=1.5) +
  scale_fill_manual("legend", values=c(fg, bg)) +
  labs(
    title = "Smoking rate in Germany (2017)", x = "Sex", y = "Smoking rate [%]",
    caption = "Data Source: Journal of Health Monitoring S2/2019, Robert Koch Institute") +
  theme(legend.position="none")

ggsave("smoker.png", width=12, height=9, dpi=96)
