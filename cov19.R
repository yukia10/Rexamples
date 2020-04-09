# cov19.R - COVID-19 deaths in Germany

# Data source:
#   https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland
#   https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Archiv.html

### Scraping
Sys.setlocale("LC_TIME", "de_DE.UTF-8")

library(rvest)
url <- "https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland"
# url <- "COVID-19-Pandemie in Deutschland – Wikipedia.html"
tbl <- url %>% read_html %>% html_nodes("table")

conv_tbl <- function(tbl) {
  t <- html_table(tbl, dec=".")[[1]]
  t <- t[-nrow(t), 1:17]
  t <- data.frame(lapply(t, function(x) gsub("[(].*[)]", "", x)), stringsAsFactors=FALSE)
  t <- data.frame(lapply(t, function(x) gsub("[.]", "", x)), stringsAsFactors=FALSE)
  t <- data.frame(lapply(t, function(x) gsub("[\u2013\u2014]", "0", x)), stringsAsFactors=FALSE)
  t[, 1] <- gsub(".*♠", "", t[, 1])
  data.frame(
    lapply(t[, -1], as.numeric),
    row.names=as.character(strptime(gsub("\u00a0", " ", t[, 1]), "%d %b %Y")))
}

t1 <- conv_tbl(tbl[4]) # Cumulative infections
t2 <- conv_tbl(tbl[5]) # Cumulative infections (per 100,000 inhabitants)
t3 <- conv_tbl(tbl[6]) # Cumulative deaths

Sys.setlocale("LC_TIME", "en_US.UTF-8")
stopifnot(
  row.names(t1)[nrow(t1)] == row.names(t2)[nrow(t2)],
  all(colnames(t1) == colnames(t2)), all(colnames(t1) == colnames(t3)))

### Analysis
library(growthcurver)

pop <- as.numeric(100000 * t1[nrow(t1), ] / t2[nrow(t2), ]) # larger number will make less erroneous
t4 <- data.frame(t(t(as.matrix(t3)) / pop * 100000)) # Cumulative deaths (per 100,000 inhabitants)

L <- lapply(t4, function(x) SummarizeGrowth(seq_along(rownames(t4)), x, bg_correct="none"))
t5 <- data.frame(sapply(L, function(x) predict(x$model))) # Fitted
row.names(t5) <- row.names(t4)

### Plot

# Federal States of Germany
state <- c(
  "BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV",
  "NI", "NW", "RP", "SL", "SN", "ST", "SH", "TH")

ew <- c(
  "w", "w", "b", "e", "w", "w", "w", "e",
  "w", "w", "w", "w", "e", "e", "w", "e")
names(ew) <- state

col <-c(e="red", w="blue", b="purple")
fg <- col[ew]
names(fg) <- state

x <- as.POSIXct(row.names(t4))
matplot(
  x, t4, xaxt="n", type="l", lty=1, col=fg,
  main="COVID-19 deaths in Germany",
  xlab="Date", ylab="Deaths per 100,000 inhabitants")
axis.POSIXct(
  side=1, at=seq(x[1], x[length(x)], "day"), las=2, format="%m/%d")
text(x[length(x)], t4[nrow(t4), ], names(t4), col=fg, adj=-0.2, cex=0.8)
legend(
  "topleft", c("Former East Germany", "Former West Germany", "Berlin"),
  text.col=col)

library(tidyverse)
library(geofacet)

# library(cowplot)
# theme_set(theme_cowplot())

td <- as_tibble(t4, rownames="date") %>%
  pivot_longer(-date, names_to="state", values_to="death")
tf <- as_tibble(t5, rownames="date") %>%
  pivot_longer(-date, names_to="state", values_to="fitted")
td <- inner_join(td, tf, by=c("date", "state"))

bg <- c(e="#FFCCCC", w="#CCCCFF", b="#FFCCFF")[ew]
names(bg) <- state

ggplot(td, aes(as.Date(date), death)) +
  facet_geo(~ state, grid="de_states_grid1") +
  geom_rect(aes(fill=state), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape=1) +
  geom_line(aes(as.Date(date), fitted), color="red") +
  scale_x_date(date_labels="%b %d", date_breaks="2 weeks") +
  scale_fill_manual("legend", values=bg) +
  labs(
    title = sprintf(
      "COVID-19 deaths in Germany (As of %s)",
      rownames(t4)[nrow(t4)]),
    caption = "Data Source: Wikipedia based on Robert Koch Institut",
    x = "Date",
    y = "Deaths per 100,000 inhabitants") +
  geom_text(
    x=-Inf, y=Inf, aes(label=sprintf("%.2f", t4[nrow(t4), state])),
    vjust=1.2, hjust=-0.2, size=5, check_overlap=TRUE) +
  theme(legend.position="none")

ggsave("cov19.png", dpi=100)
