# cov19.R - COVID-19 deaths in Germany

# Data source:
#   https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland
#   https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Archiv.html

# Federal States of Germany
state <- c(
  "BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV",
  "NI", "NW", "RP", "SL", "SN", "ST", "SH", "TH")

# Former East Germany (GDR), Former West Germany (BRD), Berlin
ew <- c(
  "w", "w", "b", "e", "w", "w", "w", "e",
  "w", "w", "w", "w", "e", "e", "w", "e")
names(ew) <- state

# Month name in Germany
month.de <- c(
  "Jan", "Feb", "M\u00e4r", "Apr", "Mai", "Jun", 
  "Jul", "Aug", "Sep", "Okt", "Nov", "Dez")

### Scraping
library(rvest)

# url <- "https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland"
url <- "COVID-19-Pandemie in Deutschland – Wikipedia.html"
tbl <- url %>% read_html %>% html_nodes("table")

conv_tbl <- function(tbl) {
  t <- html_table(tbl, dec=".")[[1]]
  t <- t[-nrow(t), 1:17]
  t <- data.frame(lapply(
    t, function(x) gsub("[(].*[)]", "", x)), stringsAsFactors=FALSE)
  t <- data.frame(lapply(
    t, function(x) gsub("[.]", "", x)), stringsAsFactors=FALSE)
  t <- data.frame(lapply(
    t, function(x) gsub("[\u2013\u2014]", "0", x)), stringsAsFactors=FALSE) # dash to 0
  t[, 1] <- gsub(".*\u2660", "", t[, 1])
  t[, 1] <- gsub("\u00a0", " ", t[, 1]) # replace space
  t[, 1] <- sapply(strsplit(t[, 1], " "), function(x) 
    sprintf("%s-%02d-%s", x[3], match(x[2], month.de), x[1]))
  data.frame(lapply(t[, -1], as.numeric), row.names=t[, 1])
}

t1 <- conv_tbl(tbl[4]) # Cum. infections
t2 <- conv_tbl(tbl[5]) # Cum. infections / 100,000 pop.
t3 <- conv_tbl(tbl[6]) # Cum. deaths

stopifnot(all(colnames(t2) == colnames(t3)), all(colnames(t1) == colnames(t3)))

k <- tail(intersect(row.names(t1), row.names(t2)), 1) # latest common date
pop <- as.numeric(100000 * t1[k, ] / t2[k, ]) # larger number will make less erroneous
t4 <- data.frame(t(t(as.matrix(t3)) / pop * 100000)) # Cum. deaths / 100,000 pop.

### Analysis
library(growthcurver)

L <- lapply(t4, function(x) SummarizeGrowth(seq_along(rownames(t4)), x, bg_correct="none"))
t5 <- data.frame(sapply(L, function(x) predict(x$model))) # Fitted
row.names(t5) <- row.names(t4)

### Plot
library(ggplot2)
library(cowplot)

ggplot(stack(t4), aes(rep(as.Date(row.names(t4)), ncol(t4)), y=values, group=ind)) +
  geom_line(aes(color=ew[ind])) + 
  geom_text(
    data=stack(t4[nrow(t4),]), 
    aes(x=as.Date(tail(row.names(t4), 1)), y=values, label=ind, color=ew[ind]), 
    hjust=-0.2) +
  scale_color_manual(values=c(e="red", w="blue", b="purple")) +
  scale_x_date(date_labels="%m/%d") +
  labs(
    title = sprintf(
      "COVID-19 deaths in Germany (As of %s)", rownames(t4)[nrow(t4)]),
    caption = "Data Source: Wikipedia based on Robert Koch Institut",
    x = "Date", y = "Deaths / 100,000 pop.") +
  theme_cowplot() +
  theme(legend.position="none")


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
  scale_x_date(date_labels="%m/%d", date_breaks="2 weeks") +
  scale_fill_manual("legend", values=bg) +
  labs(
    title = sprintf(
      "COVID-19 deaths in Germany (As of %s)",
      rownames(t4)[nrow(t4)]),
    caption = "Data Source: Wikipedia based on Robert Koch Institut",
    x = "Date",
    y = "Deaths / 100,000 pop.") +
  geom_text(
    x=-Inf, y=Inf, aes(label=sprintf("%.2f", t4[nrow(t4), state])),
    vjust=1.2, hjust=-0.2, size=5, check_overlap=TRUE) +
  theme(legend.position="none")

ggsave("cov19.png", dpi=100)
