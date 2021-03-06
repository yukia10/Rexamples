# scrape_cov19.R - Scraping COVID-19 deaths in Germany

# Data source:
#   https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland
#   https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Archiv.html

# Month name in Germany
month.de <- c(
  "Jan", "Feb", "M\u00e4r", "Apr", "Mai", "Juni", 
  "Jul", "Aug", "Sep", "Okt", "Nov", "Dez")

# Federal States of Germany
state <- c(
  "BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV",
  "NI", "NW", "RP", "SL", "SN", "ST", "SH", "TH")

# Former East Germany (GDR), Former West Germany (BRD), Berlin
ew <- c(
  "w", "w", "b", "e", "w", "w", "w", "e",
  "w", "w", "w", "w", "e", "e", "w", "e")
names(ew) <- state

library(rvest)

# url <- "https://de.wikipedia.org/wiki/COVID-19-Pandemie/Statistik"
url <- "doc/COVID-19-Pandemie_Statistik.html"
tbl <- url %>% read_html %>% html_nodes("table")

nc <- sapply(seq_along(tbl), function(i) {
  t <- try(html_table(tbl[i], dec=".")[[1]], silent=TRUE)
  ifelse(inherits(t, "try-error"), 0, ncol(t))})
tbl <- tbl[18 <= nc]
stopifnot(length(tbl) == 3)

conv_tbl <- function(tbl) {
  t <- html_table(tbl, dec=".")[[1]]
  t <- t[-nrow(t), 1:17]
  t <- data.frame(lapply(
    t, function(x) gsub("[(].*[)]", "", x)), stringsAsFactors=FALSE)
  t <- data.frame(lapply(
    t, function(x) gsub("\\[.*\\]", "", x)), stringsAsFactors=FALSE)
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

tbl <- lapply(seq_along(tbl), function(i) conv_tbl(tbl[i]))
# tbl[[1]]: Cum. infections
# tbl[[2]]: Cum. infections / 100,000 pop.
# tbl[[3]]: Cum. deaths

for (i in seq_along(tbl))
  if (colnames(tbl[[i]])[10] == "NRW")
    colnames(tbl[[i]])[10] <- "NW"

stopifnot(
  all(state == colnames(tbl[[1]])), 
  all(state == colnames(tbl[[2]])), 
  all(state == colnames(tbl[[3]])))

retract <- function(x) { # still assuming carefully curated cumulative data
  for (k in which(c(0, diff(x)) < 0)) {
    y <- x[seq_len(k - 1)] # x[seq_len(k)] may be another choice
    l <- which(y > x[k])
    x[l] <- max(0, y[-l])
  }
  x
}

tbl <- lapply(tbl, function(t) 
  as.data.frame(lapply(t, retract), row.names=row.names(t)))

# tbl[[3]]["2020-04-22", "TH"] <- 61 # 71 -> 61, force monotone, something was wrong

k <- tail(intersect(row.names(tbl[[1]]), row.names(tbl[[2]])), 1) # latest common date
pop <- as.integer(100000 * tbl[[1]][k, ] / tbl[[2]][k, ]) # larger number will make less erroneous
names(pop) <- state

tbl[[2]] <- data.frame(t(t(tbl[[1]]) / pop * 100000)) # Cum. infections / 100,000 pop.
tbl[[4]] <- data.frame(t(t(tbl[[3]]) / pop * 100000)) # Cum. deaths / 100,000 pop.

save(tbl, pop, state, ew, file="cov19de.RData")
