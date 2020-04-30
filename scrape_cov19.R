# scrape_cov19.R - Scraping COVID-19 deaths in Germany

# Data source:
#   https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland
#   https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/Archiv.html

library(rvest)

# Month name in Germany
month.de <- c(
  "Jan", "Feb", "M\u00e4r", "Apr", "Mai", "Jun", 
  "Jul", "Aug", "Sep", "Okt", "Nov", "Dez")

# url <- "https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland"
url <- "doc/COVID-19-Pandemie in Deutschland – Wikipedia.html"
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

stopifnot(
  all(colnames(tbl[[1]]) == colnames(tbl[[2]])), 
  all(colnames(tbl[[1]]) == colnames(tbl[[3]])))

tbl[[3]]["2020-04-22", "TH"] <- 61 # 71 -> 61, force monotone, something was wrong

k <- tail(intersect(row.names(tbl[[1]]), row.names(tbl[[2]])), 1) # latest common date
pop <- as.integer(100000 * tbl[[1]][k, ] / tbl[[2]][k, ]) # larger number will make less erroneous

tbl[[2]] <- data.frame(t(t(tbl[[1]]) / pop * 100000)) # Cum. infections / 100,000 pop.
tbl[[4]] <- data.frame(t(t(tbl[[3]]) / pop * 100000)) # Cum. deaths / 100,000 pop.

save(tbl, file="cov19de.RData")
