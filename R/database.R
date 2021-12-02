# file to store all databasing functions

library(DBI)

pre_to_sql <- function(clean_alleles){
  m <- dbDriver("SQLite")
  conn <- dbConnect(m,dbname = 'chim_test.db')
  dbWriteTable(conn, "peak_report", clean_alleles, append = TRUE)
  dbDisconnect(conn)
}

# post_to_sql <- function(analysis_table){
#   m <- dbDriver("SQLite")
#   conn <- dbConnect(m,dbname = 'chim_test.db')
#   dbWriteTable(conn, )
# }