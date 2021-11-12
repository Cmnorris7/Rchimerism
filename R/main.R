setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('get_alleles.R')
source('Rchimerism_v1.1.R')
source('locSD.R')
source('chiSD.R')
source('locDD.R')
source('chiDD.R')

# install.packages('svDialogs')
# install.packages('sqldf')
# install.packages('data.table')


library(svDialogs)

# garbage collection to clear up memory
gc()

# markers used in chimerism assay. These are static.
markers = c('D3S1358','TH01','D21S11','D18S51','Penta E','D5S818','D13S317','D7S820','D16S539','CSF1PO','Penta D','vWA','D8S1179','TPOX','FGA');

# Prompt the user to select workflow
analysis <- dlg_list(c('Donor Analysis (Pre)','Recipient Analysis (Pre)','Single Donor (Post)','Multidonor (Post)'),preselect = NULL, multiple = FALSE, rstudio=FALSE, title = "Chimerism Analysis",gui = .GUI)$res


pre_analysis <- function(){
  # input = donor peak file
  # output = data table peak calls. will send to sql
} #2nd


single_donor <- function(){
  
  ddata <- clean_pre_file()
  rdata <- clean_pre_file()

  sdata <- get_informative_marks_sd(ddata, rdata, )

  loci <- locSD(ddata, rdata, markers)

  profile <- loci[[2]]
  rt <- loci[[3]]
  dt <- loci[[4]]
  dm <- loci[[5]]
  rm <- loci[[6]]
  d <- loci[[7]]
  r <- loci[[8]]
  
  percents <- chiSD(sdata,markers,profile,rt,dt,d,r)
  
  
  results <- percents[[1]]
  print(results)
  return(results)
} 


multi_donor <- function(){
  
  d1data <- clean_pre_file()
  d2data <- clean_pre_file()
  rdata <- clean_pre_file()
  sdata <- get_informative_marks_dd(d1data, d2data, rdata, )
  
  loc_dd_output <- locDD(d1data, d2data, rdata, markers)
  
  profile <- loc_dd_output[[2]]
  ru <- loc_dd_output[[3]]
  rt <- loc_dd_output[[4]]
  rnn <- loc_dd_output[[5]]
  d1nn <- loc_dd_output[[6]]
  d2nn <- loc_dd_output[[7]]
  d1u <- loc_dd_output[[8]]
  d2u <- loc_dd_output[[9]]
  d1t <- loc_dd_output[[10]]
  d2t <- loc_dd_output[[11]]
  r <- loc_dd_output[[12]]
  d1m <- loc_dd_output[[13]]
  d2m <- loc_dd_output[[14]]
  rm <- loc_dd_output[[15]]
  
  chi_dd_output <-  chiDD(sdata,markers,profile,
                          ru,rt,rnn,d1nn,d2nn,d1u,d2u,d1t,d2t,r)
  
  results <- chi_dd_output[[1]]

  print(results[,c(1:3,5:7,9:11)])
  return(results)
} #1st


pre_to_sql <- function(){
  # input = data table from donor/recip analysis
  # output = updated sql file
} #3rd


fetch_pre_sql <- function(){
  # input = database file
  # output = temp pre files
} #4th




if (analysis == "Single Donor (Post)"){
  single_donor()
} else if (analysis == "Multidonor (Post)"){
  multi_donor()
} else if (analysis == "Donor Analysis (Pre)"){
  print("Donor analysis workflow")
} else if (analysis == "Recipient Analysis (Pre)"){
  print("Recipient analysis workflow")
}



