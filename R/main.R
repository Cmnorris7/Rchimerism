source('P:\\Rchimerism_github\\R\\get_alleles.R')
source('P:\\Rchimerism_github\\R\\Rchimerism_v1.1.R')
source('P:\\Rchimerism_github\\R\\locSD.R')
source('P:\\Rchimerism_github\\R\\chiSD.R')
library(svDialogs)
getwd()
markers = c('D3S1358','TH01','D21S11','D18S51','Penta E','D5S818','D13S317','D7S820','D16S539','CSF1PO','Penta D','vWA','D8S1179','TPOX','FGA');

analysis <- dlg_list(c('Donor Analysis (Pre)','Recipient Analysis (Pre)','Single Donor (Post)','Multidonor (Post)'),preselect = NULL, multiple = FALSE, title = "Chimerism Analysis",gui = .GUI)$res


single_donor <- function(){
  
  # get donor file
  dlg_open(
    'S:\\UHTL\\3130\\Molecular Lab Data\\Chimerism\\*',
    'Select Donor Peak Report',
    multiple = FALSE,
    filters = dlg_filters["All",],
    gui = .GUI)$donor_file
  
  # get recipient file
  dlg_open(
    'S:\\UHTL\\3130\\Molecular Lab Data\\Chimerism\\*',
    'Select Recipient Peak Report',
    multiple = FALSE,
    filters = dlg_filters["All",],
    gui = .GUI)$recipient_file
  
  print(donor_file)
  ddata <- clean_pre_file(donor_file)
  rdata <- clean_pre_file(recipient_file)
  sdata <- get_informative_marks_sd(ddata, rdata, )
  
  test <- locSD(ddata, rdata, markers)
  
  profile <- test[[2]]
  rt <- test[[3]]
  dt <- test[[4]]
  dm <- test[[5]]
  rm <- test[[6]]
  d <- test[[7]]
  r <- test[[8]]
  
  test2 <- chiSD(sdata,markers,profile,rt,dt,d,r)
  print(test2)
  
  results <- test2[[1]]
  print(results[,1:3])
  return(results)
}


multi_donor <- function(){
  
  d1data <- clean_pre_file()
  d2data <- clean_pre_file()
  rdata <- clean_pre_file()
  sdata <- get_informative_marks_dd(d1data, d2data, rdata, )
  
  test <- locSD(ddata, rdata, markers)
  
  profile <- test[[2]]
  rt <- test[[3]]
  dt <- test[[4]]
  dm <- test[[5]]
  rm <- test[[6]]
  d <- test[[7]]
  r <- test[[8]]
  
  test2 <- chiSD(sdata,markers,profile,rt,dt,d,r)
  print(test2)
  
  results <- test2[[1]]
  print(results[,1:3])
  return(results)
}

if (analysis == "Single Donor (Post)"){
  single_donor()
}


