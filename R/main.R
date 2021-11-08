source('P:\\Rchimerism_github\\R\\get_alleles.R')
source('P:\\Rchimerism_github\\R\\Rchimerism_v1.1.R')
source('P:\\Rchimerism_github\\R\\locSD.R')
source('P:\\Rchimerism_github\\R\\chiSD.R')

getwd()
markers = c('D3S1358','TH01','D21S11','D18S51','Penta E','D5S818','D13S317','D7S820','D16S539','CSF1PO','Penta D','vWA','D8S1179','TPOX','FGA');


ddata <- clean_pre_file()
rdata <- clean_pre_file()
sdata <- get_informative_marks(ddata, rdata, )





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


