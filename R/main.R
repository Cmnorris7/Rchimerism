source('P:\\Chimerism\\UHTL_build\\get_alleles.R')
source('P:\\Chimerism\\UHTL_build\\Rchimerism.R')


markers = c('D3S1358','TH01','D21S11','D18S51','Penta E','D5S818','D13S317','D7S820','D16S539','CSF1PO','Penta D','vWA','D8S1179','TPOX','FGA');


ddata <- clean_pre_file()
rdata <- clean_pre_file()
sdata <- get_informative_marks(ddata, rdata, )

print(ddata)
print(rdata)
print(sdata)



test <- locSD1(ddata, rdata, markers)
print(test)
profile <- test[[2]]
rt <- test[[3]]
dt <- test[[4]]
dm <- test[[5]]
rm <- test[[6]]
d <- test[[7]]
r <- test[[8]]

# test2 <- chiSD(sdata,markers,profile,rt,dt,d,r)

# results <- test2[[1]]
# print(results[,1:3])

dDataTEST <- read.delim('S:\\UHTL\\3130\\Molecular Lab Data\\Chimerism\\2021\\3137GS101821_TEST_CAMERON\\5814021489-HaidetKenneth-Pos-BM_1.txt')
# print(dData)
