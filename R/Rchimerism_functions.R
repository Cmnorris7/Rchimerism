# file to store all Rchimerism functions


#' Rchimerism locSD
#'
#' An internal function to determine informative loci for a single donor
#'
#'
#' @param ddata Donor data input text file
#' @param rdata Recipient data input text file
#' @param markers List of locus markers
#'
#' @return Internal variables used by chiSD.R

locSD <- function(ddata,rdata,markers) {
  #markers = c('D3S1358','TH01','D21S11','D18S51','Penta E','D5S818','D13S317','D7S820','D16S539','CSF1PO','Penta D','vWA','D8S1179','TPOX','FGA');
  
  #rData = read.delim('rdata.txt');
  #dData = read.delim('ddata.txt');
  
  # dData <- read.delim(ddata$datapath)
  # rData <- read.delim(rdata$datapath)
  
  # dData <- read.delim(ddata,sep = '\t')
  # rData <- read.delim(rdata,sep = '\t')
  
  dData <- data.frame(ddata)
  rData <- data.frame(rdata)
  
  
  #Function to handle invalid data text files
  coherent_input <- function(any_input) {
    if(ncol(any_input)<7 || nrow(any_input)<1) {
      return(paste("Cannot read ",deparse(substitute(any_input)),sep=""))
    }
    return(NULL)
  }
  
  ci <- c(coherent_input(dData),coherent_input(rData))
  ci_r <- ci[which(!is.null(ci))]
  
  if(!is.null(ci_r)) {
    return(ci_r)
  }
  
  
  #clean up the raw data (OL, X, '')
  # r = rData[grep("[^[:alnum:]]",rData[,4]),c(3:4,7)];
  # d = dData[grep("[^[:alnum:]]",dData[,4]),c(3:4,7)];
  
  r <- rData[,c(3:4,7)]
  d <- dData[,c(3:4,7)]
  
  
  rr = droplevels(r); #fantom levels can cause problem, rr, dd are tem variables
  dd = droplevels(d);
  
  
  # get rid of the noise (less than half of the maximum area)
  # maxR = tapply(rr[,3],rr[,1],max)/2;
  # r = rr[rr[,3]>maxR[rr[,1]],];
  # 
  # maxD = tapply(dd[,3],dd[,1],max)/2;
  # d = dd[dd[,3]>maxD[dd[,1]],];
  
  #Allele matrix and calculation;
  r$V4 = 'r';
  d$V4 = 'd';
  rd = rbind(r,d);
  
  
  #Compare rd with markers, end program if user defined markers not in input
  xtra_in_markers <- setdiff(markers,rd[,1])
  if (length(xtra_in_markers) != 0) {
    return(paste("'",xtra_in_markers,"'"," from markers not found in input data",
                 sep = ""))
  }
  
  trd = table(rd[,c(1,2,4)]);
  rt = trd[,,'r'];
  rt = rt[markers,];
  rt = rt[,sort(colnames(rt))];
  dt = trd[,,'d'];
  dt = dt[markers,];
  dt = dt[,sort(colnames(dt))];
  sum = dt + rt;
  diff = dt - rt;
  
  ###########
  #classify loci according to the matrix manipulation
  ###########
  
  profile = diff[,1];
  profile[apply(diff,1,any)] = 1; #informative locus marked as "1"
  profile[!apply(diff,1,any)] = 0; #non-informative locus marked as "0"
  
  # classify locus (ssum/sdiff: D allele number + R allele number): 4/0: 2+2; 3/1: 2+1; 3/-1: 1+2; 2/0: 1+1
  ssum = apply(sum,1,sum);
  sdiff = apply(diff,1,sum);
  
  
  # Further identify locus for three situations that use unique formulas
  # Label them uniquely
  
  l221 = apply(sum==1,1,any)& apply(sum==2,1,any)& ssum==4 & sdiff==0; #221: 2(donor)2(receipient)1(one shared)
  profile[l221] = 221;
  
  l121 = ssum==3 & sdiff==-1 & apply(sum==2,1,any); #121 1(donor)2(recipient)1(one shared)
  profile[l121] = 121;
  
  l211 = ssum==3 & sdiff==1 & apply(sum==2,1,any);
  profile[l211] = 211;
  profile=profile[markers];
  
  dm = cbind(dt,apply(dt,1,sum),profile);
  colnames(dm)[length(colnames(dm))-1] = 'Sum';
  colnames(dm)[length(colnames(dm))] = 'Profile';
  rm = cbind(rt,apply(rt,1,sum));
  colnames(rm)[length(colnames(rm))] = 'Sum';
  
  #save.image(file='data/locusSD.RData');
  
  
  # print("Donor Allele Matrix", quote = F);
  # print(dm);
  # print("Recipient Allele Matrix", quote = F);
  # print(rm);
  return(list(markers,profile,rt,dt,dm,rm,d,r))
}




#' Rchimerism locDD
#'
#' An internal function to determine informative loci for double donor cases
#'
#'
#' @param donor1_data Donor 1 data input text file
#' @param donor2_data Donor 2 data input text file
#' @param recipient_data Recipient data input text file
#' @param markers List of locus markers
#'
#'
#' @return Internal variables used by chiDD.R

# Output: profile,ru,rt,rnn,d1nn,d2nn,d1u,d2u,d1t,d2t,r,d1m,d2m,rm
#############################################



#markers = c('D3S1358','TH01','D21S11','D18S51','Penta E','D5S818','D13S317','D7S820','D16S539','CSF1PO','Penta D','vWA','D8S1179','TPOX','FGA');
locDD <- function(donor1_data, donor2_data,recipient_data,markers) {
  #rData = read.delim('rdata.txt');
  #d1Data = read.delim('d1data.txt');
  #d2Data = read.delim('d2data.txt');
  # ddata <- donor1_data
  # d2data <- donor2_data
  # rdata <- recipient_data
  #
  # d1Data <- read.delim(ddata$datapath)
  # d2Data <- read.delim(d2data$datapath)
  # rData <- read.delim(rdata$datapath)
  
  # d1Data <- read.delim(donor1_data, sep='\t')
  # d2Data <- read.delim(donor2_data, sep='\t')
  # rData <- read.delim(recipient_data, sep= '\t')
  
  d1Data <- data.frame(donor1_data)
  d2Data <- data.frame(donor2_data)
  rData <- data.frame(recipient_data)
  
  #Function to handle invalid data text files
  coherent_input <- function(any_input) {
    if(ncol(any_input)<7 || nrow(any_input)<1) {
      return(paste("Cannot read ",deparse(substitute(any_input)),sep=""))
    }
    return(NULL)
  }
  
  ci <- c(coherent_input(d1Data),coherent_input(d2Data),coherent_input(rData))
  ci_r <- ci[which(!is.null(ci))]
  
  if(!is.null(ci_r)) {
    return(ci_r)
  }
  
  
  #clean up the raw data (OL, X, '')
  r = rData[grep("[^[:alpha:]]",rData[,4]),c(3:4,6)];
  d1 = d1Data[grep("[^[:alpha:]]",d1Data[,4]),c(3:4,6)];
  d2 = d2Data[grep("[^[:alpha:]]",d2Data[,4]),c(3:4,6)];
  
  rr = droplevels(r); #fantom levels can cause problem, rr, dd are tem variables
  dd1 = droplevels(d1);
  dd2 = droplevels(d2);
  
  #get rid of the noise (less than half of the maximum height)
  maxR = tapply(rr[,3],rr[,1],max)/2;
  r = rr[rr[,3]>maxR[rr[,1]],];
  
  maxD1 = tapply(dd1[,3],dd1[,1],max)/2;
  d1 = dd1[dd1[,3]>maxD1[dd1[,1]],];
  
  maxD2 = tapply(dd2[,3],dd2[,1],max)/2;
  d2 = dd2[dd2[,3]>maxD2[dd2[,1]],];
  
  #Allele matrix (the matrix) and calculation;
  r[,4] = 'r';
  d1[,4] = 'd1';
  d2[,4] = 'd2';
  rd = rbind(r,d1,d2);
  
  #Compare rd with markers, end program if user defined markers not in input
  xtra_in_markers <- setdiff(markers,rd[,1])
  defic_in_markers <- setdiff(rd[,1],markers)
  if (length(xtra_in_markers) != 0) {
    return(paste("'",xtra_in_markers,"'"," from markers not found in input data",
                 sep = ""))
  }
  
  trd = table(rd[,c(1,2,4)]);
  rt = trd[,,'r'];
  rt = rt[markers,];
  rt = rt[,sort(colnames(rt))];
  d1t = trd[,,'d1'];
  d1t = d1t[markers,];
  d1t = d1t[,sort(colnames(d1t))];
  d2t = trd[,,'d2'];
  d2t = d2t[markers,];
  d2t = d2t[,sort(colnames(d2t))];
  
  ###########
  # Find the unique alleles in r, d1, and d2 (will be using them in chiDD.R to calculate percentage of r, d1, and d2
  ###########
  
  ru = rt-d1t-d2t;
  ru[ru!=1] = 0;
  run = apply(ru,1,sum);
  d1u = d1t-d2t-rt;
  d1u[d1u!=1] = 0;
  d1un = apply(d1u,1,sum);
  d2u = d2t-d1t-rt;
  d2u[d2u!=1] = 0;
  d2un = apply(d2u,1,sum);
  
  rn = apply(rt,1,sum);
  d1n = apply(d1t,1,sum);
  d2n = apply(d2t,1,sum);
  
  rur = sum(run)/sum(rn);
  d1ur = sum(d1un)/sum(d1n);
  d2ur = sum(d2un)/sum(d2n);
  
  # Find the heterozygous loci with one unique allele (will do a x 2 in chiDD.R)
  rnn = rn;
  d1nn = d1n;
  d2nn = d2n;
  
  rnn[!((rn==2)&(run==1))]=1;
  d1nn[!((d1n==2)&(d1un==1))]=1;
  d2nn[!((d2n==2)&(d2un==1))]=1;
  
  d1m = cbind(d1t,d1n,d1un,d1nn);
  colnames(d1m)[length(colnames(d1m))-2] = 'Sum';
  colnames(d1m)[length(colnames(d1m))-1] = 'Unique';
  colnames(d1m)[length(colnames(d1m))] = 'Factor';
  d2m = cbind(d2t,d2n,d2un,d2nn);
  colnames(d2m) = colnames(d1m);
  rm = cbind(rt,rn,run,rnn);
  colnames(rm) = colnames(d1m);
  
  #save.image(file='locusDD.RData');
  
  # print("Donor 1 Allele Matrix", quote=F);
  # print(d1m);
  # print("Donor 2 Allele Matrix", quote=F);
  # print(d2m);
  # print("Receipient Allele Matrix", quote = F);
  # print(rm);
  return(list(markers,profile,ru,rt,rnn,d1nn,d2nn,d1u,d2u,d1t,d2t,r,d1m,d2m,rm))
}



#' Rchimerism chiSD
#'
#' An internal function to determine donor percentages for a single donor
#'
#'
#' @param sdata Sample data input text file
#' @param markers List of locus markers
#' @param profile,rt,dt,d,r Internal variables for matrix operations
#'
#' @return A data frame with donor percentage results, and sample allele matrix
#'
#'


#sdata,markers,rt,dm,rm
#
chiSD <- function(sdata,markers,profile,rt,dt,d,r) {
  #load('data/locusSD.RData'); #load the previous workplace (r, d, markers, profile)
  
  
  #sData <- read.delim(sdata$datapath);
  # sData <- read.delim(sdata, sep= '\t')
  sData <- sdata
  
  #Function to handle invalid data text files
  coherent_input <- function(any_input) {
    if(ncol(any_input)<7 || nrow(any_input)<1) {
      return(paste("Cannot read ",deparse(substitute(any_input)),sep=""))
    }
    return(NULL)
  }
  
  ci <- c(coherent_input(sData))
  ci_r <- ci[which(!is.null(ci))]
  
  if(!is.null(ci_r)) {
    return(ci_r)
  }
  
  s = sData[grep("[^[:alpha:]]",sData[,4]),c(3:4,7)]; #clean up the raw data
  #s = s[s[,1]!='DYS391',] #Remove DYS391 locus on Y chromosome
  s = droplevels(s);
  s$Allele = as.factor(s$Allele);
  st =  rt;
  st[,] = 0
  
  #check point: any locus present in sample but not in recipient/donor
  # print("The following alleles are possible false calls by ABI GeneAnalyzer",quote = F)
  # print(s[!(s[,2] %in% colnames(st)),]);
  
  #check point: remove loci not specified in user markers csv file
  s = s[(s[,1] %in% rownames(st)),];
  
  #Handles noisy sample data, outputs table of possible false calls
  if (nrow(s[!(s[,2] %in% colnames(st)),]) != 0) {
    false_calls <- s[!(s[,2] %in% colnames(st)),]
    return(false_calls)
  }
  
  #check point: donor and receipient matrix
  # print("Donor Allele Matrix", quote = F);
  # print(dm);
  # print("Receipient Allele Matrix", quote = F);
  # print(rm);
  
  
  st[as.matrix(s[,1:2])] = 1;
  st[(dt+rt)==0 & st!=0]= 999; # the alleles not in donor or recipient
  sa = rt;
  sa[,] = 0;
  sa[as.matrix(s[,1:2])] = s[,3];
  
  # profile for sample, receipient, and donor
  #proS = rep(profile,table(s[,1])[markers]);
  #proR = rep(profile,table(r[,1])[markers]);
  #proD = rep(profile,table(d[,1])[markers]);
  
  C = profile; #chimerism matrix (donor percentage)
  
  # calculate percent donor chimerism based on general formula and specific formula
  for (m in markers){
    if (profile[m]==211){
      Ad = s[s[,1]==m & s[,2]==setdiff(d[d[,1]==m,2],r[r[,1]==m,2]),3];
      A = s[s[,1]==m & s[,2]==intersect(d[d[,1]==m,2],r[r[,1]==m,2]),3];
      
      if (length(Ad)==0){
        Ad = 0;
      }
      C[m]=2*Ad/(Ad+A);
    }
    
    if (profile[m]==221){
      Ad = s[s[,1]==m & s[,2]==setdiff(d[d[,1]==m,2],r[r[,1]==m,2]),3];
      Ar = s[s[,1]==m & s[,2]==setdiff(r[r[,1]==m,2],d[d[,1]==m,2]),3];
      
      if (length(Ar)==0){
        Ar = 0;
      }
      
      if (length(Ad)==0){
        Ad = 0;
      }
      C[m]=Ad/(Ad+Ar);
    }
    
    if (profile[m]==121){
      Ar = s[s[,1]==m & s[,2]==setdiff(r[r[,1]==m,2],d[d[,1]==m,2]),3];
      A = s[s[,1]==m & s[,2]==intersect(r[r[,1]==m,2],d[d[,1]==m,2]),3];
      
      if (length(Ar)==0){
        Ar = 0;
      }
      C[m]=1-(2*Ar/(Ar+A));
    }
    
    if (profile[m]==1){
      Ad = sum(s[s[,1]==m & s[,2] %in% setdiff(d[d[,1]==m,2],r[r[,1]==m,2]),3]);
      A = sum(s[s[,1]==m,3]);
      
      if (length(Ad)==0){
        Ad = 0;
      }
      C[m]=Ad/A
    }
    
    if (profile[m]==0){
      C[m]=NA
    }
    
  }
  
  #Generate final result data frame
  results = cbind(profile,C,NA,NA,NA,NA);
  results[,3] = NA;
  me = mean(results[,2],na.rm=T);
  sd = sd(results[,2],na.rm=T);
  l = !((abs(results[,2]-me) > 2*sd)| is.na(results[,2])); #non-informative locus or 2 SD locus are excluded
  results[l,3]=mean(results[l,2]);
  results[l,4]=sd(results[l,2]);
  results[l,5]=results[l,4]/results[l,3];
  results[l,6]=1-results[l,3];
  colnames(results)[1] = 'Profile';
  colnames(results)[2] = 'Donor%';
  colnames(results)[3] = 'Donor%_Mean';
  colnames(results)[4] = 'Donor%_SD';
  colnames(results)[5] = 'Donor%_CV';
  colnames(results)[6] = 'Recipient%_Mean';
  
  sm = cbind(st,apply(st,1,sum));
  colnames(sm)[length(colnames(sm))] = 'Sum';
  
  
  # print("Sample Allele Matrix", quote = F);
  # print(sm);
  # print("Sample Allele Area Matrix", quote = F);
  # print(sa);
  # print("Final Results", quote = F);
  # print(results);
  
  return(list(results,sm))
}



#' Rchimerism chiDD
#'
#' An internal function to determine donor percentages for a single donor
#'
#'
#' @param sdata Sample data input text file
#' @param markers List of locus markers
#' @param profile,ru,rt,rnn,d1nn,d2nn,d1u,d2u,d1t,d2t,r Internal variables for matrix operations
#'
#' @return A data frame with donor percentage results, and sample allele matrix
#'
#'


#load('locusDD.RData'); #load the previous workplace (r, d, markers, profile)
chiDD <- function(sdata,markers,profile,ru,rt,rnn,d1nn,d2nn,d1u,d2u,d1t,d2t,r) {
  print(getwd(),quote=F);
  
  #sData <- read.delim(sdata$datapath);
  # sData <- read.delim(sdata,sep = '\t')
  sData <- sdata
  
  #Function to handle invalid data text files
  coherent_input <- function(any_input) {
    if(ncol(any_input)<7 || nrow(any_input)<1) {
      return(paste("Cannot read ",deparse(substitute(any_input)),sep=""))
    }
    return(NULL)
  }
  
  ci <- c(coherent_input(sData))
  ci_r <- ci[which(!is.null(ci))]
  
  if(!is.null(ci_r)) {
    return(ci_r)
  }
  
  
  s = sData[grep("[^[:alpha:]]",sData[,4]),c(3:4,7)]; #clean up the raw data
  s = s[s[,1]!='DYS391',] #Remove DYS391 locus on Y chromosome
  s = droplevels(s);
  s$Allele = as.factor(s$Allele);
  
  # Fill in the matrix with alleles
  st = ru;
  st[,] = 0;
  
  
  #check point: any locus present in sample but not in recipient/donors
  #print("The following alleles are possible false calls by ABI GeneAnalyzer",quote = F)
  #print(s[!(s[,2] %in% colnames(st)),]);
  #browser()
  
  #check point: remove loci not specified in user markers csv file
  s = s[(s[,1] %in% rownames(st)),];
  
  #Handles noisy sample data, outputs table of possible false calls
  if (nrow(s[!(s[,2] %in% colnames(st)),]) != 0) {
    false_calls <- s[!(s[,2] %in% colnames(st)),]
    return(false_calls)
  }
  #chck point
  #print("Donor 1 Allele Matrix", quote=F);
  #print(d1m);
  #print("Donor 2 Allele Matrix", quote=F);
  #print(d2m);
  #print("Receipient Allele Matrix", quote = F);
  #print(rm);
  
  st[as.matrix(s[,1:2])] = 1;
  st[(d1t + d2t + rt) == 0 & st!=0] = 999; # the alleles not in donor or recipient
  sn = apply(st,1,sum);
  
  # Fill in the matrix with area values
  sa = ru;
  sa[,] = 0;
  sa[as.matrix(s[,1:2])] = s[,3];
  
  # Only keep the unique allele area values
  saru = sa;
  saru[ru==0] = 0;
  sad1u = sa;
  sad1u[d1u==0] = 0;
  sad2u = sa;
  sad2u[d2u==0] = 0;
  
  # Put NA in non-informative loci
  ss = apply(sa,1,sum);
  srs = apply(saru,1,sum);
  srs[apply(ru,1,sum)==0] = NA; #put NA in non-informative loci
  sd1s = apply(sad1u,1,sum);
  sd1s[apply(d1u,1,sum)==0] = NA; #put NA in non-informative loci
  sd2s = apply(sad2u,1,sum);
  sd2s[apply(d2u,1,sum)==0] = NA; #put NA in non-informative loci
  
  # Ratio/Percentage results
  sr = srs * rnn / ss;
  me = mean(sr,na.rm=T);
  sd = sd(sr,na.rm=T);
  l = !((abs(sr-me) > 2*sd)| is.na(sr)); #non-informative locus or 2 SD locus are excluded
  srMean = sr;
  srMean[!l] = NA;
  srMean[l] = mean(srMean,na.rm=T);
  
  srSD = sr;
  srSD[!l] = NA;
  srSD[l] = sd(sr,na.rm=T);
  srCV = srSD/srMean;
  
  sd1 = sd1s * d1nn / ss;
  me = mean(sd1,na.rm=T);
  sd = sd(sd1,na.rm=T);
  l = !((abs(sd1-me) > 2*sd)| is.na(sd1)); #non-informative locus or 2 SD locus are excluded
  sd1Mean = sd1;
  sd1Mean[!l] = NA;
  sd1Mean[l] = mean(sd1,na.rm=T);
  sd1SD = sd1;
  sd1SD[!l] = NA;
  sd1SD[l] = sd(sd1,na.rm=T);
  sd1CV = sd1SD/sd1Mean;
  
  
  sd2 = sd2s * d2nn / ss;
  me = mean(sd2,na.rm=T);
  sd = sd(sd2,na.rm=T);
  l = !((abs(sd2-me) > 2*sd)| is.na(sd2)); #non-informative locus or 2 SD locus are excluded
  sd2Mean = sd2;
  sd2Mean[!l] = NA;
  sd2Mean[l] = mean(sd2Mean,na.rm=T);
  sd2SD = sd2;
  sd2SD[!l] = NA;
  sd2SD[l] = sd(sd2,na.rm=T);
  sd2CV = sd2SD/sd2Mean;
  
  
  results = cbind(sd1,sd1Mean,sd1SD,sd1CV,sd2,sd2Mean,sd2SD,sd2CV,sr,srMean,srSD,srCV,sd1+sd2+sr);
  results = rbind(results, apply(!is.na(results),2,sum));
  rownames(results)[nrow(results)] = "Info#";
  colnames(results) = c('Donor_1%', 'Donor_1%_Mean','Donor_1%_SD','Donor_1%_CV','Donor_2%',
                        'Donor_2%_Mean','Donor_2%_SD','Donor_2%_CV','Recipient%', 'Recipient%_Mean',
                        'Recipient%_SD','Recipient%_CV',"Sum");
  sm = cbind(st,sn);
  colnames(sm)[length(colnames(sm))] = 'Sum';
  
  #write.table(results,file='results.xls',sep="\t",col.names=FALSE);
  
  
  
  # print("Sample Allele Matrix", quote = F);
  # print(sm);
  # print("Sample Allele Area", quote = F);
  # print(sa);
  # print("Final Results", quote = F);
  # print(results);
  return(list(results,sm))
}
