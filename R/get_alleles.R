library(sqldf)
library(data.table)
library(svDialogs)



# Take in peak sizing file (donor_pre or recipient_pre) and return true allele calls
clean_pre_file <- function(file){
  if(missing(file)){
    raw_file <- dlg_open(
      'S:\\UHTL\\3130\\Molecular Lab Data\\Chimerism\\*',
      'Select Donor Peak Report',
      multiple = FALSE,
      filters = dlg_filters["All",],
      gui = .GUI)
    peak_table <- fread(raw_file$res, sep = '\t', header=TRUE, na.strings=c("","NA"))
  } else{
    peak_table <- fread(file, sep = '\t', header= TRUE, na.strings=c("", "NA"))
  }
  clean_peak_table <- peak_table[!is.na(peak_table$Marker),]
  markers <- unique(clean_peak_table$Marker)
  
  # Remove OL alleles and select 2 entries for each marker that have the highest peaks
  clean_peak_table <- clean_peak_table[Allele != "OL"]
  clean_peak_table <- clean_peak_table[Allele != "?"]
  allele_table <- clean_peak_table[order(Marker,-Height)]
  allele_table <- data.table(allele_table, key = "Marker")
  allele_table <- allele_table[ , head(.SD, 2), by="Marker"]
  setcolorder(allele_table,c("Dye/Sample Peak", "Sample File Name", "Marker",
                                 "Allele", "Size", "Height", "Area", "Data Point",
                                 "V9"))
  
  # remove 2nd peak if the sample is homozygous for that marker. cutoff = 0.7
  for (marker in markers){
    alleles <- allele_table[Marker == marker]
    if(length(alleles$Height)>1 & alleles[2,6]/alleles[1,6] < 0.7){
      allele_table <- allele_table[!alleles[2,], on=.(Marker,Allele,Height)]
    }
  }
  
  # remove V9 column. Not sure where this comes from?
  # allele_table <- allele_table[,V9:=NULL]
  
  write.table(allele_table, file=paste(raw_file$res,"_TEMP.txt"), quote=FALSE, sep='\t', row.names = FALSE)
  print(is.data.frame(allele_table))
  return(allele_table)
  # return(paste(raw_file$res,"_TEMP.txt"))
}


# Take in 3 files and remove unwanted peaks from sample file
get_informative_marks_sd <- function(donor_tab, recipient_tab, sample_file){
  
  
  # if(missing(donor_tab)){
  #   donor_tab <- file.choose(new = FALSE)
  # }
  # if(missing(recipient_tab)){
  #   recipient_tab <- file.choose(new = FALSE)
  # }
  if(missing(sample_file)){
    sample_file <- file.choose(new = FALSE)
  }
  
  # read files into tables
  # donor_table <- fread(donor_tab, sep = '\t', header= TRUE)
  donor_table <- data.table(donor_tab)
  # recipient_table <- fread(recipient_tab, sep = '\t', header= TRUE)
  recipient_table <- data.table(recipient_tab)
  sample_table <- fread(sample_file, sep = '\t', header= TRUE, na.strings=c("", "NA"))
  sample_table <- sample_table[!is.na(sample_table$Marker)]
  
  # create single table with all markers/alleles from donor & recipient files
  columns <- c('Marker', 'Allele')
  donor_alleles <- donor_table[, columns, with=FALSE]
  recipient_alleles <- recipient_table[, columns, with=FALSE]
  alleles_of_interest <- rbind(donor_alleles, recipient_alleles)
  alleles_of_interest <- alleles_of_interest[!duplicated(alleles_of_interest)]
  
  # remove all markers/alleles from sample file that don't exist in pre donor/recip files
  final_table <- data.table()
  count <- 0
  for (row1 in 1:nrow(sample_table)){
    marker <- sample_table[row1,]$Marker
    allele <- sample_table[row1,]$Allele
    for (row2 in 1:nrow(alleles_of_interest)){
      if((marker == alleles_of_interest[row2,]$Marker & 
          allele == alleles_of_interest[row2,]$Allele)){
        final_table <- rbind(final_table,sample_table[Marker==marker & Allele == allele])
        break
      }
    }
  }
  # print(final_table)
  # write.table(final_table, file=paste(sample_file,"_TEMP.txt"), quote=FALSE, sep='\t', row.names = FALSE)
  return(data.frame(final_table))
  # return(paste(sample_file,"_TEMP.txt"))
}


