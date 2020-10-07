indentified <-  read.csv("[pos] Identification Result.csv")
database <- read.csv("wuxiaozhu-RXQ-POS.csv")
colnames(database)[1] <- "Compound.name"


#  tobe indentified
msDa <- database
msDa$Compound.name <- as.vector(msDa$Compound.name)


#  database file importing
moDa <- indentified
moDa$Compound.name <- as.vector(moDa$Compound.name)

#  need to be indentified
msDaNA <- msDa[is.na(msDa$Compound.name),]

#  already indentified
msDaID <- msDa[!is.na(msDa$Compound.name),]

#  match mz and rt with database
for (i in 1:nrow(msDaNA)) {
  
  msDaNA %>% slice(i) %>% select(mz) %>% as.numeric() -> mz.candi
  msDaNA %>% slice(i) %>% select(rt) %>% as.numeric() -> rt.candi
  
  moDa %>% filter(near(mz, mz.candi, tol = 0.01) &      # set the mz tolerance             
                  near(rt, rt.candi, tol = 10))  %>%    # set the rt tolerance
    mutate(deltmz = abs(mz - mz.candi),
           deltrt = abs(rt - rt.candi)) %>% 
    arrange(.,deltmz) %>%                 # choose  smaller deltmz
    slice(1) %>% select(Compound.name, MS2.score) -> 
    candidata
  
  cat(i," : ",candidata[1,1] %>% as.character(), 
      "  MS2.score =", candidata[1,2] %>% as.numeric(), "\n")
  msDaNA$Compound.name [i] <-  candidata[1,1] %>% as.character() 
  msDaNA$MS2.score [i] <-  candidata[1,2] %>% as.numeric() 
  
}

msDaID %>% mutate(database = "original") %>%
  select(database,everything()) -> msOriginal
msDaNA[!is.na(msDaNA$Compound.name),] %>% mutate(database = "updated")%>%
  select(database,everything()) -> msNew
msDaNA[is.na(msDaNA$Compound.name),] %>% mutate(database = "NA")%>%
  select(database,everything()) -> msNA

nummsNew <- nrow(msNew)
nummsOriginal <- nrow(msOriginal)
nummsNA <- nrow(msNA)
numALL = nummsNew + nummsOriginal

cat("\n","Totally, ",numALL, " metabolites annotated", "\n")
cat("\n","There are new ",nummsNew, " metabolites annotated", "\n")
cat("\n","There are still ",nummsNA, " metabolites unknown", "\n")
Bindresult <- rbind(msOriginal,msNew,msNA)
write.csv(Bindresult,"POS.Bindresult.csv", row.names = F)
