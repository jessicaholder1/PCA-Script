#new tutorial to find other principle components
#install.packages("ISLR")
#library("ISLR")
#provides data for "An introduction to statistical learning with applications in R" book so I shouldn't need it to run PCA
numerical_data <- df_prostatedata_all_genes[,4:60486]
#the input for prcomp is numerical data (df_prostate_all_genes had columns like case ID that weren't numerical so I had to remove them)
#tried inputting into prcomp and got an error :( 
#pc <- prcomp(numerical_data, scale = TRUE)
#Error in prcomp.default(numerical_data, scale = TRUE) : cannot rescale a constant/zero column to unit variance

#someone said you can just get rid of those values so i am going to try that 

##identify constant or zero columns
df_PD_all_genes_const_0 <- sapply(df_prostatedata_all_genes, function(x) is.atomic(x) && length(unique(x)) == 1)
######df_PD_all_genes_const_0_test <- sapply (df_prostatedata_all_genes,function(x) length(unique(x)) == 1)
######returns the same thing so the is.atomic is not necessary probably just checking to make sure that the data is an atomic type, but the function(x) is a part of hte sapply command 
#sapply applies a function over and over to a list,vector, or data frame (kind of like looping in other languages) and outputs a vector or matrix
#https://www.guru99.com/r-apply-sapply-tapply.html 
#length(unique(x) == 1) find the values that are unique and give me the ones that have only 1 unique value

#outputs a table with true (this column has constant/zero values) and false (there is more than one value in this column) statements
#so now we can remove the genes where the value is true to fix the error

df_PD_all_genes <- df_prostatedata_all_genes[,!df_PD_all_genes_const_0] ##process called subsetting
#instead of 60486 genes we now have 56968 genes

#again we need just numeric values so only keep the columns with the genes (get rid of patientID, caseID, and genotype)
df_PD_all_genes <- df_PD_all_genes[,4:56968]

#now we can run prcomp
pc <- prcomp(df_PD_all_genes, scale = TRUE)
summary(pc)
#gives us the principle components with their associated standard deviation and proportion of variance values

##############---now let's try plotting it---#################
#trying to add color that reflects genotype
#the df_prostatedata_genotype df has genotype 1 which we aren't interested in right now so I got rid of genotype 1 
genotype_0_2 <- df_prostatedata_genotype$genotype!=1 #returns a TRUE or FALSE table (TRUE if the genotype is not 1, FALSE if the genotype is 1)
#this table basically gives back the samples with genotype 0 and 2 
df_prostatedata_genotype <- df_prostatedata_genotype[genotype_0_2,] #keep the rows where the genotype is 0 or 2 aka a TRUE row in genotype_0_2
#419 samples with NA values

#remove all NA values
df_prostatedata_genotype <- na.omit(df_prostatedata_genotype)
#352 samples (2 columns: patient Id and genotype)

##creating a variable with all the genotypes
labs_genotype <- df_prostatedata_genotype$genotype

#declare function to map genotype to a different color
Cols=function(x){
  Cols=rainbow (length (unique(x)))
  return(Cols[as.numeric (as.factor(x))])}
#Cols=rainbow assigns the color
##as.factor is calling the genotype of each patient (either 0 or 2)
#as.numeric is assigning each genotype (0 or 2) to a number either 1 (=genotype 0) or 2 (= genotype 2)

######---generate plots---######
#PC1 vs. PC2 
plot(pc$x[,1:2], col = Cols(labs_genotype), pch =19, xlab ="PC1",ylab="PC2")
#pc variable is where the prcomp data is stored (all of the PC) calling the 

#plotting different principle components
plot(pc$x[,3:4], col = Cols(labs_genotype), pch =19, xlab ="PC3",ylab="PC4")
plot(pc$x[,5:6], col = Cols(labs_genotype), pch =19, xlab ="P5",ylab="PC6")
plot(pc$x[,7:8], col = Cols(labs_genotype), pch =19, xlab ="PC7",ylab="PC8")
plot(pc$x[,9:10], col = Cols(labs_genotype), pch =19, xlab ="PC9",ylab="PC10")

#how many of each genotype do we have
summary(as.factor(labs_genotype))
#matching color to genotype
#genotype 0 is red (less samples) and genotype 2 is blue (more samples)

#############---downloading clinical data from val data2/han_lab/mitra/clinical.tsv---###########
clinical_subset <- read.table("clinical_subset.tsv", header = TRUE, sep = "\t")
#just the case_submitter_id and self-reported race columns
##for some reason there are two of each sample (duplicates) (978 rows/samples) which will mess up my data later in DESeq2 so I will get rid of those duplicates now

##applying the unique command to the dataframe will remove any duplicates (it  returns an object of the same class as x (typically a vector-like, data-frame-like, or array-like object) but with duplicate elements/rows removed.)
clinical_subset <- unique(clinical_subset)
#now the number of rows is 489

#creating a variable with all of the races
labs_race <- clinical_subset$race
#there are 5 races (american Indian/Alaska native,Asian,black/African american, white or not reported)
plot(pc$x[,1:2], col = Cols(labs_race), pch =19, xlab ="PC1",ylab="PC2")
plot(pc$x[,3:4], col = Cols(labs_race), pch =19, xlab ="PC3",ylab="PC4")
plot(pc$x[,5:6], col = Cols(labs_race), pch =19, xlab ="PC5",ylab="PC6")
plot(pc$x[,7:8], col = Cols(labs_race), pch =19, xlab ="PC5",ylab="PC6")
plot(pc$x[,9:10], col = Cols(labs_race), pch =19, xlab ="PC5",ylab="PC6")
plot(pc$x[,11:12], col = Cols(labs_race), pch =19, xlab ="PC5",ylab="PC6")
plot(pc$x[,13:14], col = Cols(labs_race), pch =19, xlab ="PC5",ylab="PC6")

#how many of each race do we have
summary(as.factor(labs_race))
#trying to match color to race
#white = purple, american Indian/Alaska native = red, black/African american = green, have 12 asian and not reported so I am not sure what colors they are associated with



########conclusion: there's not much clustering so we can't conclude much about the influence of genotype and race
####there's not much influence within these variables


