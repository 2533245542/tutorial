# This is the tests for the filter module of mistudio.
# It consists of 4 tests. The first three ones are individual tests, and the last one is a combined tests.
# Unlike an universial test, this test file is not automated. It relies on manually checking if the plot 
# these tests generate are the same as the ones generated in mistudio.

# Each test is self-contained and can be run individually.

# Each test loads a phyobj called phyloseq object called barb. The otu table in this barb has sample names as row names

# The first one test the aggregation function. The variable in this test is "ta5". It should generate two plots and user should compared with the one generated in mistudio
# also, in mistudio, the user should set other variables inrellevent. That means min count should be 0, and zero fraction should be 1

# The second test has the variable 0.9, the rest guide lines follow the first test

# the third test has the variable "minsum = 60000", the rest guide lines follow the first test

# The forth one is combined. It is important to see that it is sequential and the order matters. 
# This test is designed in this order because at the moment this file is created, the order of mistudio is aggregate->min count->zero fraction. Changing the order will changed the result.

########################################end of usage########################################

make_plot = function(otu){
  qplot(rowSums(otu), geom = "histogram", 
        xlab = "Adding each sample's otu count together", 
        ylab = "The number of sample with that amount of count"
  ) %>% print()
  otu_zerofraction = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x))))       
  
  qplot(otu_zerofraction, geom = "histogram", 
        xlab = "Proportion of zero in a otu", 
        ylab = "Number of sample with otu with that proportion"
  ) %>% print
  return(0)
}

library(phyloseq)
library(miLineage)
load("newbarb.RData")
phyobj = barb
otu = t(otu_table(phyobj)@.Data) 
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
	col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

uniq_ta5 = unique(tax[,"ta5"]) 
tax_filter = matrix(character(), nrow = length(uniq_ta5), ncol = 5)
i = 1
for(x in uniq_ta5){
  tax_filter[i,] = tax[tax[,"ta5"] == uniq_ta5[i], , drop = FALSE][1, 1:5, drop = FALSE]
  i = i + 1
}

otu_filter = t(aggregate(t(otu), list(tax[,"ta5"]), FUN = sum)[-1])
colnames(otu_filter) = paste0("ta5", uniq_ta5)
rownames(tax_filter) = paste0("ta5", uniq_ta5)
phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(tax_filter),
  sample_data(phyobj)
)

make_plot(otu_filter)
#######################################end of mistudio_filter_otu_aggregate#######################################

make_plot = function(otu){
  qplot(rowSums(otu), geom = "histogram", 
        xlab = "Adding each sample's otu count together", 
        ylab = "The number of sample with that amount of count"
  ) %>% print()
  otu_zerofraction = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x))))       
  
  qplot(otu_zerofraction, geom = "histogram", 
        xlab = "Proportion of zero in a otu", 
        ylab = "Number of sample with otu with that proportion"
  ) %>% print
  return(0)
}

library(phyloseq)
library(miLineage)
load("newbarb.RData")
phyobj = barb
otu = t(otu_table(phyobj)@.Data) 
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
	col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

criteria_remove = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x)) > 0.9)) # filter otu by 0.9 of zero 
otu_filter = otu[, !criteria_remove, drop = FALSE]
tax_filter = tax[!criteria_remove, , drop = FALSE]

phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(tax_filter),
  sample_data(phyobj)
)
make_plot(otu_filter)

#######################################end of mistudio_filter_otu_zerofraction#######################################
make_plot = function(otu){
  qplot(rowSums(otu), geom = "histogram", 
        xlab = "Adding each sample's otu count together", 
        ylab = "The number of sample with that amount of count"
  ) %>% print()
  otu_zerofraction = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x))))       
  
  qplot(otu_zerofraction, geom = "histogram", 
        xlab = "Proportion of zero in a otu", 
        ylab = "Number of sample with otu with that proportion"
  ) %>% print
  return(0)
}

library(phyloseq)
library(miLineage)
load("newbarb.RData")
phyobj = barb
minsum = 60000
otu = t(otu_table(phyobj)@.Data) 
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
	col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

criteria = apply(otu, MARGIN = 1, function(x) (sum(x) > minsum)) # find all rows that have otu sum > minsum		
otu_filter = otu[criteria, , drop = FALSE]
sample_filter = sample[criteria, , drop = FALSE]

phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(phyobj),
  sample_data(sample_data)
)
make_plot(otu_filter)
#######################################end of mistudio_filter_otu_minsum#######################################
library(phyloseq)
library(miLineage)
load("newbarb.RData")
phyobj = barb
otu = t(otu_table(phyobj)@.Data) 
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
                       col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

uniq_ta5 = unique(tax[,"ta5"]) 
tax_filter = matrix(character(), nrow = length(uniq_ta5), ncol = 5)
i = 1
for(x in uniq_ta5){
  tax_filter[i,] = tax[tax[,"ta5"] == uniq_ta5[i], , drop = FALSE][1, 1:5, drop = FALSE]
  i = i + 1
}

otu_filter = t(aggregate(t(otu), list(tax[,"ta5"]), FUN = sum)[-1])
colnames(otu_filter) = paste0("ta5", uniq_ta5)
rownames(tax_filter) = paste0("ta5", uniq_ta5)
phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(tax_filter),
  sample_data(phyobj)
)

make_plot(otu_filter)
#######################################end of mistudio_filter_otu_aggregate#######################################
make_plot = function(otu){
  qplot(rowSums(otu), geom = "histogram", 
        xlab = "Adding each sample's otu count together", 
        ylab = "The number of sample with that amount of count"
  ) %>% print()
  otu_zerofraction = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x))))       
  
  qplot(otu_zerofraction, geom = "histogram", 
        xlab = "Proportion of zero in a otu", 
        ylab = "Number of sample with otu with that proportion"
  ) %>% print
  return(0)
}

library(phyloseq)
library(miLineage)
load("newbarb.RData")
phyobj = barb


###
otu = t(otu_table(phyobj, taxa_are_rows = TRUE)@.Data)
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
                       col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

uniq_ta5 = unique(tax[,"ta5"]) 
tax_filter = matrix(character(), nrow = length(uniq_ta5), ncol = 5)
i = 1
for(x in uniq_ta5){
  tax_filter[i,] = tax[tax[,"ta5"] == uniq_ta5[i], , drop = FALSE][1, 1:5, drop = FALSE]
  i = i + 1
}

otu_filter = t(aggregate(t(otu), list(tax[,"ta5"]), FUN = sum)[-1])
colnames(otu_filter) = paste0("ta5", uniq_ta5)
rownames(tax_filter) = paste0("ta5", uniq_ta5)
phyobj = phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(tax_filter),
  sample_data(phyobj)
)



#####
minsum = 60000
otu = otu_table(phyobj)@.Data
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
                       col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

criteria = apply(otu, MARGIN = 1, function(x) (sum(x) > minsum)) # find all rows that have otu sum > minsum		
otu_filter = otu[criteria, , drop = FALSE]
sample_filter = sample[criteria, , drop = FALSE]

phyobj = phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(phyobj),
  sample_data(sample)
)


####
otu = otu_table(phyobj)@.Data
tax = tax_table(phyobj)@.Data 		
sample = as.data.frame(sample_data(phyobj)@.Data, stringsAsFactors = FALSE,  # sample_data gives list of data lists with no list names
                       col.names = sample_data(phyobj)@names, row.names = sample_data(phyobj)@row.names) # manually reassign column and row names		

criteria_remove = apply(otu, MARGIN = 2, function(x) return((sum(x == 0) / length(x)) > 0.9)) # filter otu by 0.9 of zero 
otu_filter = otu[, !criteria_remove, drop = FALSE]
tax_filter = tax[!criteria_remove, , drop = FALSE]

phyobj = phyloseq(
  otu_table(otu_filter, taxa_are_rows = FALSE),
  tax_table(tax_filter),
  sample_data(phyobj)
)
make_plot(otu_filter)

