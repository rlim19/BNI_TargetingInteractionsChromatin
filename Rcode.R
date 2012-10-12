#################################################################################
# Reproduce BNI of Targeting Interactions in Chromatin                          #
# Paper :Bayesian Network Analysis of Targeting Interaction in Chromatin (2010) #
# Author:Bas van Steensel, Ulrich Braunschweig, Guillaume J.Filion,             #
# Menzies Chen, Joke G.van Bemmel, and Trey Ideker.                             #  
#################################################################################     

#################################################################################
# Author : Ricky Lim                                                            #
# Course : Interdisciplinary seminar Biological Networks at UU, 2012            # 
#################################################################################

#################################################################################

# library required
library(reshape) # for shaping the data structure
library(igraph)  # for the analysis of network topology
library(ggplot2) # for plotting
library("grid")  # for layout of the plots

# Get the edges (interactions) of chromatin using bootstraap-BNI algorithm
#################################################################################

# load the dataset (DAMid map of 43 chromatin proteins) 
# this dataset was provided from suppl.figure 2 in the paper
data <- read.csv("~/Dropbox/Projects/BNI_chromatin/RawData/data.csv")
n.data <- data[,-1] # data without id label
# check the data 
head(n.data[1:5,1:5])
# each col corresponds to each chromatin component
# each row corresponds to the locus (the gene location in the genome)

# binarize the dataset (based on 95% percentile threshold)
# < 95% is assigned 0 (nontarget) otherwise 1(target)
discrete.data <- apply(n.data, 2, function(x) 
	                     ifelse(x < quantile(x, prob=95/100),0,1))
              # quantile for each chromatin (in the col of the matrix)

# check the binarized data
head(discrete.data[1:5,1:5])

# create the input data for the BNI algorithm and store into the input folder 
# the input data is created by 
# resampling 1000x, 4380 rows (genes) with replacement for each col(chromatin)
# ! note: the input data that were generated stored in input folder !
for (i in 1:1000){
	resampled.data <- discrete.data[sample(4380, replace=T), ]
	write.table(resampled.data, file=paste("./input/sample",i,".txt", sep=""),
							row.names=F, col.names=F)	
}

# running banjo in R using the input data (resampling data)
# for each input data, the best network is formed from the Banjo-BNI algorithm
# banjo is run for 1000 times using 1000 input files stored in the input directory
# ! note: the output data that were generated stored in output folder !
# ! in this output folder the files were already parsed !
for (i in 1:1000){
	# get the sample name of the input data 
	sample.name <- paste("sample",i, ".txt/",sep="")
	# update the observation sample in line 33, each iteration with the new sample
	# the new samples are ordered by numbers 1,2,3.....1000
	obs.name <- paste("33s/sample.txt/", sample.name, sep="")
	# update the settings file iteratively and store as settings.txt
	# note this file is being updated in each iteration
	settingcommand <- paste("sed",shQuote(obs.name), "RL.settings.txt > settings.txt")
	settingfile <- system(settingcommand)
	# note in the settings file, 
	# each output from banjo'd have different name based on the timestamp 
	system("java -jar banjo.jar settingsFile=settings.txt")
}
# in the output file, 43 chromatin components are described with numbers from 0-42

# parsing the output file
# combine the files as a list
files.output <- as.character(list.files(path="output"))
length(files.output)

# parse the output collecting only the edges from the influence scores matrix
for (i in seq_along(files.output)){
	file_name <- paste("output/", files.output[[i]], sep="")
	# get the lines of the infl_score matrix
	infl_score <- paste("gsed --i -n", shQuote('/->/p'), file_name)
	infl_score_file <- system(infl_score)	
	# get only the nodes within the parentheses
	node_regex<-shQuote("s/^.*(\\([0-9]\\+\\),\\([0-9]\\+\\))[ tab]*->[ tab]*(\\([0-9]\\+\\),\\([0-9]\\+\\)).*/\\1 \\3/")
	nodes <- paste ("gsed --i -e", node_regex, file_name)
	nodes_file <- system(nodes)
}

# combine all the parsed output in the output folder
all.the.files <- list.files("output",full=TRUE)
all.the.data <- lapply( all.the.files,  read.csv, header=F)
DATA <- do.call("rbind", all.the.data)
countTable<- as.data.frame(table(DATA))										 
colnames(countTable) <- c("parent.child", "bootstrap.freq")


Bootstrap.table.rl <- data.frame(countTable, colsplit(countTable$parent.child, 
																							 split = " ", 
																							 names =c("parent", "child")))
head(Bootstrap.table.rl)

													 
# get the chromatin names and store these as var.names' object
var.names<- colnames(n.data)
length(var.names) # check if there are 43 chromatin proteins

## make the names of the var.names similar to 
## the names used in the pubmed search term (paper)
# change Su.var3.7 -> Su(var)3-7
var.names[which(var.names=="Su.var.3.7")] <- "Su(var)3-7"
# change Su.var3.9 -> Su(var)3-9
var.names[which(var.names=="Su.var.3.9")] <- "Su(var)3-9"
# change MBD.like -> MBD-like
var.names[which(var.names=="MBD.like")] <- "MBD-like"
# change Su.UR. -> Su(UR)
var.names[which(var.names=="Su.UR.")] <- "Su(UR)"
# change Su.Hw. -> Su(Hw)
var.names[which(var.names=="Su.Hw.")] <- "Su(Hw)"
# change BEAF.32 -> BEAF-32
var.names[which(var.names=="BEAF.32")] <- "BEAF-32"
# check the 43 chromatin names
var.names

## change the chromatin component's number with its name
## 1. adjust the numbering of the chromatin component
## note: banjo output starts with 0 while in R numbering starts with 1
Bootstrap.table.rl$parent <- (Bootstrap.table.rl$parent + 1)
Bootstrap.table.rl$child <- (Bootstrap.table.rl$child + 1)
head(Bootstrap.table.rl)													 
nrow(Bootstrap.table.rl)													 

# function to get the chromatin's names from chromatin's number
map_chrom.comp <- function(x, chrom.names=var.names){
	comp.no <- x[1]
	names <- chrom.names[as.numeric(comp.no)]
	return (names)
}

# map chromatin component (Parent)
Bootstrap.table.rl$parent <- sapply(Bootstrap.table.rl$parent, map_chrom.comp)
# map chromatin component (Child)
Bootstrap.table.rl$child <- sapply(Bootstrap.table.rl$child, map_chrom.comp)

# get the freq.parent2child
Bootstrap.table.rl$bootstrap.freq.parent2child <- Bootstrap.table.rl$bootstrap.freq 

# get the freq. child2parent
Bootstrap.table.rl$bootstrap.freq.child2parent <- 0
for (i in seq_along(Bootstrap.table.rl$bootstrap.freq.child2parent)){
# to select the array index in which the corresponding component (parent) acts as child
	sub.array <- subset(Bootstrap.table.rl, parent == Bootstrap.table.rl$child[i] & 
		                child == Bootstrap.table.rl$parent[i])
	freq <- sub.array$bootstrap.freq
	ifelse (freq== 0, 
					Bootstrap.table.rl$bootstrap.freq.child2parent[i] <- 0, 
					Bootstrap.table.rl$bootstrap.freq.child2parent[i] <- freq)
}

# get the freq.combined
Bootstrap.table.rl$bootstrap.freq.combined<- apply(Bootstrap.table.rl[,c(5,6)], 1, sum)

# extract only cols with labels 
BootstraapScore.rl <- Bootstrap.table.rl[,3:7]
head(BootstraapScore.rl)

# order the table based on the alphabet of the parent's name
BootstraapScore.rl<-BootstraapScore.rl[order(BootstraapScore.rl$parent),]
rownames(BootstraapScore.rl) <- NULL
head(BootstraapScore.rl)
nrow(BootstraapScore.rl) # there are 1668 edges found (from the BNI)

# write the bootstraap result 													 
write.table(BootstraapScore.rl,quote=F, 
						file = "BootstraapScoreRL.txt", row.names=F)

# select pairs with Bootstraap Confidence Score(BCS) > 80% for plotting
#################################################################################

# select pairs (re-analysis) with bootstraap score > 800 out 1000 (80%) 
Selected.pairs.rl<-subset(BootstraapScore.rl, 
													bootstrap.freq.parent2child > 800 | 
													bootstrap.freq.child2parent> 800 | 
													bootstrap.freq.combined >800)
head(Selected.pairs.rl)
nrow(Selected.pairs.rl) 

# select pairs (paper) with bootstraap score > 800 out 1000 (80%)
# load the data from the paper containing pairs derived from the bootstrapping-BNI 
BootstraapScore.paper<- read.delim("~/Dropbox/Projects/BNI_chromatin/RawData/BootstraapScore.vanSteensel098822_SuppData2.txt")
Selected.pairs.paper<-subset(BootstraapScore.paper, 
													  bootstrap.freq.parent2child > 800 | 
														bootstrap.freq.child2parent> 800 | 
														bootstrap.freq.combined >800)
head(Selected.pairs.paper)
nrow(Selected.pairs.paper) # 52 pairs of targeting interactions

# collect edges list for plotting (reanalysis)
edges.list.rl <- Selected.pairs.rl[,1:3]
edges.list.rl[,3] <-  round(edges.list.rl[,3]/100) # scale down the bootstrap.freq 0-10
# exclude pairs with edges confidence score < 0
edges.list.rl <-subset(edges.list.rl, bootstrap.freq.parent2child>0) 
head(edges.list.rl)
nrow(edges.list.rl)
# collect edges list for plotting (paper)
edges.list.paper <- Selected.pairs.paper[,1:3]
edges.list.paper[,3] <-  round(edges.list.paper[,3]/100) # scale down the bootstrap.freq 0-10
# exclude pairs with edges confidence score <0
edges.list.paper <-subset(edges.list.paper, bootstrap.freq.parent2child>0) 
head(edges.list.paper)
nrow(edges.list.paper)

# input file for cytoscape of the reinferred network (reanalysis)
write.table(edges.list.rl, file="input_graph/EdgesList_RL",row.names=F, col.names=F, quote=F, sep="\t")
# input file for cytoscape of the reinferred network (paper)
write.table(edges.list.paper, file="input_graph/EdgesList_paper",row.names=F, col.names=F, quote=F, sep="\t")


# Computational Validation 
#################################################################################

#######################################################
# Function to get the pairs of chromatin components   #
# with a certain Bootstrap Confidence Score (BCS)     #
# for the computational validation                    #
#######################################################
# input   : pairs obtained from the applied BNI algorithm
#  - optional: threshold by default is set to be 800 from 1000 (80%) according to the paper
#  - threshold: bootstrap confidence score for the interaction of pairs
# output : selected pairs according to the chosen bootstrap confidence score (BCS)

pairs.validation <- function (allpairs, threshold =800){
	# select pairs with a certain level of bootstrap score
	Selected.pairs<- subset(allpairs, 
													bootstrap.freq.parent2child > threshold | 
													bootstrap.freq.child2parent > threshold |  
													bootstrap.freq.combined >threshold)
	
	# remove all the duplicated pairs
	# since each component in all pairs in the dataset is present as parent and child
	# this leads to the duplication of pairs.
	for (i in 1:nrow(Selected.pairs)){
		duplo_idx <-which(arr.ind=T,
											Selected.pairs$parent==as.character(Selected.pairs$child[i])
											& Selected.pairs$child==as.character(Selected.pairs$parent[i]))
		Selected.pairs[duplo_idx,] <- NA # NAs are given to the found duplicated pairs
	}
	
	# remove the NAs (NAs in this case corresponds to the duplicated pairs)
	Selected.pairs<-na.omit(Selected.pairs)
	# get only the pairs' names in col 1 and 2
	Selected.pairs<- Selected.pairs[,1:2]
	return (Selected.pairs)
}

##########################################################
# Function to map the PubMed search term from the paper  #
# for the corresponding chromatin components             #
##########################################################
# input : - x, chromatin component
#         - SearchTerm file, default file is provided from the paper (Suppl Dataset 5)
# load the pubmed search term for the 43 chromatin components 
PubMedSearchTerm <- read.delim("~/Dropbox/Projects/BNI_chromatin/RawData/PubMedSearchTerm.txt")
head(PubMedSearchTerm)
Map_PubMedTerm <- function(x, SearchTerm=PubMedSearchTerm){	
	idx_term <- which(SearchTerm$ChromatinComponent== x[1])
	PubMedTerm <- as.character(SearchTerm$PubMedTerm[idx_term])
	return (PubMedTerm)
}

# get the pairs of chromatin from the reanalysis
pairs.validation.rl <- pairs.validation(BootstraapScore.rl)
head(pairs.validation.rl)
nrow(pairs.validation.rl) 
# there are 50 chromatin pairs with 80% BCS


# get the pairs of chromatin from the paper
pairs.validation.paper <- pairs.validation(BootstraapScore.paper)
head(pairs.validation.paper)
nrow(pairs.validation.paper) 
# there are 52 chromatin pairs with 80% BCS

# create the dataframe with the search terms for the components of the pairs(reanalysis)
pairs.term.validation.rl <- data.frame(parent.term=
	                          sapply(pairs.validation.rl$parent, Map_PubMedTerm), 
																					child.term=
												    sapply(pairs.validation.rl$child, Map_PubMedTerm))

head(pairs.term.validation.rl)
# create the dataframe with the search terms for the components of the pairs(paper)
pairs.term.validation.paper <- data.frame(parent.term=
	                            sapply(pairs.validation.paper$parent, Map_PubMedTerm), 
											  			      			child.term=
								  						sapply(pairs.validation.paper$child, Map_PubMedTerm))
head(pairs.term.validation.paper)

# write this dataframe as an input for the validation (reanalysis)
write.table(pairs.term.validation.rl, file = "input_validation/pairs.rl.txt", 
						row.names=F, col.names=F, quote=F, sep ="\t")

# write this dataframe as an input for the validation (paper)
write.table(pairs.term.validation.paper, file = "input_validation/pairs.paper.txt", 
						row.names=F, col.names=F, quote=F, sep ="\t")


## Get all possible combination of pairs of 43 chromatin components

## combinne all the pairs of these 43 chromatin components 
all43pairs <- combn(var.names,2)
all43pairs <- t(all43pairs)
all43pairs <- as.data.frame(all43pairs)
colnames(all43pairs) <- c("parent", "child")
head(all43pairs)
nrow(all43pairs) # there are 903 possible pairs of components

# create the dataframe with the search terms for all 43 pairs
ALLpairs.term.validation <- data.frame(parent.term=
	                          sapply(as.character(all43pairs$parent), Map_PubMedTerm), 
																			 child.term=
										   			sapply(as.character(all43pairs$child), Map_PubMedTerm))
head(ALLpairs.term.validation)
nrow(ALLpairs.term.validation)
# write this dataframe as an input for the validation 
write.table(ALLpairs.term.validation, file = "input_validation/ALLpairs.txt", 
						row.names=F, col.names=F, quote=F, sep ="\t")


# data frame of the output validation
pubmed_val <- data.frame(Prediction=c("BN80_08", "AllComb_08", "BN80_12","BN80_12*", "AllComb_12"),
												 PubMedYear = as.factor(c("2008", "2008", "2012","2012","2012")),
												 Recovery = c("48%", "13%","81%","84%","61%"))
pubmed_val

# plot the output

# reorder the prediction based on the PubMedYear
pubmed_val<- transform(pubmed_val, Prediction = reorder(Prediction, order(PubMedYear)))
ggplot(pubmed_val, aes(x = Prediction, y=Recovery, fill = PubMedYear, order=PubMedYear))+
	geom_bar()+
	opts(title="Retrieval of Targeting Interactions in PubMed Abstracts")



# Network Analysis
#################################################################################

# note that connection in this analysis refers to out-degree

# get the edgelist of BNI's paper
el.paper <- as.matrix(edges.list.paper[,1:2], colnames=F)
graph.el.paper <- graph.edgelist(el.paper) # make the graph object
# compute the degree-out of the graph (number of neighboring edges directed out from a node)
outdegree.node.paper<-data.frame(Chromatin=get.vertex.attribute(graph.el.paper, 'name'),
					                       Connection = degree(graph.el.paper, mode="out"))
head(outdegree.node.paper)
mean(outdegree.node.paper$Connection) # mean of degree out is 2.5

# draw the distribution of out-degrees of previous BN80
degree.plot.paper<- ggplot(outdegree.node.paper, aes(x=Connection))+
	                  geom_histogram(binwidth = 1,colour="black",
	                  							 aes(fill = ..count..), origin = -0.5)+
	                  geom_vline(aes(xintercept=mean(Connection)),
	                  					 color="blue", linetype="dashed", size=1)+
						 	      scale_fill_gradient("Count", low = "green", high = "red")+
						      	xlab("out-degree")+
						 	      opts(title="Out-Degree Distribution of the previous BN80")

# write the output file for the previous BN80
write.table(outdegree.node.paper, file="input_graph/OutDegree_paper",
						row.names=F, col.names=F, quote=F, sep="\t")


# repeat the above analysis for the reconstructed network
el.rl <- as.matrix(edges.list.rl[,1:2], colnames=F)
graph.el.rl <- graph.edgelist(el.rl)
outdegree.node.rl<-data.frame(Chromatin=get.vertex.attribute(graph.el.rl, 'name'),
					                  Connection = degree(graph.el.rl, mode="out"))

# draw the distribution of out-degrees of reconstructed BN80
degree.plot.rl<-ggplot(outdegree.node.rl, aes(x=Connection))+
	              geom_histogram(binwidth = 1,colour="black",
	              							 aes(fill = ..count..), origin = -0.5)+
                geom_vline(aes(xintercept=mean(Connection)),
                					 color="blue", linetype="dashed", size=1)+
						    scale_fill_gradient("Count", low = "green", high = "red")+
		  	        xlab("out-degree")+
				        opts(title="Out-Degree Distribution of the reconstructed BN80")

# write the output of outdegree dist for reconstructed BN80
write.table(outdegree.node.rl, file="input_graph/OutDegree_rl",
						row.names=F, col.names=F, quote=F, sep="\t")

# plot the outdegree dist of previous and reconstructed BN80
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(degree.plot.rl, vp = vplayout(1, 1))
print(degree.plot.paper, vp = vplayout(1, 2))
