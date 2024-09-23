library(RSiena)
library(sna)
setwd(" ")

### Some fucntions
  #This adds a function to model geodesic distances, not included by default (as opposed to the IndegreeDistribution and OutdegreeDistribution functions).
  # simply copy this function in its entirety)
GeodesicDistribution <- function (i, data, sims, period, groupName,
                                  varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative)
  { gdi <- sapply(levls, function(i){ sum(a<=i) }) }
  else
  { gdi <- sapply(levls, function(i){ sum(a==i) }) }
  names(gdi) <- as.character(levls)
  gdi
}

plot_GOFs <- function(Model.results, network_name) {
  gofi <- sienaGOF(Model.results, IndegreeDistribution, verbose=TRUE, join=TRUE, varName=network_name)
  gofo <- sienaGOF(Model.results, OutdegreeDistribution, verbose=TRUE, join=TRUE, varName=network_name)
  gofgeo <- sienaGOF(Model.results, GeodesicDistribution, levls=1:12, verbose=TRUE, join=TRUE, varName=network_name)
  list(plot(gofi), plot(gofo), plot(gofgeo))
}


#### STEP 1: data pre paration
#first tell the script the number of panels and the threshold to binarize networks
panel_n <- 3
binarize_threshold <- 1
panel_n_range <- seq(1, panel_n, by=1)

#load nodelist data, and prepare attribute data
nodelist_raw <- read.table('nodelist_eg.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
#recode some values in the nodelist
nodelist <- nodelist_raw

nodelist[nodelist == 'missing'] <- ""
nodelist[nodelist == 'high'] <- 3
nodelist[nodelist == 'medium'] <- 2
nodelist[nodelist == 'low'] <- 1
nodelist[nodelist == 'below'] <- 0
nodelist[nodelist == 'above'] <- 1


nodelist$CCengageHigh <- ifelse(nodelist$relevance == 3, 1, 0)
nodelist$CCengageMed <- ifelse(nodelist$relevance == 2, 1, 0)
nodelist$CCengageLow <- ifelse(nodelist$relevance == 1, 1, 0)

#nodelist$CCinterestHigh <- ifelse(nodelist$interest == 3, 1, 0)
#nodelist$CCinterestMed <- ifelse(nodelist$interest == 2, 1, 0)
#nodelist$CCinterestLow <- ifelse(nodelist$interest == 1, 1, 0)


nodelist$verification[nodelist$verification == 'no'] <- 0
nodelist$verification[nodelist$verification == 'org'] <- 1
nodelist$verification[nodelist$verification == 'indiv'] <- 1


nodelist$environment[nodelist$environment == 'no'] <- 0
nodelist$environment[nodelist$environment == 'yes'] <- 1
nodelist$environment[nodelist$environment == ''] <- 0

nodelist$science[nodelist$science == 'no'] <- 0
nodelist$science[nodelist$science == 'yes'] <- 1
nodelist$science[nodelist$science == ''] <- 0

nodelist$economy[nodelist$economy == 'no'] <- 0
nodelist$economy[nodelist$economy == 'yes'] <- 1
nodelist$economy[nodelist$economy == ''] <- 0

nodelist$energy[nodelist$energy == 'no'] <- 0
nodelist$energy[nodelist$energy == 'yes'] <- 1
nodelist$energy[nodelist$energy == ''] <- 0

nodelist$expertise <- rowSums(sapply(nodelist[, c('environment', 'energy', 'science')], function(x) as.numeric(as.character(x))))

nodelist$laypublic[nodelist$expertise >= 1] <- 0
nodelist$laypublic[nodelist$expertise == 0] <- 1

nodelist$opleader[nodelist$expertise >= 1] <- 1
nodelist$opleader[nodelist$expertise == 0] <- 0



#prepare variable data: developmentalism_panel0 environmentalism_panel0.....
T_variables <- c('policy')
all_T_v <- list()
for (k in seq(length(T_variables))) {
  i_T_v <- T_variables[k]
  v_cols_ls <- c()
  for (n in panel_n_range)
  {
    v_cols <- paste( 'P', n-1, '_', i_T_v, sep="")
    v_cols_ls <- c(v_cols_ls, v_cols)
  }
  all_T_v[[k]] <- v_cols_ls
}

Tpolicy <- as.matrix(nodelist[, all_T_v[[4]]])
dimnames(Tpolicy) <- NULL
Tpolicy <- matrix(as.numeric(Tpolicy), ncol = ncol(Tpolicy))


#variable: relevance, gov, media ...
#CCengagement <- as.matrix(as.numeric(nodelist[, 'relevance']))
#dimnames(CCengagement) <- NULL

CCinterest <- as.matrix(as.numeric(nodelist[, 'interest']))
dimnames(CCinterest) <- NULL

#CCengagement_Median <- as.matrix(as.numeric(nodelist[, 'CCengagementMedian']))
#dimnames(CCengagement_Median) <- NULL

CCengagement_High <- as.matrix(as.numeric(nodelist[, 'CCengageHigh']))
dimnames(CCengagement_High) <- NULL
CCengagement_Med <- as.matrix(as.numeric(nodelist[, 'CCengageMed']))
dimnames(CCengagement_Med) <- NULL
CCengagement_Low <- as.matrix(as.numeric(nodelist[, 'CCengageLow']))
dimnames(CCengagement_Low) <- NULL

#CCinterest_High <- as.matrix(as.numeric(nodelist[, 'CCinterestHigh']))
#dimnames(CCinterest_High) <- NULL
#CCinterest_Med <- as.matrix(as.numeric(nodelist[, 'CCinterestMed']))
#dimnames(CCinterest_Med) <- NULL
#CCinterest_Low <- as.matrix(as.numeric(nodelist[, 'CCinterestLow']))
#dimnames(CCinterest_Low) <- NULL


verify <- as.matrix(as.numeric(nodelist[, 'verification']))
dimnames(verify) <- NULL
fans <- as.matrix(as.numeric(nodelist[, 'fans']))
dimnames(fans) <- NULL

gov <- as.matrix(as.numeric(nodelist[, 'gov']))
dimnames(gov) <- NULL
media <- as.matrix(as.numeric(nodelist[, 'media']))
dimnames(media) <- NULL
civil <- as.matrix(as.numeric(nodelist[, 'civil']))
dimnames(civil) <- NULL
business <- as.matrix(as.numeric(nodelist[, 'business']))
dimnames(business) <- NULL
#edu <- as.matrix(as.numeric(nodelist[, 'edu']))
#dimnames(edu) <- NULL

enviro <- as.matrix(as.numeric(nodelist[, 'environment']))
dimnames(enviro) <- NULL
sci <- as.matrix(as.numeric(nodelist[, 'science']))
dimnames(sci) <- NULL
economy <- as.matrix(as.numeric(nodelist[, 'economy']))
dimnames(economy) <- NULL
energy <- as.matrix(as.numeric(nodelist[, 'energy']))
dimnames(energy) <- NULL
laypublic <- as.matrix(as.numeric(nodelist[, 'laypublic']))
dimnames(laypublic) <- NULL
opleader <- as.matrix(as.numeric(nodelist[, 'opleader']))
dimnames(laypublic) <- NULL
expertise <- as.matrix(as.numeric(nodelist[, 'expertise']))
dimnames(expertise) <- NULL


#load in network matrices data and binarize relations
panel_array <- c()
panel_matrix <-  c()
for (n in panel_n_range)
{
  filename <- paste("mx_panel_", n-1, ".csv", sep = "")
  pn <- as.matrix(read.table(filename, row.names = 1, header = TRUE, sep =',', stringsAsFactors = TRUE))
  dimnames(pn) <- NULL
  pn[pn < binarize_threshold] <- 0
  pn[pn >= binarize_threshold] <- 1 #binarize network
  panel_array <- c(panel_array, pn)
  panel_matrix[[n]] <- pn
}



#plot networks, #color by Tparticipation #coordinate by all three nets
net1 <- network(panel_matrix[[1]], matrix.type="adjacency", directed=TRUE, multiple=FALSE)
net1 %v% 'Tpolicy'<- Tpolicy[, 1]
net2 <- network(panel_matrix[[2]], matrix.type="adjacency", directed=TRUE, multiple=FALSE)
net2 %v% 'Tpolicy'<- Tpolicy[, 2]
net3 <- network(panel_matrix[[3]], matrix.type="adjacency", directed=TRUE, multiple=FALSE)
net3 %v% 'Tpolicy'<- Tpolicy[, 3]
#net4 <- network(panel_matrix[[4]], matrix.type="adjacency", directed=TRUE, multiple=FALSE)
#net4 %v% 'Tpolicy'<- Tpolicy[, 4]


coordinates <-  plot(net1 + net2 + net3)
par( mfrow = c( 1, 3 ) )
par(mar = c(1, 0, 1, 0))
pal <- colorRampPalette(c("white", "black"))
color <- pal(4) #number of attribute value levels
gplot(net1, displaylabels=F, label = nodelist$username, label.cex=.5, vertex.col=color[Tpolicy[, 1]], edge.col = 'grey',
         sub = paste('panel', 1, sep=""), coord = coordinates, arrowhead.cex = 0.8, label.pos=1)
gplot(net2, displaylabels=F, label = nodelist$username, label.cex=.5, vertex.col=color[Tpolicy[, 2]], edge.col = 'grey',
         sub = paste('panel', 2, sep=""), coord = coordinates, arrowhead.cex = 0.8, label.pos=1)
gplot(net3, displaylabels=F, label = nodelist$username, label.cex=.5, vertex.col=color[Tpolicy[, 3]], edge.col = 'grey',
         sub = paste('panel', 3, sep=""), coord = coordinates, arrowhead.cex = 0.8, label.pos=1)
#gplot(net4, displaylabels=TRUE, label = nodelist$username, label.cex=.5, vertex.col=color[Tpolicy[, 4]], edge.col = 'grey',
#         main = paste('panel', 4, sep=""), coord = coordinates, arrowhead.cex = 0.5, label.pos=1)
#legend("bottomright",fill=color,legend=paste("Tparticipation",0:4),cex=0.75)


# Moran's autocorrelation for outgoing and incoming ties:
nacf(net1, Tpolicy[, 1], type="moran", neighborhood.type='total')[2]
nacf(net2, Tpolicy[, 2], type="moran", neighborhood.type='total')[2]
nacf(net3, Tpolicy[, 3], type="moran", neighborhood.type='total')[2]


# is there sufficient change in alcohol use?
# examine correlations and discrete changes
cor(Tpolicy)
table(Tpolicy[, 1], Tpolicy[, 2], useNA='always')
table(Tpolicy[, 2], Tpolicy[, 3], useNA='always')
#table(Tpolicy[, 3], Tpolicy[, 4], useNA='always')
# & total:
( totalChange <- table(Tpolicy[, 1], Tpolicy[, 3], useNA='always') )
sum(diag(totalChange)) / sum(totalChange)

t.test(Tpolicy[,1], Tpolicy[,2])
t.test(Tpolicy[,2], Tpolicy[,3])

### STEP 2
# creates a dependent network object called "communication" which specifies the ordering of the networks
# and their dimensions (n x n people, n timepoints)
n_actor <- dim(nodelist)[1]
communication <- sienaDependent(array(panel_array, dim=c(n_actor, n_actor, panel_n)))


#create constant attributes
Gov <- coCovar(gov[,1])
Media <- coCovar(media[,1])
Civil <- coCovar(civil[,1])
Business <- coCovar(business[,1])
#Edu <- coCovar(edu[,1])

Verify <- coCovar(verify[,1])
Fans <- coCovar(fans[,1])
#CCEngagement <- coCovar(CCengagement[, 1])
CCInterest <- coCovar(CCinterest[, 1])

CCengagement_High <- coCovar(CCengagement_High[, 1])
CCengagement_Med <- coCovar(CCengagement_Med[, 1])
CCengagement_Low <- coCovar(CCengagement_Low[, 1])

#CCinterest_High <- coCovar(CCinterest_High[, 1])
#CCinterest_Med <- coCovar(CCinterest_Med[, 1])
#CCinterest_Low <- coCovar(CCinterest_Low[, 1])

#CCengagement_Median <- coCovar(CCengagement_Median[, 1])

Enviro <- coCovar(enviro[,1])
Sci <- coCovar(sci[,1])
Economy <- coCovar(economy[,1])
Energy <- coCovar(energy[,1])
Laypublic <- coCovar(laypublic[,1])
Opleader <- coCovar(opleader[,1])
Expertise <- coCovar(expertise[,1])


#create changing covariates which serves as both an independent variable and as a dependent variable
# the dual role of the vairable is specified by the sienaDependent, type="behavior" command,
# which also indicates to RSiena that the model is one for the co-evolution of networks and behavior
T_Policy <- sienaDependent(Tpolicy, type="behavior")
#T_Crisis <- sienaDependent(Tcrisis, type="behavior")
#T_Proximity <- sienaDependent(Tproximity, type="behavior")
#T_Environmentalism <- sienaDependent(Tenvironmentalism, type="behavior")
#T_Policy <- varCovar(Tpolicy)


### STEP 3: create MyData and initial MyEff
#creates the full RSiena object(MyData) by adding the networks and attribute objects


MyData <- sienaDataCreate(communication,
                          T_Policy,
                          CCInterest,
                          Fans, Verify,
                          Gov, Media, Civil,
                          Opleader)


#create an initial effects object which includes outdegree and reciprocity
# (as well as the rate parameters and shape parameters)
#additional effects will be added to this object later
MyEff <- getEffects(MyData)
#prints the initial report from the data object and given the effects object
print01Report(MyData, modelname = 'init')


#Model4: co-evolution
#MyEff <- getEffects(MyData)
#MyEff <- includeEffects(MyEff, gwespFF, inPop, isolateNet, inIsDegree, isolatePop, outAct, inAct)
#MyEff <- setEffect(MyEff, outTrunc, parameter=2)

MyEff <- getEffects(MyData)
MyEff <- includeEffects(MyEff, transTrip, inPopSqrt, outActSqrt, inActSqrt)


MyEff<- includeEffects(MyEff,altX, interaction1="Verify")
MyEff<- includeEffects(MyEff, altX, interaction1="Fans")

MyEff<- includeEffects(MyEff,altX, sameX, interaction1="Gov")
MyEff<- includeEffects(MyEff,altX, interaction1="Media")
MyEff<- includeEffects(MyEff,altX,  interaction1="Civil")
#MyEff<- includeEffects(MyEff,altX,  interaction1="Business")


#MyEff<- includeEffects(MyEff,egoX,  interaction1="Enviro")
#MyEff<- includeEffects(MyEff,egoX, interaction1="Sci")
#MyEff<- includeEffects(MyEff,egoX,  interaction1="Economy")

MyEff<- includeEffects(MyEff,altX,  interaction1="Opleader")

MyEff<- includeEffects(MyEff,egoX, altX, simX, interaction1="T_Policy") #higher interest in devlp, more incoming/outgoing ties


#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "Gov")
#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "Media")
#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCEngagement")

#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCengagement_High")
#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCengagement_Med")
#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCengagement_Low")

#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCinterest_High")
#MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCinterest_Med")

MyEff <- includeEffects(MyEff, name = "T_Policy", effFrom, interaction1 = "CCInterest")

MyEff <- includeEffects(MyEff, name = "T_Policy", avSim, interaction1 = "communication")

MyEff <- includeEffects(MyEff, name = "T_Policy", indeg, interaction1 = "communication")
MyEff <- includeEffects(MyEff, name = "T_Policy", outdeg, interaction1 = "communication")

MyEff <- includeTimeDummy(MyEff, name="T_Policy", linear, timeDummy="all")


Model4 <- sienaModelCreate(useStdInits = TRUE, projname = 'Model4-results')
Model4.results <- siena07(Model4, data=MyData, effects=MyEff, batch=FALSE, verbose=FALSE, returnDeps=TRUE, useCluster=TRUE)
Model4.results
siena.table(Model4.results, type="html", sig=TRUE)
plot_GOFs(Model4.results, "communication")


print(p)

## Conduct the score type test to assess whether heterogeneity is present.
tt <- sienaTimeTest(Model4.results)
plot(tt, effects=1:2)
summary(tt)


#if need prevAns
Model4 <- sienaModelCreate(useStdInits = TRUE, projname = 'Model4-results')
Model4.results <- siena07(Model4, data=MyData, effects=MyEff, batch=FALSE, verbose=FALSE,
                          returnDeps=TRUE, useCluster=TRUE, prevAns = Model4.results)
Model4.results
siena.table(Model4.results, type="html", sig=TRUE)
plot_GOFs(Model4.results, "communication")


