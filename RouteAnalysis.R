#### project created to analyse trajectories and shared days
#### among 9/11 hijackers
#### First load data
#### then create analysis tables
#### See https://github.com/rafaelprietocuriel/RoutesSeptember11 for a description
#### and to understand the source of files
####
#### Created by Rafael Prieto Curiel and Olivier Walther. 2022/08/24

require(igraph); require(viridis); require(RColorBrewer); require(fossil)
windowsFonts("Arial" = windowsFont("Arial")); windowsFonts("Helvetica" = windowsFont("Helvetica"))

#### load data
{
TR <- read.csv("Routes.csv")
AT <- read.csv("Attributes.csv")
Locations <- read.csv("Locations.csv")
Contacts <- read.csv("SocialContactMatrix.csv")
Contacts <- as.matrix(Contacts[, 2:20])
Contacts <- pmax(Contacts, t(Contacts))
Contacts[lower.tri(Contacts)] <- 0
}

#### analyse data, sim, index and traj and net
{
#### create data matrices
{
VecR <- function(AcDate, POI, startd =35750, endd = 37145){
    u <- order(AcDate, decreasing = F)
    AcDate <- AcDate[u]; POI <- POI[u]
    ToD <- seq(from = startd, to = endd)
    locs  <- rep(NA, endd - startd+1)
    for(k in 1:length(u)){
        f <- which(ToD == AcDate[k])[1]
        locs[f:length(locs)] <- POI[k]
    }
    return(locs)
}

jaccard <- function(user1, user2) {
    mean(user1 == user2, na.rm = T)
}

#### obtain patterns
{
CleanPOI <- function(POI){
    jump <- rep(NA, 50)
    a <- POI[1]
    l <- a
    for(k in 2:length(POI)){
        if(POI[k] != l){a <- c(a, POI[k]); l <- POI[k]}
    }
    jump[1:length(a)] <- a
    return(jump)
}
}

#### create JUMPS CLEAN
for(u in 1:19){
    filt <- TR$Actor == u
    AcPOI <- TR$POI[filt]
    if(u == 1){
        AcJUMP <- CleanPOI(AcPOI)}else{
            AcJUMP <- cbind(AcJUMP, CleanPOI(AcPOI))}
}

#### create JUMPS LAST YEAR ONLY
for(u in 1:19){
    filt <- TR$Actor == u
    tFilt <- TR$date > (37145-365)
    AcPOILastYear <- TR$POI[filt & tFilt]
    if(u == 1){
        AcJUMPLastYear <- CleanPOI(AcPOILastYear)}else{
            AcJUMPLastYear <- cbind(AcJUMPLastYear, CleanPOI(AcPOILastYear))}
}

#### jaccard trajectories
jaccardTraj <- function(user1, user2){
    M1 <- c(); M2 <- c()
    user1 <- user1[!is.na(user1)]; user2 <- user2[!is.na(user2)]
    for(k in 1:(length(user1)-1)){M1 <- c(M1, paste(user1[k+1]," - ",user1[k+1]))}
    for(k in 1:(length(user2)-1)){M2 <- c(M2, paste(user2[k+1]," - ",user2[k+1]))}
    return( (sum(M1 %in% M2)/length(M1) + sum(M2 %in% M1)/length(M2))/2)
}

#### create matrix with POI by day
for(u in 1:19){
    filt <- TR$Actor == u
    AcDate <- TR$date[filt]
    AcPOI <- TR$POI[filt]
    if(u == 1){
        AcSEQ <- VecR(AcDate, AcPOI)}else{
        AcSEQ <- cbind(AcSEQ, VecR(AcDate, AcPOI))}
}

#### create similarity by days, all days
{
simM <- matrix(rep(0, 19*19), ncol = 19)
for(k in 1:19){for(j in 1:19){
    simM[k, j] <- jaccard(AcSEQ[, k], AcSEQ[, j])
}}
}

#### last year only
{
AcSEQLastY <- AcSEQ[(1396-365):1396,]
simMLastY <- matrix(rep(0, 19*19), ncol = 19)
for(k in 1:19){for(j in 1:19){
    simMLastY [k, j] <- jaccard(AcSEQLastY[, k], AcSEQLastY[, j])
}}
}

#### create trajectories all days
{
simTraj <- matrix(rep(0, 19*19), ncol = 19)
for(k in 1:19){for(j in 1:19){
    simTraj[k, j] <- jaccardTraj(AcJUMP[, k], AcJUMP[, j])
}}
}
    
#### create trajectories last year
{
simTrajLastY <- matrix(rep(0, 19*19), ncol = 19)
for(k in 1:19){for(j in 1:19){
    simTrajLastY[k, j] <- jaccardTraj(AcJUMPLastYear[, k], AcJUMPLastYear[, j])
}}
}

#### intercells sim All days and simulate 10K times
{
diag(simM) <- NA
InterGAll <- c()
for(k in 1:4){
    filt <- which(AT$cell == k)
    tsiM <- simM[filt, ]
    tsiM <- tsiM[, filt]
    InterGAll <- c(InterGAll, mean(tsiM, na.rm = T))
}

SimGAll <- c()
for(k in 1:10000){
    q <- sample.int(19, size = 5)   
    tsiM <- simM[q, ]
    tsiM <- tsiM[, q]
    SimGAll <- c(SimGAll, mean(tsiM, na.rm = T))
}
}

#### InterGroup sim Last year and simulate 10K times
{
    diag(simMLastY) <- NA
    InterGLastY <- c()
    for(k in 1:4){
        filt <- which(AT$cell == k)
        tsiM <- simMLastY[filt, ]
        tsiM <- tsiM[, filt]
        InterGLastY <- c(InterGLastY, mean(tsiM, na.rm = T))
    }
    
    SimGLastYear <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simMLastY[q, ]
        tsiM <- tsiM[, q]
        SimGLastYear <- c(SimGLastYear, mean(tsiM, na.rm = T))
    }
}

#### InterGroup TRAJ All days
{
    diag(simTraj) <- NA
    InterTrajAll <- c()
    for(k in 1:4){
        filt <- which(AT$cell == k)
        tsiM <- simTraj[filt, ]
        tsiM <- tsiM[, filt]
        InterTrajAll <- c(InterTrajAll, mean(tsiM, na.rm = T))
    }

    SimTrajAll <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simTraj[q, ]
        tsiM <- tsiM[, q]
        SimTrajAll <- c(SimTrajAll, mean(tsiM, na.rm = T))
    }
}

#### InterGroup TRAJ Last Year
{
    diag(simTrajLastY) <- NA
    InterTrajLastYear <- c()
    for(k in 1:4){
        filt <- which(AT$cell == k)
        tsiM <- simTrajLastY[filt, ]
        tsiM <- tsiM[, filt]
        InterTrajLastYear <- c(InterTrajLastYear, mean(tsiM, na.rm = T))
    }
    
    SimTrajLastYear <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simTraj[q, ]
        tsiM <- tsiM[, q]
        SimTrajLastYear <- c(SimTrajLastYear, mean(tsiM, na.rm = T))
    }
}

#### pilot vs muscle sim All days and simulate 10K times
{
    jobs <- c("pilot", "muscle")
    diag(simM) <- NA
    InterGAllPM <- c()
    for(k in 1:2){
        filt <- which(AT$job == jobs[k])
        tsiM <- simM[filt, ]
        tsiM <- tsiM[, filt]
        InterGAllPM <- c(InterGAllPM, mean(tsiM, na.rm = T))
    }
    
    SimGAllPM <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simM[q, ]
        tsiM <- tsiM[, q]
        SimGAllPM <- c(SimGAllPM, mean(tsiM, na.rm = T))
    }
}

#### pilot vs muscle sim Last year and simulate 10K times
{
    diag(simMLastY) <- NA
    InterGLastYPM <- c()
    for(k in 1:2){
        filt <- which(AT$job == jobs[k])
        tsiM <- simMLastY[filt, ]
        tsiM <- tsiM[, filt]
        InterGLastYPM <- c(InterGLastYPM, mean(tsiM, na.rm = T))
    }
    
    SimGLastYearPM <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simMLastY[q, ]
        tsiM <- tsiM[, q]
        SimGLastYearPM <- c(SimGLastYearPM, mean(tsiM, na.rm = T))
    }
}

#### pilot vs muscle TRAJ All days
{
    diag(simTraj) <- NA
    InterTrajAllPM <- c()
    for(k in 1:2){
        filt <- which(AT$job == jobs[k])
        tsiM <- simTraj[filt, ]
        tsiM <- tsiM[, filt]
        InterTrajAllPM <- c(InterTrajAllPM, mean(tsiM, na.rm = T))
    }
    
    SimTrajAllPM <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simTraj[q, ]
        tsiM <- tsiM[, q]
        SimTrajAllPM <- c(SimTrajAllPM, mean(tsiM, na.rm = T))
    }
}

#### pilot vs muscle TRAJ Last Year
{
    diag(simTrajLastY) <- NA
    InterTrajLastYearPM <- c()
    for(k in 1:2){
        filt <- which(AT$job == jobs[k])
        tsiM <- simTrajLastY[filt, ]
        tsiM <- tsiM[, filt]
        InterTrajLastYearPM <- c(InterTrajLastYearPM, mean(tsiM, na.rm = T))
    }
    
    SimTrajLastYearPM <- c()
    for(k in 1:10000){
        q <- sample.int(19, size = 5)   
        tsiM <- simTraj[q, ]
        tsiM <- tsiM[, q]
        SimTrajLastYearPM <- c(SimTrajLastYearPM, mean(tsiM, na.rm = T))
    }
}
}

#### igraph network
{
J <- data.frame(from = TR$POI[1],
                to = TR$POI[2],
                Actor = TR$Actor[2],
                Cell = TR$Cell[2])
for(k in 3:dim(TR)[1]){
    if(TR$Actor[k] == TR$Actor[k-1] & TR$POI[k] != TR$POI[k-1]){
        J <- rbind(J, data.frame(from = TR$POI[k-1],
                        to = TR$POI[k],
                        Actor = TR$Actor[k],
                        Cell = TR$Cell[k]))
    }
}

G <- graph_from_data_frame(J,
                           directed = T,
                           vertices = Locations)
G1 <- subgraph.edges(G, 
                     which(J$Cell == 1),
                     delete.vertices = F)
G2 <- subgraph.edges(G, 
                     which(J$Cell == 2),
                     delete.vertices = F)
G3 <- subgraph.edges(G, 
                     which(J$Cell == 3),
                     delete.vertices = F)
G4 <- subgraph.edges(G, 
                     which(J$Cell == 4),
                     delete.vertices = F)
GP <- subgraph.edges(G, 
                     which(J$Actor %in% which(AT$job == "pilot")),
                     delete.vertices = F)
GM <- subgraph.edges(G, 
                     which(J$Actor %in% which(AT$job == "muscle")),
                     delete.vertices = F)

E(G)$w <- count.multiple(G)
E(G1)$w <- count.multiple(G1)
E(G2)$w <- count.multiple(G2)
E(G3)$w <- count.multiple(G3)
E(G4)$w <- count.multiple(G4)
E(GP)$w <- count.multiple(GP)
E(GM)$w <- count.multiple(GM)

AllEd <- rbind(as_data_frame(G1, what = "edges"),
               as_data_frame(G2, what = "edges"),
               as_data_frame(G3, what = "edges"),
               as_data_frame(G4, what = "edges"))

ContactM <- graph_from_adjacency_matrix(Contacts)
for(k in 1:dim(AT)[2]){
    ContactM <- set_vertex_attr(ContactM,
                          name = names(AT)[k],
                          value = AT[,k])
}
ContactM <- simplify(ContactM)
#write.csv(AllEd, file = "WeightedEdges.csv", row.names = F)
}

#### relational graph
{
RG <- graph_from_adjacency_matrix(mode = "undirected",
                                  #simM,
                                  #simMLastY,
                                  #simTrajLastY,
                                  simTraj,
                                  weighted = T,
                                  diag = F)

for(k in 1:dim(AT)[2]){
RG <- set_vertex_attr(RG,
                       name = names(AT)[k],
                       value = AT[,k])
}
cl <- cluster_louvain(RG, weights = E(RG)$weight)
}

#### compare the index
{
rand.index(AT$cell, cl$membership)
adj.rand.index(AT$cell, cl$membership)
rand.index(AT$job == "pilot", cl$membership)
adj.rand.index(AT$job == "pilot", cl$membership)
}

#### simulate 10000 times
{
simRandInd <- c()
for(k in 1:10000){
x <- AT$cell[sample(1:19)] 
y <- AT$cell[sample(1:19)]  
simRandInd <- c(simRandInd, 
                adj.rand.index(x,y))
}
simRandInd <- simRandInd[!is.na(simRandInd)]
}

#### how similar are two social matrices
{
ToSimN<- simTraj
diag(ToSimN) <- 0
Res <- data.frame(tresh = c(), sim = c())
CM <- as.matrix(as_adj(ContactM, names = F, sparse = F))
CM <- CM + t(CM)
for(k in 1:1000){
    x <- runif(1) 
    tempM <- (ToSimN > x)*1
    TS <- sum(abs(CM - tempM))/(sum(tempM) + sum(CM))
    Res <- rbind(Res, data.frame(tresh = x, sim = TS))
}

vv = order(Res$tresh)
plot(Res$tresh[vv], 
     ylab = "Differences",
     xlab = "Treshold M",
     Res$sim[vv],
     lwd = 4,
     type = "l", col = "tomato")

u <- which.min(Res$sim)
Res[u,]
sum(abs(CM - (1-CM)))/(sum(CM) + sum(1-CM))
#### based on the whole trajectory
### with a threshold of 0.7556848 `3/4 we obtain optimum
#### based on the trajectory last year
### with a threshold of 0.7327808 we obtain optimum
#### so optimum threshold for both: 3/4
}
}
    
#### FIGURES
{
#### igraph configs for 48 locations
{
    u <- sample.int(48)
    L <- L[u, ]
}

##### 19 * 19 figures
{
    #### PDF plot known days index
    {
        pdf("SharedDaysByCellAllTime.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,5,10,14,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*simM[k, j])+1],
                   cex = 2*simM[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot last year index
    {
        pdf("SharedDaysByCellLastYear.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,5,10,14,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*simMLastY[k, j])+1],
                   cex = 2*simMLastY[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot known days trajectories
    {
        pdf("TrajByCellAllTime.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,5,10,14,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*simTraj[k, j])+1],
                   cex = 2*simTraj[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot last year  trajectories
    {
        pdf("TrajByCellLastYear.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,5,10,14,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*simTrajLastY[k, j])+1],
                   cex = 2*simTrajLastY[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot known days index PM
    {
        u <- which(AT$job == "pilot"); v <- which(AT$job == "muscle")
        PilotMuscle <- simM
        PilotMuscle <- PilotMuscle[c(u,v), ]
        PilotMuscle <- PilotMuscle[, c(u,v)]
        pdf("SharedDaysPilotMuscleAllTime.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,4,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot last year index PM
    {
        u <- which(AT$job == "pilot"); v <- which(AT$job == "muscle")
        PilotMuscle <- simMLastY
        PilotMuscle <- PilotMuscle[c(u,v), ]
        PilotMuscle <- PilotMuscle[, c(u,v)]
        pdf("SharedDaysPilotMuscleLastYear.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,4,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot known days trajectories PM
    {
        u <- which(AT$job == "pilot"); v <- which(AT$job == "muscle")
        PilotMuscle <- simTraj
        PilotMuscle <- PilotMuscle[c(u,v), ]
        PilotMuscle <- PilotMuscle[, c(u,v)]
        pdf("TrajPilotMuscleAllTime.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,4,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot last year  trajectories PM
    {
        u <- which(AT$job == "pilot"); v <- which(AT$job == "muscle")
        PilotMuscle <- simTrajLastY
        PilotMuscle <- PilotMuscle[c(u,v), ]
        PilotMuscle <- PilotMuscle[, c(u,v)]
        pdf("TrajPilotMuscleLastYear.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,4,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    
    #### PDF plot known days index brothers
    {
        u <- c(3,5,9,10,12,19, #brothers
               1,2,4,6,7,8,11, #notb
               13,14,15,16,17,18)
        PilotMuscle <- simM
        PilotMuscle <- PilotMuscle[u, ]
        PilotMuscle <- PilotMuscle[, u]
        pdf("SharedDaysBrothersAllTime.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,2, 4, 6,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot last year index brothers
    {
        u <- c(3,5,9,10,12,19, #brothers
               1,2,4,6,7,8,11, #notb
               13,14,15,16,17,18)
        PilotMuscle <- simM
        PilotMuscle <- PilotMuscle[u, ]
        PilotMuscle <- PilotMuscle[, u]
        pdf("SharedDaysBrothersLastYear.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,2, 4, 6,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot known days trajectories brothers
    {
        u <- c(3,5,9,10,12,19, #brothers
               1,2,4,6,7,8,11, #notb
               13,14,15,16,17,18)
        PilotMuscle <- simTraj
        PilotMuscle <- PilotMuscle[u, ]
        PilotMuscle <- PilotMuscle[, u]
        pdf("TrajBrothersAllTime.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,2, 4, 6,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }
    
    #### PDF plot last year  trajectories brothers
    {
        u <- c(3,5,9,10,12,19, #brothers
               1,2,4,6,7,8,11, #notb
               13,14,15,16,17,18)
        PilotMuscle <- simTrajLastY
        PilotMuscle <- PilotMuscle[u, ]
        PilotMuscle <- PilotMuscle[, u]
        pdf("TrajBrothersLastYear.pdf", width = 6, height = 6)
        par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
        plot(1, col = NA,
             xaxt = "n", yaxt = "n",
             asp = 1,
             xlim = c(-1,20), ylim = c(0,21))
        ed <- c(0,2, 4, 6,19)+0.5
        for(k in 1:length(ed)){
            points(c(0.5,19.5), 
                   20-c(ed[k], ed[k]),
                   type = "l")
            points(c(ed[k], ed[k]),
                   c(0.5,19.5),
                   type = "l")
        }
        for(k in 2:length(ed)){
            polygon(c(ed[k], ed[k], ed[k-1], ed[k-1]), 
                    20-c(ed[k], ed[k-1], ed[k-1], ed[k]),
                    col = "gray90")
        }
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(100))
        for(k in 1:19){for(j in 1:19){
            points(k, 20-j,
                   pch = 22,
                   bg = cols[floor(99*PilotMuscle[k, j])+1],
                   cex = 2*PilotMuscle[k, j]+.1)
        }}
        for(k in 1:19){
            points(k, 20-k,
                   pch = 21,
                   bg = "gray80",
                   cex = 2)
            text(0, 20-k, k, adj = 1)
            text(k, 20, k)
        }
        dev.off()
    }    
}

#### igraph plot all
{
    pdf(file = "NetW.pdf",
        height = 6, width = 6)
    par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
    #### PILOT
    {
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        text(0,1.1, "all journeys", cex = 1)
        plot(G, add = T, 
             vertex.frame.color="gray10",
             vertex.frame.width=.1,
             vertex.color = "skyblue",
             edge.color = gray.colors(5)[E(G)$w],
             layout= L,
             edge.width = E(G)$w,
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.4, 
             vertex.label.dist=0.,
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label = NA,
             vertex.size = 1.*sqrt(degree(G))+3)
        for(k in 1:48){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
    }
    dev.off()
}

#### igraph plot muscle pilots
{
    pdf(file = "PilotMuscleNetW.pdf",
        height = 6, width = 12)
    par(mar = c(0,0,0,0)+.1, mfrow = c(1,2))
    #### PILOT
    {
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        text(0,1.1, "pilots", cex = 1)
        plot(GP, add = T, 
             vertex.frame.color="gray10",
             vertex.frame.width=.1,
             vertex.color = "skyblue",
             edge.color = gray.colors(5)[E(GP)$w],
             layout= L,
             edge.width = E(GP)$w,
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.4, 
             vertex.label.dist=0.,
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label = NA,
             vertex.size = 1.*sqrt(degree(GP))+3)
        for(k in 1:20){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
    }
    
    #### MUSCLE
    {
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        text(0,1.1, "muscle", cex = 1)
        plot(GM, 
             vertex.size = 1.*sqrt(degree(GM))+3,
             edge.width = E(GM)$w,
             vertex.frame.color="gray10",
             vertex.frame.width=.1,
             vertex.color = "skyblue",
             edge.color = gray.colors(100),
             layout= L,
             add = T, 
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.4, 
             vertex.label.dist=0,
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label = NA)
        for(k in 1:20){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
        }
    dev.off()
}

#### igraph plot 4 cells
{
pdf(file = "FourCellsNetW.pdf",
    height = 12, width =12)
par(mar = c(0,0,0,0), mfrow = c(2,2))
#### FIRST CELL
{
plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
         asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    text(0,1.1, "cell 1", cex = 1)
    plot(G1, 
         vertex.size = 1.*sqrt(degree(G1))+3,
         edge.width = E(G1)$w,
         vertex.frame.color="gray10",
         vertex.frame.width=.1,
         vertex.color = "skyblue",
         edge.color = gray.colors(100),
         layout= L,
         add = T, 
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.label.font=2,
         vertex.label.cex=0.4, 
         vertex.label.dist=0,
         edge.arrow.size = 0.01,
         edge.curved = 0.2,
         vertex.label = NA)
    for(k in 1:20){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
}

#### 2 CELL
{
    plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
         asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    text(0,1.1, "cell 2", cex = 1)
    plot(G2, 
         vertex.size = 1.*sqrt(degree(G2))+3,
         edge.width = E(G2)$w,
         vertex.frame.color="gray10",
         vertex.frame.width=.1,
         vertex.color = "skyblue",
         edge.color = gray.colors(100),
         layout= L,
         add = T, 
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.label.font=2,
         vertex.label.cex=0.4, 
         vertex.label.dist=0,
         edge.arrow.size = 0.01,
         edge.curved = 0.2,
         vertex.label = NA)
    for(k in 1:20){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
}

#### 3 CELL
{
    plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
         asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    text(0,1.1, "cell 3", cex = 1)
    plot(G3, 
         vertex.size = 1.*sqrt(degree(G3))+3,
         edge.width = E(G3)$w,
         vertex.frame.color="gray10",
         vertex.frame.width=.1,
         vertex.color = "skyblue",
         edge.color = gray.colors(100),
         layout= L,
         add = T, 
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.label.font=2,
         vertex.label.cex=0.4, 
         vertex.label.dist=0,
         edge.arrow.size = 0.01,
         edge.curved = 0.2,
         vertex.label = NA)
    for(k in 1:20){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
}

#### 4 CELL
{
plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
         asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    text(0,1.1, "cell 4", cex = 1)
    plot(G4, 
         vertex.size = 1.*sqrt(degree(G4))+3,
         edge.width = E(G4)$w,
         vertex.frame.color="gray10",
         vertex.frame.width=.1,
         vertex.color = "skyblue",
         edge.color = gray.colors(100),
         layout= L,
         add = T, 
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.label.font=2,
         vertex.label.cex=0.4, 
         vertex.label.dist=0,
         edge.arrow.size = 0.01,
         edge.curved = 0.2,
         vertex.label = NA)
    for(k in 1:20){text(L[k, 1]+0.02, L[k, 2], Locations$POI[k], cex = 0.4, adj = 0)}
}
dev.off()
}

#### igraph plot color 19 routes
{
    pdf(file = "19RoutesNetW.pdf",
        height = 7, width =10)
    par(mar = c(0,0,0,0), mfrow = c(4,5))
    cols <- rainbow(19)
    labs <- c("C1 - pilot","C1 - muscle","C1 - muscle","C1 - muscle","C1 - muscle",
              "C2 - pilot","C2 - muscle","C2 - muscle","C2 - muscle","C2 - muscle",
              "C3 - pilot","C3 - muscle","C3 - muscle","C3 - muscle","C3 - muscle",
              "C4 - pilot","C4 - muscle","C4 - muscle","C4 - muscle","C4 - muscle")
    toK <- c(1:14,0, 15:19) 
    for(u in 1:20){
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        if(toK[u]>0){
        text(0,1.1, labs[u], cex = 1)
        plot(G, 
             vertex.size = 1.*sqrt(degree(G))+3,
             edge.width = E(G)$w,
             vertex.frame.color="gray80",
             vertex.frame.width=.1,
             vertex.color = "skyblue",
             edge.color = "gray90",
             layout= L,
             add = T, 
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.4, 
             vertex.label.dist=0,
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label = NA)
        GT <- subgraph.edges(G, 
                             which(J$Actor == toK[u]),
                                 delete.vertices = F)
        E(GT)$w <- count.multiple(GT)    
            plot(GT, add = T, 
                 vertex.frame.color="gray10",
                 vertex.frame.width=.1,
                 vertex.color = "skyblue",
                 edge.color = cols[toK[u]],
                 layout= L,
                 edge.width = E(GT)$w,
                 vertex.label.family="Helvetica",
                 vertex.label.color=NA,
                 vertex.label.font=2,
                 vertex.label.cex=0.4, 
                 vertex.label.dist=0.5,
                 edge.arrow.size = 0.01,
                 edge.curved = 0.2,
                 vertex.label = NA,
                 vertex.size = 1.9*sqrt(degree(G))+1.)
        }
            }
    dev.off()
}

#### igraph social relation for plotting relation matrix only
{
    simTrajPlot <- simTraj
    diag(simTrajPlot) <- 1
    simTrajPlot <- simTrajPlot*(simTrajPlot>0.5)
    RGPlot <- graph_from_adjacency_matrix(mode = "undirected",
                                          simTrajPlot,
                                          weighted = T,
                                          diag = F)
L <- layout_nicely(RGPlot)
pdf(file = "SocialN.pdf",
    height = 6, width = 6)
par(mar = c(0,0,0,0)+.1, mfrow = c(1,1))
plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(4))

for(k in 1:dim(AT)[2]){
    RGPlot <- set_vertex_attr(RGPlot,
                              name = names(AT)[k],
                              value = AT[,k])
}
par(mar = c(0,0,0,0))
plot(RGPlot,
     add = T,
     layout = L,
     vertex.frame.color="gray10",
     vertex.frame.width=.1,
     vertex.color =cols[V(RGPlot)$cell],
     edge.color = gray.colors(100),
     edge.width = E(RGPlot)$weight*5,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.label.font=2,
     vertex.label.cex=0.8 , 
     edge.arrow.size = 0.01,
     edge.curved = 0.2,
     vertex.label.dist=1.5,
     vertex.label.degree=pi/2, 
     vertex.label = V(RGPlot)$name_short,
     vertex.size = 5+ (V(RGPlot)$job == "pilot")*10)
for(k in 1:4){
    text(0.95, 1.05-k/10, paste("cell", k), col = cols[k], adj = 1, cex = 2)
}
points(0.55, 0.55, cex = 5, pch = 21, bg = "gray80")
text(0.95, 0.55, "pilots", col = 1, adj = 1, cex = 2)
dev.off()
}

#### igraph social relation and contacts
{
    simTrajPlot <- simTraj
    diag(simTrajPlot) <- 1
    simTrajPlot <- simTrajPlot*(simTrajPlot>0.5)
    RGPlot <- graph_from_adjacency_matrix(mode = "undirected",
                                          simTrajPlot,
                                          weighted = T,
                                          diag = F)
    L <- layout_nicely(RGPlot)
    pdf(file = "SocialNContactN.pdf",
        height = 6, width = 12)
    par(mar = c(0,0,0,0)+.1, mfrow = c(1,2))
    {
    plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
         asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(4))
    
    for(k in 1:dim(AT)[2]){
        RGPlot <- set_vertex_attr(RGPlot,
                                  name = names(AT)[k],
                                  value = AT[,k])
    }
    par(mar = c(0,0,0,0))
    plot(RGPlot,
         add = T,
         layout = L,
         vertex.frame.color="gray10",
         vertex.frame.width=.1,
         vertex.color =cols[V(RGPlot)$cell],
         edge.color = gray.colors(100),
         #edge.width = E(RGPlot)$weight*5,
         vertex.label.family="Helvetica",
         vertex.label.color="black",
         vertex.label.font=2,
         vertex.label.cex=0.8 , 
         edge.arrow.size = 0.01,
         edge.curved = 0.2,
         vertex.label.dist=1.5,
         vertex.label.degree=pi/2, 
         vertex.label = V(RGPlot)$name_short,
         vertex.size = 5+ (V(RGPlot)$job == "pilot")*10)
    for(k in 1:4){
        text(0.95, 1.05-k/10, paste("cell", k), col = cols[k], adj = 1, cex = 2)
    }
    points(0.55, 0.55, cex = 5, pch = 21, bg = "gray80")
    text(0.95, 0.55, "pilots", col = 1, adj = 1, cex = 2)
    }
    
    {
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(4))
        plot(ContactM,
             add = T,
             layout = L,
             vertex.frame.color="gray10",
             vertex.frame.width=.1,
             vertex.color =cols[V(ContactM)$cell],
             edge.color = gray.colors(100),
             #edge.width = E(RGPlot)$weight*5,
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.8 , 
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label.dist=1.5,
             vertex.label.degree=pi/2, 
             vertex.label = V(ContactM)$name_short,
             vertex.size = 5+ (V(ContactM)$job == "pilot")*10)
        for(k in 1:4){
            text(0.95, 1.05-k/10, paste("cell", k), col = cols[k], adj = 1, cex = 2)
        }
        points(0.55, 0.55, cex = 5, pch = 21, bg = "gray80")
        text(0.95, 0.55, "pilots", col = 1, adj = 1, cex = 2)
    }
    
    
    dev.off()
}

#### igraph social relation and contacts graph after threshold is 3/4
{
    simTrajPlot <- simTraj
    diag(simTrajPlot) <- 1
    simTrajPlot <- simTrajPlot*(simTrajPlot>3/4)
    RGPlot <- graph_from_adjacency_matrix(mode = "undirected",
                                          simTrajPlot,
                                          weighted = T,
                                          diag = F)
    RGPlot <- simplify(RGPlot)
    L <- layout_nicely(ContactM)
    pdf(file = "SocialNContactN.pdf",
        height = 6, width = 12)
    par(mar = c(0,0,0,0)+.1, mfrow = c(1,2))
    {
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(4))
        
        for(k in 1:dim(AT)[2]){
            RGPlot <- set_vertex_attr(RGPlot,
                                      name = names(AT)[k],
                                      value = AT[,k])
        }
        par(mar = c(0,0,0,0))
        plot(RGPlot,
             add = T,
             layout = L,
             vertex.frame.color="gray10",
             vertex.frame.width=.1,
             vertex.color =cols[V(RGPlot)$cell],
             edge.color = gray.colors(100),
             #edge.width = E(RGPlot)$weight*5,
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.8 , 
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label.dist=1.5,
             vertex.label.degree=pi/2, 
             vertex.label = V(RGPlot)$name_short,
             vertex.size = 5+ (V(RGPlot)$job == "pilot")*10)
        for(k in 1:4){
            text(0.95, 1.05-k/10, paste("cell", k), col = cols[k], adj = 1, cex = 2)
        }
        points(0.55, 0.55, cex = 5, pch = 21, bg = "gray80")
        text(0.95, 0.55, "pilots", col = 1, adj = 1, cex = 2)
    }
    
    {
        plot(0, col = NA, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1),
             asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        cols = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(4))
        plot(ContactM,
             add = T,
             layout = L,
             vertex.frame.color="gray10",
             vertex.frame.width=.1,
             vertex.color =cols[V(ContactM)$cell],
             edge.color = gray.colors(100),
             #edge.width = E(RGPlot)$weight*5,
             vertex.label.family="Helvetica",
             vertex.label.color="black",
             vertex.label.font=2,
             vertex.label.cex=0.8 , 
             edge.arrow.size = 0.01,
             edge.curved = 0.2,
             vertex.label.dist=1.5,
             vertex.label.degree=pi/2, 
             vertex.label = V(ContactM)$name_short,
             vertex.size = 5+ (V(ContactM)$job == "pilot")*10)
        for(k in 1:4){
            text(0.95, 1.05-k/10, paste("cell", k), col = cols[k], adj = 1, cex = 2)
        }
        points(0.55, 0.55, cex = 5, pch = 21, bg = "gray80")
        text(0.95, 0.55, "pilots", col = 1, adj = 1, cex = 2)
    }
    
    
    dev.off()
}

#### histograms of simulated vs observed
{
hist(SimGAll, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
points(c(InterGAll[k], InterGAll[k]),
       c(0,10000),
       type = "l")
}

hist(SimGAllPM, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterGAllPM[k], InterGAllPM[k]),
           c(0,10000),
           type = "l")
}

hist(SimGLastYear, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterGLastY[k], InterGLastY[k]),
           c(0,10000),
           type = "l")
}

hist(SimGLastYearPM, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterGLastYPM[k], InterGLastYPM[k]),
           c(0,10000),
           type = "l")
}

hist(SimTrajAll, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterTrajAll[k], InterTrajAll[k]),
           c(0,10000),
           type = "l")
}

hist(SimTrajAllPM, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterTrajAllPM[k], InterTrajAllPM[k]),
           c(0,10000),
           type = "l")
}

hist(SimTrajLastYear, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterTrajLastYear[k], InterTrajLastYear[k]),
           c(0,10000),
           type = "l")
}

hist(SimTrajLastYearPM, breaks = c(0:100)/100, col = 2)
for(k in 1:4){
    points(c(InterTrajLastYearPM[k], InterTrajLastYearPM[k]),
           c(0,10000),
           type = "l")
}

hist(as.numeric(SimulD), 
     main = "Metric of difference between networks",
     n = 100,
     col = "tomato")
points(c(120, 120),
       c(0,90),
       type = "l", lwd = 3, col = "cyan3")

hist(SimulB, 
     n = 100,
     main = "Metric of difference between networks",
     col = "tomato")
points(c(391.5, 391.5),
       c(0,90),
       type = "l", lwd = 3, col = "cyan3")

}
}
