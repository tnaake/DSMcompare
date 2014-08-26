## a collection to store manual measurements
## a collection to store model values
## setClass("manualMeasure", )

## functions to read in shp data (vectorized, contains path)
.readManualMeasure <- .readModel <- function(files) {
    .len <- length(files)
    ans <- list()
    for (i in 1:.len)
        ans[[i]] <- readShapePoints(files)
    return(ans)
}

# .cutData <- function(uncutfiles, 
#     columns = list("x_coord" = 1, "y_coord" = 2, "z_coord" = 3, "class" = 4)) {
#     
#     
# }

#' @title
#' @aliases
#' @description
#' @param manualMeasure a list of manual measurements, each list entry is of 
#' class \code{data.frame} which contains x-coordinates (first column), 
#' y-coordinates (second column), z-coordinates (third column) and
#' (optionally) class (fourth column)
#' @export

mM <- list("/Users/thomasnaake/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra1.shp",
    "/Users/thomasnaake/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra2.shp",
    "/Users/thomasnaake/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra3.shp")

.meanManualMeasure <- function(manualMeasure) {
    if (!is.list(manualMeasure))
        stop("argument is not a list")
    if (length(manualMeasure) == 0)
        stop ("argument is of length 0")
    .len <- length(manualMeasure)
    .meta <- manualMeasure[[1]]
    .mean <- lapply(manualMeasure, "[", 3)
    .mean <- as.data.frame(.mean)
    .mean <- apply(.mean, MARGIN = 1, mean)
    
    .sub <- rep(NA, dim(.meta)[1])
    ans <- data.frame("x_coord" = .sub, "y_coord" = .sub, 
                      "mean" = .sub, "class" = .sub)
    ## load column 1 with x-coordinates
    ans[, 1] <- .meta[, 1]
    ## load column 2 with y-coordinates
    ans[, 2] <- .meta[, 2]
    ## load column 3 with mean heights 
    ans[, 3] <- .mean
    ## load column 4 with class
    ans[, 4] <- .meta[, 3]
    
    return(ans)
}



.calculateModelValues <- function(coordinates, model, method = c("2D", "IDW"),
                                    idw = list("p" = 2, "m" = 5, "rad" = 5)) {
    .lenCoord <- dim(coordinates)[1]
    .mat <- matrix(NA, ncol = 2, nrow = .lenCoord)
    
    .lenModel <- dim(model)[1]
    
    ## START: 2D ##
    if (method == "2D") {
        .zMinDist <- vector(mode = "numeric", length = .lenCoord)
        .minDist <- vector(mode = "numeric", length = .lenCoord)
        for (i in 1:.lenCoord) {
            .dist <- sqrt((coordinates[i, 1] - model[, 1])^2 + 
                                (coordinates[i, 2] - model[, 2])^2)
            .minDist[i] <- min(.dist)
            .zMinDist[i] <- model[which.min(.dist), 3] 
        }
        .mat[, 1] <- .zMinDist
        .mat[, 2] <- .minDist
        ## END: 2D ##
    } else {
        ## START: IDW ## 
        ## exponent according to Shepard, 1968
        .p <- idw$p 
        ## number of points which are used
        .m <- iwd$m
        ## use points which are iwd$rad m far from the coordinate point
        .rad <- iwd$rad
        
        ## vector which will contain height
        .zIDW <- vector("numeric", length = .lenCoord)
        
        for (i in 1:.lenCoord) {
            .dist <- sqrt((coordinates[i, 1] - model[, 1])^2 + 
                             (coordinates[i, 2] - model[, 2])^2)
            .ind <- which(.dist <= .rad)
            .Cp <- .dist[.ind] ## Cp = {Di | d <= r}, r = rad
            ## ordering of the Di by increasing distance from P
            .Cpn <- sort(.Cp) 
            if (length(.Cp) < .m)
                .Cpn <- .Cpn[1:length(.Cpn)]
            else
                .Cpn <- .Cpn[1:m]
            ## vector for weighting factors
            .method <- vector("numeric", length(.Cpn))
            if (length(.Cpn) != 0 && !is.na(.Cpn[1])) {
                for (j in 1:length(.Cpn)) {
                    if (.Cpn[j] > 0 && .Cpn[j] <= .rad/3)
                        .method <- "m1"
                    if (.Cpn[j] > .rad/3 && .Cpn[j] <= .rad)
                        .method <- "m2"
                    if (.Cpn[j] > .rad)
                        .method[j] <- "m3"
                }
                ## matrix whose diagonal axis bears indices of points .Cpn i 
                .matInd <- matrix(NA, length(.Cpn), length(.Cpn))
                for (j in 1:length(.Cpn)) {
                    if (q == 1) {
                        .vec <- which(.Cpn[j] == .dist)[1]
                        .vec <- c(rep(NA, j - 1), 
                                  .vec, 
                                  rep(NA, 10 - length(.vec)))
                        .matInd[j, j] <- .vec[j]
                        .n <- 2
                    }
                    if (j > 1) {
                        .vec <- which(.Cpn[j] == .dist)
                        if (vec[n - 1] == .matInd[j - 1, j - 1]) {
                            .matInd[j, j] <- vec[j]
                            n <- n + 1
                            if (n == length(.vec) + 1)
                                n <- 2 
                        } else {
                            .matInd[j, j] <- vec[1]
                            n <- 2
                        }
                    }
                }
                ## vector which contains indices of .matInd
                .ind <- vector("numeric", dim(.matInd)[1])
                for (j in 1:length(.ind))
                    .ind[j] <- .matInd[j, j]
                ## vector with weighting factors according to Shepard (1968)
                .wf <- vector("numeric", length(.Cpn))
                .n <- 1
                for (j in .ind) {
                    .n <- .n 
                    if (.method[n] == "m1")
                        .wf[.n] <- 1/sqrt((.coordinates[i, 1] - model[j, 1])^2 +
                            (.coordinates[i, 2] - model[j, 2])^2)^.p
                    if (.method[.n] == "m2")
                        .wf[.n] <- (27/(4*.rad) * (sqrt(
                            (.coordinates[i, 1] - model[j, 1])^2 + 
                                (.coordinates[i, 2] - model[j, 2])^2) / rad - 1))^.p
                    if (.method[.n] == "m3")
                        .wf[.n] <- 0
                    n <- n + 1
                }
                ## if distance of next point == 0 to reference point (identical)
                ## take this height
                if (min(.Cpn) == 0 && !is.na(.Cpn[1]))
                    .zIDW[i] <- model[which(min(.Cpn) == .dist), 3]
                else
                    .zIDW[i] <- sum(.wf * model[.ind, 3]) / sum(.wf)
                
                ## write IDW height and longest distance to .mat
                if (!is.na(.Cpn[1]))
                    .mat[i, 1] <- .zIDW[i]
                else
                    .mat[i, 1] <- NA
                
                if (!is.na(.Cpn[1]))
                    .mat[i, 2] <- max(.dist[which(.dist <= max(.Cpn))[1:length(.Cpn)]])
                else 
                    .mat[i, 2] <- NA
            }
                
        }
        ## END: IDW ## 
}

.errorModel <- function(manual, model) {
    ans <- manual - model
}

## a function to create and concatenate heights from different models 
.calcModelHeights <- function(coordinates, model, method = c("2D", "IDW"),
                              idw = list("p" = 2, "m" = 5, "rad" = 5)) {
    
    .lenCoord <- dim(coordinates)[1]
    .mat <- (NA, nrow = .lenCoord, ncol = 2*length(model))
    for (i in 1:length(model)) {
        .values <- .meanManualMeasure(coordinates = coordinates, model = model, 
                        method = method, idw = list(idw$p, idw$m, idw$rad))
        .mat[, i] <- .values[, 1]
        .mat[, i + 1] <- .values[, 2]
    }
    return(.mat)
}


#'@name Kolmogorov-Smirnow and plot
#'@author Thomas Naake <thomasnaake@@gmx.de>
#'@param manual a list of length 1 containing values to compare
#'modell a list of modells to compare, each list entry is of class
#' \code{data.frame} which contains x-coordinates (first column), 
#' y-coordinates (second column), z-coordinates (third column) and
#' (optionally) class (fourth column)
#'@aliases
# .ksPlot <- function(manual, model, class = FALSE) {
#     if (!is.list(modell))
#         stop("argument is not a list")
#     if (length(modell) == 0)
#         stop("argument is of length 0")
#     
#     modell <- modell[[1]]
#     .error <- modell - mean(modell)
#     xyplot()
#     
# }