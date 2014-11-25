## to do
## a collection to store manual measurements
## a collection to store model values
## setClass("manualMeasure", )

#'@title Read manual measurements
#'@name readManualMeasure
#'@author Thomas Naake <naake@@stud.uni-heidelberg.de>
#'@usage readManualMeasure(files)
#'@examples \dontrun{
#'readManualMeasure(files)
#'}
#'@description Read .shp files using maptools::readShapePoints() of a 
#'vector with given file paths pointing to the favoured files. The function 
#'will read the files and create an list object which contains the files.
#'@param files a vector with characters giving the file paths.
#'@return a list which contains the individual models 
#'@export
readManualMeasure <- .readModel <- function(files) {
    .len <- length(files)
    ans <- list()
    for (i in 1:.len)
        ans[[i]] <- readShapePoints(files[i])
    return(ans)
}

#'@title cutData
#'@author Thomas Naake <naake@@stud.uni-heidelberg.de>
#'@usage cutData(uncutfiles, 
#'          columns = list("x_coord" = 1, "y_coord" = 2, 
#'                              "z_coord" = 3, "class" = 4),
#'          omit.class = FALSE)
#'@examples \dontrun{
#'cutData(uncutfiles, 
#'          columns = list("x_coord" = 1, "y_coord" = 2, 
#'                              "z_coord" = 3, "class" = 4),
#'          omit.class = FALSE)
#'}
#'@description a function to extract x, y, z and class from raw data.
#'@param uncutfiles a list of raw data e.g. uncut manual measures
#'columns a list given the column of x-coordinates, y-coordinates, z-coordinates 
#'and class
#'omit.class logical, if \code{TRUE} a fourth column will be added containing 
#'the assigned class
#'@return a list of data frames with x-coordinates, y-coordinates, z-coordinates
#'and class (optional).
#'@export
## function to extract x, y, z and class from raw data
## returns a list of data frames
cutData <- function(uncutfiles, 
    columns = list("x_coord" = 1, "y_coord" = 2, "z_coord" = 3, "class" = 4), omit.class = FALSE) {
    .len <- length(uncutfiles)
    ans <- list()
    for (i in 1:.len) {
        .data <- uncutfiles[[i]]
       ## .newAns <- data.frame(row.names = 1:dim(.data)[1])
        .newAns <- matrix(nrow = dim(.data)[1], ncol = ifelse(omit.class, 3, 4))
        .newAns <- as.data.frame(.newAns)
        .newAns[, 1] <- as.data.frame(.data[, columns$x_coord])[, 1]
        .newAns[, 2] <- as.data.frame(.data[, columns$y_coord])[, 1]
        .newAns[, 3] <- as.data.frame(.data[, columns$z_coord])[, 1]
        if (!omit.class)
            .newAns[, 4] <- as.data.frame(.data[, columns$class])[, 1]
        ans[[i]] <- .newAns 
    }  
    return(ans)
}

#'@title meanManualMeasure
#'@author Thomas Naake <naake@@stud.uni-heidelberg.de>
#'@description Calculate mean values of manual measurements. 
#'@param manualMeasure a list of manual measurements, each list entry is of 
#' class \code{data.frame} which contains x-coordinates (first column), 
#' y-coordinates (second column), z-coordinates (third column) and
#' (optionally) class (fourth column)
#'@return a data.frame containing x-coordinates, y-coordinates, z-coordinates and
#' class
#'@export
meanManualMeasure <- function(manualMeasure) {
    if (!is.list(manualMeasure))
        stop("argument is not a list")
    if (length(manualMeasure) == 0)
        stop ("argument is of length 0")
    
    .len <- length(manualMeasure)
    .meta <- manualMeasure[[1]]
    
    .sub <- rep(NA, dim(.meta)[1])
    ans <- data.frame("x_coord" = .sub, "y_coord" = .sub, 
                      "mean" = .sub, "class" = .sub)
    ## load column 1 with x-coordinates
    ans[, 1] <- .meta[, 1]
    ## load column 2 with y-coordinates
    ans[, 2] <- .meta[, 2]
    ## load column 4 with class
    ans[, 4] <- .meta[, 4]
    
    .mean <- lapply(manualMeasure, "[", 3)
    .mean <- as.data.frame(.mean)
    .mean <- apply(.mean, MARGIN = 1, mean)
    ## load column 3 with mean heights 
    ans[, 3] <- .mean
    return(ans)
}


## function to calculate model calues at given coordinates
## model takes x/y/z coords of a model
## coordinates: data frame of coordinates (1st column x, 2nd column y values)
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
        colnames(.mat) <- c("z_mod_2D", "min_2Ddist_z_mod")
        ## END: 2D ##
    } else {
        ## START: IDW ## 
        ## exponent according to Shepard, 1968
        .p <- idw$p 
        ## number of points which are used
        .m <- idw$m
        ## use points which are iwd$rad [m] far from the coordinate point
        .rad <- idw$rad
        
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
                .Cpn <- .Cpn[1:length(.Cp)]
            else 
                .Cpn <- .Cpn[1:.m]
            ## vector for weighting factors
            .method <- vector("numeric", length(.Cpn))
            if (length(.Cpn) != 0 && !is.na(.Cpn[1])) {
                for (j in 1:length(.Cpn)) {
                    if (.Cpn[j] > 0 && .Cpn[j] <= .rad/3)
                        .method[j] <- "m1"
                    if (.Cpn[j] > .rad/3 && .Cpn[j] <= .rad)
                        .method[j] <- "m2"
                    if (.Cpn[j] > .rad)
                        .method[j] <- "m3"
                }
                ## matrix whose diagonal axis bears indices of points .Cpn i 
                .matInd <- matrix(NA, length(.Cpn), length(.Cpn))
                for (j in 1:length(.Cpn)) {
                    if (j == 1) {
                        .vec <- which(.Cpn[j] == .dist)[1]
                        .vec <- c(rep(NA, j - 1), 
                                  .vec, 
                                  rep(NA, 10 - length(.vec)))
                        .matInd[j, j] <- .vec[j]
                        n <- 2 ## was 2 ??
                    }
                    if (j > 1) {
                        .vec <- which(.Cpn[j] == .dist)
                        if (.vec[n - 1] == .matInd[j - 1, j - 1]) {
                            .matInd[j, j] <- .vec[n]
                            n <- n + 1
                            if (n == length(.vec) + 1)
                                n <- 2 
                        } else 
                            .matInd[j, j] <- .vec[1]; n <- 2
                    
                    }
                } ## end for loop
                ## vector which contains indices of .matInd
                .ind <- vector("numeric", dim(.matInd)[1])
                for (j in 1:length(.ind))
                    .ind[j] <- .matInd[j, j]
                ## vector with weighting factors according to Shepard (1968)
                .wf <- vector("numeric", length(.Cpn))
                .n <- 1
                for (j in .ind) {
                    .n <- .n 
                    if (.method[.n] == "m1")
                        .wf[.n] <- 1/ sqrt( 
                            (coordinates[i, 1] - model[j, 1])^2 + (coordinates[i, 2] - model[j, 2])^2
                        ) ^ .p
                    if (.method[.n] == "m2")
                        .wf[.n] <- (27/(4*.rad) * (sqrt(
                            (coordinates[i, 1] - model[j, 1])^2 + 
                                (coordinates[i, 2] - model[j, 2])^2) / .rad - 1))^.p
                    if (.method[.n] == "m3")
                        .wf[.n] <- 0
                    .n <- .n + 1
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
        colnames(.mat) <- c("z_mod_IDW", "max_2Ddist_z_mod")
    } ## END: IDW ## 
    return(.mat)
}

#'@title calcModelHeights
#'@author Thomas Naake <naake@@stud.uni-heidelberg.de>
#'@usage calcModelHeights(coordinates, model, method = c("2D", "IDW"),
#'                          idw = list("p" = 2, "m", "rad" = 5))
#'@examples \dontrun{
#'calcModelHeights(coordinates, model, method = c("2D", "IDW"),
#'                          idw = list("p" = 2, "m", "rad" = 5))
#'}
#'@description a function to interpolate z-coordinates at given x- and y-coordinates
#'@param coordinates a matrix or data.frame containing x- and y-coordinates
#'model a list of \code{data.frame}s with model values with x-coordinates,
#'y-coordinates, z-coordinates and classes
#'method character, use "2D" or "IDW"
#'idw a list with parameters for IDW interpolarion
#'@return a list of data frames with x-coordinates, y-coordinates, z-coordinates
#'and class.
#'@details see Shepard for IDW.
#'@export
## a function to create and concatenate heights from different models 
## models: 1st column x, 2nd column y, 3rd column z
calcModelHeights <- function(coordinates, model, method = c("2D", "IDW"),
                              idw = list("p" = 2, "m" = 5, "rad" = 5)) {
    
    method <- match.arg(method)
    if (!(is.data.frame(coordinates) || is.matrix(coordinates)))
        stop("argument coordinates is neither a data frame nor a matrix")
    if (!is.list(model))
        stop("argument model is not a list")
    if (!min(nchar(names(model))))
        stop("model is not a fully named list")
    
    .lenMod <- length(model)

    for (i in 1:.lenMod) {
        .values <- .calculateModelValues(coordinates = coordinates, 
                                        model = model[[i]], 
                                        method = method, 
                                        idw = list("p" = idw$p, 
                                                   "m" = idw$m, 
                                                   "rad" = idw$rad))
        model[[i]] <- .values
        print(paste0("model ", names(model)[i], ": ", i, " / ", .lenMod))
    }
   return(model)
}

#'@title Calculate error values 
#'@usage errorModel(manual, model)
#'@author Thomas Naake <naake@@stud.uni-heidelberg.de>
#'@examples \dontrun{
#'errorModel(manual, model)
#'}
#'@description a function to calculate the errors of model and manual 
#'measurement values
#'@param manual a matrix or data.frame containing x-, y-, z-coordinates and classes
#'model a list containing model values
#'@return a list of data frames with x-coordinates, y-coordinates, z-coordinates
#'and class.
#'@details see Shepard for IDW.
#'@export
errorModel <- function(manual, model) {
    
    if (!(is.data.frame(manual) || is.matrix(manual)))
        stop("argument manual is neither a data frame nor a matrix")
    if (!is.list(model))
        stop("argument model is not a list")
    
    .points <- dim(manual)[1]
    .model <- lapply(cModh, "[", 1:.points) 
    .manual <- manual[, 3]
    
    ans <- list()
    for (i in 1:length(.model)) {
        ans[[i]] <- manual
        ans[[i]][, 3] <- .manual - .model[[i]]
        names(ans)[i] <- names(.model)[i]
    }
    return(ans)
}

#'@name errorNormalTest
#'@title Kolmogorov-Smirnow and plot
#'@author Thomas Naake <naake@@stud.uni-heidelberg.de>
#'@param error a list of errors to compare, each list entry is of class
#' \code{data.frame} which contains x-coordinates (first column), 
#' y-coordinates (second column), z-coordinates (third column) and
#' (optionally) classes (fourth column)
#' @export
errorNormalTest <- function(error, 
        hist = TRUE, 
        ksTest = FALSE, 
        qq = FALSE, classes = FALSE) {
    
    if (!is.list(error))
        stop("argument is not a list")
    if (length(error) == 0)
        stop("argument is of length 0")
    
    .lenError <- length(error)
    
    for (i in 1:.lenError) {
        .nameE <- names(error[i])
        print(.nameE)
        .error <- error[[i]][, 3]
        .bins <- seq(floor(min(.error)), ceiling(max(.error)), by = 1)  
        set.seed(1)
        ## model normal distribution
        .r <- rnorm(n = 1000, mean = mean(.error), sd = sd(.error))
        
        if (hist) {
            hist(.error, breaks = .bins, 
                col = rgb(0, 0, 1, 1/4), 
                freq = F, 
                main = paste0("Histogram of e", .nameE))
            points(density(.r), type = "p")
            if (qq || (!qq && i < .lenError))
                readline(prompt = "Press <Enter> to continue... ")
        }
        
        if (qq) {
            qqplot(x = .r, y = .error, ylab = "Sample Quantiles", 
                xlab = "Theoretical Quantiles", main = paste0("Normal Q-Q plot of e", .nameE))
            points(sort(.r), sort(.r), col = "red", type = "l", lty = 2, lwd = 2)
            if (ksTest || (!ksTest && i < .lenError))
                readline(prompt = "Press <Enter> to continue... ")
        }
        
        if (ksTest) {
            .ks <- ks.test(.error, .r) ## null hypothesis: same distribution
            .ks$data.name <- paste0("error of ", .nameE, 
                                " and modelled values of normal distribution")
            print(.ks)          
            if (i < .lenError)
                readline(prompt = "Press <Enter> to continue... ")
        }     
    } ## end for loop
}

#'@name stat
#'@title Non-parametric statistics for non-normally distributed errors
#'@return A list. Each list entry comprises a data frame which bears the 
#'statistical parameters
#'@param error A list containing error values of models, each list entry is of class
#' \code{data.frame} which contains x-coordinates (first column), 
#' y-coordinates (second column), z-coordinates (third column) and
#' (optionally) classes (fourth column)
#'@param cfi logical, should confidence intervall be returned
#'@param classes logical, indicating if classes will be used concerning calculation 
#'of the statistics
#'@export
stat <- function(error, cfi = TRUE, classes = FALSE) {
    
    .lenError <- length(error)
    ans <- error
    .df <- as.data.frame(matrix(ncol = ifelse(cfi, 3, 1), nrow = 5))
    colnames(.df) <- ifelse(cfi, c("value", "cfi_l", "cfi_r"), c("value"))
    rownames(.df) <- c("68%quantile", "95%quantile", "median", "NMAD", "max|h|")

    for (i in 1:.lenError) {
        
        ## index each list entry
        .error <- error[[i]]
        .n <- dim(.error)[1]
        ## 68% quantile
        .q68 <- vector("numeric", length = 2000)
        .q68[2000] <- quantile(abs(.error[, 3]), 0.683)
        .df[1, 1] <- .q68[2000]
        ## 95% quantile
        .q95 <- vector("numeric", length = 2000)
        .q95[2000] <- quantile(abs(.error[, 3]), 0.95)
        .df[2, 1] <- .q95[2000]
        ## median
        .med <- vector("numeric", length = 2000)
        .med[2000] <- median(.error[, 3])
        .df[3, 1] <- .med[2000]
        ## NMAD
        .nmad <- vector("numeric", length = 2000)
        .nmad[2000] <- 1.4826 * median(abs(.error[, 3] - median(.error[, 3])))
        .df[4, 1] <- .nmad[2000]
        ## hmax
        .df[5, 1] <- max(abs(.error[, 3]))
        
        if (cfi) {
            for (j in 1:1999) {
                ## calculate via bootstrapping 95% confidence intervall for parameter
                .sample <- sample(.error[, 3], size = .n, replace = T)
                ## 68.3 % quantile
                .q68[j] <- quantile(abs(.sample), 0.683)
                ## 95% quantile
                .q95[j] <- quantile(abs(.sample), 0.95)
                ## median
                .med[j] <- median(.sample)
                ## NMAD
                .nmad[j] <- 1.4826 * median(abs(.sample - median(.sample)))
            }
            ## 68% quantile
            .df[1, 2] <- quantile(.q68, 0.025)
            .df[1, 3] <- quantile(.q68, 0.975)
            ## 95% quantile            
            .df[2, 2] <- quantile(.q95, 0.025)
            .df[2, 3] <- quantile(.q95, 0.975)
            ## median    
            .df[3, 2] <- quantile(.med,0.025)
            .df[3, 3] <- quantile(.med,0.0975)
            ## NMAD
            .df[4, 2] <- quantile(.nmad, 0.025)
            .df[4, 3] <- quantile(.nmad, 0.975)
        }
        
        ans[[i]] <- .df
        
    }
    return(ans)
    
}

#'@name plotStats
#'@title Visualize nonparametric parameters
#'@export
plotStats <- function(stats, param = c("median", "NMAD", "max|h|")) {
   
    param <- match.arg(arg = param, 
        choices = c("median", "NMAD", "max|h|"), several.ok = TRUE)
    
    .valq68 <- unlist(lapply(stats, "[[", 1, 1))
    .valq95 <- unlist(lapply(stats, "[[", 2, 1))
    .valmed <- unlist(lapply(stats, "[[", 3, 1))
    .valnmad <- unlist(lapply(stats, "[[", 4, 1))
    .valhmax <- unlist(lapply(stats, "[[", 5, 1))
    
    .maxval <- max(.valq95)
    if ("max|h|" %in% param)
        .maxval <- c(.maxval, .valhmax)
    .maxval <- max(.maxval)
        
    barplot(.valq95, border = FALSE, beside = TRUE, ylab=c("|dh|"), 
        col="lightgrey", axes = T, xpd = FALSE, 
        ylim=c(0, .maxval))
    barplot(.valq68, border = FALSE, beside = TRUE, add = T, col = "darkgrey",
        axes = FALSE, xpd = FALSE)
    
    .lenError <- length(stats)
    x0 <- seq(0.2, .lenError * 3, 1.2)
    x1 <- seq(1.2, .lenError * 3, 1.2)
    
    
    for (i in 1:.lenError) {
        segments(x0 = x0[i], x1 = x1[i], y0 = .valq68[i], lwd = 1, col = "lightgrey")
        segments(x0 = x0[i], x1 = x1[i], y0 = .valq95[i], lwd = 1, col = "darkgrey")
        if ("max|h|" %in% param)
            segments(x0 = x0[i], x1 = x1[i], y0 = .valhmax[i], lwd = 1, col = "black")
        if ("median" %in% param)
            segments(x0 = x0[i], x1 = x1[i], y0 = .valmed[i], lwd = 2, col = "chocolate1")
        if ("NMAD" %in% param)
            segments(x0 = x0[i], x1 = x1[i], y0 = .valnmad[i], lwd = 2, col = "chartreuse3")
    }    
}