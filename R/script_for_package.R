
mM <- c("/home/thomas/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra1.shp",
        "/home/thomas/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra2.shp",
        "/home/thomas/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra3.shp")
##test <- readShapePoints(mM[[1]]) ## read manual Measurements
##.newAns[, 1] <- as.data.frame(.data[, 9])[, 1] ## cut data
##.newAns[, 2] <- as.data.frame(.data[, 10])[, 1]
##.newAns[, 3] <- as.data.frame(.data[, 6])[, 1]
##.newAns[, 4] <- as.character(as.data.frame(.data[, 3])[, 1])
##.newAns <- list(.newAns) ## for calculation of means


## work flow
## read manual measurements
lmM <- .readManualMeasure(mM)

## cut data
cmM <- .cutData(uncutfiles =  list(lmM[[1]], lmM[[2]]), ## Messpunkte_Petra3 has different allocation
         columns = list("x_coord" = 9, "y_coord" = 10, "z_coord" = 6, "class" = 3), 
         omit.class= FALSE)

## calculate means of cut manual Measure
mcmM <- .meanManualMeasure(cmM)

## modelValues
    coordinatesModel <- mcmM[, 1:2]
    ## read model
    MTb1 <- read.table(file = "/home/thomas/Documents/University/Bachelor/FVA/Daten/Models/MT_balanced1_buf5m.txt", sep = ";")
    
    ## not run
        model2D <- .calculateModelValues(coordinates = coordinatesModel, model = MTb1, method = "2D")
    ## not run 
        modelIDW <- .calculateModelValues(coordinates = coordinatesModel, model = MTb1, method = "IDW", 
            idw = list("p" = 2, "m" = 5, "rad" = 5))
    
    .calcModelHeights(coordinates = coordinatesModel, 
                  model = list(MTb1), 
                  method = "2D", 
                  idw = list("p" = 2, "m" = 5, "rad" = 5))
