library(maptools)
mM <- c("/home/thomas/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra1.shp",
        "/home/thomas/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra2.shp",
        "/home/thomas/Documents/University/Bachelor/FVA/Daten/Messpunkte/Messpunkte_Petra3.shp")

## work flow
## read manual measurements
lmM <- readManualMeasure(files = mM)

## cut data
cmM <- cutData(uncutfiles =  list(lmM[[1]], lmM[[2]]), ## Messpunkte_Petra3 has different allocation
         columns = list("x_coord" = 9, "y_coord" = 10, "z_coord" = 6, "class" = 3), 
         omit.class= FALSE)

## calculate means of cut manual Measure
mcmM <- meanManualMeasure(manualMeasure = cmM)

## modelValues
    coordinatesModel <- mcmM[, 1:2]
    ## read model
    MTb1 <- read.table(file = "/home/thomas/Documents/University/Bachelor/FVA/Daten/Models/MT_balanced1_buf5m.txt", sep = ";")
    MTb2 <- read.table(file = "/home/thomas/Documents/University/Bachelor/FVA/Daten/Models/MT_balanced2_buf5m.txt", sep = ";")
    MTb3 <- read.table(file = "/home/thomas/Documents/University/Bachelor/FVA/Daten/Models/MT_balanced3_buf5m.txt", sep = ";")
    
    ## not run
        model2D <- .calculateModelValues(coordinates = coordinatesModel, model = MTb1, method = "2D")
    ## not run 
        modelIDW <- .calculateModelValues(coordinates = coordinatesModel, model = MTb1, method = "IDW", 
            idw = list("p" = 2, "m" = 5, "rad" = 5))
    
    cModh <- calcModelHeights(coordinates = coordinatesModel, 
                  model = list(MTb1 = MTb1, MTb2 = MTb2, MTb3 = MTb3), 
                  method = "2D", 
                  idw = list("p" = 2, "m" = 5, "rad" = 5))

## calculate error between manual measure and model values
errorMod <- errorModel(manual = mcmM, model = cModh)

## test if the errors are normally distributed
errorNormalTest(errorMod, hist = TRUE, ksTest = FALSE, 
                qq = FALSE, classes = FALSE)
errorNormalTest(errorMod, hist = TRUE, ksTest = FALSE,
                qq = TRUE, classes = FALSE)
errorNormalTest(errorMod, hist = FALSE, ksTest = TRUE,
                qq = FALSE, classes = FALSE)
