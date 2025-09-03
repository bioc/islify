#The point of this function is to exclude regions with much higher than
#expected yields. We talk about "background" below, but that is a historic
#remnant, from a time when it was the background that was calculated. 
medianChunkFilter <- function(locFileRaggedIsles){
    origDims <- dim(locFileRaggedIsles)
    
    #Now, which side is longer? 
    fileRows <- nrow(locFileRaggedIsles)
    fileCols <- ncol(locFileRaggedIsles)
    rowColRatio <- fileRows/fileCols
    
    if(rowColRatio < 1){
        nChunksRow <- 4
        nChunksCol <- round(4/rowColRatio)
    } else {
        nChunksRow <- round(4*rowColRatio)
        nChunksCol <- 4
    }
    rowNumList <- split(seq(1, fileRows), 
                        ceiling(seq(1, fileRows)/(fileRows/nChunksRow)))
    colNumList <- split(seq(1, fileCols), 
                        ceiling(seq(1, fileCols)/(fileCols/nChunksCol)))
    #Now, we calculate the background for each of the squares
    backgroundMat <- matrix(NA, nChunksRow, nChunksCol)
    ragged01 <- as.matrix(locFileRaggedIsles)
    ragged01[which(ragged01 != 0)] <- 1
    for(i in seq_along(rowNumList)){
        for(j in seq_along(colNumList)){
            rows <- rowNumList[[i]]
            cols <- colNumList[[j]]
            backgroundMat[i,j] <- mean(ragged01[rows,cols])
        }
    }
  
    backgroundAll <- mean(backgroundMat)
    backgroundSd <- sd(backgroundMat)
    locFileDark <- locFileRaggedIsles
    for(i in seq_along(rowNumList)){
        for(j in seq_along(colNumList)){
            rows <- rowNumList[[i]]
            cols <- colNumList[[j]]
            locBackground <- backgroundMat[i,j]
            if(locBackground > backgroundAll+(2*backgroundSd)){
                locFileDark[rows,cols] <- 0
            }
        }
    }
    
    bigIslesPixels <- islandPixels(locFileRaggedIsles)
    darkIslesPixels <- islandPixels(locFileDark)
    if(length(bigIslesPixels) > length(darkIslesPixels)){
        overlappingIsles <- which(names(bigIslesPixels) %in%
                                      names(darkIslesPixels))
        missingIsles <- 
            names(bigIslesPixels)[-overlappingIsles]
        bigIslePixelsRed <- bigIslesPixels[overlappingIsles]
        if(any(unlist(bigIslePixelsRed) > 
               unlist(darkIslesPixels))){
          smallerIsles <- names(bigIslePixelsRed)[which(unlist(bigIslePixelsRed) > 
                                                          unlist(darkIslesPixels))]
          darkIsles <- c(missingIsles, smallerIsles)
        } else {
          darkIsles <- missingIsles
        }
    } else {
        darkIsles <- names(bigIslesPixels)[which(unlist(bigIslesPixels) > 
                                              unlist(darkIslesPixels))]
    }
    
    if(length(darkIsles) > 0){
      
        #Now,  we create the new, reduced file. 
        sm <- as.data.frame(Matrix::summary(locFileRaggedIsles))
        colnames(sm) <- c("row", "column", "value")
        sm$value[which(sm$value %in% as.numeric(darkIsles))] <- 0
        # Here, we remove all new zeros.
        smSmall <- sm[which(sm$value > 0), ]
        locSM <- sparseMatrix(
            i = smSmall$row,
            j = smSmall$column,
            x = smSmall$value,
            dims = origDims
        )
    } else {
        locFileRaggedIsles
    }
}