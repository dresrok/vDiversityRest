getIvi <- function(url){
  library(jsonlite);
  data <- fromJSON(url);
  forest <- data.frame(table(data$especie), stringsAsFactors=FALSE, row.names = NULL);
  names(forest) <- c('Species', 'absDen');
  forest$Species <- as.character(forest$Species);
  forest$relDen <- (forest$absDen/sum(forest$absDen))*100;
  speciesFreq <- data.frame(table(data$transecto, data$subparcela, data$especie));
  names(speciesFreq) <- c('transect', 'quadrant', 'Species', 'absFreq');
  dtFreq <- data.frame(table(speciesFreq$Species[speciesFreq$absFreq > 0]));
  forest$absFreq <- dtFreq$Freq;
  forest$relFreq <- (forest$absFreq/sum(forest$absFreq))*100;
  data$dap[is.na(data$dap)] <- 0;
  data$dA <- (pi/40000) * (data$dap^2);
  dtDbh <- data.frame(tapply(data$dA, data$especie, sum), row.names=NULL);
  names(dtDbh) <- c('sumDbh');
  forest$absDom <-  dtDbh$sumDbh;
  forest$relDom <- (forest$absDom/sum(forest$absDom))*100;
  forest$IVI <- forest$relDen + forest$relFreq + forest$relDom;
  total <- data.frame(
    Species = "Total", absDen = sum(forest$absDen), relDen = sum(forest$relDen),
    absFreq = sum(forest$absFreq), relFreq = sum(forest$relFreq),
    absDom = sum(forest$absDom), relDom = sum(forest$relDom), IVI = sum(forest$IVI)
  );
  forest <- rbind(forest, total);
  return(forest);
}

getAlphadiversity <- function(url){
  library(jsonlite);
  data <- fromJSON(url);
  forest <- data.frame(table(data$especie), row.names = NULL);
  names(forest) <- c('Species', 'absDen');
  forest$relDen <- (forest$absDen/sum(forest$absDen));
  s <- length(unique(data$especie));
  n <- nrow(data);
  index <- NULL;
  index$Margalef <- (s-1)/log(n);
  index$Menhinick <- s/sqrt(n);
  index$ShannonH <- -sum((forest$relDen * log(forest$relDen)));
  index$ShannonE <- index$ShannonH / log(s);
  index$Simpson <- sum(forest$relDen^2);
  index$BergerParker <- max(forest$absDen) / n;
  return(index);
}
