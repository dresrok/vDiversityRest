## Función para calcular el IVI
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
## Función para calcular alfadiversidad
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
## Función para calcular betadiversidad
getBetadiversity <- function(urlEcosystemA, urlEcosystemB){
  library(jsonlite);
  ecosystemA <- fromJSON(urlEcosystemA);
  ecosystemB <- fromJSON(urlEcosystemB);
  speciesA <- unique(ecosystemA$especie);
  speciesB <- unique(ecosystemB$especie);

  AB <- intersect(tolower(speciesA), tolower(speciesB));

  a <- length(unique(ecosystemA$especie));
  b <- length(unique(ecosystemB$especie));
  j <- length(AB);

  index <- NULL;
  index$JaccardCj <- j/(a+b-j);
  index$SorensenCs <- 2*j/(a+b);

  aN <- nrow(ecosystemA);
  bN <- nrow(ecosystemB);

  match <- integer(0);
  ani <- integer(0);
  bni <- integer(0);
  ani_bni <- integer(0);
  Xij_Xik <- integer(0);
  maxS <- integer(0);

  for(i in 1:length(AB)){
    sA <- subset(ecosystemA, tolower(especie) == AB[i]);
    sB <- subset(ecosystemB, tolower(especie) == AB[i]);

    if(nrow(sA) < nrow(sB)){
      match[i] <- nrow(sA);
      maxS[i] <- nrow(sB);
    } else{
      match[i] <- nrow(sB);
      maxS[i] <- nrow(sA);
    }
    ani_bni[i] <- nrow(sA)*nrow(sB);
    Xij_Xik[i] <- ( nrow(sA) - nrow(sB) )^2
  }

  jN <- sum(match);

  index$SorensenCn <- 2*jN/(aN + bN);
  S_ani_bni <- sum(ani_bni);

  freqSpeciesA <- as.data.frame(table(ecosystemA$especie));
  freqSpeciesB <- as.data.frame(table(ecosystemB$especie));

  ani <- sum(freqSpeciesA$Freq^2);
  bni <- sum(freqSpeciesB$Freq^2);

  da <- ani / (aN^2);
  db <- bni / (bN^2);

  index$MorisitaHornCmH <- (2*S_ani_bni) / ((da+db)*(aN*bN));

  tmpXij <- subset(ecosystemA, !( tolower(especie) %in% AB ) );
  Xij <- as.data.frame(table(factor(tmpXij$especie)));
  tmpXik <- subset(ecosystemB, !( tolower(especie) %in% AB ) );
  Xik <- as.data.frame(table(factor(tmpXik$especie)));

  S_Xij_Xik <- sum(Xij$Freq^2) + sum(Xik$Freq^2) + sum(Xij_Xik);

  index$DJK <- sqrt(S_Xij_Xik);
  n <- (a+b)-j;
  index$djk <- sqrt( (index$DJK)^2 / n );

  PS <- 200 * ( jN / (nrow(ecosystemA) + nrow(ecosystemB)) );

  index$PD <- 100 - PS;

  RI <- 100 * ( jN / (sum(Xij$Freq) + sum(Xik$Freq) + sum(maxS)) );

  index$PR <- 100 - RI;

  return(index)
}
