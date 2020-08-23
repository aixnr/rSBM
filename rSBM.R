# Perform QC on hemagglutinin sequence.
# Original code by Dr. Chris Anderson.
# Refactored by A.K.E. (Aizan), Topham Lab 2019-2020
# This function accepts protein FASTA file, properly delimited.
# Requires seqinr package.

HemagglutininQC <- function(fasta, AA_size) {
  # Read .fasta file
  sequence <- seqinr::read.alignment(file = fasta, format = 'fasta')
  fastafile <- seqinr::read.fasta(file = fasta, seqtype = 'AA', as.string = T, set.attributes = F)
  
  # Local variable for expected length, use with grep()
  #length <- c("565|566")
  length <- c(AA_size)
  
  # Catch for matched length
  PullIndex <- grep(length, sapply(sequence$seq, nchar))
  PulledNames <- names(fastafile)[PullIndex]
  
  # Pull sequences associated with PulledNames, using Catch (logical array) as intermediate
  Catch <- names(fastafile) %in% PulledNames
  QCSequence <- fastafile[ c(which(Catch)) ]
  
  # Find for missing amino acids
  #   start by creating an empty list
  namList <- list()
  
  # Use for loop to check each and every strain
  for (i in 1:length(QCSequence)) {
    myName <- names(QCSequence[i])    # access name of the strain
    aSeq <- QCSequence[[i]]           # access sequence of the name
    
    # Split elements of character vector into substrings with strsplit(x, "")
    #  then locate the occurrence of X with which()
    hasX <- which(strsplit(aSeq, "")[[1]] == "X")
    
    # When found, dump strain name into the pre-made list
    if (length(hasX) > 0) {
      namList <- c(namList, myName)
    }
  }
  
  # Exclude sequence with reference to namList, call this new object as QCSequence2
  QCSequence2 <- QCSequence[ c(which(!(names(QCSequence) %in% namList))) ]
  
  
  # finally, return the QC'ed sequence
  return(QCSequence2)
}

#####################################################################

# Create matrix and calculate Hamming distance
# Original code by Dr. Chris Anderson.
# Refactored by A.K.E. (Aizan), Topham Lab 2019-2020
# ALIGN YOUR FASTA FIRST!

matrix_generator <- function(aligned_fasta) {
  aligned <- seqinr::read.alignment(file = aligned_fasta, format = "fasta")
  
  # Copy the dataframe
  virAln <- aligned
  
  # Loop
  for (i in 1:length(virAln$nam)) {
    name <- virAln$nam[[i]]
    seq <- virAln$seq[[i]]
    
    # Create key:value dictionary in memory for each seq identifier and sequence
    assign(name, seq) 
    
    # The epitope
    a <- substr(seq, 1,579)
    A <- paste0(a)
  
    if (i == 1) {
      namList <- name
    }
    else {
      namList <- c(namList, name)
    }
  }
  
  # Initialize an empty matrix table
  table_matrix <- matrix(data = NA, nrow = length(namList), ncol = length(namList))
  
  # Name the columns and rows
  colnames(table_matrix) <- namList
  rownames(table_matrix) <- namList
  
  # Perform Hamming distance calculation
  for (i in 1:length(namList)) {
    name <- paste0(namList[i], "_vec")
    v1 <- get(namList[i])
    for (j in 1:length(namList)) {
      v2 <- get(namList[j])
      dis <- stringdist(v1, v2, method = "h")
      if (j == 1) {vec <- dis} else {vec <- c(vec, dis)}
      assign(name, vec)
    }
    if (j == length(namList)) {
      table_matrix[, i] <- vec
    }
  }
  
  return(table_matrix)
}

#####################################################################

# Create matrix and calculate Hamming distance
# Original code by Dr. Chris Anderson.
# Refactored by A.K.E. (Aizan), Topham Lab 2019-2020

plot_distance <- function(table_matrix) {
  # Perform pricipal compenent analysis on distance matrix
  res.pca <- PCA(table_matrix)
  
  # Pull out numbers from PCA()-generated object (3 out 5 dimensions pulled, but why?)
  PC1 <- res.pca$ind$coord[,1]
  PC2 <- res.pca$ind$coord[,2]
  PC3 <- res.pca$ind$coord[,3]
  
  # Rest of the code (copied verbatim, commented code produced error)
  labs <- rownames(res.pca$ind$coord)
  PCs <- data.frame(cbind(PC1,PC2,PC3))
  rnPC <- rownames(PCs)
  #PC1var <- round(res.pca$eig[3][1, 1], digits = 0)
  PC1var <- round(res.pca$eig[3][1], digits = 0)
  myxLab <- paste0("PC1","(",PC1var,"%",")")
  # PC2var<- round((res.pca$eig[3][2, 1]) - (res.pca$eig[3][1, 1]), digits=0)
  PC2var <- round((res.pca$eig[3][2]) - (res.pca$eig[3][1]), digits=0)
  myyLab <- paste0("PC2","(",PC2var,"%",")")
  rownames(PCs) <- labs 
  PCs$label = rownames(PCs)
  plot(PCs$PC1, PCs$PC2, pch = 19, col = "black", ylim=c(min(PCs$PC2),max(PCs$PC2)), xlim=c(min(PCs$PC1),max(PCs$PC1)), 
       ylab=paste("PC2", paste0("(",PC2var,"%",")")), xlab=paste("PC1", paste0("(",PC1var,"%",")")))
  pointLabel(PCs$PC1, PCs$PC2, PCs$label, cex = 0.7, col="grey48")
}
