read_docs <- function (filename="mult.dat") 
{
  # Reads documents in LDA-C format. This 
  # function considers that vocab-id starts at 0
  # 
  # Arguments: 
  #   filename - file in the LDAC format  
  # 
  # Reference: 
  #    'lda' R package 
  # 
  
  one <- scan(filename, what = "", sep = "\n")
  two <- chartr(":", " ", one)
  three <- strsplit(two, " ", fixed = TRUE)
  
  docs <- lapply(three, function(x) matrix(as.integer(x[-1]), nrow = 2));
  
  docs 
}



vectorize_docs <- function (docs) 
{
  # Convert each document to to Gibbs sampling 
  # format (e.g. David Newman LDA implementation). 
  # 
  # Arguments: 
  #   docs - a list of documents  
  # 
  
  D <- length(docs)
  did <- c();
  wid <- c();
  
  for (d in 1:D){
    doc <- docs[[d]]; 
    u <- dim(doc)[2]; # the number of unique words in document d 
    if (u > 0){
      doc.n <- sum(doc[2,]);
      did <- rbind(did, array(1, dim=c(doc.n, 1)) * (d-1)); # document instances
      
      for (i in 1:u){
        wid <- rbind( wid, array(1, dim=c(doc[2,i], 1)) * doc[1,i] ); 
      }
    }
  }
  
  # We assume that both vocab-id and doc-id starts at 0
  list(did=did, wid=wid);
}


calc_doc_lengths <- function(docs)
{
  # Computes the number of words in each document
  # 
  # Arguments: 
  #   docs - a list of documents  
  # 
  
  D <- length(docs);
  doc.N <- array(0, dim=c(D, 1));
  
  for (d in 1:D){
    doc <- docs[[d]];
    u <- dim(doc)[2];
    if (u > 0) { 
      doc.N[d] <- sum(doc[2,]); 
    }
    
  }
  
  doc.N 
}