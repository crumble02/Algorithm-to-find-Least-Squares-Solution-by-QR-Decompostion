## QR DECOMPOSITION
## Roll: bs2016, Debabrata Sarkar

#We define a function to normalize a vector

unit = function(v){
  
  #the formula we use: unit(v) = v/||v||, where ||v|| is the norm of v
  #since unit vector is not defined for v=0, we exclude this case
  #we use the fact that ||v||=0 iff v=0
  if( sum(v*v)!=0 ){
    return(v/sqrt(sum(v*v)));
  }
  else{
    return(0);
  }
}

#We define a function to multiply a vector and the householder matrix

hmult = function(u,x){
  
  #we directly calculate the product of a vector x with a householder matrix (I -2u*u')
  #it is given by: (I - 2u*u')*x = x - 2u(u'*x)
  return(x - 2*sum(u*x)*u);
}

#We define a shaver function

shaver = function(x){
  
  #we shall shave a vector, or subtract its norm from its first component
  x[1] = x[1] - sqrt(sum(x*x));
  return(unit(x));
}

# We define a function to resize vectors to make them of equal length by adding 0's

size = function(x, len){
  
  z = rep(0,len);
  l = length(x);
  if(l==len){
    return(x);
  }
  else{
    for(i in 1:l){
      z[len - l + i] = x[i];
    }
    return(z)
  }
}
  
#We define a function to solve for the least squares solution of the system Ax=b

least.sq = function(A, b){
  
  R = A;
  m = nrow(A);
  n = ncol(A);
  
  
  #the motive behind creating this list would be discussed later
  vec = list();
  
  
  #if length of b is not equal to length of Ax (i.e number of columns of A), the system is inconsistent
  #so we address that case in the following code block
  stopifnot( length(b) == m )
  
  for(i in 1:min(m,n) ){
    
    #starting from the i-th column of A from its i-th element onwards
    #we subtract its norm from the first component (shave)
    k = A[,i][i:m];
    
    s = shaver(k);
    
    #the shaving process alters the length of the vector s
    #to make length of vectors = number of rows of A, without altering their norm
    #we add 0's to fill up the remaining components 
    
    k = size(s, m);
    
    vec[[i]] = k;
    
    
    for(j in 1:n){
      
      A[,j] <- hmult(k, A[,j]);
    }
  }
  
  # due to errors in calculation
  # the S step may not make the swept column 0
  # in such case, the residual left would be of very small order
  # so we round that off to 0
  for(i in 1 : m){
    for(j in 1 : n){
      
      if(abs(A[i, j]) <= 1e-13)
        A[i, j] = 0;
    }
  }
  
  print("The upper triangular matrix R: ");
  print(A);
  
  
  # We must stop if the matrix is not of full column rank
  #in A=QR, Q has full rank and
  #rank of R (upper triangular) = number of non-zero rows
  
  rank = 0
  for(i in 1:m){
    
    if( sum(A[i,]*A[i,]) ){
      
      #for each non-zero row, rank increases by 1
      #at the end, 'rank' stores the value of the rank of A
      rank <- rank+1;
    }
     
  }
  if(rank != n){
    
    stop("Matrix A is not of full column rank.");
  }
  else{
    for(i in 1:min(m,n)){
      
      u = vec[[i]];
      b = hmult(u,b);
    }
  }
    
    # We solve Ax = b
    
    x = rep(0,n);

  
  #we made sure that A has full column rank
  #we now partition the upper triangular matrix R into R1 and R2
  #R1 contains the non-zero rows and R2 = O
  R1 = matrix(0 , rank, n);
  for(i in 1:rank){
    
    for(j in 1:n){
      
      R1[i,j] = A[i,j];
    }
  }
  
  #We now solve the system by back substitution
  for(i in seq(n, 1, by = -1)){
    
    if(i != n){
      for(j in (i+1):n){
        
        b[i] <- b[i] - x[j]*R1[i,j];
      }
    }
    x[i] <- b[i]/R1[i,i];
  }
  
  print("The least squares solutions of the system Ax=b: ");
  print(x);
}