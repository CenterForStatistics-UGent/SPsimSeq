#Test matrix for positive definiteness
isPosDefMat <- function (x, tol = 1e-08) 
{
  if (!is.square.matrix2(x)) 
    stop("argument x is not a square matrix")
  if (!is.symmetric.matrix2(x)) 
    stop("argument x is not a symmetric matrix")
  if (!is.numeric(x)) 
    stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues <= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

is.symmetric.matrix2 <- function( x )
{
  ###
  ### this function determines if the matrix is symmetric
  ###
  ### argument
  ### x = a numeric matrix object
  ###
  if ( !is.matrix( x ) ) {
    stop( "argument x is not a matrix" )
  }
  if ( !is.numeric( x ) ) {
    stop( "argument x is not a numeric matrix" )
  }    
  if ( !is.square.matrix2( x ) )
    stop( "argument x is not a square numeric matrix" )
  return( sum( x == t(x) ) == ( nrow(x) ^ 2 ) )
} 
is.square.matrix2 <- function( x )
{
  ###
  ### determines if the given matrix is a square matrix
  ###
  ### arguments
  ### x = a matrix object
  ###
  if ( !is.matrix( x ) )
    stop( "argument x is not a matrix" )
  return( nrow(x) == ncol(x) )
}
  
# D <- matrix( c( -2, 1, -2, 1, -2, 1, -2, 1, -2 ), nrow=3, byrow=TRUE )
# isPosDefMat( D )
