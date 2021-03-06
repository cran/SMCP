\name{SMCP}
\alias{SMCP}
\alias{sp}
\title{Smoothed minimax concave penalization (SMCP) method for genome-wide association studies.}
\description{
  Fit coefficients paths for a linear model with smoothed minimax concave penalty (SMCP). The model is marginalized so that missing values can be accommendated.
}
\usage{
SMCP(x,y,alpha,lambda,gamma,eps=1e-20,n.iter=100)
sp(x,y,alpha,n.lambda,lambda.min=ifelse(n>p,.001,.05),gamma)
}

\arguments{
  \item{x}{The design matrix which can include missing values.}
  \item{y}{The response variable.}
  \item{alpha}{The proportion of the tuning parameter for MCP and QA parts. alpha*lambda is the tuning parameter for MCP. (1-alpha)*lambda is the tuning parameter for QA part.}
  \item{lambda}{The overall tuning parameter.}
  \item{gamma}{The MCP tuning parameter which affects the magnitude of the shrinkage.}
  \item{eps}{Covergence criterion. The iteration will stop if the relative change is smaller than eps}
  \item{n.iter}{The maximum number of iterations}
  \item{n.lambda}{The number of lambdas to be solved for the solution paths}
  \item{lambda.min}{The minimum lambda to find solution paths}
}

\details{
  The function minimizes 1/(2n)MLS+SMCP, where MLS is the marginalized least squares, SMCP is a smoothed minimax concave penalty.
}
\value{
  An object with S3 class '"SMCP"' containing:
  beta: The fitted values of coefficients. 
  w2 : The weights for QA parts computed as adjacent correlation.

  An object with S3 class '"sp"' containing:
  sp : The solution paths for the lambda from the maximum to the minimum. The number of columns equals the number of coefficients, and the number of rows is equal to 'n.lambda'. The first row is the coefficients for the maximum lambda, and the last row is for the minimum lambda.

}
\references{Liu, J., Wang, K. , Ma, S. and Huang, J. (2010)  Accouting for linkage disequilibrium in genome-wide association studies:A penalized regression method. Technical report #404, Department of Statistics & Actuarial Science, University of Iowa}
\author{Jin Liu <jin-liu@uiowa.edu>}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit=SMCP(x,y,0.3,0.3,3)
#solution path for this data
s.p=sp(x,y,0.3,10,0.05,3)
}
\keyword{ model fitting }
