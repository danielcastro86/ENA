\name{ENA}
\alias{ENA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected Network Activity
}
\description{
Estimates the Expected Network Activity following @@ proposal.}
\usage{
ENA(data, factors, activity, symptoms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{mgm network model}
  \item{factors}{Node number of the network nodes to change}
  \item{activity}{Value for each of the network nodes to change}
  |item{symptoms}{Remaining nodes of the network}
}
\details{
the mgm network model to use should be mgmModel$results$nodemodels}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Daniel Castro}
\note{dcastro@umaia.pt}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
