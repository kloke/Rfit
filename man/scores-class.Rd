\name{scores-class}
\Rdversion{1.1}
\docType{class}
\alias{scores-class}
\alias{scores}
%\alias{getScores}
%\alias{getScoresDeriv}

\title{Class "scores"}
\description{
A score function and it's corresponding derivative is required for rank-based estimation.  
This object puts them together.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("scores", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{phi}:}{Object of class \code{"function"} the score function }
    \item{\code{Dphi}:}{Object of class \code{"function"} the first derivative of the score function }
    \item{\code{param}:}{Object of class \code{"param"} }
  }
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{ John Kloke }


\seealso{
\code{\linkS4class{param}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("scores")
}
\keyword{classes}
