\name{kwayr}
\alias{kwayr}
\alias{cellx}
\alias{khmat}
\alias{pasteColsRfit}
\alias{redmod}
\alias{subsets}
\title{
Internal Functions for K-Way analysis of variance
}
\description{
These are internal functions used to construct the robust anova table.
The function \code{raov} is the main program.
}
\usage{
kwayr(levs, data)
cellx(X)
khmat(levsind,permh)
pasteColsRfit(x,sep="")
redmod(xmat,amat)
subsets(k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{levs}{
vector of levels corresponding to each of the factors
}
  \item{data}{
data matrix in the form y, factor 1,..., factor k
}
\item{X}{
n x k matrix where the columns represent the levels of the k factors.
}
\item{levsind}{
Internal parameter.
}
\item{permh}{
Internal parameter.
}
\item{x}{
n x k matrix where the columns represent the levels of the k factors.
}
\item{xmat}{
n x p full model design matrix 
}
\item{amat}{
Internal parameter.
}
\item{k}{
Internal parameter.
}
\item{sep}{
Seperator used in pasteColsRfit 
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Hocking, R. R. (1985), \emph{The Analysis of Linear Models}, Monterey, California: Brooks/Cole.
}
\author{
Joseph McKean, John Kloke
}
\note{
  Renamed pasteCols of library plotrix written by Jim Lemon et. al. June 2011 under GPL 2
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{raov}}
}
