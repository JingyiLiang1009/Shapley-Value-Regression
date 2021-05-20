\name{ShapleyValue}
\alias{ShapleyValue}
\title{ShapleyValueRegression -- to calculate the relative importance of attributes in linear regression}
\usage{
shapleyvalue(y = the dependent variable, x = independent variables)
}
\description{
Prints a data table with two rows: the unstandardized and standardized relative importance of each attributes using shapley value regression method.
}
\arguments{
y            A coloumn or data set of the dependent variable

x            A matrix or data set of the independent variables
}
\examples{
library(MASS)
library(tidyverse)
data <- Boston
y <- data$medv
x <- as.data.frame(data[,5:8])
shapleyvalue(y,x)
}
