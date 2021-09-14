inflect <- function(x, threshold = 1){
  up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  a    <- cbind(x,up,down)
  list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
}

path="C:/Users/Alish Chelackal/Desktop/master thesis/split code"
setwd(path)

# install.packages("BiocManager")
# BiocManager::install("EBImage")
library("EBImage")
library(pdftools)
library(pracma)
library(zoo)
library(ggplot2)
instimg_orig = readImage("rTLC_demopicture_4_01.jpg")
# Rows und Col are interchanged
display(img_orig, method="raster")
# Infos of picture
print(img_orig)
img<-img_orig
# To grayscale
colorMode(img) = Grayscale
# get out the Matrix of greyscales
mat_img <- img[,,1]
# Calculating the mean greyscale 
meanrow<-apply(mat_img,1,mean)

# Mean greyscale of x
plot(meanrow)

# Making the fast fourier Transformation
# Division to Normalize

ffft<-fft(meanrow)/sqrt(length(meanrow));
plot(meanrow)

# Absolutvalues of the Amplitudes
ffftabs<-abs(fft(meanrow))
head(ffftabs)
tail(ffftabs)
# because the values are freal all is double
# Taking one have is enough
# first is the Constant so we omit it

ffta<-ffftabs[2:ceiling(length(ffftabs)/2)]

# Amplitudes vs. k

plot(1:length(ffta),ffta,typ="l")

kmax<-which.max(ffta)
plot((kmax-5):(kmax+5),ffta[(kmax-5):(kmax+5)],typ="l")

# Taking Maximal Am^plitude +- 8 
kde<-8
indicetest<-(kmax-kde):(kmax+kde)
indicetesta<-indicetest+1
indicetestb<-length(ffft)-indicetest+1

# upper and lower part 
# should be the same with minus imaginary part
ffft[indicetesta]
ffft[indicetestb]


# only take a short part of the spektrum

ffftinv<-ffft*0
ffftinv[indicetesta]<-ffft[indicetesta]
ffftinv[indicetestb]<-ffft[indicetestb]

# inverse fourier Transform
fftinv<-fft(ffftinv,inverse = T)/sqrt(length(meanrow));

# function with only a small part of the Spektrum
# the essentiell part

plot(Re(fftinv), type = "l")

n<-2

# finds local Minima with the function inflect
bottoms <- lapply(1:n, function(x) inflect(Re(fftinv), threshold = x)$minima)
cf.1 <- grDevices::colorRampPalette(c("pink","red"))

# Are these the minimas for the inverse fourier ?
plot(Re(fftinv), type = "l")
for(i in 1:n){
  points(bottoms[[i]], Re(fftinv)[bottoms[[i]]], pch = 16, col = cf.1(n)[i], cex = i/1.5)
}

# original Graph
# with the position of the minimas

plot(meanrow, type = "l")
for(i in 1:n){
  points(bottoms[[i]], meanrow[bottoms[[i]]], pch = 16, col = cf.1(n)[i], cex = i/1.5)
}

# one is missing

# both graphs

plot(meanrow, type = "l", ylim = c(-0.1,0.5))
lines(Re(fftinv), type = "l", col="red")
for(i in 1:n){
  points(bottoms[[i]], meanrow[bottoms[[i]]], pch = 16, col = cf.1(n)[i], cex = i/1.5)
}

# which is the mean distance beetween 2 mimina ?
deltamean<-mean(diff(bottoms[[1]]))
bottomsneu<-bottoms[[1]]
# adding the missing point by hand
if (bottomsneu[1]-deltamean>0){
  bottomsneu<-c(bottomsneu[1]-deltamean,bottomsneu)
}
if (tail(bottomsneu,1)+deltamean<length(meanrow)){
  bottomsneu<-c(bottomsneu,tail(bottomsneu,1)+deltamean)
}


plot(meanrow, type = "l")
for(i in 1:length(bottomsneu)){
  points(bottomsneu[i], meanrow[bottomsneu[i]], pch = 16, col = "red", cex = 2/1.5)
}


# write the lanes as png


lanes<-list()
k<-1
for (k in 1:(length(bottomsneu)-1)){
  img_crop = img[bottomsneu[k]:bottomsneu[k+1],,1]  
  lanes[[k]]<-img_crop
  filename<-paste("Lane",k,".png",sep="")
  display(img_crop)
  writeImage(img_crop, filename,"png")
}








