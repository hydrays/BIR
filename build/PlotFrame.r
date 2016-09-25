require(raster)
Ndim = 73
Mdim = 86
r <- raster(nrow=Ndim, ncol=Mdim, xmn=0, xmx=Mdim, ymn=0, ymx=Ndim)
A = matrix(0, nrow=Ndim*Mdim, ncol=1)
B = matrix(0, nrow=Ndim*Mdim, ncol=1)
phi = matrix(0, nrow=Ndim*Mdim, ncol=1)
phix = matrix(0, nrow=Ndim*Mdim, ncol=1)
phiy = matrix(0, nrow=Ndim*Mdim, ncol=1)
AA = matrix(0, nrow=Ndim*Mdim, ncol=1)
BB = matrix(0, nrow=Ndim*Mdim, ncol=1)
phiphi = matrix(0, nrow=Ndim*Mdim, ncol=1)
Nsample <- 100
index <- seq(Ndim*Mdim)
X <- (index-1)%%Mdim + 1
Y <- Ndim - floor((index-1)/Mdim)
for ( i in seq(Nsample)){
    res <- as.matrix(read.table(paste('result', i, '.txt', sep='')))
    A <- res[,1]
    B <- res[,1] + res[,2]
    phi <- acos(cos(res[,3]))
    phix <- phix + B*cos(res[,3])
    phiy <- phiy + B*sin(res[,3])
    AA <- AA + A
    BB <- BB + B
    phiphi <- phiphi + phi
    png(paste('septin', i, '.png', sep=''), width=800, heigh=800)
    plot(0, 0, type='n', xlim=c(0, Mdim+1), ylim=c(0, Ndim+1), asp=1)
    lcoef = 0.25*B
    s <- A / B
    xstart <- X-s*cos(phi)*lcoef
    xend <- X+s*cos(phi)*lcoef
    ystart <- Y-s*sin(phi)*lcoef
    yend <- Y+s*sin(phi)*lcoef
    arrows(xstart, ystart, xend, yend, length=0.0)
    dev.off()
}
AA <- AA/Nsample
BB <- BB/Nsample
phix <- phix/Nsample
phiy <- phiy/Nsample
phiphi <- atan2(phiy, phix)
##phiphi <- phiphi/Nsample


plot(0, 0, type='n', xlim=c(0, Mdim+1), ylim=c(0, Ndim+1), asp=1)
lcoef = 0.5*BB
s <- AA/BB
xstart <- X-s*cos(phiphi)*lcoef
xend <- X+s*cos(phiphi)*lcoef
ystart <- Y-s*sin(phiphi)*lcoef
yend <- Y+s*sin(phiphi)*lcoef
arrows(xstart, ystart, xend, yend, length=0.0)
