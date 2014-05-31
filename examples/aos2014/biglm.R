## big glm example
require(ffbase)
require(LaF)
require(ETLUtils)

# download.file("http://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/BSAPUFS/Downloads/2010_Carrier_PUF.zip", "2010_Carrier_PUF.zip")
unzip(zipfile="2010_Carrier_PUF.zip")

## the LaF package is great if you are working with fixed-width format files but equally good for csv files and laf_to_ffdf does what it
## has to do: get the data in an ffdf
dat <- laf_open_csv(filename = "2010_BSA_Carrier_PUF.csv",
                    column_types = c("integer", "integer", "categorical", "categorical", "categorical", "integer", "integer", "categorical", "integer", "integer", "integer"),
                    column_names = c("sex", "age", "diagnose", "healthcare.procedure",
                                     "typeofservice", "service.count", "provider.type", "servicesprocessed",
                                     "place.served", "payment", "carrierline.count"), 
                    skip = 1)
x <- laf_to_ffdf(laf = dat)


##
## Data Profiling using table.ff
##
table.ff(x$age)
table.ff(x$sex)
table.ff(x$typeofservice)
barplot(table.ff(x$age), col = "lightblue")
barplot(table.ff(x$sex), col = "lightblue")
barplot(table.ff(x$typeofservice), col = "lightblue")

print(x[1,])

##
## Make a linear model using biglm
##
require(biglm)
print(object.size(x), units="Mb")
mymodel <- bigglm(payment ~ sex + age + place.served, data = x)
print(summary(mymodel))
# This will overflow your RAM as it will get your data from ff into RAM
#summary(glm(payment ~ sex + age + place.served, data = x[,c("payment","sex","age","place.served")]))

##
## Ok, we were working only on +/- 2.8Mio records which is not big, let's explode the data by 100 to get 280Mio records 
##
explode.data <- function(x) {
  x$id <- ffseq_len(nrow(x))
  print("Exploding to 200M records")
  xexploded <- expand.ffgrid(x$id, ff(1:100)) # Had to wait 3 minutes on my computer
  colnames(xexploded) <- c("id","explosion.nr")
  xexploded <- merge(xexploded, x, by.x="id", by.y="id", all.x=TRUE, all.y=FALSE) ## this uses merge.ffdf, might take 30 minutes
  dim(xexploded) ## hopsa, 280 Mio records and 13.5Gb created
  sum(.rambytes[vmode(xexploded)]) * (nrow(xexploded) * 9.31322575 * 10^(-10))
  ## And build the linear model again on the whole dataset
  mymodel <- bigglm(payment ~ sex + age + place.served, data = xexploded)
  summary(mymodel)  
}

