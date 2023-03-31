root <- "D:/OneDrive - CGIAR/Data/HERplus"
options(timeout=36000)
getSalinity <- function(url, path){
  temp <- tempfile()
  download.file(url,temp ,mode = "wb")
  unzip(temp, exdir = path)
}
#https://figshare.com/articles/dataset/Predicting_Long-term_Dynamics_of_Soil_Salinity_and_Sodicity_on_a_Global_Scale/13295918/1?file=25617230
#https://www.pnas.org/doi/10.1073/pnas.2013771117
if(!dir.exists(paste0(root, "/Salinity/"))){
  getSalinity("https://figshare.com/ndownloader/files/25617323", root)
  getSalinity("https://figshare.com/ndownloader/files/25616990", root)
  getSalinity("https://figshare.com/ndownloader/files/25617050", root)
  getSalinity("https://figshare.com/ndownloader/files/25617098", root)
}
library(readstata13)
hh <- read.dta13(paste0(root, '/TMRI baseline GPS_upazila identifiers.dta'), convert.factors=T) 
names(hh)
library(terra)
shp <- vect(hh, geom=c("longitude", "latitude"))
library(geodata)
bgd4 <- gadm(country="BGD", level=4, path=root)
bgd0 <- gadm(country="BGD", level=0, path=root)
crs(shp) <- crs(bgd0)
x11()
plot(shp, col='red', cex=1.1)
dim(shp)
plot(bgd0, add=TRUE)
#remove coordinates falling outside the country
shp <- crop(shp, bgd0)
dim(shp)
plot(bgd0, cex=1.1)
plot(shp, col='red', add=TRUE)

#' Extract spatial salinity metrics by aggregating salinity in a given year in each admin 4 boundary.


filename <- paste0(root,"/Salinity/Bangladesh_spatial_temporal_salinity_chirps.csv")

files <- list.files(paste0(root,"/Salinity/"),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
knitr::kable(head(files), caption = 'Salinity raster files')
df <- data.frame(NAME_4=bgd4$NAME_4)
for(i in 1:length(files)){
  #temp <- extract(rast(files[i]), shp, ID=F) #point level
  temp <- extract(rast(files[i]), bgd4, fun=mean, na.rm=T)#Level 4 admin boundaries
  df <- cbind(df, temp)
}
colnames(df)[-1] <- gsub("ECe_","",names(df)[-1])
#write.csv(df, filename)  
                   
#' CHIRPS DOWNLOAD


library(chirps)
dates <- c("2013-5-16","2013-5-20")
#lonlat <- data.frame(lon = c(shp@bbox[1,1],shp@bbox[1,2]), #lat = c(shp@bbox[2,1], shp@bbox[2,2]))
data <- get_chirps(shp, dates, server = "CHC", as.raster = TRUE)
data



### Spatially aggregated rainfall

dff <- extract(data[["chirps-v2.0.2013.05.16"]], shp, bind=T)
dff <- merge(df, as.data.frame(dff), by="a01")
write.csv(dff, filename)

# library(geodata)
# ken <- gadm(country="KEN", level=0, path=root)
# kenjson <- geojsonio::geojson_json(sf::st_as_sf(ken)) 
# geojsonio::geojson_write(kenjson, file = paste0(root, "/Kenya.geojson") )
#aa=geojsonR::FROM_GeoJson("D:/OneDrive - CGIAR/Data/HERplus/Kenya.geojson")

