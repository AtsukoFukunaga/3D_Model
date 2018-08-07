library(raster)
library(rgeos)
library(ggplot2)

ras <- raster("GitHub/DEM_processing/FFS_4185_DEM_1cm.tif")

extent <- aggregate(ras, fac = 128, fun = mean, expand = FALSE, na.rm = FALSE)
extent_uniform <- extent > -Inf
extent_polygon <- rasterToPolygons(extent_uniform, dissolve = TRUE)

ras_clip <- mask(crop(ras, extent(extent_polygon)), extent_polygon)

dat64 <- data.frame(fac = c(1, 2, 4, 8, 16, 32, 64), s_area = NA, area = NA)
dat128 <- data.frame(fac = c(1, 2, 4, 8, 16, 32, 64, 128), s_area = NA, area = NA)

dat <- list(dat64, dat128)

for (j in 1:length(dat)) {
  for (i in 1:nrow(dat[[j]])) {
    if (dat[[j]]$fac[i] == 1) {
        reef <- ras_clip
    } else {
        reef <- aggregate(ras_clip, fac = dat[[j]]$fac[i], fun = mean, expand = FALSE, na.rm = FALSE)
    }
    reef_g <- as(reef, "SpatialGridDataFrame")
    dat[[j]]$s_area[i] <- surfaceArea(reef_g)
    dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * 0.01) ^ 2
  }
}

p64 <- ggplot(dat[[1]], aes(x = log(fac * 0.01), y = log(s_area))) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

p128 <- ggplot(dat[[2]], aes(x = log(fac * 0.01), y = log(s_area))) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

p64; p128

d64 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[1]])
slope64 <- coef(d64)[[2]]
fd64 <- 2 - slope64

d128 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[2]])
slope128 <- coef(d128)[[2]]
fd128 <- 2 - slope128

fd64; fd128

rugosity <- dat[[1]]$s_area[1]/dat[[1]]$area[1]
rugosity
