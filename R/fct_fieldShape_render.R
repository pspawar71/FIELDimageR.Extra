#' fieldShape_render 
#' 
#' @title Building the plot \code{fieldshape} file using Zoom Visualization
#' 
#' @description The user should select the four experimental field corners and the shape file with plots will be automatically 
#' built using a grid with the number of ranges and rows. Attention: The clicking sequence must be (1) upper left, (2) upper right, (3) lower right, and (4) lower left.
#' 
#' @param mosaic object of class 'rast' or 'stars' obtained from function \code{\link{terra::rast()}}.
#' @param ncols number of columns.
#' @param nrows number of rows.
#' @param fieldData data frame with plot ID and all attributes of each plot (Traits as columns and genotypes as rows).
#' @param fieldMap matrix with plots ID identified by rows and ranges, please use first the funsction \code{\link{fieldMap}}.
#' @param plotID name of plot ID in the fieldData file to combine with fieldShape.
#' @param buffer negative values should be used to remove boundaries from neighbor plot.
#' @param plot_size specific plot shape size c(x,y). For example, mosaic with pixels in meter \code{\link{plot_size=c(0.5,1.5)}} means 0.5m x 1.5cm the plot shape.   
#' @param r red layer in the mosaic (for RGB image normally is 1). If NULL the first single layer will be plotted. 
#' @param g green layer in the mosaic (for RGB image normally is 2). If NULL the first single layer will be plotted.
#' @param b blue layer in the mosaic (for RGB image normally is 3). If NULL the first single layer will be plotted.
#' @param color_options single layer coloring options. Check more information at \code{\link{color_options}} from \code{\link{leafem}} package.
#' @param max_pixels maximun pixels allowed before down sampling. Reducing size to accelerate analysis. Default = 100000000.
#' @param downsample  numeric downsample reduction factor. Default = 5.
#'  
#' @importFrom sf st_crs st_bbox st_transform st_is_longlat st_crop st_make_grid st_cast st_coordinates st_buffer st_sf st_centroid
#' @importFrom terra crop nlyr rast
#' @importFrom stars write_stars st_warp st_as_stars 
#' @importFrom mapview mapview
#' @importFrom mapedit editFeatures editMap
#' @importFrom leafem addGeoRaster
#' @importFrom dplyr mutate
#'
#' @return A field plot shape file class "sf" & "data.frame".
#'
#' @export
fieldShape_render <- function(mosaic,
                              ncols, 
                              nrows, 
                              fieldData=NULL,
                              fieldMap=NULL,
                              PlotID=NULL,
                              buffer=NULL,
                              plot_size=NULL,
                              r=1,
                              g=2,
                              b=3,
                              color_options=viridisLite::viridis,
                              max_pixels=100000000,
                              downsample=5
) {
  print("Starting analysis ...")
  if (is.null(mosaic)) {
    stop("The input 'mosaic' object is NULL.")
  }
  
  if (class(mosaic) %in% c("RasterStack", "RasterLayer", "RasterBrick")) {
    mosaic <- terra::rast(mosaic)
  }
  pixels <- prod(dim(mosaic))
  if (pixels > max_pixels) {
    print("Your 'mosaic' is too large and downsampling is being applied.")
  }
  if (pixels < max_pixels) {
    stars_object <- mosaic
    if (!inherits(stars_object, "stars")) {
      stars_object <- st_as_stars(mosaic)
      if (!st_is_longlat(stars_object) && nlyr(mosaic) > 2) {
        stars_object <- st_warp(stars_object, crs = 4326)
      }
    }
  } else {
    stars_object <- mosaic
    if (!inherits(stars_object, "stars")) {
      stars_object <- st_as_stars(mosaic, proxy = TRUE)
    }
  }
  
  print("Use 'Draw Marker' to select 4 points at the corners of the field and press 'DONE'. Attention is very important; start clicking from left to the right and top to bottom.")
  if (nlyr(mosaic) > 2 && pixels < max_pixels) {
    stars_object[is.na(stars_object)] <- 0
    four_point <- mapview() %>%
      leafem:::addRGB(
        x = stars_object, r = r, g = g, b = b,
        fieldData = path_csv_file
      ) %>%
      editMap("mosaic", editor = "leafpm")
  } else { 
    if (nlyr(mosaic) > 2 && pixels > max_pixels) {
      starsRGB <- read_stars(stars_object[[1]], proxy = TRUE)
      starsRGB <- st_downsample(starsRGB, n = downsample)
      starsRGB[is.na(starsRGB)] <- 0
      four_point <- mapview() %>%
        leafem:::addRGB(
          x = starsRGB, r = r, g = g, b = b,
          fieldData = path_csv_file
        ) %>%
        editMap("mosaic", editor = "leafpm")
    } else {
      if (nlyr(mosaic) == 1 && pixels > max_pixels) {
        stars_object[is.na(stars_object)] <- NA
        four_point <- mapview() %>%
          leafem:::addGeotiff(
            stars_object[[1]], colorOptions = leafem:::colorOptions(palette = color_options, na.color = "transparent"),
            fieldData = path_csv_file
          ) %>%
          editMap("mosaic", editor = "leafpm")  
      } else {
        stars_object[is.na(stars_object)] <- NA
        four_point <- mapview() %>%
          leafem:::addGeoRaster(
            stars_object, colorOptions = leafem:::colorOptions(palette = color_options, na.color = "transparent"),
            fieldData = path_csv_file
          ) %>%
          editMap("mosaic", editor = "leafpm")
      }
    }
  }
  if (length(four_point$finished$geometry) == 4) {
    grids <- st_make_grid(four_point$finished$geometry, n = c(ncols, nrows)) %>% st_transform(st_crs(mosaic))
    point_shp <- st_cast(st_make_grid(four_point$finished$geometry, n = c(1, 1)), "POINT")
    sourcexy <- rev(point_shp[1:4]) %>% st_transform(st_crs(mosaic))
    Targetxy <- four_point$finished$geometry %>% st_transform(st_crs(mosaic))
    controlpoints <- as.data.frame(cbind(st_coordinates(sourcexy), st_coordinates(Targetxy)))
    linMod <- lm(formula = cbind(controlpoints[, 3], controlpoints[, 4]) ~ controlpoints[, 1] + controlpoints[, 2], data = controlpoints)
    parameters <- matrix(linMod$coefficients[2:3, ], ncol = 2)
    intercept <- matrix(linMod$coefficients[1, ], ncol = 2)
    geometry <- grids * parameters + intercept
    grid_shapefile <- st_sf(geometry, crs = st_crs(mosaic)) %>% mutate(ID = seq(1:length(geometry)))
    
    rect_around_point <- function(x, xsize, ysize) {
      bbox <- st_bbox(x)
      bbox <- bbox + c(xsize / 2, ysize / 2, -xsize / 2, -ysize / 2)
      return(st_as_sfc(st_bbox(bbox)))
    }
    
    if (!is.null(plot_size)) {
      if (length(plot_size) == 1) {
        cat("\033[1;31mError:\033[0m Please provide x and y distance. e.g., plot_size=c(0.5,2.5)\n")
      } else {
        if (st_is_longlat(grid_shapefile)) {
          grid_shapefile <- st_transform(grid_shapefile, crs = 3857)
          cen <- suppressWarnings(st_centroid(grid_shapefile))
          bbox_list <- lapply(st_geometry(cen), st_bbox)
          points_list <- lapply(bbox_list, st_as_sfc)
          rectangles <- lapply(points_list, function(pt) rect_around_point(pt, plot_size[1], plot_size[2]))
          points <- rectangles[[1]]
          for (i in 2:length(rectangles)) {
            points <- c(points, rectangles[[i]])
          }
          st_crs(points) <- st_crs(cen)
          grid <- st_as_sf(points)
          grid<-st_transform(grid, st_crs('EPSG:4326'))
          b<-st_transform(grid_shapefile, crs = 4326)
          ga = st_geometry(grid)
          cga = st_centroid(ga)
          grid_shapefile = (ga-cga) *parameters+cga
          if(!is.null(mosaic_layer)){
            st_crs(grid_shapefile) <- st_crs(mosaic)
            grid_shapefile<-st_as_sf(grid_shapefile)
          }
          if(is.null(mosaic_layer)){
            st_crs(grid_shapefile) <- st_crs(points_layer)
            grid_shapefile<-st_as_sf(grid_shapefile)
          }
        
        }else
        {
          cen <- suppressWarnings(st_centroid(grid_shapefile))
          
          bbox_list <- lapply(st_geometry(cen), st_bbox)
          points_list <- lapply(bbox_list, st_as_sfc)
          
          rectangles <- lapply(points_list, function(pt) rect_around_point(pt, plot_size[1], plot_size[2]))
          
          points <- rectangles[[1]]
          for (i in 2:length(rectangles)) {
            points <- c(points, rectangles[[i]])
          }
          st_crs(points) <- st_crs(cen)
          grid <- st_as_sf(points)
          if(!is.null(mosaic_layer)){
            st_crs(grid) <- st_crs(mosaic)
          }
          if(is.null(mosaic_layer)){
            st_crs(grid) <- st_crs(points_layer)
          }
          b<-st_transform(grid_shapefile, crs = 4326)
          ga = st_geometry(grid)
          cga = st_centroid(ga)
          grid_shapefile = (ga-cga) *parameters+cga
          if(!is.null(mosaic_layer)){
            st_crs(grid_shapefile) <- st_crs(mosaic)
            grid_shapefile<-st_as_sf(grid_shapefile)
          }
          if(is.null(mosaic_layer)){
            st_crs(grid_shapefile) <- st_crs(points_layer)
            grid_shapefile<-st_as_sf(grid_shapefile)
          }
          
        }
      }
    }
    if (!is.null(buffer)) {
      if (st_is_longlat(grid_shapefile)) {
        grid_shapefile <- st_transform(grid_shapefile, crs = 3857)
        grid_shapefile <- st_buffer(grid_shapefile, dist = buffer)
        grid_shapefile <- st_transform(grid_shapefile, st_crs(mosaic))
      } else {
        grid_shapefile <- st_buffer(grid_shapefile, dist = buffer)
        grid_shapefile <- st_transform(grid_shapefile, st_crs(mosaic))
      }
    }
    grid_shapefile$PlotID <- seq(1, dim(grid_shapefile)[1])                          
    print("Almost there ...")
    if (!is.null(fieldMap)) {
      id <- NULL
      # for(i in 1:dim(fieldMap)[1]){
      #   id<-c(id,rev(fieldMap[i,]))
      # }
      for (i in dim(fieldMap)[1]:1) {
        id <- c(id, fieldMap[i,])
      }
      grid_shapefile$PlotID <- as.character(id)
    }
    
    if (!is.null(fieldData)) {
      if (is.null(fieldMap)) {
        cat("\033[31m", "Error: fieldMap is necessary", "\033[0m", "\n")
      }
      fieldData <- as.data.frame(fieldData)
      fieldData$PlotID <- as.character(fieldData[, colnames(fieldData) %in% c(PlotID)])
      plots <- merge(grid_shapefile, fieldData, by = "PlotID")
    } else {
      if (!is.null(grid_shapefile)) {
        # Plot
        plots <- grid_shapefile
      } 
    }
    
    return(plots)
  } else {
    cat("\033[31m", "Error: Select four points only. Points must be set at the corners of the field of interest under the plots space", "\033[0m", "\n")
  }
}
