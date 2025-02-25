#'
#' @title Create an interpolated yield object from raw data
#' @description  Create an interpolated yield object from
#'   raw data
#' @name pa_yield
#' @rdname pa_yield
#' @param input an sf object containing the raw yield
#'   monitor data
#' @param data.columns When algorithm is \sQuote{simple},
#'   this argument should be a vector of length 2 or 3
#'   (depends on whether the user wants to adjust for time
#'   lag) indicating which column contains the yield data
#'  , a column containing moisture
#'   information, and a column indicating the time between
#'   readings. When algorithm is \sQuote{ritas}, an optional
#'   named vector with the column names for the variables
#'   \sQuote{mass, flow, moisture, interval, angle, swath,
#'   distance}. If a an unnamed vector is supplied, the
#'   vector is assumed to be in this order. The default is
#'   NULL, in which case the function attempts to guess the
#'   columns by using a dictionary of possible guesses.
#' @param data.units When algorithm is \sQuote{simple},
#'   should be a vector of length two, indicating the units
#'   of the yield column and the moisture column. Common
#'   values would be \sQuote{c('bu/ac', '\%')}. When
#'   algorithm is \sQuote{ritas}, an optional named vector
#'   with strings representing units for the variables
#'   \sQuote{mass, flow, moisture, interval, angle, swath,
#'   distance}. If a an unnamed vector is supplied, the
#'   vector is assumed to be in this order. A typical value
#'   for this argument would be \sQuote{c(flow = 'lb/s',
#'   moisture = '\%', interval = 's', angle = 'degreeN',
#'   width = 'ft', distance = 'ft')}. Please see
#'   \link[units]{valid_udunits} for help with specifying
#'   units. The default is NULL, in which case the function
#'   attempts to guess the units according to the values of
#'   the variable.
#' @param grid an sf or pa_trial object containing the prediction grid.
#'   If the user is processing yield data coming from a
#'   research trial (i.e. follows a trial design), the user
#'   can pass the sf object containing the trial design
#'   information to this argument. If the argument
#'   \sQuote{formula} contains any predictions, the
#'   predictor should be included in the sf object supplied
#'   to this argument. polygons for which the predictions
#'   generated.
#' @param algorithm algorithm used to generate the yield
#'   object.
#' @param formula formula defining the relationship between
#'   the dependent and independent variables. If the
#'   dependent variable is a linear function of the
#'   coordinates, the formula can be \sQuote{z ~ X + Y}. If
#'   the dependent variable is modeled only as a function of
#'   the mean spatial process, the formula can be \sQuote{z
#'   ~ 1}. If no formula is supplied, it defaults to
#'   \sQuote{z ~ 1}.
#' @param overlap.threshold a fraction threshold to remove
#'   observations when there is overlap between the
#'   vehicular polygons. A value of 0 does not remove any
#'   observations. A value of 1 removes all observations
#'   that overlap even minimally with neighboring
#'   observations.
#' @param var.label optional string to name the final
#'   product. Defaults to \sQuote{yield}.
#' @param boundary optional sf object representing the
#'   field's outer boundary. If it not supplied, the
#'   function attempts to generate a boundary from the
#'   observed points.
#' @param clean whether to clean the raw data based on
#'   distance from the field edge and global standard
#'   deviation.
#' @param clean.sd standard deviation above which the
#'   cleaning step will remove data. Defaults to 3.
#' @param clean.edge.distance distance (m) from the field
#'   edge above which the cleaning step will remove data.
#'   Defaults to 0.
#' @param smooth.method the smoothing method to be used. If
#'   \sQuote{none}, no smoothing will be conducted. If
#'   \sQuote{idw}, inverse distance weighted interpolation
#'   will be conducted. If \sQuote{krige}, kriging will be
#'   conducted.
#' @param fun a function used to transform the data.
#'   Currently, the option are \sQuote{none} and
#'   \sQuote{log}. If none, data operations are carried out
#'   in the data scale. If log, the function will
#'   use\link[gstat]{krigeTg} to perform kriging in the log
#'   scale. For now, only relevant when \sQuote{method} is
#'   krige. the log scale and back transform predictions to
#'   the data scale. When TRUE, \sQuote{fomula} should be
#'   \sQuote{z ~ 1}.
#' @param lbs.per.bushel a numeric value representing the
#'   number of pounds in a bushel (e.g., 60 for soybean and
#'   56 for corn). This argument can be ommitted when the
#'   input and output units are in the metric system. It is
#'   necessary otherwise.
#' @param moisture.adj an optional numeric value to set the
#'   moisture value to which the yield map predictions
#'   should be adjusted (e.g., 15.5 for corn, and 13.0 for
#'   soybean). If NULL, the function will adjust the
#'   moisture to the average moisture of the field.
#' @param lag.adj an optional numeric value used to account
#'   for the time lag between the crop being cut by the
#'   combine and the time at which the combine records a
#'   data point.
#' @param unit.system a string representing the unit system
#'   to be used in the function output. If
#'   \sQuote{standard}, the function output will be in
#'   bushel/acre. Alternatively, if \sQuote{metric}, outputs
#'   will be in metric tonnes/hectare.
#' @param remove.crossed.polygons logical, whether to remove
#'   vehicle polygons that crossed different experimental
#'   units of the grid. This is intented to prevent from
#'   diluting the treatment effects. When this argument is
#'   TRUE, the argument \sQuote{grid} must be supplied.
#' @param cores the number of cores used in the operation
#' @param steps EXPERIMENTAL - whether to return the intermediate steps 
#' of the yield processing algorithm
#' @param verbose whether to print function progress.
#'   \sQuote{FALSE or 0} will suppress details. \sQuote{TRUE
#'   or 1} will print a progress bar. \sQuote{>1} will print
#'   step by step messages.
#' @param ... additional arguments to be passed
#'   \link[gstat]{krige} and \link[gstat]{idw}
#' @details This function will follow the steps in the
#'   selected algorithm to produce a yield map from the raw
#'   data.
#' @return an object of class yield
#' @author Caio dos Santos and Fernando Miguez
#' @export
#' @examples
#' \dontrun{
#' extd.dir <- system.file("extdata", package = "pacu")
#' raw.yield <- sf::read_sf(file.path(extd.dir, '2012-basswood.shp'),
#'                          quiet = TRUE)
#' ## the simple algorithm
#' pa_yield(input = raw.yield,
#' algorithm = 'simple',
#' unit.system = 'metric',
#' lbs.per.bushel = 56) ## 56 lb/bushel of maize
#'
#' ## the ritas algorithm
#' pa_yield(input = raw.yield,
#' algorithm = 'ritas',
#' unit.system = 'metric',
#' lbs.per.bushel = 56)
#' }
#'
pa_yield <- function(input,
                     data.columns = NULL,
                     data.units = NULL,
                     grid = NULL,
                     algorithm = c('none', 'simple', 'ritas'),
                     formula = NULL,
                     overlap.threshold = 0.5,
                     var.label = 'yield',
                     boundary = NULL,
                     clean = FALSE,
                     clean.sd = 3,
                     clean.edge.distance = 0,
                     smooth.method = c('none', 'krige', 'idw'),
                     fun = c('none', 'log'),
                     lbs.per.bushel = NULL,
                     moisture.adj = NULL,
                     lag.adj = 0,
                     unit.system = c('none', 'metric', 'standard'),
                     remove.crossed.polygons = FALSE,
                     steps = FALSE,
                     cores = 1L,
                     verbose = TRUE,
                     ...) {
  
  algorithm <- match.arg(algorithm)
  smooth.method <- match.arg(smooth.method)
  unit.system <- match.arg(unit.system)
  fun <- match.arg(fun)
  pb <- ifelse(verbose == 1, TRUE, FALSE)
  verbose <- ifelse(verbose > 1, 1, 0)
  
  
  if (algorithm == 'none')
    stop('Please choose between the simple and ritas algorithms')
  
  if (unit.system == 'none')
    stop('Please choose between the metric and standard unit systems')
  
  if (any(!(data.columns[!is.na(data.columns)] %in% names(input))))
    stop('One or more of the data.columns supplied does not match any columns in the input.')
  
  if ((is.null(lbs.per.bushel) && unit.system == 'standard') ||
      (any(grepl('bu', data.units)) && is.null(lbs.per.bushel))){
    stop('"lbs.per.bushel" argument is needed to convert units to/from US standard unit system.')
  }
  
  if(verbose) cat("Starting... \n")
  
  crt.crs <- sf::st_crs(input)
  if(is.na(crt.crs)) {
    if (verbose) cat("No CRS found. Defaulting to EPSG:4326 \n")
    sf::st_crs(input) <- 'epsg:4326'
  }
  
  input <- pa_2utm(input, verbose) ## This step appears to be quick
  
  if(!is.null(boundary)){
    if (sf::st_crs(boundary) != sf::st_crs(input)) {
      boundary <- sf::st_transform(boundary, sf::st_crs(input))
    }
  }
  
  
  
  if(!is.null(grid)){
    
    if (inherits(grid, 'trial')){
      grid <- grid[['trial']]
    }
    
    if (!inherits(grid, 'sf')) {
      grid <- sf::st_as_sf(grid)
    }
    
    if(is.na(sf::st_crs(grid))) {
      if (verbose) cat("No CRS found for grid. Defaulting to EPSG:4326 \n")
      sf::st_crs(grid) <- 'epsg:4326'
    }
    
    if (sf::st_crs(grid) != sf::st_crs(input)) {
      grid <- sf::st_transform(grid, sf::st_crs(input))
    }
    
    if (!all(c('X', 'Y') %in% names(grid))) {
      grid <- cbind(grid, suppressWarnings(sf::st_coordinates(sf::st_centroid(grid))))
    }
    
  }
  
  
  if (is.null(formula)) {
    form <- formula(z ~ 1)
  }else{
    form <- formula(formula)
  }
  
  
  if(!is.null(formula) && smooth.method != 'krige') {
    stop('formula should only be used when smooth.method = krige')
  }
  
  
  exp.vars <-  all.vars(form)
  exp.vars <- exp.vars[exp.vars != 'z']
  
  if (length(exp.vars) > 0) {
    if (is.null(grid))
      stop('When formula contains explanatory variables, grid must be supplied.')
    
    if (!all(exp.vars %in% names(grid)))
      stop('One or more of the explanatory variables are not present in the grid.')
  }
  
  
  if (!is.null(grid) && !is.null(boundary)) {
    grid <- suppressWarnings(sf::st_intersection(grid, boundary))
  }
  
  if (algorithm == 'simple') {
    
    if(pb) {
      progress.bar <- utils::txtProgressBar(min = 0, max = 5, style = 3, initial = -1)
      on.exit(close(progress.bar))
    }
    
    exp.order <- c('yield', 'moisture', 'interval')
    if (is.null(data.columns)) data.columns <- rep(NA, 3)
    if (is.null(data.units)) data.units <- rep(NA, 3)
    if(length(data.columns) < 3) length(data.columns) <- 3
    if(length(data.units) < 3) length(data.units) <- 3
    if(is.null(names(data.columns))) names(data.columns) <- exp.order
    if(is.null(names(data.units))) names(data.units) <- exp.order
    data.units <- data.units[exp.order]
    data.columns <- data.columns[exp.order]
    
    if(lag.adj > 0) {
      interval <- .pa_get_variable(input, 'interval', data.units['interval'], data.columns['interval'], verbose)
      
      if (is.null(interval)) {
        time <- try(input[[.pa_get_variable_columns(input, 'time', verbose)]], silent = TRUE)
        if(inherits(time, 'try-error'))
          stop('interval when "lag.adj" > 0. Could not find the interval column or derive it from a timestamp')
        interval <- .pa_time2interval(time)
        interval <- .pa_enforce_units(interval, 'time')
        input$interval <- interval
      }
      
      input <- .pa_adjust_lag(input, lag.adj, 'interval')
    }
    
    
    moisture <- .pa_get_variable(input, 'moisture', data.units['moisture'], data.columns['moisture'], verbose)
    yield <- .pa_get_variable(input, 'yield', data.units['yield'], data.columns['yield'], verbose)
    
    not.found <- sapply(list(yield, moisture), is.null)
    if(any(not.found)) {
      not.found.i <- which(not.found == TRUE)
      stop('unable to find column(s): ', paste(exp.order[not.found.i], collapse = ', '))
    }
    
    if(pb)
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    
    if (units(yield)[[1]] == 'bushel') {
      data.units[1] <- 'bushel/acre'
      attributes(yield) <- NULL
    }
    
    if(pb)
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    
    tgt <- data.frame(yield = yield)
    tgt <- cbind(tgt, sf::st_geometry(input))
    tgt <- st_as_sf(tgt)
    
    attributes(moisture) <- NULL
    if(is.null(moisture.adj)) {moisture.adj <- mean(moisture, na.rm = TRUE)}
    tgt[[1]] <- .pa_moisture(tgt[[1]], crt = moisture, to = 0, verbose = FALSE)
    
    
    initial.geometries <- sf::st_geometry(tgt)
    if(clean){
      if (is.null(boundary)){
        if(verbose)
          cat('The argument boundary is needed to clean the data from the field edges. Estimating it internally.')
        
        boundary <- .pa_field_boundary(sf::st_geometry(input))
      }
      tgt <- .pa_clean_yield_monitor(tgt, 'yield', boundary, clean.edge.distance, clean.sd, verbose)
      
    }
    
    final.geometries <- sf::st_geometry(tgt)
    
    if (!is.null(grid)) {
      if (!is.null(boundary)) {
        grid <- suppressWarnings(sf::st_intersection(grid, boundary))
      }
      f.grid <- stats::aggregate(tgt, grid, FUN = function(x) mean(x, na.rm = TRUE))
      f.grid <- stats::na.omit(f.grid)
      tgt <- f.grid
    }
    
    app.pols <- tgt
    names(app.pols) <- c('mass', 'geometry')
    
    if (tolower(data.units[1]) %in% c('bu/ac', 'bushel/acre')) {
      
      if(is.null(lbs.per.bushel)) {
        stop('When the units of the yield data are in the U.S. standard system of units, lbs.per.bushel must be supplied.')
      }
      app.pols$mass <- .pa_bushel_metric(app.pols$mass, lbs.per.bushel, 'metric', 1)
    }else{
      units(app.pols$mass) <- units::as_units('g/m2')
      attributes(app.pols$mass) <- NULL
    }
    
    if(pb)
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    
  sbs <- list(initial.geometries, final.geometries)
  steps.names <- c('initial.geometries', 
                   'final.geometries',
                   'grid')
  }
  
  
  if (algorithm == 'ritas') {
    
    sbs <- list()
    ## handling units and column names
    exp.order <- c('mass', 'flow', 'moisture', 'interval', 'angle', 'width', 'distance')
    if (is.null(data.columns)) data.columns <- rep(NA, 7)
    if (is.null(data.units)) data.units <- rep(NA, 7)
    if(is.null(names(data.columns))) names(data.columns) <- exp.order
    if(is.null(names(data.units))) names(data.units) <- exp.order
    data.units <- data.units[exp.order]
    data.columns <- data.columns[exp.order]
    
    interval <- .pa_get_variable(input, 'interval', data.units['interval'], data.columns['interval'], verbose)
    if (is.null(interval)) {
      time.col <- .pa_get_variable_columns(input, 'time', verbose)
      if (!is.null(time.col)){
        time <- input[[time.col]]
        interval <- .pa_time2interval(time)
        interval <- .pa_enforce_units(interval, 'time')
        input$interval <- interval
      }
    }
    
    if(lag.adj > 0) {
      int.names <- .pa_get_variable_names('interval')
      if(is.null(int.names)){
        stop('Interval is needed when adjusting for lag.')
      }
      int.col <-  which(tolower(names(input)) %in% int.names)
      input <- .pa_adjust_lag(input, lag.adj, int.col)
      interval <- .pa_get_variable(input, 'interval', data.units['interval'], data.columns['interval'], verbose)
    }
    
    ## keeping track of the units. this is intend this to prevent mistakes.
    mass <- .pa_get_variable(input, 'mass', data.units['mass'], data.columns['mass'], verbose)
    moisture <- .pa_get_variable(input, 'moisture', data.units['moisture'], data.columns['moisture'], verbose)
    
    if(is.null(mass)){
      flow <- .pa_get_variable(input, 'flow', data.units['flow'], data.columns['flow'], verbose)
      if (is.null(interval) || is.null(flow))
        stop('Either mass or, simutaneously, flow and interval are necessary to continue. Those conditions were not met.')
      mass <- .pa_flow2mass(flow, interval, moisture)
    }else{
      units(mass) <- NULL
      units(moisture) <- NULL
      flow <- 'mass is present'
      mass <- .pa_moisture(mass, moisture, 0, verbose)
    }
    
    angle <- .pa_get_variable(input, 'angle', data.units['angle'], data.columns['angle'], verbose)
    if(is.null(angle)) {
      if(verbose) cat('Trajectory angle not found. estimating it from geographical coordinates.\n')
      angle <- .pa_estimate_angle(sf::st_geometry(input))
    }
    swath <- .pa_get_variable(input, 'width', data.units['width'], data.columns['width'], verbose)
    distance <- .pa_get_variable(input, 'distance', data.units['distance'], data.columns['distance'], verbose)
    
    ## checking that all necessary variables were found
    not.found <- sapply(list(mass, flow, moisture, interval, angle, swath, distance), is.null)
    if(any(not.found)) {
      not.found.i <- which(not.found == TRUE)
      stop('unable to find column(s): ', paste(exp.order[not.found.i], collapse = ', '))
    }
    
    ### might need to change this after testing
    if(pb) {
      progress.bar <- utils::txtProgressBar(min = 0, max = 7, style = 3)
      on.exit(close(progress.bar))
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    }
    
    
    units(moisture) <- NULL
    if(is.null(moisture.adj)) {moisture.adj <- mean(moisture, na.rm = TRUE)}
    
    sbs[[length(sbs) + 1]] <- sf::st_geometry(input)
    
    ## now, we can drop the units because we know which units are
    ## in and out of each operation
    swath <- units::drop_units(swath)
    distance <- units::drop_units(distance)
    angle <- units::drop_units(angle)
    
    if(pb)
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    
    v.pols <- pa_make_vehicle_polygons(sf::st_geometry(input),
                                       swath,
                                       distance,
                                       angle,
                                       cores = cores,
                                       verbose = verbose)
    
    sbs[[length(sbs) + 1]] <- sf::st_geometry(v.pols)
    
    if(pb)
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    
    adj.pols <- pa_adjust_obs_effective_area(v.pols,
                                             mass,
                                             var.label = 'mass',
                                             overlap.threshold = overlap.threshold,
                                             cores = cores,
                                             verbose = verbose)
    
    sbs[[length(sbs) + 1]] <- sf::st_geometry(adj.pols)
    
    if(remove.crossed.polygons) {
      
      if(is.null(grid)){
        stop('when remove.crossed.polygons is true, grid needs to be supplied')
      }
      
      grid <- sf::st_transform(grid, sf::st_crs(input))
      cpi <- sf::st_covered_by(adj.pols, grid)
      cpi <- as.numeric(cpi)
      adj.pols <- adj.pols[!is.na(cpi), ]
      if (verbose) {cat('removing ', sum(is.na(cpi)), 'vehicle polygons that crossed experimental units from the grid \n')}
    }
    

    if(clean){
      if (is.null(boundary)){
        if(verbose)
          cat('The argument boundary is needed to clean the data from the field edges. Estimating it internally.\n')
        boundary <- .pa_field_boundary(sf::st_geometry(input))
      }
      adj.pols <- .pa_clean_yield_monitor(adj.pols, 'adj.mass', boundary, clean.edge.distance, clean.sd, verbose)
    }
    
    sbs[[length(sbs) + 1]] <- sf::st_geometry(adj.pols)
    
    if(pb)
      utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
    
    app.pols <- pa_apportion_mass(polygons =  sf::st_geometry(adj.pols),
                                  mass.vector = adj.pols$adj.mass,
                                  cores = cores,
                                  verbose = verbose)
    
    sbs[[length(sbs) + 1]] <- sf::st_geometry(app.pols)
    
    
  steps.names <- c('initial.points', 'harvest.polygons', 'adjusted.polygons', 
                  'cleaned.polygons', 'apportioned.polygons', 'grid')
  }
  
  ## the following steps are the same regardless of the algorithm
  
  if (!is.null(grid)){
    app.pols <- suppressWarnings(sf::st_join(app.pols, grid, join = sf::st_intersects, left = TRUE, largest = TRUE))
  }
  
  
  
  if(pb)
    utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
  
  if (smooth.method == 'krige'){
    app.pols$mass <- .pa_moisture(app.pols$mass, 0, moisture.adj, verbose)
    app.pols$mass[app.pols$mass == 0] <- 1e-3
    app.pols$z <- app.pols$mass
    app.pols <- cbind(app.pols, suppressWarnings(sf::st_coordinates(sf::st_centroid(app.pols))))
    app.pols <- stats::na.omit(app.pols)
    preds <- .pa_predict(formula = form,
                         smooth.method = smooth.method,
                         df = app.pols,
                         new.df = grid,
                         cores = cores,
                         fun = fun,
                         verbose = verbose,
                         ...)
    
    variogram.model <- preds[[2]]
    variogram <- preds[[3]]
    preds <- preds[[1]]
    predicted.var <- data.frame(preds$var1.pred, preds$var1.var)
    names(predicted.var) <- c(var.label, paste0(var.label,'.var'))
    preds <- cbind(preds, predicted.var)
    preds <- preds[c(var.label,  paste0(var.label,'.var'), 'geometry')]
    sf::st_geometry(preds) <- 'geometry'
    
    if (!is.null(grid)){
      if(length(exp.vars) > 0)
        preds <- cbind(preds, as.data.frame(grid)[exp.vars])
    }
    
    preds[[1]] <- .pa_unit_system(preds[[1]], unit.system, lbs.per.bushel)
    preds[[2]] <- .pa_unit_system(preds[[2]], unit.system, lbs.per.bushel, 2)
    attr(preds[[2]], 'units') <- paste(unlist(units(preds[[2]])[1:2]), collapse = '/')
    attributes(preds[[2]]) <- NULL
  }
  
  if (smooth.method == 'idw'){
    
    if (form != formula(z ~ 1)){
      stop('The IDW smoothing method does not allow for predictors in the formula. The "formula" argument should be: z ~ 1')
    }
    
    app.pols$z <- app.pols$mass
    app.pols <- subset(app.pols, !is.na(mass))
    preds <- .pa_predict(formula = form,
                         smooth.method = smooth.method,
                         df = app.pols,
                         new.df = grid,
                         cores = cores,
                         verbose = verbose,
                         ...)
    
    variogram.model <- NULL
    variogram <- NULL
    preds <- preds[[1]]
    
    
    predicted.var <- data.frame(preds$var1.pred)
    names(predicted.var) <- var.label
    preds <- cbind(preds, predicted.var)
    preds <- preds[c(var.label, 'geometry')]
    sf::st_geometry(preds) <- 'geometry'
    preds[[1]] <- .pa_moisture(preds[[1]], 0, moisture.adj, verbose)
    preds[[1]] <- .pa_unit_system(preds[[1]], unit.system, lbs.per.bushel)
  }
  if (smooth.method == 'none'){
    preds <- app.pols['mass']
    preds <- stats::na.omit(preds)
    preds[['mass']] <- .pa_moisture(preds[['mass']], 0, moisture.adj, verbose)
    preds[['mass']] <- .pa_unit_system(preds[['mass']], unit.system, lbs.per.bushel)
    names(preds) <- c(var.label, 'geometry')
    sf::st_geometry(preds) <- 'geometry'
    preds <- preds[c(var.label, 'geometry')]
    variogram.model <- NULL
    variogram <- NULL
  }
  
  if(pb)
    utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
  
  attr(preds, 'moisture') <- moisture.adj
  attr(preds, 'units') <- paste(unlist(units(preds[[var.label]])[1:2]), collapse = '/')
  attr(preds, 'algorithm') <- algorithm
  attr(preds, 'smooth.method') <- smooth.method
  attr(preds, 'formula') <- form
  attr(preds, 'resp') <- var.label
  if(!is.null(lbs.per.bushel))
    attr(preds, 'lbs.per.bushel') <- lbs.per.bushel
  
  
  attributes(preds[[var.label]]) <- NULL
  
  sbs[[length(sbs) + 1]] <- sf::st_geometry(preds)
  names(sbs) <- steps.names
  
  res <- list(yield = preds,
              variogram = variogram,
              variogram.model = variogram.model,
              steps = NULL)
  
  
  if (steps){
  res[[length(res)]] <- sbs
  }
  
  if(pb)
    utils::setTxtProgressBar(progress.bar, utils::getTxtProgressBar(progress.bar) + 1)
  
  if (verbose)
    cat('Processing complete!\n')
  
  class(res) <- c('yield', class(res))
  return(res)
}








