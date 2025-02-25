#' @title Create a plot from a pacu object
#' @name pa_plot
#' @rdname pa_plot
#' @param x object to be plotted
#' @return No return value, called for side effects
#' @export
pa_plot <- function(x, ...){
  UseMethod('pa_plot', x)
}

#' @param ... additional arguments. None used currently.
#' @param plot.type type of plot to be produced Defaults to
#'   trial.
#' @param palette a string representing a color palette from
#'   \link[grDevices]{hcl.pals}. Defaults to \sQuote{Temps}.
#' @param main a main title for the plot
#' @param plot.var the name of the column to be plotted.
#'   Defaults to \sQuote{yield}
#' @param style style applied to the colors
#' @param interactive logical. Whether to produce
#'   interactive plots.
#' @param border.col color of the border for the polygons
#'   plotted in the yield map
#' @param scale a numerical value indicating the
#'   magnification of the graph. A value of 1 produces a
#'   plot using the default magnification. Greater values
#'   will produce zoomed in plots.
#' @param frame logical. Whether to draw the frame around
#'   the plotting area.
#' @param legend.outside logical. Whether to place the legend outside of the graph.
#' @param nbreaks numerical value indicating the number of breaks for the color scale.
#' @param breaks a vector indicating numerical breaks for the color scale.
#' @rdname pa_plot
#' @export
#'
pa_plot.trial <- function(x,
                          ...,
                          plot.type = c('trial'),
                          palette = 'Temps',
                          main = '',
                          plot.var = NULL,
                          interactive = FALSE,
                          border.col = 'black',
                          style =  c("quantile", "pretty", 'equal'),
                          scale = 1,
                          nbreaks = 5,
                          breaks = NULL,
                          frame = TRUE,
                          legend.outside = FALSE){
  
  plot.type <- match.arg(plot.type)
  style <- match.arg(style)
  
  s.wrns <-  get("suppress.warnings", envir = pacu.options)
  s.msgs <-  get("suppress.messages", envir = pacu.options)
  
  if(is.null(plot.var))
    plot.var <- attr(x$trial, 'resp')
  
  if(plot.type == 'trial'){
    
    tmap::tmap_options(overlays = NULL, basemaps = NULL)
    ## controlling the colors
    cols <- function(n) {hcl.colors(n, palette, rev = TRUE)}
    ## setting the tmap mode
    if(interactive){suppressMessages(tmap::tmap_mode("view"))} else {suppressMessages(tmap::tmap_mode('plot'))}
    ## the basic plot
    p <- tmap::tm_shape(x$trial)
    
    
    if (sf::st_geometry_type(x$trial[1, ]) %in% c("POLYGON", 'MULTIPOLYGON')){
      p <- p +
        tmap::tm_borders(col = border.col,
                         lwd = 0.5) +
        tmap::tm_fill(plot.var,
                      palette = cols(nbreaks), title = attr(x$trial, 'units'),
                      style = style,
                      style.args = list(n = nbreaks))
      
    }
    
    if (sf::st_geometry_type(x$trial[1, ]) %in% c("POINT", 'MULTIPOINT')){
      p <- p + tmap::tm_dots(col = plot.var,
                             palette = cols(nbreaks), title = attr(x$trial, 'units'),
                             style = style,
                             style.args = list(n = nbreaks))
      
    }
    ## adjusting the layout
    p <- p + tmap::tm_layout(main.title = main,
                             scale = scale,
                             frame = frame,
                             legend.outside = legend.outside,
                             title.size = 1)
    
    print(p)
  }
}

#' @param ... additional arguments. None used currently.
#' @param plot.type type of plot to be produced Defaults to
#'   yieldmap.
#' @param palette a string representing a color palette from
#'   \link[grDevices]{hcl.pals}. Defaults to \sQuote{Temps}.
#' @param main a main title for the plot
#' @param plot.var the name of the column to be plotted.
#'   Defaults to \sQuote{yield}
#' @param style style applied to the colors
#' @param interactive logical. Whether to produce
#'   interactive plots.
#' @param border.col color of the border for the polygons
#'   plotted in the yield map
#' @param scale a numerical value indicating the
#'   magnification of the graph. A value of 1 produces a
#'   plot using the default magnification. Greater values
#'   will produce zoomed in plots.
#' @param frame logical. Whether to draw the frame around
#'   the plotting area.
#' @param extent a bbox object indicating the geographical area to be plotted
#' @param legend.outside logical. Whether to place the legend outside of the graph.
#' @param nbreaks numerical value indicating the number of breaks for the color scale.
#' @param breaks a vector indicating numerical breaks for the color scale.
#' @param ask whether to ask for user before starting a new page of output. If FALSE, plots are 
#' arranged using \link[patchwork]{wrap_plots}
#' @rdname pa_plot
#' @export
#'
pa_plot.yield <- function(x,
                          ...,
                          plot.type = c('yieldmap', 'variogram', 'steps'),
                          palette = 'Temps',
                          main = '',
                          plot.var = NULL,
                          interactive = FALSE,
                          border.col = 'black',
                          style =  c("quantile", "pretty", 'equal'),
                          scale = 1,
                          nbreaks = 5,
                          breaks = NULL,
                          frame = TRUE,
                          extent = sf::st_bbox(x[['yield']]),
                          legend.outside = FALSE,
                          ask = TRUE){
  
  plot.type <- match.arg(plot.type)
  style <- match.arg(style)
  
  s.wrns <-  get("suppress.warnings", envir = pacu.options)
  s.msgs <-  get("suppress.messages", envir = pacu.options)
  
  if(is.null(plot.var))
    plot.var <- attr(x$yield, 'resp')
  
  if(plot.type == 'yieldmap'){
    
    tmap::tmap_options(overlays = NULL, basemaps = NULL)
    ## controlling the colors
    cols <- function(n) {hcl.colors(n, palette, rev = TRUE)}
    ## setting the tmap mode
    if(interactive){suppressMessages(tmap::tmap_mode("view"))} else {suppressMessages(tmap::tmap_mode('plot'))}
    ## the basic plot
    p <- tmap::tm_shape(x$yield)
    
    
    if (sf::st_geometry_type(x$yield[1, ]) %in% c("POLYGON", 'MULTIPOLYGON')){
      p <- p +
        tmap::tm_borders(col = border.col,
                         lwd = 0.5) +
        tmap::tm_fill(plot.var,
                      palette = cols(nbreaks), title = attr(x$yield, 'units'),
                      style = style,
                      style.args = list(n = nbreaks))
      
    }
    
    if (sf::st_geometry_type(x$yield[1, ]) %in% c("POINT", 'MULTIPOINT')){
      p <- p + tmap::tm_dots(col = plot.var,
                             palette = cols(nbreaks), title = attr(x$yield, 'units'),
                             style = style,
                             style.args = list(n = nbreaks))
      
    }
    ## adjusting the layout
    p <- p + tmap::tm_layout(main.title = main,
                             scale = scale,
                             frame = frame,
                             legend.outside = legend.outside,
                             title.size = 1)
    
    print(p)
  }
  
  if(plot.type == 'variogram'){
    if(is.null(x$variogram)){
      stop('Could not find a variogram to plot.')
    }
    vob <- x$variogram
    vf <- x$variogram.model
    invisible(print(plot(vob, vf)))
  }
  
  
  if(plot.type == 'steps'){
    if(is.null(x$steps))
      stop('No intermediate steps found to plot. Did you set steps = TRUE in pa_yield?')
    
    alg <- attr(x$yield, 'algorithm')
    
    if(alg == 'simple'){
      sts <- x$steps
      
      colmap <- c('#1A85FF', '#DC3220')
      names(colmap) <- c('initial points','removed points')
      
      fillmap <-  c('initial points' = 'black',
                    'removed points' = 'transparent')
      
      point.size = 0.2 * scale
      
      p1 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$initial.geometries,
                         cex = point.size,
                         ggplot2::aes(col = 'initial points',
                                      fill = 'initial points'))+
        ggplot2::labs(col = '', fill = '', title = 'Initial points')
      
      removed.indices <- !(sf::st_geometry(sts$initial.geometries) %in% 
                             sf::st_geometry(sts$final.geometries))
      
      p2 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$initial.geometries[removed.indices, ],
                         linewidth = point.size,
                         ggplot2::aes(col = 'removed points',
                                      fill = 'removed points'))+
        ggplot2::labs(col = '', fill = '', title = 'Removed points')
      
      if (!any(removed.indices)){
        
        if (!s.wrns)
          warning('No points removed. Plot showing removed points will not be produced.')
        p2 <- NULL 
      }
      
      plot_list <- list(p1, p2)
      plot_list <- plot_list[!sapply(plot_list, is.null)] 
      plot_list <- lapply(plot_list, function(x){
        x +    
          ggplot2::scale_color_manual(values = colmap)+
          ggplot2::scale_fill_manual(values = fillmap)+
          ggplot2::coord_sf(xlim = extent[c(1,3)], 
                            ylim = extent[c(2,4)],
                            expand = TRUE)+
          ggplot2::theme_bw()+
          ggplot2::theme(axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank())
      })
      
      if (!ask){
        pp <- patchwork::wrap_plots(plot_list)+
          patchwork::plot_layout(guides = 'collect')
        print(pp)
      }
      
      if (ask){
        grDevices::devAskNewPage(ask = TRUE)
        invisible(lapply(plot_list, print))
        grDevices::devAskNewPage(ask = FALSE)
      }
    }
    
    if (alg == 'ritas'){
      sts <- x$steps
      point.size = 0.2 * scale
      colmap <- c('#E66100', '#D41159', '#1A85FF', '#DC3220', '#40B0A6', 'black')
      names(colmap) <- c('initial points','harvest polygons', 'tessellated polygons',
                         'removed polygons', 'apportioned polygons', 'grid')
      
      fillmap <-  c('initial points' = 'black',
                    'harvest polygons' = 'transparent',
                    'tessellated polygons' = 'transparent',
                    'removed polygons' = 'transparent',
                    'apportioned polygons' = 'transparent',
                    'grid' = 'transparent')
      
      p1 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$initial.points,
                         cex = point.size,
                         ggplot2::aes(col = 'initial points',
                                      fill = 'initial points'))+
        ggplot2::labs(col = '', fill = '', title = 'Initial points')
      
      p2 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$initial.points,
                         cex = point.size,
                         ggplot2::aes(col = 'initial points',
                                      fill = 'initial points'))+
        ggplot2::geom_sf(data = sts$harvest.polygons,
                         cex = point.size,
                         ggplot2::aes(fill = 'harvest polygons',
                                      col = 'harvest polygons'))+
        ggplot2::labs(col = '', fill = '', title = 'Harvest polygons')
      
      
      p3 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$harvest.polygons,
                         linewidth = point.size * 0.5,
                         ggplot2::aes(fill = 'harvest polygons',
                                      col = 'harvest polygons'))+
        ggplot2::geom_sf(data = sts$adjusted.polygons,
                         linewidth = point.size,
                         ggplot2::aes(col = 'tessellated polygons',
                                      fill = 'tessellated polygons'))+
        ggplot2::labs(col = '', fill = '', title = 'Intersected polygons')
      
      
      
      removed.indices <- !(sf::st_geometry(sts$adjusted.polygons) %in% 
                             sf::st_geometry(sts$cleaned.polygons))
      
      p4 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$adjusted.polygons[removed.indices, ],
                         linewidth = point.size,
                         ggplot2::aes(col = 'removed polygons',
                                      fill = 'removed polygons'))+
        ggplot2::labs(col = '', fill = '', title = 'Removed polygons')
      
      if (length(removed.indices) < 1){
        if (!s.wrns)
          warning('No points removed. Plot showing removed points will not be produced.')
        p4 <- NULL 
      }
      
      
      p5 <- ggplot2::ggplot()+
        ggplot2::geom_sf(data = sts$apportioned.polygons,
                         linewidth = point.size,
                         ggplot2::aes(col = 'apportioned polygons',
                                      fill = 'apportioned polygons'))+
        ggplot2::labs(col = '', fill = '', title = 'Apportioned polygons')
      
      
      plot_list <- list(p1, p2, p3, p4, p5)
      plot_list <- plot_list[!sapply(plot_list, is.null)] 
      plot_list <- lapply(plot_list, function(x){
        x +    
          ggplot2::scale_color_manual(values = colmap)+
          ggplot2::scale_fill_manual(values = fillmap)+
          ggplot2::coord_sf(xlim = extent[c(1,3)], 
                            ylim = extent[c(2,4)],
                            expand = TRUE)+
          ggplot2::theme_bw()+
          ggplot2::theme(axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank())
      })
      
      if (!ask){
        pp <- patchwork::wrap_plots(plot_list)+
          patchwork::plot_layout(guides = 'collect')
        print(pp)
      }
      
      if (ask){
        grDevices::devAskNewPage(ask = TRUE)
        invisible(lapply(plot_list, print))
        grDevices::devAskNewPage(ask = FALSE)
      }
    }
    
  }
}





#' @rdname pa_plot
#' @param by a string or vector of strings used to group the data 
#' when plotting. Defaults to 'year' 
#' @param xlab a string used as label for x axis
#' @param ylab a string used as label for y axis
#' @param legend.title a string used as title for the legend
#' @param pch an integer indicating which shape to use for points
#' @export
pa_plot.veg.index <- function(x,
                              ...,
                              palette = ifelse(plot.type == 'timeseries', 'Dark 2','Temps'),
                              plot.type= c('spatial', 'timeseries'),
                              main = '',
                              plot.var = NULL,
                              by = 'year',
                              xlab = NULL,
                              ylab = NULL,
                              style =  c("quantile", "pretty", 'equal'),
                              nbreaks = 5,
                              border.col = 'black',
                              frame = TRUE,
                              legend.outside = FALSE,
                              legend.title = NULL,
                              pch = 16) {

  plot.type <- match.arg(plot.type)
  
  
  
  if (is.null(plot.var))
    plot.var <- attr(x, 'vegetation.index')
  
  if (plot.type == 'spatial'){
  
  initial.options <- tmap::tmap_options()
  tmap::tmap_options(show.warnings = FALSE)
  on.exit(suppressMessages(tmap::tmap_options(initial.options)))

  is.raster <- inherits(try(sf::st_geometry_type(x), silent = TRUE), 'try-error')
  style <- match.arg(style)

  ## controlling the colors
  cols <- function(n) {hcl.colors(n, palette, rev = TRUE)}
  ndates <- stars::st_get_dimension_values(x, 'time')
  
  if(is.raster){
    p <- tmap::tm_shape(x)
    p <- p + tmap::tm_raster(palette = cols(nbreaks),
                             style = style,
                             style.args = list(n = nbreaks))
  }


  if(!is.raster){
    p <- tmap::tm_shape(x[plot.var])
    p <- p +
      tmap::tm_borders(col = border.col,
                       lwd = 0.5) +
      tmap::tm_fill(plot.var,
                    palette = cols(nbreaks),
                  style = style,
                  style.args = list(n = nbreaks))

  }


    p <- p + tmap::tm_layout(legend.outside = T)

  print(suppressWarnings(p))
  }
  
  if (plot.type == 'timeseries'){
    .pa_plot_ts(x = x, 
                plot.var = plot.var,
                by = by,
                time = 'time',
                ylab = ifelse(is.null(ylab), plot.var, ylab),
                xlab = ifelse(is.null(xlab), "Day of the year", xlab),
                legend.outside = legend.outside,
                palette = palette,
                legend.title = legend.title,
                main = main,
                pch = pch)
  }
  
  
}


#' @param saturation numeric. Controls the image saturation. 0 maps to grayscale. 1 maps to the default value.  See \link[tmap]{tm_rgb} for details.
#' @param alpha numeric between 0 and 1. See \link[tmap]{tm_rgb} for details.
#' @param interpolate logical. Whether the raster image should be interpolated. See \link[tmap]{tm_rgb} for details.
#' @rdname pa_plot
#' @export
pa_plot.rgb <- function(x,
                        ...,
                        interactive = FALSE,
                        saturation = 1,
                        alpha = 1,
                        interpolate = FALSE){
  
  if(interactive){suppressMessages(tmap::tmap_mode("view"))} else {suppressMessages(tmap::tmap_mode('plot'))}
  
  sx <- structure(x)
  sm <- sapply(sx, max, na.rm = TRUE)
  
  time.points <- stars::st_get_dimension_values(x, 'time')
  
  
  p <- list()
  for (t in 1:length(time.points)) {
    one.img <- x[, , , t]
    one.img <-  stars::st_redimension(one.img,
                                      new_dims = dim(one.img)[1:2])
    
    sx <- structure(one.img)
    sm <- sapply(sx, max, na.rm = TRUE)
    g <- tmap::tm_shape(one.img)+
      tmap::tm_rgb(max.value = max(sm),
                   saturation = saturation,
                   alpha = alpha,
                   interpolate = interpolate)+
      tmap::tm_layout(main.title = as.character(time.points[t]))
    p[[length(p) + 1]] <-  g
    
    
    p <- tmap::tmap_arrange(p)
    
    
  }
  
  
  print(p)
}



#'
#' @rdname pa_plot
#' @param unit.system unit system to be used: international
#'   (metric) or stanrdard (imperial)
#' @param start day of the year to start computing the
#'   climate normals. Defaults to 1.
#' @param end day of the year to finish computing the
#'   climate normals. Defaults to 365.
#' @param vars which variables to include in the summary
#'   plot
#' @param tgt.year which year to focus and compare to the
#'   historical mean. Defaults to the last year in the data
#'   set.
#' @param months a numerical vector indicating which months
#'   to produce a plot for in the case of monthly
#'   distribution plots. Defaults to 1:12.
#' @author Caio dos Santos and Fernando Miguez
#' @export
#' @examples
#' \dontrun{
#' ## for examples, please see the pacu vignette
#' }
#'
pa_plot.met <- function(x,
                        ...,
                        plot.type = c('climate_normals', 'monthly_distributions'),
                        unit.system = c('international', 'standard'),
                        start = 1,
                        end = 365,
                        months = 1:12,
                        vars = c('maxt', 'mint', 'crain', 'cradn'),
                        tgt.year = 'last') {



  req.namespaces <- c('ggplot2', 'patchwork')
  for (ns in req.namespaces) {
    if(!requireNamespace(ns, quietly = TRUE)){
      warning('The ', ns, ' package is required for this function')
      return(NULL)
    }
  }

  weather.data <- x
  plot.type <- match.arg(plot.type)

  if(!inherits(weather.data, 'met'))
    stop('weather.data must be a met object created with pa_get_weather_shp')

  if(tgt.year != 'last' && (!is.numeric(tgt.year) || !any(grepl(tgt.year, weather.data$year))))
    stop('tgt.year has to be last or a year in the data set')

  if(tgt.year != 'last'){crt.year <- tgt.year} else{ crt.year <- max(unique(weather.data$year))}

  unit.system <- match.arg(unit.system)
  vars <- match.arg(vars, several.ok = TRUE)
  weather.data <- as.data.frame(weather.data)
  weather.data <- subset(weather.data, day >= start & day <= end)
  plt.units <- c('\u00B0C', '\u00B0C', 'MJ/m2', 'mm')

  if(unit.system == 'standard'){
    weather.data <- .pa_convert_met_to_standard(weather.data)
    plt.units <- c('\u00B0F', '\u00B0F', 'MJ/m2', 'in')

  }


  if (plot.type == 'climate_normals'){

    weather.data$crain <- with(weather.data, stats::ave(rain, year, FUN = cumsum))
    weather.data$cradn <- with(weather.data, stats::ave(radn, year, FUN = cumsum))
    weather.summary <-  do.call(data.frame, stats::aggregate(weather.data[c('maxt', 'mint', 'radn', 'rain', 'crain', 'cradn')],
                                                             by = weather.data['day'],
                                                             function(x) c(mean = mean(x), sd = stats::sd(x), max = max(x), min = min(x))))

    weather.summary$category <- 'historical'
    crt.weather <- subset(weather.data, year == crt.year)
    crt.weather$category <- crt.year

    weather.summary$date <- as.Date(weather.summary$day, '%j', origin = as.Date('2019-12-31'))
    crt.weather$date <- as.Date(crt.weather$day, '%j', origin = as.Date('2019-12-31'))


    cmaps <- c('historical normal' = 'darkolivegreen3',
               'record maximum' = 'tomato1',
               'record minimum' = 'steelblue3',
               'current year' = 'navyblue')
    names(cmaps)[4] <- crt.year

    plt.list <- list(
      maxt =  ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = weather.summary, ggplot2::aes(x = date, ymin = .data[['maxt.mean']] - .data[['maxt.sd']], ymax = .data[['maxt.mean']] + .data[['maxt.sd']], fill =  names(cmaps)[1]), alpha = 0.5)+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['maxt.max']], col = names(cmaps)[2]))+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['maxt.min']], col = names(cmaps)[3]))+
        ggplot2::geom_line(data = crt.weather, ggplot2::aes(x = date, y = .data[['maxt']], col = names(cmaps)[4])) +
        ggplot2::labs(y = paste0('Maximum\ntemperature, ', plt.units[1])),

      mint = ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = weather.summary, ggplot2::aes(x = date, ymin = .data[['mint.mean']] - .data[['mint.sd']], ymax = .data[['mint.mean']] + .data[['mint.sd']], fill =  names(cmaps)[1]), alpha = 0.5)+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['mint.max']], col = names(cmaps)[2]))+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['mint.min']], col = names(cmaps)[3]))+
        ggplot2::geom_line(data = crt.weather, ggplot2::aes(x = date, y = .data[['mint']], col = names(cmaps)[4])) +
        ggplot2::labs(y = paste0('Minimum\ntemperature, ', plt.units[2])),

      cradn = ggplot2::ggplot() +
        ggplot2::geom_ribbon(data = weather.summary, ggplot2::aes(x = date, ymin = .data[['cradn.mean']] - .data[['cradn.sd']], ymax = .data[['cradn.mean']] + .data[['cradn.sd']], fill =  names(cmaps)[1]), alpha = 0.5)+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['cradn.max']], col = names(cmaps)[2]))+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['cradn.min']], col = names(cmaps)[3]))+
        ggplot2::geom_line(data = crt.weather, ggplot2::aes(x = date, y = .data[['cradn']], col = names(cmaps)[4])) +
        ggplot2::labs(y = paste0('Cumulative\nradiation, ', plt.units[3])),

      crain =  ggplot2::ggplot()+
        ggplot2::geom_ribbon(data = weather.summary, ggplot2::aes(x = date, ymin = .data[['crain.mean']] - .data[['crain.sd']], ymax = .data[['crain.mean']] + .data[['crain.sd']], fill =  names(cmaps)[1]), alpha = 0.5)+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['crain.max']], col = names(cmaps)[2]))+
        ggplot2::geom_line(data = weather.summary, ggplot2::aes(x = date, y = .data[['crain.min']], col = names(cmaps)[3]))+
        ggplot2::geom_line(data = crt.weather, ggplot2::aes(x = date, y = .data[['crain']], col = names(cmaps)[4])) +
        ggplot2::labs(y = paste0('Cumulative\nrain, ', plt.units[4])) )

    out <- patchwork::wrap_plots(plt.list[names(plt.list) %in% vars])+
      patchwork::plot_layout(ncol = 1, guides = 'collect') &
      ggplot2::scale_x_date(date_labels = '%b-%d')&
      ggplot2::scale_color_manual(values = cmaps) &
      ggplot2::scale_fill_manual(values = cmaps) &
      ggplot2::theme_bw() &
      ggplot2::theme(legend.position = 'top',
                     legend.direction = 'horizontal',
                     plot.margin = ggplot2::unit(c(0, 1, 0, 1), 'line')) &
      ggplot2::labs(x = 'Date', col = '', fill = '')
    return(out)
  }

  if (plot.type == 'monthly_distributions'){

    year <- .data <- month <- day <- NULL

    ## this plot assumes that the latest year in the data set is the current year
    crt.year <- max(unique(weather.data$year))

    plt.units <- c('\u00B0C', '\u00B0C', 'MJ/m2', 'mm')
    if(unit.system == 'standard'){
      weather.data <- .pa_convert_met_to_standard(weather.data)
      plt.units <- c('\u00B0F', '\u00B0F', 'MJ/m2', 'in')
    }

    weather.data$date <- as.Date(weather.data$day, '%j', origin = '2019-12-31')
    weather.data$month <- as.numeric(strftime(weather.data$date, '%m'))
    weather.data$month.abb <- factor(strftime(weather.data$date, '%b'),
                                     ordered = TRUE, levels = month.abb)
    weather.data <- subset(weather.data, month %in% months)


    cols.to.agg.by <- c('year', 'month', 'month.abb')
    weather.summary <- do.call(data.frame,
                               stats::aggregate(weather.data[c('maxt', 'mint', 'rain', 'radn')],
                                         weather.data[cols.to.agg.by],
                                         function(x) c(mean = mean(x), sum = sum(x))))

    crt.weather <- subset(weather.data, year == crt.year)
    crt.weather <- do.call(data.frame,
                           stats::aggregate(crt.weather[c( 'maxt', 'mint', 'rain', 'radn')],
                                     crt.weather[cols.to.agg.by],
                                     function(x) c(mean = mean(x), sum = sum(x))))


    cmaps <- c('maximum temperature' = 'tomato1',
               'minimum temperature' = 'steelblue3',
               'current year' = 'navyblue',
               'cumulative rain' = 'darkolivegreen3',
               'cumulative radiation' = 'orange',
               'historical mean' = 'black')
    names(cmaps)[3] <- crt.year


    plt.list <- list(

      maxt = ggplot2::ggplot()+
        ggplot2::geom_density(data = weather.summary, ggplot2::aes( x = .data[['maxt.mean']], fill = names(cmaps)[1]), alpha = 0.5)+
        ggplot2::geom_vline(data = crt.weather, ggplot2::aes( xintercept = .data[['maxt.mean']], col = names(cmaps)[3]), linetype = 2)+
        ggplot2::geom_text(data = crt.weather, ggplot2::aes(x = .data[['maxt.mean']], y = Inf, label = crt.year), vjust = 1, hjust = 1)+
        ggplot2::scale_fill_manual(values = cmaps)+
        ggplot2::scale_color_manual(values = cmaps)+
        ggplot2::labs(col = '', fill = '', x = paste0('Max.\ntemperature, ', plt.units[1]), y = 'Density')+
        ggplot2::theme_bw()+
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_blank())+
        ggplot2::facet_wrap(~month.abb,  ncol = 1,  strip.position = 'right'),

      mint = ggplot2::ggplot()+
        ggplot2::geom_density(data = weather.summary, ggplot2::aes( x = .data[['mint.mean']], fill = names(cmaps)[2]), alpha = 0.5)+
        ggplot2::geom_vline(data = crt.weather, ggplot2::aes( xintercept = .data[['mint.mean']], col = names(cmaps)[3]), linetype = 2)+
        ggplot2::geom_text(data = crt.weather, ggplot2::aes(x = .data[['mint.mean']], y = Inf, label = crt.year), vjust = 1, hjust = 1)+
        ggplot2::scale_fill_manual(values = cmaps)+
        ggplot2::scale_color_manual(values = cmaps)+
        ggplot2::labs(col = '', fill = '', x = paste0('Min \ntemperature, ', plt.units[2]), y = 'Density')+
        ggplot2::theme_bw()+
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_blank())+
        ggplot2::facet_wrap(~month.abb,  ncol = 1,  strip.position = 'right'),

      cradn = ggplot2::ggplot()+
        ggplot2::geom_density(data = weather.summary, ggplot2::aes( x = .data[['radn.sum']], fill = names(cmaps)[5]), alpha = 0.5)+
        ggplot2::geom_vline(data = crt.weather, ggplot2::aes( xintercept = .data[['radn.sum']], col = names(cmaps)[3]), linetype = 2)+
        ggplot2::geom_text(data = crt.weather, ggplot2::aes(x = .data[['radn.sum']], y = Inf, label = crt.year), vjust = 1, hjust = 1)+
        ggplot2::scale_fill_manual(values = cmaps)+
        ggplot2::scale_color_manual(values = cmaps)+
        ggplot2::labs(col = '', fill = '', x = paste0('Cumulative\nradiation, ', plt.units[3]), y = 'Density')+
        ggplot2::theme_bw()+
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_blank())+
        ggplot2::facet_wrap(~month.abb,  ncol = 1,  strip.position = 'right'),

      crain = ggplot2::ggplot()+
        ggplot2::geom_density(data = weather.summary, ggplot2::aes( x = .data[['rain.sum']], fill = names(cmaps)[4]), alpha = 0.5)+
        ggplot2::geom_vline(data = crt.weather, ggplot2::aes( xintercept = .data[['rain.sum']], col = names(cmaps)[3]), linetype = 2)+
        ggplot2::geom_text(data = crt.weather, ggplot2::aes(x = .data[['rain.sum']], y = Inf, label = crt.year), vjust = 1, hjust = 1)+
        ggplot2::scale_fill_manual(values = cmaps)+
        ggplot2::scale_color_manual(values = cmaps)+
        ggplot2::labs(col = '', fill = '', x = paste0('Cumulative\nrain, ', plt.units[4]), y = 'Density')+
        ggplot2::theme_bw()+
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_blank())+
        ggplot2::facet_wrap(~month.abb,  ncol = 1,  strip.position = 'right')
    )

    plt.list <- plt.list[names(plt.list) %in% vars]

    plt.list[[length(plt.list)]] <- plt.list[[length(plt.list)]] +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = 'transparent', colour = 'transparent'),
                     strip.text = ggplot2::element_text())

    out <- patchwork::wrap_plots(plt.list)+
      patchwork::plot_layout(guides = 'collect', nrow = 1) &
      ggplot2::theme(panel.spacing = ggplot2::unit(0, 'line'),
                     strip.background = ggplot2::element_rect(fill = 'transparent', colour = 'transparent'),
                     legend.position = 'top',
                     legend.direction = 'horizontal',
                     plot.margin = ggplot2::unit(c(1, 0, 1, 0), 'line'),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
    return(out)

  }

}
