#' Format input data to plot with ipath
#' @param data input data.frame
#' @param id_column name of the column containing ids to plot on the ipath map. Must be one of the supported data types and in the correct format (see: https://pathways.embl.de/help.cgi)
#' @param color_column name of the column containing values to vary colors in the ipath plot
#' @param width_column name of the column containing values to vary width in the ipath plot
#' @param opacity_column name of the column containing values to vary opacity in the ipath plot
#' @param min_path_width minimum width of paths with a match
#' @param max_path_width maximum width of paths with a match. Default width for when no wwidth_column is provided.
#' @param color_type type of color palette to map values to (discrete or continous)
#' @param color_cutoff value which seperates values into higher or lower for discrete color mapping
#' @param palette RColorBrewer pallete name
#' @return data.frame formatted for sending to ipath
#' @export
#' @author Beth Signal
create_ipath_data = function(data,
                             id_column,
                             color_column=NA,
                             width_column=NA,
                             opacity_column=NA,
                             min_path_width = 3,
                             max_path_width = 10,
                             color_type = "discrete",
                             color_cutoff = 0,
                             pallete = "RdBu"){
  data = as.data.frame(data)
  output_data = as.data.frame(data[,which(colnames(data) == id_column)])

  if(!is.na(color_column)){
    output_data = cbind(output_data,
                        color = map_values_to_colors(as.numeric(data[,which(colnames(data) == color_column)]),
                                                     color_type = color_type,color_cutoff=color_cutoff,pallete = pallete))
  }else{
    output_data = cbind(output_data,
                        color = map_values_to_colors(rep(1, nrow(output_data)),
                                                     color_type = color_type,color_cutoff=color_cutoff,pallete = pallete))
  }
  if(!is.na(width_column)){
    output_data$width = map_values_to_range(as.numeric(data[,which(colnames(data) == width_column)]),
                                            min_val = 5, max_val = 10)
    output_data$width = paste0("W", round(output_data$width, 2))
  }else{
    output_data$width = paste0("W", min_path_width)
  }
  if(!is.na(opacity_column)){
    output_data = cbind(output_data,
                        opacity = map_values_to_range(data[,which(colnames(data) == opacity_column)],
                                                    min_val = 0, max_val = 1))
  }else{
    output_data = cbind(output_data, opacity = 1)
  }

  colnames(output_data)[1] = id_column
  return(output_data)
}

#' Map a vector of values to hex colors
#' @param values numeric vector of values
#' @param color_type type of color palette to map values to (discrete or continous)
#' @param color_cutoff value which seperates values into higher or lower for discrete color mapping
#' @param palette RColorBrewer pallete name, viridis palette name (if continous), or two colors to create a colorRampPalette from.
#' @return vector of hex colors
#' @keywords internal
#' @importFrom RColorBrewer brewer.pal
#' @author Beth Signal
#' @examples
#' map_values_to_colors(seq(-4, 5))
#' map_values_to_colors(seq(-4, 5), color_type = "continous")
#' map_values_to_colors(seq(-4, 5), color_type = "continous", pallete = "magma")
#' map_values_to_colors(seq(-4, 5), color_type = "continous", pallete = c("lightblue", "blue"))

map_values_to_colors = function(values,
                                color_type = "discrete",
                                color_cutoff = 0,
                                pallete = "RdBu"){

  if(color_type == "discrete"){
    # make RColorBrewer pallete from pallete name
    # uses n=5 to get good insense extremes
    colors = RColorBrewer::brewer.pal(5, pallete)[c(1,5)]
    value_cols = ifelse(values > color_cutoff, colors[1], colors[2])
    return(value_cols)
  }else{

    ## Use n equally spaced breaks to assign each value to n-1 equal sized bins
    ii = cut(values, breaks = seq(min(values), max(values), len = 100),
             include.lowest = TRUE)
    if(pallete[1] %in% c("A", "B","C","D","E", "magma","inferno", "plasma","viridis", "cividis")){
      value_cols = viridis(100, option = pallete[1])

      if(all(nchar(value_cols) == 9)){
        value_cols = gsub("FF$", "", value_cols)
      }
    }else{
      if(length(pallete) == 2){
        value_cols = colorRampPalette(pallete)(99)
      }else{
        colors = RColorBrewer::brewer.pal(5, pallete)
        value_cols = colorRampPalette(colors)(99)
      }
    }

    ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
    value_cols = value_cols[ii]

  }

}

#' Map a vector of values to a prespecified range
#' @param values numeric vector of values
#' @param min_val value to map minimum input value to
#' @param max_val value to map maximum input value to
#' @keywords internal
#' @return vector of values
#' @author Beth Signal
#' @examples
#' map_values_to_range(seq(0, 10), min_val = 2, max_val=4)
#'
map_values_to_range = function(values, min_val= 3, max_val=10){

  new_values = (values - min(values))/(max(values) - min(values))

  new_values = (new_values*(max_val - min_val)) + min_val

  return(new_values)
}

#' Format ipath data.frame as parameters
#' @param ipath_data data.frame formatted for sending to ipath
#' @param default_color Hex color value. Non-selected parts of the map will use this color by default.
#' @param default_width Default path/edge width (px) for the non-selected parts of the map.
#' @param default_radius Compound radius (px) for the non-selected parts of the map.
#' @param default_opacity Opacity setting for the non-selected parts of the map. Zero represents full transparency, while one is fully opaque.
#' @param keep_colors Keep original colors? When TRUE, non-selected parts of the map will not change color to the default one.
#' @param whole_modules Select whole modules? When TRUE, any KEGG module with at least one matching edge or compound will be highlighted.
#' @param whole_pathways Select whole pathways? If TRUE, any pathway with at least one matching edge or compound will be highlighted.
#' @param query_reactions Query reaction compounds. If TRUE, compound presence within each edges reactions will also be checked
#' @param tax_filter string containing either NCBI tax ID or KEGG 3 letter species code(s). Only pathways present in selected species will be included in the map.
#' @param export_dpi numeric value for DPI of SVG output.
#' @keywords internal
#' @return list of parameters for ipath
#' @author Beth Signal
to_parameters = function(ipath_data,
                         default_color = '#bdbdbd',
                         default_width = 3,
                         default_radius = 7,
                         default_opacity = 1,
                         keep_colors = FALSE,
                         whole_modules = FALSE,
                         whole_pathways = FALSE,
                         query_reactions = FALSE,
                         tax_filter = '',
                         export_dpi = 1200){

  selection  = paste0(apply(ipath_data, 1, function(x) paste(x, collapse =" ")),
                      collapse = "\n")

  ipath_parameters = list(selection = selection,
            export_type = 'svg',
            keep_colors = ifelse(keep_colors, 1, 0),
            include_metabolic = 1,
            include_secondary = 0,
            include_antibiotic = 0,
            include_microbial = 0,
            whole_modules = ifelse(whole_modules, 1, 0),
            default_opacity = ifelse(default_opacity, 1, 0),
            whole_pathways = ifelse(whole_pathways, 1, 0),
            default_width = default_width,
            default_color = default_color,
            default_radius = default_radius,
            query_reactions = ifelse(query_reactions, 1, 0),
            tax_filter = tax_filter,
            export_dpi = export_dpi
  )
  return(ipath_parameters)
}
