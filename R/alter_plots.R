#' Change color of boxes surrounding major pathway names
#' @param pathway_map text containing svg source code
#' @param box_color hex value for box color
#' @param text_color hex value for text (in box) color
#' @return pathway map with altered box colors
#' @keywords internal
#' @import stringr
#' @author Beth Signal
recolor_textbox = function(pw_map, box_color = "#969696", text_color = '#FFFFFF'){

  # recolor rectangles around text
  if(!is.null(box_color)){
    g_splits = unlist(stringr::str_split(pw_map, "<g>"))
    rects = which(stringr::str_sub(g_splits,3,6) == "rect")
    g_splits[rects] = gsub("#[A-G,0-9]{6}", box_color, g_splits[rects])
  }
  #recolor text
  if(!is.null(text_color)){
    layer_locs = find_layer_locs(g_splits)
    rect_text = g_splits[layer_locs$indices[5]:layer_locs$indices[6]]
    g_splits[rect_text] = gsub("#[A-G,0-9]{6}", text_color, g_splits[rects])
  }
  pw_map2 = paste0(g_splits, collapse = "<g>")
  return(pw_map2)
}

#' Get indices for the start of each layer in svg text
#' @param g_splits vector of svg text split by the delimiter '<g>'
#' @return data.frame of layer numbers and their starting index
#' @keywords internal
#' @import stringr
#' @author Beth Signal
find_layer_locs = function(g_splits){

  indices = grep("layer", g_splits)
  layer_text = g_splits[indices]
  layer_number = unlist(lapply(stringr::str_split(lapply(stringr::str_split(layer_text, "layer"),"[[", 2), "'"), "[[",1))

  layer_number = c(layer_number, "end")
  indices = c(indices, length(g_splits))

  return(data.frame(layer_number, indices))

}

#' Reorder paths so selected (colored) paths are on top
#' @param svg_text text containing svg source code
#' @param ipath_data data.frame formatted for sending to ipath. Needs a 'color' column
#' @return pathway map (svg formatted text) with colored paths on top
#' @keywords internal
#' @import stringr
#' @importFrom spatialEco insert.values
#' @author Beth Signal
colour_on_top = function(svg_text, ipath_data){
  g_splits = unlist(stringr::str_split(svg_text, "<g>"))
  layer_locs = find_layer_locs(g_splits)
  mapped_cols = unique(ipath_data$color)

  g_splits_by_layer = list()
  g_splits_layer = g_splits[layer_locs$indices[-nrow(layer_locs)]]
  for(j in 1:(nrow(layer_locs)-1)){

    if(j == (nrow(layer_locs)-1)){
      g_splits_by_layer[[j]] = g_splits[(layer_locs$indices[j]+1):(layer_locs$indices[j+1])]
    }else{
      g_splits_by_layer[[j]] = g_splits[(layer_locs$indices[j]+1):(layer_locs$indices[j+1]-1)]
    }
  }

  for(j in 1:(nrow(layer_locs)-1)){

    indices = unlist(lapply(mapped_cols, function(x) grep(x, g_splits_by_layer[[j]])))

    if(length(indices) > 0){
      g_split_move = g_splits_by_layer[[j]][indices]
      g_split_rm = g_splits_by_layer[[j]][-indices]
      g_split_ins = c(g_split_rm, g_split_move)

      g_splits_by_layer[[j]] = c(g_split_rm, g_split_move)
    }


  }

  g_splits_v2 = vector()
  for(j in 1:(nrow(layer_locs)-1)){
    g_splits_v2 = c(g_splits_v2, g_splits_layer[j], g_splits_by_layer[[j]])
  }

  svg_text2 = paste0(g_splits_v2, collapse = "<g>")
  return(svg_text2)

}
