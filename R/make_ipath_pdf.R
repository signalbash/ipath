#' Write an iPath3 metabolic pathways map to pdf
#' @param ipath_data data.frame formatted for sending to ipath
#' @param file_name pdf file name to export plot.
#' @param default_color Hex color value. Non-selected parts of the map will use this color by default.
#' @param default_width Default path/edge width (px) for the non-selected parts of the map.
#' @param default_radius Compound radius (px) for the non-selected parts of the map.
#' @param textbox_textcol hex value for text (in box) color
#' @param textbox_col hex value for box color
#' @export
#' @import httr
#' @importFrom rsvg rsvg_pdf
#' @author Beth Signal
#' @examples
print_ipath_pdf = function(ipath_data, file_name=NULL, default_width = 2,
                default_color = '#bdbdbd',default_radius = 5,
                textbox_textcol = NULL, textbox_col = NULL){

  params = to_parameters(ipath_data,
                         default_width = default_width,
                         default_color = default_color,
                         default_radius = default_radius)

  ipath_result = httr::POST("https://pathways.embl.de/mapping.cgi?map=metabolic",
                             body = params, encode = "form")

  pathway_map = httr::content(ipath_result, as="text", encoding = "UTF-8")
  pathway_map = colour_on_top(pathway_map, ipath_data)

  if(!is.null(textbox_textcol) | !is.null(textbox_col)){
    pathway_map = recolor_textbox(pathway_map, box_color = textbox_col, text_color = textbox_textcol)
  }

  #rsvg::rsvg_pdf(charToRaw(pathway_map), file = file_name)
  # scale the svg so it's easily viewble
  pathway_map = gsub("height='2250' width='3774' viewBox=\"-10 -10 3794 2270\"", "height='16cm' width='26cm' viewBox=\"0 0 26.0cm 16.0cm\"", pathway_map)
  pathway_map = gsub("<style type='text/css'>", "<g transform=\"translate(0, 0) scale(0.265) \"><style type='text/css'>", pathway_map)
  pathway_map = gsub("</g>\n</svg>", "</g>\n</g>\n</svg>", pathway_map)

  # remove crap
  #pathway_map = gsub("[\n][<]text(.*?)Kanehisa(.*?)text[>]", "", pathway_map)
  #pathway_map = gsub("[\n][<]text(.*?)5/11/17(.*?)text[>]", "", pathway_map)

  all_newlines = str_locate_all(pathway_map, "\n")

  rm_loc = str_locate_all(pathway_map, "5/11/17")
  diffs = all_newlines[[1]][,1] - rm_loc[[1]][1]
  str_sub(pathway_map, all_newlines[[1]][which(diffs == max(diffs[diffs < 0])),1],all_newlines[[1]][which(diffs == min(diffs[diffs > 0])),1]-1) = ""

  rm_loc = str_locate_all(pathway_map, "Kanehisa")
  diffs = all_newlines[[1]][,1] - rm_loc[[1]][1]
  str_sub(pathway_map, all_newlines[[1]][which(diffs == max(diffs[diffs < 0])),1],all_newlines[[1]][which(diffs == min(diffs[diffs > 0])),1]-1) = ""

  if(!is.null(file_name)){
    write.table(pathway_map, file = file_name, row.names = F, col.names = F, quote = F, sep="\n")
    #svg_xml = data.table::fread(file_name,  data.table=FALSE, sep='\n', header = F)[,1]
  }

  svg_xml = unlist(str_split(pathway_map, "\n"))

  return(svg_xml)

}

#' Convert svg data to tables for ggplot
#' @param svg_xml data.frame containing svg data
#' @export
#' @import stringr
#' @author Beth Signal
#' @examples
svg_to_table = function(svg_xml){

  #svg_xml = data.table::fread(svg_file, data.table = FALSE, header = FALSE)

  layer = c(1:5)
  layer_loc = vector()
  for(i in layer){
    layer_loc[i] = grep(paste0("layer", i), svg_xml)
  }

  layer_locations = data.frame(layer,layer_loc)

  ## Layer 1: paths ##
  layer_1_paths = grep("path", svg_xml[layer_locations$layer_loc[1]:layer_locations$layer_loc[2]], value = TRUE)
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_1_paths, "[ ]+"), "[[",3))))
  stroke = gsub("'", "",gsub('"',"",(gsub('stroke=',"",unlist(lapply(stringr::str_split(layer_1_paths, "[ ]+"), "[[",5))))))
  stroke_width = as.numeric(gsub("'", "",gsub('"',"",(gsub('stroke-width=',"",unlist(lapply(stringr::str_split(layer_1_paths, "[ ]+"), "[[",7)))))))

  path_locs = gsub("'", "",gsub('"',"",gsub('/>',"",gsub('d=',"",unlist(lapply(stringr::str_split(lapply(stringr::str_split(layer_1_paths, "d="), "[[",2), ">"), "[[",1))))))

  layer_1 = data.frame(opacity, stroke, stroke_width, id = path_locs)

  path_locs_CZ = path_locs[grepl("C", path_locs) | grepl("Z", path_locs)]
  path_locs_L = path_locs[!(grepl("C", path_locs) | grepl("Z", path_locs))]

  coords_CZ = do.call("rbind", lapply(unique(path_locs_CZ), function(x) cbind(id = x, make_svg_df(x))))
  layer_1_curves = layer_1[layer_1$id %in% path_locs_CZ,]
  m = match(coords_CZ$id, layer_1_curves$id)
  layer_1_curves = cbind(layer_1_curves[m,], coords_CZ[,-which(colnames(coords_CZ) == "id")])
  color_order = as.data.frame(table(layer_1_curves$stroke))
  base_col = color_order$Var1[which.max(color_order$Freq)]
  layer_1_curves = layer_1_curves[c(which(layer_1_curves$stroke == base_col),which(layer_1_curves$stroke != base_col)),]


  coords_L = cbind(id = path_locs_L, make_svg_df_straight(path_locs_L))
  layer_1_lines = dplyr::full_join(layer_1[layer_1$id %in% path_locs_L,], coords_L, by="id")
  layer_1_lines = layer_1_lines[c(which(layer_1_lines$stroke == base_col),which(layer_1_lines$stroke != base_col)),]


  ## Layer 2: points ##
  layer_2_points = gsub("'","",grep("ellipse", svg_xml[layer_locations$layer_loc[2]:layer_locations$layer_loc[3]], value = TRUE))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",3))))
  rx = as.numeric(gsub("'", "",gsub('"',"",(gsub('rx=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",4)))))))
  ry = as.numeric(gsub("'", "",gsub('"',"",(gsub('ry=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",5)))))))
  cx = as.numeric(gsub("'", "",gsub('"',"",(gsub('cx=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",7)))))))
  cy = as.numeric(gsub("'", "",gsub(">(.*)","",gsub('"',"",(gsub('cy=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",8))))))))
  fill = gsub("'", "",gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",6))))))
  layer_2 = data.frame(type = "ellipse", opacity, rx,ry,cx,cy,fill)

  ## Layer 3: Path text ##
  layer_3_text = gsub("'","",grep("text", svg_xml[layer_locations$layer_loc[3]:layer_locations$layer_loc[4]], value = TRUE))
  font_size = as.numeric(gsub('px;opacity:',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",3))))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",4))))
  fill = gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",5)))))
  x = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('x=',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",6)))))))
  y = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('y=',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",7)))))))
  text = gsub("[ ]$","",gsub("</text", "",unlist(lapply(stringr::str_split(layer_3_text, ">"), "[[", 2))))
  layer_3 = data.frame(type="text", font_size, opacity, fill, x,y,text)

  ## layer 4: textrects ##
  layer_4_rect = gsub("'","",grep("rect", svg_xml[layer_locations$layer_loc[4]:layer_locations$layer_loc[5]], value = TRUE))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",3))))
  fill = gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",4)))))
  stroke_width = as.numeric(gsub('"',"",(gsub('stroke-width=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",5))))))
  x = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('x=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",6)))))))
  y = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('y=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",7)))))))
  height = as.numeric(gsub('"',"",(gsub('height=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",8))))))
  width = as.numeric(gsub('"',"",(gsub('width=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",9))))))
  rx = as.numeric(gsub('"',"",(gsub('rx=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",10))))))
  ry = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('ry=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",11)))))))
  layer_4 = data.frame(type="rect", opacity, fill, stroke_width, x, y, height, width, rx,ry)

  ## Layer 5: Big text ##
  layer_5_text = gsub("'","",grep("text", svg_xml[layer_locations$layer_loc[5]:length(svg_xml)], value = TRUE))
  layer_5_text = layer_5_text[which(stringr::str_sub(layer_5_text, 2,5) == "text")]
  font_size = as.numeric(gsub('px;opacity:',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",3))))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",4))))
  fill = gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",5)))))
  x = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('x=',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",6)))))))
  y = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('y=',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",7)))))))
  text = gsub("[ ]$","",gsub("<(.*)","", gsub(">(.*)", "",unlist(lapply(stringr::str_split(layer_5_text, ">"), "[[", 2)))))
  layer_5 = data.frame(type="text", font_size, opacity, fill, x,y,text)

  all_layers = list()
  all_layers[['layer_1_curves']] = layer_1_curves
  all_layers[['layer_1_lines']] = layer_1_lines
  all_layers[['layer_2']] = layer_2
  all_layers[['layer_3']] = layer_3
  all_layers[['layer_4']] = layer_4
  all_layers[['layer_5']] = layer_5

  return(all_layers)

}

#' Convert a svg path specification to data.frame compatible with ggplot
#' Used for paths with curves (C)
#' @param input string of svg path
#' @param rebase set the minimum path value to 0?
#' @export
#' @import stringr
#' @importFrom plyr arrange
#' @importFrom reshape2 melt
#' @importFrom plyr rbind.fill
#' @author Beth Signal
#' @examples
#'
make_svg_df = function(input, rebase = F){
  # SVG path vector types
  types = unlist(str_split(gsub("[0-9,. ]*","",input), ""))
  # locations of each path type indicator in input
  locs = unlist(lapply(str_locate_all(input, unique(types)), function(x) x[,1]))
  # make data.frame of each vector type in input, and its location within the input string
  types = str_sub(input, locs, locs)
  type_locs = rbind(data.frame(types, locs), data.frame(types = "end",locs = nchar(input)))
  type_locs = plyr::arrange(type_locs, locs)
  type_locs = type_locs[!duplicated(with(type_locs, paste(types,locs))),]

  # get SVG coordniates for each path type
  all_coords=NULL
  for(i in 1:(nrow(type_locs)-1)){
    all_coords = rbind(all_coords,
                       cbind(type = type_locs$types[i],
                             as.data.frame(matrix(as.numeric(unlist(str_split(
                               gsub(" ","",str_sub(input, type_locs$locs[i]+1, type_locs$locs[i+1]-1)),
                               ","))), ncol=6, byrow = TRUE))
                       )
    )
  }

  all_coords$type = as.character(all_coords$type)
  all_coords[which(all_coords$type != "C"),c(4:7)] = NA

  if(all_coords$type[nrow(all_coords)] == "Z"){
    all_coords[nrow(all_coords),-1] = all_coords[1,-1]
  }

  if(any(all_coords$type == "C")){
    all_coords[all_coords$type=="C",] = all_coords[all_coords$type=="C",c(1,6,7,2,3,4,5)]
  }

  # if TRUE start coordinates at 0
  if(rebase == TRUE){
    min_x = min(all_coords$V1)
    min_y = min(all_coords$V2)
    all_coords = all_coords
    all_coords$V1 = all_coords$V1 - min_x
    all_coords$V2 = all_coords$V2 - min_y
    all_coords$V3 = all_coords$V3 - min_x
    all_coords$V4 = all_coords$V4 - min_y
    all_coords$V5 = all_coords$V5 - min_x
    all_coords$V6 = all_coords$V6 - min_y
  }

  # convert from x-y points to x1,y1,x2,y2 segments (or curves)
  all_coords_plot = all_coords[-1,]
  colnames(all_coords_plot)[-1] = c("x2","y2", "xb1","yb1","xb2","yb2")
  all_coords_plot = cbind(x1=all_coords$V1[-nrow(all_coords)],y1=all_coords$V2[-nrow(all_coords)],all_coords_plot)
  all_coords_plot = all_coords_plot[,c(3,1,2,4:9)]
  # fix for Z
  all_coords_plot$type = ifelse((all_coords_plot$type == "C"), "curve","segment")

  if(any(all_coords_plot$type == "curve")){
    all_coords_plot_circles = all_coords_plot[all_coords_plot$type == "curve",]
    all_coords_plot_circles$bez_id = paste0(input,".b.curve", c(1:nrow(all_coords_plot_circles)))
    all_coords_plot_bez = cbind(reshape2::melt(all_coords_plot_circles[,c(2,4,6,8,10)], id.vars = c("bez_id")),
                                reshape2::melt(all_coords_plot_circles[,c(3,5,7,9,10)], id.vars = c("bez_id"))[,3])
    colnames(all_coords_plot_bez) = c("bez_id", "b_type", "x1","y1")
    all_coords_plot_bez$b_type = factor(all_coords_plot_bez$b_type, levels =c("x1","xb1","xb2","x2"))
    all_coords_plot_bez = plyr::arrange(all_coords_plot_bez, bez_id, b_type)
    all_coords_plot_bez$type="curve"
    all_coords_plot = all_coords_plot[all_coords_plot$type == "segment", c(1:5)]

    all_coords_plot = plyr::rbind.fill(all_coords_plot, all_coords_plot_bez)
  }else{
    all_coords_plot = all_coords_plot[, c(1:5)]
    all_coords_plot$bez_id = NA
    all_coords_plot$b_type = NA
  }
  all_coords_plot$id = input

  return(all_coords_plot)
}

#' Convert a svg stright line path specification to data.frame compatible with ggplot
#' Used for paths without curves (C)
#' @param input string of svg path
#' @export
#' @import stringr
#' @import ggplot2
#' @import ggforce
#' @author Beth Signal
#' @examples
#'
make_svg_df_straight = function(input){

  type = "segment"
  x1 = as.numeric(str_sub(unlist(lapply(str_split(input, ","),"[[",1)),2,-1))
  y1 = as.numeric(unlist(lapply(str_split(lapply(str_split(input, ","),"[[",2), "L"),"[[",1)))
  x2 = as.numeric(unlist(lapply(str_split(lapply(str_split(input, ","),"[[",2), "L"),"[[",2)))
  y2 = as.numeric(unlist(lapply(str_split(input, ","),"[[",3)))

  all_coords_plot = data.frame(type,x1,y1,x2,y2,curvature=NA)
  return(all_coords_plot)
}

#' make a ipath ggplot
#' @param plot_layers list of svg layer data produced by svg_to_table()
#' @export
#' @import stringr
#' @import ggplot2
#' @import ggforce
#' @author Beth Signal
#' @examples
#'
make_ipath_ggplot = function(plot_layers){

  layer_1_curves = plot_layers[['layer_1_curves']]
  layer_1_lines = plot_layers[['layer_1_lines']]
  layer_2 = plot_layers[['layer_2']]
  layer_3 = plot_layers[['layer_3']]
  layer_4 = plot_layers[['layer_4']]
  layer_5 = plot_layers[['layer_5']]

  bg_stroke = table(c(as.character(layer_1_curves$stroke), as.character(layer_1_lines$stroke), as.character(layer_2$fill)))
  bg_stroke_col = names(which.max(bg_stroke))


  layer_1 = plyr::rbind.fill(layer_1_curves, layer_1_lines)
  stroke_cols = as.character(sort(unique(layer_1$stroke)))
  layer_1_bottom = layer_1[layer_1$stroke == bg_stroke_col,]
  layer_1_top = layer_1[layer_1$stroke != bg_stroke_col,]

  layer_1_bottom$stroke = factor(layer_1_bottom$stroke, levels = stroke_cols)
  layer_1_top$stroke = factor(layer_1_top$stroke, levels = stroke_cols)

  fill_cols = as.character(sort(unique(layer_4$fill)))
  layer_4$fill = factor(layer_4$fill, levels=fill_cols)

  rects_left = layer_4
  rects_left$x = rects_left$x
  rects_left$width = rects_left$height

  rects_right = layer_4
  rects_right$x = rects_right$x + rects_right$width - (rects_right$rx*2)
  rects_right$width = rects_right$height
  rects_rounded = rbind(rects_left, rects_right)

  p =
    ggplot() +
    geom_bezier(aes(x = x1, y = y1, group=bez_id,col=stroke, size=stroke_width),
                data = layer_1_bottom[layer_1_bottom$type == "curve",]) +
    geom_segment(data=layer_1_bottom[layer_1_bottom$type == "segment", ],
                 aes(x=x1, y=y1, xend = x2, yend=y2, col=stroke, size=stroke_width)) +
    geom_bezier(aes(x = x1, y = y1, group=bez_id,col=stroke, size=stroke_width),
                  data = layer_1_top[layer_1_top$type == "curve",]) +
    geom_segment(data=layer_1_top[layer_1_top$type == "segment", ],
                   aes(x=x1, y=y1, xend = x2, yend=y2, col=stroke, size=stroke_width)) +

    geom_point(data=layer_2, aes(x=cx,y=cy), size= 0.5, col=layer_2$fill[1]) +
    geom_text(data=layer_3, aes(x=x,y=y, label = text), size = 1.6, hjust = 0, vjust=0) +
    #geom_rect(data=layer_4, aes(xmin=x, ymin = y, xmax =x+width, ymax=y+height, fill=fill)) +
    geom_ellipse(data = rects_rounded, aes(x0=x+(width/2),y0=y+(height/2), a=width/2, b=height/2, angle=0, m1=2, m2=2, fill=fill), col=NA) +
    geom_rect(data = layer_4, aes(xmin=x+rx, ymin=y,ymax=y+height, xmax=x+width-rx, fill=fill), col=NA) +
    geom_text(data=layer_5, aes(x=x,y=y, label = text), size = 2.2, hjust = 0, vjust=0, fontface = "bold", col="white") +
    scale_y_continuous(trans="reverse") +
    theme_void() +
    theme(legend.position = "none") +
    scale_color_manual(values = stroke_cols)+
    scale_fill_manual(values = fill_cols) +
    scale_size_continuous(range = c(0.5,2), limits = c(2,10))

  return(p)

}

