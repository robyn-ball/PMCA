#' Generates a k-partite graph from the mapped results
#'
#' @param mapped.results a list of doubles containing the FPRs for the desired data
#' @param focus a string specifying the column of interest
#' @param threshold a float for the number to threshold on
#' @param debug boolean to show debug statements
#'
#' @return a plotly figure
#'
#'
#' @export
#'
#' @examples see Visualization_README
make.graph <- function(dataframe_datalist, focus, default_threshold=0.001, apply_threshold=NULL, debug=FALSE) {

  # Obtain a list of edgelists from the dataframe list; specify node of interest
  edgelist_datalist <- kp$get_edgelists(dataframe_datalist, focus, debug=debug)

  # Obtain a list of nodes from the set of edgelists
  nodelist <- kp$get_nodes(edgelist_datalist)

  # Concatenate all the edgelists together and position the edges/nodes
  positioned_edges <- kp$get_edge_positions(edgelist_datalist, nodelist, debug=debug)

  # Plot the visualization
  fig <- plot_ly(height=1000)

  # Add traces for positioned edges
  for (i in 1:nrow(positioned_edges)) {
    e <- positioned_edges[i, ]

    anti_value <- e[["anti"]]

    fig <- fig %>%
      add_trace(
        x = c(e[['x1']], e[['x2']]),
        y = c(e[['y1']], e[['y2']]),
        type = "scatter",
        mode = "lines",
        line = list(
          width = e[['weight']] * 3,
          color = ifelse(anti_value, "#f54542", "#4560ba"),
          dash = ifelse(!e[['direct']], "dash", "solid")
        ),
        opacity = e[['weight']]
      )
  }

  # Add trace for nodelist markers and text
  fig <- fig %>%
    add_trace(
      x = nodelist[['x']],
      y = nodelist[['y']],
      type = "scatter",
      mode = "markers+text",
      marker = list(
        color = "#333",
        size = 10
      ),
      name = NULL
      # text = c(ifelse(grepl("_", n), paste(strsplit(n, "_")[[1]][-1], collapse = " "), n) for n in nodelist[['node']])
    )

  # Add annotations for nodelist nodes
  annotations <- list()
  for (i in 1:nrow(nodelist)) {
    node <- nodelist[i, ]

    annotation <- list(
      x = node[['x']],
      y = node[['y']],
      text = kp$format_text(node[['node']]),
      showarrow = FALSE,
      yshift = 10,
      font = list(
        color = "#333",
        size = 12
      ),
      bgcolor = "#fff"
    )

    annotations <- append(annotations, list(annotation))
  }

  # Update layout settings
  fig <- fig %>%
    layout(
      template = "plotly_white",
      showlegend = FALSE,
      xaxis = list(
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        showticklabels = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
      ),
      annotations = annotations
    )

  return(fig)

}
