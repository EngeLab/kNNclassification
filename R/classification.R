
#' kNNclassify
#'
#' The function classifies samples in an unsupervised fashion by:
#' \enumerate{
#'   \item Running a principal component analysis.
#'   \item Uses Horn's technique to evaluate components to retain via
#'   \code{\link[paran]{paran}}.
#'   \item Finds k nearest neighbors in PCA space.
#'   \item Calculates the Euclidean distance between samples in PCA space.
#'   \item Constructs a weighted graph where each sample is connected to the k
#'   nearest neighbors with an edge weight = 1 - Euclidean distance.
#'   \item Uses the Louvain community detection algorithm to classify the
#'   samples.
#' }
#'
#' @name kNNclassify
#' @rdname kNNclassify
#' @author Jason T. Serviss
#' @param cpm matrix; Counts per million.
#' @param geneIdx Integer; Indices of genes to include in PCA.
#' @param PCiter Integer; Length 1 vector indicating the number of iterations to
#' perform when determining the numer of retained principal components.
#' @param k Integer; Length 1 vector indicating the number of nearest neighbors
#' for each sample.
#' @param pca Matrix; Optional pre-computed PCA. If NULL, PCA will be computed within
#' the function.
#' @param quietly Logical; indicates if function should be verbose.
#' @return Returns a tibble with two columns; the first indicating the sample
#' name and the second indicating the classification.
#' @examples
#'
#' #setup input data
#' s <- stringr::str_detect(colnames(testCounts), "^s")
#' e <- stringr::str_detect(rownames(testCounts), "^ERCC\\-[0-9]*$")
#' c <- testCounts[!e, s]
#' cpm <- t(t(c) / colSums(c) * 10^6)
#'
#' #pre-run PCA
#' pca <- gmodels::fast.prcomp(t(cpm), scale. = TRUE)$x
#'
#' #run KNN graph classification
#' kc <- kNNclassify(cpm, 1:nrow(c), 20, 15, pca = pca)
#'
#' #plot
#' pData <- merge(kc, matrix_to_tibble(pca[, 1:2], "sample"))
#' plot(pData$PC1, pData$PC2, col = rainbow(4)[pData$louvain], pch = 16)
#'
#' @export
#' @importFrom gmodels fast.prcomp
#' @importFrom paran paran
#' @importFrom dplyr left_join %>%
#' @importFrom tidygraph as_tibble

kNNclassify <- function(cpm, geneIdx, PCiter, k, pca = NULL, quietly = TRUE) {
  if(is.null(pca)) {
    pca <- gmodels::fast.prcomp(t(cpm[geneIdx, ]), scale. = TRUE)$x
  }

  pcs <- paran::paran(cpm[geneIdx, ], PCiter, quietly = quietly)$Retained
  knn <- getKNN(pca, pcs, k)
  dists <- distInPCAspace(pca, pcs)
  data <- dplyr::left_join(knn, dists, by = c("from", "to"))
  g <- constructGraph(data)
  g <- louvain(g, dist)
  tidygraph::as_tibble(g) %>%
    dplyr::select(name, louvain) %>%
    dplyr::rename(sample = name)
}

#' distInPCAspace
#'
#' Calculates Euclidean distance in PCA space.
#'
#' @name distInPCAspace
#' @keywords internal
#' @rdname distInPCAspace
#' @author Jason T. Serviss
#' @param pca Matrix; Principal components as columns, samples as rows.
#' @param pcs Integer; Length 1 vector indicating the max principal component to
#' retain.
#' @importFrom tidyr gather
#' @importFrom dplyr mutate "%>%"

distInPCAspace <- function(pca, pcs) {
  pca[, 1:pcs] %>%
    dist() %>%
    as.matrix() %>%
    `-` (1) %>%
    matrix_to_tibble("from") %>%
    tidyr::gather(to, dist, -from) %>%
    dplyr::mutate(dist = round(dist))
}

#' getKNN
#'
#' Calculates Euclidean distance in PCA space.
#'
#' @name getKNN
#' @keywords internal
#' @rdname getKNN
#' @author Jason T. Serviss
#' @param pca Matrix; Principal components as columns, samples as rows.
#' @param pcs Integer; Length 1 vector indicating the max principal component to
#' retain.
#' @param k Integer; Number of nearest neighbors.
#' @importFrom FNN get.knn
#' @importFrom tidyr gather
#' @importFrom dplyr mutate select distinct "%>%"
#' @importFrom stringr str_replace
#' @importFrom purrr map2_chr

getKNN <- function(pca, pcs, k) {
  near_data <- FNN::get.knn(pca[, 1:pcs], k = k)
  index <- near_data$nn.index
  rownames(index) <- rownames(pca)
  matrix_to_tibble(index, "from") %>%
    tidyr::gather(rank, index, -from) %>%
    dplyr::mutate(
      rank = as.numeric(stringr::str_replace(rank, "^.(.*)", "\\1"))
    ) %>%
    dplyr::mutate(to = rownames(pca)[index]) %>%
    dplyr::select(from, to, rank) %>%
    dplyr::mutate(cmb = purrr::map2_chr(from, to, function(x, y) {
      paste(sort(c(x, y)), collapse = "-")
    })) %>%
    dplyr::distinct(cmb, .keep_all = TRUE) %>%
    dplyr::select(-cmb)
}

#' constructGraph
#'
#' Constructs a graph from a data frame or tibble.
#'
#' @name constructGraph
#' @keywords internal
#' @rdname constructGraph
#' @author Jason T. Serviss
#' @param data Tibble; Expected columns are "from" and "to", indicating the
#' edges and "dist" indicating the distance between the samples.
#' @importFrom tidygraph as_tbl_graph activate .N
#' @importFrom igraph graph_from_data_frame
#' @importFrom dplyr mutate "%>%"

constructGraph <- function(data) {
  igraph::graph_from_data_frame(data, directed = FALSE) %>%
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(
      from.name = tidygraph::.N()$name[from],
      to.name = tidygraph::.N()$name[to]
    ) %>%
    tidygraph::activate(nodes)
}

#' louvain
#'
#' Runs the louvain community detection algorithm on the input graph.
#'
#' @name louvain
#' @keywords internal
#' @rdname louvain
#' @author Jason T. Serviss
#' @param g tbl_graph; A tidygraph graph.
#' @param weights.col Bare word; indicating the column name of the column
#' including the weights.
#' @importFrom tidygraph group_louvain
#' @importFrom readr parse_factor
#' @importFrom dplyr mutate enquo
#' @importFrom rlang "!!"

louvain <- function(g, weights.col) {
  dplyr::mutate(g,
    louvain = tidygraph::group_louvain(weight = !! enquo(weights.col)),
    louvain = readr::parse_factor(louvain, levels = unique(louvain))
  )
}

#' matrix_to_tibble
#'
#' Converts a matrix to a tibble without removing rownames.
#'
#' @name matrix_to_tibble
#' @rdname matrix_to_tibble
#' @author Jason T. Serviss
#' @param data matrix; The matrix to be converted.
#' @param rowname character; Length 1 vector indicating the colname that
#'  rownames should have upon tibble conversion.
#' @param drop logical; indicated if rownames should be dropped.
#'  Default = FALSE.
#' @keywords matrix_to_tibble
#' @examples
#'
#' m <- matrix(rnorm(20), ncol = 2, dimnames = list(letters[1:10], LETTERS[1:2]))
#' output <- matrix_to_tibble(m)
#'
#' @export
#' @importFrom tibble as_tibble rownames_to_column add_column
#' @importFrom rlang "!!" ":=" quo_name enquo

matrix_to_tibble <- function(data, rowname = "rowname", drop = FALSE) {
  if(!is.matrix(data)) stop("The 'data' argument is not a matrix")
  if(drop) return(as_tibble(data))
  rn.quo <- enquo(rowname)
  rn <- rownames(data)
  if(is.null(rn)) rn <- 1:nrow(data)

  rownames(data) <- NULL

  data %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tibble::add_column(!! quo_name(rn.quo) := rn, .before = 1) %>%
    as_tibble()
}

