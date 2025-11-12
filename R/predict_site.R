#' @importFrom stats plogis setNames
NULL
#' Predict m6A site probability (base-resolution)
#'
#' @param seq A single DNA/RNA sequence as character string (length ≥ 41).
#' @param pos Integer position within the sequence to score (1-based).
#'
#' @return A single numeric probability.
#'
#' @export
#'
#' @examples
#' predict_site("AGAGGACTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG", pos = 21)
predict_site <- function(seq, pos) {
  # dummy: central 21-mer GC content
  k <- 21
  half <- (k %/% 2)
  start <- pos - half
  end   <- pos + half
  if (start < 1 || end > nchar(seq)) stop("Position too close to boundary.")
  sub <- substr(seq, start, end)
  gc <- sum(strsplit(toupper(sub), "")[[1]] %in% c("G", "C")) / k
  plogis(4 * (gc - 0.5))
}
