split_vec <- function(vector, start, end, length, exp = FALSE) {
  #' Split vector into list
  #'
  #' This function splits a vector from the index start to the index end into
  #' segments of length length and complies them into a list. This function can
  #' also apply the exponential function to the vector elements if desired.
  #'
  #' @param vector The vector to be split into a list.
  #' @param start The first index of the vector vector that will be included in
  #'   the list.
  #' @param end The last index of the vector vector that will be included in the
  #'   list.
  #' @param length A value indicating the length of the segments that vector
  #'   will be split into.
  #' @param exp A logical variable indicating whether the elements of vector
  #'   should be transformed with `exp()`.
  #'
  #' @return A list of the sections of the split vector.
  #' @export
  #' @examples
  #' split_vec(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, 6, 2)
  #' split_vec(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 2, 7, 3, exp = TRUE)

  if (exp) {
    return(split(exp(vector[start:end]),
                 ceiling(seq_along(exp(vector[start:end]))/length)))
  }
  split(vector[start:end], ceiling(seq_along(vector[start:end])/length))
}
