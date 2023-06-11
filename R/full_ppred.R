#' Full forward sampling of the PPC
#'
#' @param fit An object of class `brmsfit`
#' @param newdata An optional data.frame for which to evaluate predictions. If NULL (default), the original data of the model is used.
#' @param draws An integer vector specifying the posterior draws to be used. If NULL (the default), all draws are used.
#'
#' @return A list of data.frames containing the draws.
#' @export
#'
#' @examples # Pending
full_ppred <- function(fit, newdata = NULL, draws = NULL) {
  # 1. determine term hierarchy
  resp <- response_sequence(fit)
  # 2.1. initialize dataframe using the original fit's data
  if(is.null(newdata)) newdata <- fit$data
  n <- nrow(newdata)
  # 2.3. if no draws set, range from 1 to all iters (check draws < iters)
  if(is.null(draws)) draws <- seq_len(sum(fit$fit@sim$n_save))
  # 2.4. create list to hold data
  pp_data <- list()

  for (i in draws) {
    pp_data[[i]] <- newdata
    for (vars in resp) {
      pp_data[[i]][, vars] <- array(
        brms::posterior_predict(
          fit, newdata = pp_data[[i]],
          resp = vars, draw_ids = i),
        dim = c(1, n, length(vars)))[1,,]
    }
  }
  pp_data
}

nodes_by_depth <- function(adj_matrix) {
  depth_list <- list()
  var_names <- rownames(adj_matrix)
  while(nrow(adj_matrix)) {
    pos <- which(apply(adj_matrix, 1, sum) == 0)
    depth_list <- c(depth_list, list(var_names[pos]))

    var_names <- var_names[-pos]
    adj_matrix <- adj_matrix[-pos, -pos, drop = FALSE]
  }
  depth_list
}

response_sequence <- function(x){ UseMethod("response_sequence") }
#' @method response_sequence brmsfit
#' @export
response_sequence.brmsfit <- function(x){ response_sequence(x$formula) }
#' @method response_sequence bform
#' @export
response_sequence.bform <- function(x){
  term_list <- response_sequence(brmsterms(x))
  resp_vars <- names(term_list)

  adjacency <- t(sapply(term_list, \(x)is.element(resp_vars, x)))
  attr(adjacency, "dimnames") <- list(resp_vars, resp_vars)
  nodes_by_depth(adjacency)
}
#' @method response_sequence mvbrmsterms
#' @export
response_sequence.mvbrmsterms <- function(x){
  names(x$terms) <- NULL
  sapply(x$terms, response_sequence)
}
#' @method response_sequence brmsterms
#' @export
response_sequence.brmsterms <- function(x){
  vars <- list(unique(unlist(lapply(x$dpars, response_sequence))))
  names(vars) <- all.vars(x$respform)
  vars
}
#' @method response_sequence btl
#' @export
response_sequence.btl <- function(x){ c("1", all.vars(x$formula)) }
