# functions pasted here as the fuzzy_join package is not installed in opensafely
# package repo: https://github.com/dgrtwo/fuzzyjoin

# utils functions
"%||%" <- function(x, y) if (is.null(x)) y else x


common_by <- function(by = NULL, x, y) {
  if (is.list(by)) return(by)
  
  if (!is.null(by)) {
    x <- names(by) %||% by
    y <- unname(by)
    
    # If x partially named, assume unnamed are the same in both tables
    x[x == ""] <- y[x == ""]
    
    return(list(x = x, y = y))
  }
  
  by <- intersect(dplyr::tbl_vars(x), dplyr::tbl_vars(y))
  if (length(by) == 0) {
    stop("No common variables. Please specify `by` param.", call. = FALSE)
  }
  message("Joining by: ", utils::capture.output(dput(by)))
  
  list(
    x = by,
    y = by
  )
}

# Make sure there's a distance column included in the output
ensure_distance_col <- function(ret, distance_col, mode) {
  if (!(mode %in% c("semi", "anti")) &&
      !is.null(distance_col) &&
      is.null(ret[[distance_col]])) {
    if (nrow(ret) == 0) {
      ret[[distance_col]] <- numeric(0)
    } else {
      ret[[distance_col]] <- NA
    }
  }
  ret
}

unrowwname <- function(x) {
  rownames(x) <- NULL
  x
}

# fuzzy_join function
fuzzy_join <- function (x, y, by = NULL, match_fun = NULL, multi_by = NULL, 
          multi_match_fun = NULL, index_match_fun = NULL, mode = "inner", 
          ...) 
{
  x_groups <- dplyr::groups(x)
  x <- dplyr::ungroup(x)
  regroup <- function(d) {
    if (length(x_groups) == 0) {
      return(d)
    }
    g <- purrr::map_chr(x_groups, as.character)
    missing <- !(g %in% colnames(d))
    g[missing] <- paste0(g[missing], ".x")
    dplyr::group_by_at(d, g)
  }
  mode <- match.arg(mode, c("inner", "left", "right", "full", 
                            "semi", "anti"))
  non_nulls <- (!is.null(multi_match_fun)) + (!is.null(match_fun)) + 
    (!is.null(index_match_fun))
  if (sum(non_nulls) != 1) {
    stop("Must give exactly one of match_fun, multi_match_fun, and index_match_fun")
  }
  if (!is.null(match_fun)) {
    by <- common_by(by, x, y)
    if (is.list(match_fun)) {
      match_fun <- purrr::map(match_fun, purrr::as_mapper)
    }
    else {
      match_fun <- purrr::as_mapper(match_fun)
    }
    if (length(match_fun) == 1) {
      match_fun <- rep(c(match_fun), length(by$x))
    }
    if (length(match_fun) != length(by$x)) {
      stop("Length of match_fun not equal to columns specified in 'by'.", 
           call. = FALSE)
    }
    matches <- dplyr::bind_rows(lapply(seq_along(by$x), function(i) {
      col_x <- x[[by$x[i]]]
      col_y <- y[[by$y[i]]]
      indices_x <- tibble::tibble(col = col_x, indices = seq_along(col_x)) %>% 
        dplyr::group_by(col) %>% tidyr::nest() %>% dplyr::mutate(indices = purrr::map(data, 
                                                                                      "indices"))
      indices_y <- tibble::tibble(col = col_y, indices = seq_along(col_y)) %>% 
        dplyr::group_by(col) %>% tidyr::nest() %>% dplyr::mutate(indices = purrr::map(data, 
                                                                                      "indices"))
      u_x <- indices_x$col
      u_y <- indices_y$col
      if (!is.null(names(match_fun))) {
        mf <- match_fun[[by$x[[i]]]]
      }
      else {
        mf <- match_fun[[i]]
      }
      extra_cols <- NULL
      n_x <- length(u_x)
      n_y <- length(u_y)
      m <- mf(rep(u_x, n_y), rep(u_y, each = n_x), ...)
      if (is.data.frame(m)) {
        if (ncol(m) > 1) {
          extra_cols <- m[, -1, drop = FALSE]
        }
        m <- m[[1]]
      }
      w <- which(m) - 1
      if (length(w) == 0) {
        ret <- tibble::tibble(i = numeric(0), x = numeric(0), 
                              y = numeric(0))
        return(ret)
      }
      x_indices_l <- indices_x$indices[w%%n_x + 1]
      y_indices_l <- indices_y$indices[w%/%n_x + 1]
      xls <- sapply(x_indices_l, length)
      yls <- sapply(y_indices_l, length)
      x_rep <- unlist(purrr::map2(x_indices_l, yls, function(x, 
                                                             y) rep(x, each = y)))
      y_rep <- unlist(purrr::map2(y_indices_l, xls, function(y, 
                                                             x) rep(y, x)))
      ret <- tibble::tibble(i = i, x = x_rep, y = y_rep)
      if (!is.null(extra_cols)) {
        extra_indices <- rep(w, xls * yls)
        extra_cols_rep <- extra_cols[extra_indices + 
                                       1, , drop = FALSE]
        ret <- dplyr::bind_cols(ret, extra_cols_rep)
      }
      ret
    }))
    if (length(by$x) > 1) {
      accept <- matches %>% dplyr::count(x, y) %>% dplyr::ungroup() %>% 
        dplyr::filter(n == length(by$x))
      matches <- matches %>% dplyr::semi_join(accept, by = c("x", 
                                                             "y"))
      if (ncol(matches) > 3) {
        matches <- matches %>% dplyr::semi_join(accept, 
                                                by = c("x", "y")) %>% dplyr::mutate(name = by$x[i]) %>% 
          dplyr::select(-i) %>% tidyr::gather(key, value, 
                                              -x, -y, -name) %>% tidyr::unite(newname, name, 
                                                                              key, sep = ".") %>% tidyr::spread(newname, 
                                                                                                                value)
      }
      else {
        matches <- dplyr::distinct(matches, x, y)
      }
    }
  }
  else if (!is.null(multi_match_fun)) {
    multi_match_fun <- purrr::as_mapper(multi_match_fun)
    by <- common_by(multi_by, x, y)
    number_x_rows <- nrow(x)
    number_y_rows <- nrow(y)
    indices_x <- x %>% dplyr::select_at(by$x) %>% dplyr::mutate(indices = seq_len(number_x_rows)) %>% 
      dplyr::group_by_at(dplyr::vars(-dplyr::one_of("indices"))) %>% 
      tidyr::nest() %>% dplyr::mutate(indices = purrr::map(data, 
                                                           "indices"))
    indices_y <- y %>% dplyr::select_at(by$y) %>% dplyr::mutate(indices = seq_len(number_y_rows)) %>% 
      dplyr::group_by_at(dplyr::vars(-dplyr::one_of("indices"))) %>% 
      tidyr::nest() %>% dplyr::mutate(indices = purrr::map(data, 
                                                           "indices"))
    ux <- as.matrix(indices_x[by$x])
    uy <- as.matrix(indices_y[by$y])
    pairs <- matrix(NA, nrow(ux), nrow(uy))
    ix <- row(pairs)
    iy <- col(pairs)
    ux_input <- ux[ix, ]
    uy_input <- uy[iy, ]
    m <- multi_match_fun(ux_input, uy_input)
    extra_cols <- NULL
    if (is.data.frame(m)) {
      if (ncol(m) > 1) {
        extra_cols <- m[, -1, drop = FALSE]
      }
      m <- m[[1]]
    }
    if (sum(m) == 0) {
      matches <- tibble::tibble(x = numeric(0), y = numeric(0))
    }
    else {
      x_indices_l <- indices_x$indices[ix[m]]
      y_indices_l <- indices_y$indices[iy[m]]
      xls <- purrr::map_dbl(x_indices_l, length)
      yls <- purrr::map_dbl(y_indices_l, length)
      x_rep <- unlist(purrr::map2(x_indices_l, yls, function(x, 
                                                             y) rep(x, each = y)))
      y_rep <- unlist(purrr::map2(y_indices_l, xls, function(y, 
                                                             x) rep(y, x)))
      matches <- tibble::tibble(x = x_rep, y = y_rep)
      if (!is.null(extra_cols)) {
        extra_indices <- rep(which(m), xls * yls)
        extra_cols_rep <- extra_cols[extra_indices, , 
                                     drop = FALSE]
        matches <- dplyr::bind_cols(matches, extra_cols_rep)
      }
    }
  }
  else {
    index_match_fun <- purrr::as_mapper(index_match_fun)
    by <- common_by(multi_by, x, y)
    d1 <- x[, by$x, drop = FALSE]
    d2 <- y[, by$y, drop = FALSE]
    matches <- index_match_fun(d1, d2)
  }
  matches$i <- NULL
  if (mode == "semi") {
    return(regroup(x[sort(unique(matches$x)), , drop = FALSE]))
  }
  if (mode == "anti") {
    if (nrow(matches) == 0) {
      return(regroup(x))
    }
    return(regroup(x[-sort(unique(matches$x)), , drop = FALSE]))
  }
  matches <- dplyr::arrange(matches, x, y)
  n <- intersect(colnames(x), colnames(y))
  x <- dplyr::rename_at(x, .vars = n, ~paste0(.x, ".x"))
  y <- dplyr::rename_at(y, .vars = n, ~paste0(.x, ".y"))
  if (mode == "left") {
    matches <- tibble::tibble(x = seq_len(nrow(x))) %>% dplyr::left_join(matches, 
                                                                         by = "x")
  }
  else if (mode == "right") {
    matches <- tibble::tibble(y = seq_len(nrow(y))) %>% dplyr::left_join(matches, 
                                                                         by = "y")
  }
  else if (mode == "full") {
    matches <- matches %>% dplyr::full_join(tibble::tibble(x = seq_len(nrow(x))), 
                                            by = "x") %>% dplyr::full_join(tibble::tibble(y = seq_len(nrow(y))), 
                                                                           by = "y")
  }
  ret <- dplyr::bind_cols(unrowwname(x[matches$x, , drop = FALSE]), 
                          unrowwname(y[matches$y, , drop = FALSE]))
  if (ncol(matches) > 2) {
    extra_cols <- unrowwname(matches[, -(1:2), drop = FALSE])
    ret <- dplyr::bind_cols(ret, extra_cols)
  }
  ret <- regroup(ret)
  if (!inherits(x, "tbl_df")) {
    ret <- as.data.frame(ret)
  }
  ret
}