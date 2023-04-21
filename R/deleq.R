# DELEQ functions
crossover_deleq <- function(mat, newmat, cr) {
  return((1 - cr) * mat + cr * newmat)
}

isapprox <- function(x, y, atol = 0, rtol = sqrt(.Machine$double.eps)) {
  return(x == y | (is.finite(x) & is.finite(y) & abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)))))
}

deleq <- function(fun, ..., NP, boxbounds, Emat, constr, x0 = NULL, cr = 0.7, f.param = 0.9, maxgen = 100, show.progress = TRUE) {
  if(!is.matrix(boxbounds) | ncol(boxbounds) != 2) {
    stop("'boxbounds' must be a two-column matrix.")
  }
  if(!is.matrix(Emat)) {
    stop("'Emat' must be a matrix.")
  }
  if(ncol(Emat) != nrow(boxbounds)) {
    stop("The number of rows of 'boxbounds' must be the same as the number of columns of 'Emat'.")
  }
  if(!is.matrix(constr)) {
    stop("'constr' must be a matrix.")
  }
  if(nrow(constr) != nrow(Emat)) {
    stop("Number of rows of 'constr' must be the same as that of 'Emat'.")
  }
  if(maxgen < 1 | maxgen %% 1 != 0) {
    stop("'maxgen' must be an integer greater than 0.")
  }
  if(cr < 0 | cr > 1) {
    stop("'cr' must be between 0 and 1.")
  }
  if(f.param < 0 | f.param > 1) {
    stop("'f.param' must be between 0 and 1.")
  }
  if(qr(Emat %*% t(Emat))$rank < nrow(Emat) & is.null(x0)) {
    stop("'Emat' * transpose('Emat') is singular and you have not provided an 'x0'.")
  }
  if(NP < nrow(boxbounds) * 10) {
    warning("'NP' is recommended to be at least 10 times the number of parameters.")
  }
  if(!is.null(x0)) {
    if(!is.matrix(x0) | ncol(x0) > 1) {
      stop("'x0' must be a one-column matrix.")
    }
    if(any(abs(Emat %*% x0 - constr) > sqrt(.Machine$double.eps))) {
      stop("'x0' does not satisfy the equality constraint(s).")
    }
  }
  gen <- 1
  if(is.null(x0)) {
    mat <- gen_init_pop(NP = NP, boxbounds = boxbounds, Emat = Emat, constr = constr)
  } else {
    mat <- gen_init_pop_x0(NP = NP, boxbounds = boxbounds, Emat = Emat, constr = constr, x0 = x0)
  }
  funvals <- apply(X = mat, MARGIN = 1, FUN = fun, ...)
  rbest <- which.max(funvals)
  trugen <- maxgen
  while(gen <= trugen) {
    if(show.progress) {
      print(paste0('Generation ', gen))
    }
    pbest <- rbest
    pbest_ind <- mat[pbest, ]
    pbest_val <- fun(pbest_ind, ...)
    newmat <- mutate_deleq(mat = mat, boxbounds = boxbounds, fParam = f.param)
    newmat <- crossover_deleq(mat = mat, newmat = newmat, cr = cr)
    funvals1 <- apply(X = newmat, MARGIN = 1, FUN = fun, ...)
    mat[funvals1 > funvals, ] <- newmat[funvals1 > funvals, ]
    funvals <- data.table::fifelse(funvals1 > funvals, funvals1, funvals)
    rbest <- which.max(funvals)
    if(!any(isapprox(Emat %*% matrix(mat[rbest, ], ncol = 1), constr))) {
      if(show.progress) {
        print('* Unfeasible, projecting back.')
      }
      trugen <- trugen - 1
      if(gen <= trugen) {
        projmat <- project_population(mat = mat, Emat = Emat, constr = constr)
        funvals1 <- apply(X = projmat, MARGIN = 1, FUN = fun, ...)
        rbest <- which.max(funvals1)
        if(pbest_val > fun(projmat[rbest, ], ...)) {
          projmat[rbest, ] <- pbest_ind
          rbest <- pbest
        }
        mat <- projmat
      } else {
        mat[rbest, ] <- pbest_ind
      }
    }
    gen <- gen + 1
    if(show.progress) {
      print(paste0("** Value = ", funvals[rbest]))
    }
  }
  return(list(params = mat[rbest, ], value = fun(mat[rbest, ], ...), generations = trugen))
}
