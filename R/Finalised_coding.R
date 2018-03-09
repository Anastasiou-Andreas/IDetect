options(expressions = 5e+05)
#' Derives a subset of integers from a given set
#'
#' This function finds two subsets of integers in a given interval \code{[s,e]}.
#' The routine is typically not called directly by the user; its result
#' is used in order to construct the expanding intervals, where the Isolate-Detect method
#' is going to be applied. For more details on how the Isolate-Detect methodology works, see
#' References.
#' @export
#' @param r A positive integer vector containing the set, from which the end-points
#'   of the expanding intervals are to be chosen.
#' @param l A positive integer vector containing the set, from which the start-points
#'   of the expanding intervals are to be chosen.
#' @param s A positive integer indicating the starting position, in the sense that we will
#'   choose the elements from \code{r} and \code{l} that are greater than \code{s}.
#' @param e A positive integer indicating the finishing position, in the sense that we will
#'   choose the elements from \code{r} and \code{l} that are less than \code{e}.
#' @return
#'   \code{e_points}  A vector containing the points that will be used as end-points,
#'   in order to create the left-expanding intervals. It consists of the input \code{e} and
#'   all the elements in the input vector \code{r} that are in \code{(s,e)}.
#'
#'   \code{s_points}  A vector containing the points that will be used as start-points,
#'   in order to create the left-expanding intervals. It consists of the input \code{s} and
#'   all the elements in the input vector \code{l} that are in \code{(s,e)}
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @references Anastasiou, A. and Fryzlewicz, P. (2018). Detecting multiple generalized change-points
#' by isolating single ones.
#' @examples
#' s_e_points(r = seq(10,1000,10), l = seq(991,1,-10), s=435, e = 786)
#' s_e_points(r = seq(3,100,3), l = seq(98,1,-3), s=43, e = 86)
s_e_points <- function(r, l, s, e) {
  r <- sort(r)
  l <- sort(l, decreasing = TRUE)
  if (s > e){
  stop("s should be less than or equal to e")
  }
  if (!(is.numeric(c(r, l, s, e))) | (r[1] <= 0) | (l[length(l)] <= 0) | s <= 0 | e <= 0){
    stop("The input arguments must be positive integers")
  }
  if (any(abs(r - round(r)) > .Machine$double.eps ^ 0.5)){
    warning("The input for r should be a vector of positive integers. If there is at least a positive real
            number then the integer part of that number is used.")
  }
  if (any(abs(l - round(l)) > .Machine$double.eps ^ 0.5)){
    warning("The input for l should be a vector of positive integers. If there is at least a positive real
            number then the integer part of that number is used.")
  }
  if (abs(s - round(s)) > .Machine$double.eps ^ 0.5){
    warning("The input for s should be a positive integer. If it is a positive real
            number then the integer part of that number is used.")
  }
  if (abs(e - round(e)) > .Machine$double.eps ^ 0.5){
    warning("The input for e should be a positive integer. If it is a positive real
            number then the integer part of that number is used.")
  }
  r <- as.integer(r)
  l <- as.integer(l)
  e <- as.integer(e)
  s <- as.integer(s)
  e_points <- unique(c(r[which( (r > s) & (r < e))], e))
  s_points <- unique(c(l[which( (l > s) & (l < e))], s))
  return(list(e_points = e_points, s_points = s_points))
}

cusum_function <- function(x) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data
         for which the CUSUM function will be calculated.")
  }
  n <- length(x)
  y <- cumsum(x)
  res <- sqrt( ( (n - 1):1) / n / (1:(n - 1))) * y[1:(n - 1)] - sqrt( (1:(n - 1)) / n / ( (n - 1):1)) * (y[n] - y[1:(n - 1)])
  return(res)
}

cumsum_lin <- function(x) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector.")
  }
  res <- numeric()
  n <- length(x)
  if (n <= 2) {
    res <- 0
  } else {
    b <- 2:(n - 1)
    y1 <- cumsum(x * (1:n))
    y <- cumsum(x)
    a <- sqrt(6 / ( (n - 1) * n * (n + 1) * (2 - 2 * b ^ 2 + 2 * b * n - 1 + 2 * b - n)))
    be <- sqrt( ( (n - b + 1) * (n - b)) / ( (b - 1) * b))
    res[1] <- 0
    res[b] <- a * be * ( (2 * b + n - 1) * y1[b] - (n + 1) * b * y[b]) - (a / be) * ( ( 3 * n - 2 * b + 1) * (y1[n] - y1[b]) - (n + 1) * (2 * n - b) * (y[n] - y[b]))
  }
  return(res)
}

#' Multiple change-point detection in the mean via thresholding
#'
#' This function performs the Isolate-Detect methodology (see Details for the
#' relevant literature reference) with the thresholding-based stopping rule
#' in order to detect multiple change-points in the mean of a noisy input vector
#' \code{x}, with Gaussian noise. See Details for a brief explanation of the
#' Isolate-Detect methodology, and of the thresholding-based stopping rule.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param sigma A positive real number. It is the estimate of the standard deviation
#'   of the noise in \code{x}. The default value is the median absolute deviation of \code{x}
#'   computed under the assumption that the noise is independent and identically distributed
#'   from the Gaussian distribution.
#' @param thr_const A positive real number with default value equal to 1. It is
#'   used to define the threshold; see \code{thr_fin}.
#' @param thr_fin With \code{T} the length of the data sequence, this is a positive real
#'   number with default value equal to \code{sigma * thr_const * sqrt(2 * log(T))}. It is
#'   the threshold, which is used in the detection process.
#' @param s,e Positive integers with \code{s} less than \code{e}, which indicate
#'   that you want to check for change-points in the data sequence with subscripts
#'   in \code{[s,e]}. The default values are \code{s} equal to 1 and
#'   \code{e} equal to \code{T}, with \code{T} the length of the data sequence.
#' @param points A positive integer with default value equal to 3. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively; see Details for more information.
#' @param k_l,k_r Positive integer numbers that get updated whenever the function
#'   calls itself during the detection process. They are not essential for the
#'   function to work, and we include them only to reduce the computational time.
#' @details The change-point detection algorithm that is used in \code{\link{pcm_th}} is the
#'   Isolate-Detect methodology described in ``Detecting multiple generalized
#'   change-points by isolating single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#'   The concept is simple and is split into two stages; firstly, isolation of each
#'   of the true change-points in subintervals of the data domain, and secondly their detection.
#'   ID first creates two ordered sets of \eqn{K = \lceil T/\code{points}\rceil} right- and left-expanding
#'   intervals as follows. The \eqn{j^{th}} right-expanding interval is \eqn{R_j = [1, j\times \code{points}]},
#'   while the \eqn{j^{th}} left-expanding interval is \eqn{L_j = [T - j\times \code{points} + 1, T]}.
#'   We collect these intervals in the ordered set \eqn{S_{RL} = \lbrace R_1, L_1, R_2, L_2, ... , R_K, L_K\rbrace}.
#'   For a suitably chosen contrast function, ID first identifies the point with the maximum contrast
#'   value in \eqn{R_1}. If its value exceeds a certain threshold, then it is taken as a change-point.
#'   If not, then the process tests the next interval in \eqn{S_{RL}} and repeats the above process.
#'   Upon detection, the algorithm makes a new start from estimated location.
#' @return
#'   A numeric vector with the detected change-points.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{win_pcm_th}}, \code{\link{ID_pcm}}, and \code{\link{ID}}, which employ
#' this function. In addition, see \code{\link{cplm_th}} for the case of detecting changes in
#' a continuous, piecewise-linear signal via thresholding.
#' @examples
#' single.cpt <- c(rep(4,1000),rep(0,1000))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single.th <- pcm_th(single.cpt.noise)
#'
#' three.cpt <- c(rep(4,500),rep(0,500),rep(-4,500),rep(1,500))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' cpt.three.th <- pcm_th(three.cpt.noise)
#'
#' multi.cpt <- rep(c(rep(0,50),rep(3,50)),20)
#' multi.cpt.noise <- multi.cpt + rnorm(2000)
#' cpt.multi.th <- pcm_th(multi.cpt.noise)
pcm_th <- function(x, sigma = stats::mad(diff(x) / sqrt(2)), thr_const = 1,
                   thr_fin = sigma * thr_const * sqrt(2 * log(length(x))),
                   s = 1, e = length(x), points = 3, k_l = 1, k_r = 1) {
  if (!(is.numeric(x))){
      stop("The input in `x' should be a numeric vector containing the data in
           which you would like to find change-points.")
    }
  if ( (thr_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  points <- as.integer(points)
  l <- length(x)
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s <= 1) {
      cpt <- 0
  } else {
      pos_r <- numeric()
      CUSUM_r <- numeric()
      pos_l <- numeric()
      CUSUM_l <- numeric()
      moving_points <- s_e_points(r_e_points, l_e_points, s, e)
      right_points <- moving_points[[1]]
      left_points <- moving_points[[2]]
      lur <- length(left_points)
      rur <- length(right_points)
      if (k_r < k_l) {
          while ( (chp == 0) & (k_r < min(k_l, rur))) {
              x_temp_r <- x[s:right_points[k_r]]
              ipcr <- cusum_function(x_temp_r)
              pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
              CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
              if (CUSUM_r[k_r] > thr_fin) {
                chp <- pos_r[k_r]
              } else {
                k_r <- k_r + 1
              }
          }
      }
      if (k_l < k_r) {
          while ( (chp == 0) & (k_l < min(k_r, lur))) {
              x_temp_l <- x[left_points[k_l]:e]
              ipcl <- cusum_function(x_temp_l)
              pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
              CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
              if (CUSUM_l[k_l] > thr_fin) {
                chp <- pos_l[k_l]
              } else {
                k_l <- k_l + 1
              }
          }
      }
      if (chp == 0) {
          while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
              x_temp_r <- x[s:right_points[k_r]]
              ipcr <- cusum_function(x_temp_r)
              pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
              CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
              if (CUSUM_r[k_r] > thr_fin) {
                chp <- pos_r[k_r]
              } else {
                x_temp_l <- x[left_points[k_l]:e]
                ipcl <- cusum_function(x_temp_l)
                pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
                CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
                if (CUSUM_l[k_l] > thr_fin) {
                  chp <- pos_l[k_l]
                } else {
                  k_r <- k_r + 1
                  k_l <- k_l + 1
                }
              }
          }
      }
      if (chp != 0) {
          if (chp > ( (e + s) / 2)) {
              r <- pcm_th(x, s = s, e = chp, points = points,
                          thr_fin = thr_fin, k_r = k_r, k_l = 1)
          } else {
              r <- pcm_th(x, s = chp + 1, e = e, points = points,
                          thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
          }
          cpt <- c(chp, r)
      } else {
          cpt <- chp
      }
  }
  cpt <- cpt[cpt != 0]
  return(sort(cpt))
}

#' A windows-based approach for multiple change-point detection in the mean via
#' thresholding
#'
#' This function performs the windows-based variant of the Isolate-Detect methodology
#' with the thresholding-based stopping rule in order to detect multiple change-points
#' in the mean of a noisy data sequence, with noise that is Gaussian. It is particularly
#' helpful for very long data sequences, as due to applying Isolate-Detect on moving windows,
#' the computational time is reduced. See Details for a brief explanation of this approach and
#' for the relevant literature reference.
#'
#' @export
#' @param xd A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param sigma A positive real number. It is the estimate of the standard deviation
#'   of the noise in \code{xd}. The default value is the median absolute deviation of \code{xd}
#'   computed under the assumption that the noise is independent and identically distributed
#'   from the Gaussian distribution.
#' @param thr_con A positive real number with default value equal to 1. It is
#'   used to define the threshold, which is equal to \code{sigma * thr_con * sqrt(2 * log(T))},
#'   where \code{T} is the length of the data sequence \code{xd}.
#' @param c_win A positive integer with default value equal to 3000. It is the length
#'   of each window for the data sequence in hand. Isolate-Detect will be applied
#'   in segments of the form \code{[(i-1) * c_win + 1, i * c_win]}, for \eqn{i=1,2,...,K},
#'   where \eqn{K} depends on the length \code{T} of the data sequence.
#' @param w_points A positive integer with default value equal to 3. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @param l_win A positive integer with default value equal to 12000. If the length of
#'   the data sequence is less than or equal to \code{l_win}, then the windows-based approach
#'   will not be applied and the result will be obtained by the classical Isolate-Detect
#'   methodology based on thresholding.
#' @details The method that is implemented by this function is based on splitting the given
#'   data sequence uniformly into smaller parts (windows), to which Isolate-Detect, based on the
#'   threshold stopping rule (see \code{\link{pcm_th}}), is then applied. An idea of the computational
#'   improvement that this structure offers over the classical Isolate-Detect in the case of large data
#'   sequences is given in the supplement of ``Detecting multiple generalized change-points by isolating
#'   single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A numeric vector with the detected change-points.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{pcm_th}}, which is the function that \code{\link{win_pcm_th}} is based on. Also,
#' see \code{\link{ID_pcm}} and \code{\link{ID}}, which employ \code{\link{win_pcm_th}}. In addition,
#' see \code{\link{win_cplm_th}} for the case of detecting changes in the slope of a
#' piecewise-linear and continuous signal via thresholding.
#' @examples
#' single.cpt <- c(rep(4,1000),rep(0,1000))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single.th <- win_pcm_th(single.cpt.noise)
#'
#' three.cpt <- c(rep(4,4000),rep(0,4000),rep(-4,4000),rep(1,4000))
#' three.cpt.noise <- three.cpt + rnorm(16000)
#' cpt.three.th <- win_pcm_th(three.cpt.noise)
win_pcm_th <- function(xd, sigma = stats::mad(diff(xd) / sqrt(2)), thr_con = 1,
                        c_win = 3000, w_points = 3, l_win = 12000) {
  if (!(is.numeric(xd))){
    stop("The input in `xd' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_con <= 0) || (w_points <= 0) || (c_win <= 0) || (l_win <= 0)){
    stop("The threshold constant as well as the `w_points', `c_win', `l_win' arguments should
         be positive numbers.")
  }
  if ( (abs(w_points - round(w_points)) > .Machine$double.eps ^ 0.5)
      || (abs(c_win - round(c_win)) > .Machine$double.eps ^ 0.5)
      || (abs(l_win - round(l_win)) > .Machine$double.eps ^ 0.5)){
    warning("The input values  for `w_points', `c_win', and  `l_win' should be positive integers.
            If either of them is a positive real number then the integer part of the given number
            is used to obtain the result.")
  }
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  l_win <- as.integer(l_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= l_win) {
      u <- pcm_th(x = xd, thr_const = thr_con, points = w_points)
      return(u)
  } else {
      K <- ceiling(lg / c_win)
      tsm <- list()
      u <- list()
      ufin <- numeric()
      uaddition <- list()
      tsm[[1]] <- xd[1:c_win]
      ufin <- pcm_th(tsm[[1]], thr_fin = t, points = w_points)
      uaddition[[1]] <- numeric()
      uaddition[[1]] <- pcm_th(x = xd[(max(1, c_win - (5 * w_points) + 1)):min( (c_win + (5 * w_points)), lg)], thr_fin = t, points = 2) + c_win - (5 * w_points)
      i <- 2
      while (i < K) {
          tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
          u[[i]] <- pcm_th(x = tsm[[i]], thr_fin = t, points = w_points) + (i - 1) * c_win
          uaddition[[i]] <- numeric()
          uaddition[[i]] <- pcm_th(x = xd[(max(1, i * c_win - (5 * w_points) + 1)):(min(i * c_win + (5 * w_points), lg))], thr_fin = t, points = 2) + i * c_win - (5 * w_points)
          ufin <- c(ufin, u[[i]], uaddition[[i]])
          i <- i + 1
      }
      tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
      u[[K]] <- pcm_th(tsm[[K]], thr_fin = t, points = w_points) + (K - 1) * c_win
      ufinl <- c(ufin, u[[K]])
      return(cpt = sort(unique(ufinl)))
  }
}

cusum_one <- function(x, s, e, b) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector.")
  }
  y <- cumsum(x)
  l <- numeric()
  d <- numeric()
  result <- numeric()
  if ( (length(s) != length(b)) || (length(s) != length(e)) || (length(e) != length(b))){
    stop("The vectors s, b, e, should be of the same length")
  }
  if (any(s < 1) | any(b < 1) | any(e < 1)){
    stop("The entries of the vectors s, b, e should be positive integers.")
  }
  if (any(s > b) | any(b >= e)){
    stop("The value for b should be in the interval [s,e)")
  }
  if ( (any(abs( (s - round(s))) > .Machine$double.eps ^ 0.5))
      || (any(abs( (b - round(b))) > .Machine$double.eps ^ 0.5))
      || (any(abs( (e - round(e))) > .Machine$double.eps ^ 0.5))){
    stop("The input values  for s, b, and  e should be positive integers.")
  }
  for (j in 1:length(b)) {
      l[j] <- e[j] - s[j] + 1
      d[j] <- ifelse(s[j] == 1, 0, y[s[j] - 1])
      result[j] <- abs(sqrt( (e[j] - b[j]) / (l[j] * (b[j] - s[j] + 1))) * (y[b[j]] - d[j]) - sqrt( (b[j] - s[j] + 1) / (l[j] * (e[j] - b[j]))) * (y[e[j]] - y[b[j]]))
  }
  return(result)
}

#' The solution path for the case of piecewise-constant signals
#'
#' This function starts by overestimating the number of true change-points.
#' After that, following a CUSUM-based approach, it sorts the estimated change-points
#' in a way that the estimate, which is most-likely to be correct appears first, whereas
#' the least likely to be correct, appears last. The routine is typically not called
#' directly by the user; it is employed in \code{\link{pcm_ic}}. For more information, see
#' References.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param thr_ic A positive real number with default value equal to 0.9. It is
#'   used to define the threshold. The change-points are estimated by thresholding
#'   with threshold equal to \code{sigma * thr_ic * sqrt(2 * log(T))}, where
#'   \code{T} is the length of the data sequence \code{x} and \code{sigma = mad(diff(x)/sqrt(2))}.
#'   Because we would like to overestimate the number of true change-points in \code{x}, it is
#'   suggested to keep \code{thr_ic} smaller than 1, which is the default value used as
#'   the threshold constant in the function \code{\link{pcm_th}}.
#' @param points A positive integer with default value equal to 3. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @return
#'   The solution path for the case of piecewise-constant signals.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @references Anastasiou, A. and Fryzlewicz, P. (2018). Detecting multiple generalized change-points
#' by isolating single ones.
#' @examples
#' three.cpt <- c(rep(4,4000),rep(0,4000),rep(-4,4000),rep(1,4000))
#' three.cpt.noise <- three.cpt + rnorm(16000)
#' solution.path <- sol_path_pcm(three.cpt.noise)
sol_path_pcm <- function(x, thr_ic = 0.9, points = 3) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which the solution path will be calculated.")
  }
  if ( (thr_ic <= 0) || (points <= 0)){
      stop("The threshold constant as well as the `points' argument that represents the
           magnitude of the expansion for the intervals should be positive numbers.")
    }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
      warning("The input for `points' should be a positive integer. If it is a positive real
              number then the integer part of the given number is used as the value of `points'.")
    }
  lx_ic <- length(x)
  points <- as.integer(points)
  cpt_lower <- win_pcm_th(x, thr_con = thr_ic, w_points = points)
  lcpt_ic <- length(cpt_lower)
  if ( (lcpt_ic == 1) | (lcpt_ic == 0)) {
      return(cpt_lower)
  } else {
      seb_set <- c(unique(c(1, cpt_lower)), lx_ic)
      lseb_set <- length(seb_set)
      min_C <- numeric()
      while (lseb_set >= 3) {
          Rs <- cusum_one(x, seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)],
                          seb_set[2:(lseb_set - 1)])
          min_Rs <- seb_set[2:(lseb_set - 1)][which.min(Rs)]
          min_C <- c(min_C, min_Rs)
          seb_set <- seb_set[-which(seb_set == min_Rs)]
          lseb_set <- lseb_set - 1
      }
      return(min_C[length(min_C):1])
  }
}

sic_pen <- function(n, n_param) {
    pen <- log(n)
    return(n_param * pen)
}

ssic_pen <- function(n, n_param, alpha = 1.01) {
  alpha <- as.numeric(alpha)
  pen <- log(n) ^ alpha
  return(n_param * pen)
}

#' Multiple change-point detection in the mean via minimising an information
#' criterion
#'
#' This function performs the Isolate-Detect methodology based on an information
#' criterion approach, in order to detect multiple change-points in the mean of
#' a noisy data sequence, with the noise following the Gaussian distribution.
#' More information on how this approach works as well as the relevant literature
#' reference are given in Details.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param th_const A positive real number with default value equal to 0.9. It is
#'   used to define the threshold value that will be used at the first step of the
#'   model selection based Isolate-Detect method; see Details for more information.
#' @param Kmax A positive integer with default value equal to 200. It is the
#'   maximum allowed number of estimated change-points in the solution path algorithm,
#'   described in Details below.
#' @param penalty A character vector with names of the penalty functions used.
#' @param points A positive integer with default value equal to 10. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @details The approach followed in \code{\link{pcm_ic}} in order to detect
#'   the change-points is based on identifying the set of change-points that
#'   minimise an information criterion. At first, we employ \code{\link{sol_path_pcm}},
#'   which overestimates the number of change-points using \code{th_const} in order to define the
#'   threshold, and then sorts the obtained estimates in a way that the estimate, which
#'   is most likely to be correct appears first, whereas the least likely to
#'   be correct, appears last. Let \eqn{J} be the number of estimates
#'   that this overestimation approach returns. We will obtain a vector
#'   \eqn{b = (b_1, b_2, ..., b_J)}, with the estimates ordered as explained above. We define
#'   the collection \eqn{\left\{M_j\right\}_{j = 0,1,\ldots,J}}, where \eqn{M_0} is the empty set
#'   and \eqn{M_j = \left\{b_1,b_2,...,b_j\right\}}. Among the collection of models
#'   \eqn{M_j, j=0,1,...,J}, we select the one that minimises a predefined Information
#'   Criterion. The obtained set of change-points is apparently a subset of the solution path
#'   given in \code{\link{sol_path_pcm}}. More details can be found in
#'   ``Detecting multiple generalized change-points by isolating single ones'',
#'   Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{sol_path} \tab A vector containing the solution path.
#'    \cr \code{ic_curve}   \tab A list with values of the chosen information criteria.
#'    \cr \code{cpt_ic} \tab A list with the change-points detected for each information
#'   criterion considered.
#'    \cr \code{no_cpt_ic } \tab The number of change-points detected for each information
#'   criterion considered.
#'  }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{ID_pcm}} and \code{\link{ID}}, which employ this function.
#' In addition, see \code{\link{cplm_ic}} for the case of detecting changes in
#' a continuous, piecewise-linear signal using the information criterion based approach.
#' @examples
#' single.cpt <- c(rep(4,1000),rep(0,1000))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single.ic <- pcm_ic(single.cpt.noise)
#'
#' three.cpt <- c(rep(4,500),rep(0,500),rep(-4,500),rep(1,500))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' cpt.three.ic <- pcm_ic(three.cpt.noise)
pcm_ic <- function(x, th_const = 0.9, Kmax = 200,
                       penalty = c("ssic_pen", "sic_pen"), points = 10) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you want to look for change-points.")
  }
  if ( (th_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  result <- list()
  penalty <- as.character(penalty)
  lx <- length(x)
  if (length(penalty)) {
      if (Kmax == 0 || lx <= 3) {
          result$cpt_ic[c("ssic_pen", "sic_pen")] <- c(NA, NA)
          result$no_cpt_ic[c("ssic_pen", "sic_pen")] <- as.integer(c(0, 0))
          result$ic_curve <- list()
          result$ic_curve$ssic_pen <- result$ic_curve$ssic_pen <- NA
          if (Kmax == 0) {
              stop("No change-points found, choose larger Kmax")
          } else {
              stop("Sample size is too small")
          }
      } else {
          cpt_cand <- sol_path_pcm(x, thr_ic = th_const, points = points)
          if (length(cpt_cand) == 0){
            result$cpt_ic[c("ssic_pen", "sic_pen")] <- c(0, 0)
            result$no_cpt_ic[c("ssic_pen", "sic_pen")] <- as.integer(c(0, 0))
            result$ic_curve <- list()
            result$ic_curve$ssic_pen <- result$ic_curve$ssic_pen <- NA
          }
          if (length(cpt_cand) > min(Kmax, lx - 2)){
            cpt_cand <- cpt_cand[1:min(Kmax, lx - 2)]
          }
          len_cpt <- length(cpt_cand)
          result$sol_path <- cpt_cand
          result$ic_curve <- list()
          for (j in 1:length(penalty)) {
            result$ic_curve[[penalty[j]]] <- rep(0, len_cpt + 1)
          }
          if (len_cpt){
            for (i in len_cpt:1) {
              min_log_lik <- lx / 2 * log(sum( (x - est_signal(x, cpt_cand[1:i], type = "mean")) ^ 2) / lx)
              for (j in 1:length(penalty)) {
                result$ic_curve[[penalty[j]]][i + 1] <- min_log_lik + do.call(penalty[j], list(n = lx, n_param = length(cpt_cand[1:i])))
              }
            }
          }
          for (j in 1:length(penalty)) {
            result$ic_curve[[penalty[j]]][1] <- lx / 2 * log(stats::var(x))
          }
          result$cpt_ic <- list()
          for (j in 1:length(penalty)) {
            tmp <- stats::quantile(which.min(result$ic_curve[[penalty[j]]]), 0.5, type = 3)
            if (tmp == 1) {
              result$cpt_ic[[penalty[j]]] <- NA
              result$no_cpt_ic[penalty[j]] <- as.integer(0)
            } else {
              result$cpt_ic[[penalty[j]]] <- sort(cpt_cand[1:(tmp - 1)])
              result$no_cpt_ic[penalty[j]] <- as.integer(tmp - 1)
            }
          }
      }
  }
  return(result)
}

#' Estimate the signal
#'
#' This function estimates the signal in a given data sequence \code{x} with change-points
#' at \code{cpt}. The type of the signal depends on whether the change-points represent changes
#' in a piecewise-constant or continuous, piecewise-linear signal. For more information see
#' Details below.
#'
#' @export
#' @param x A numeric vector containing the given data.
#' @param cpt A positive integer vector with the locations of the change-points.
#'   If missing, the \code{\link{ID_pcm}} or the \code{\link{ID_cplm}} function
#'   (depending on the type of the signal) is called internally to extract the
#'   change-points in \code{x}.
#' @param type A character string, which defines the type of the detected change-points.
#'   If \code{type = ``mean''}, then the change-points represent the locations of changes
#'   in the mean of a piecewise-constant signal. If \code{type = ``slope''}, then the
#'   change-points represent the locations of changes in the slope of a continuous, piecewise-linear
#'   signal.
#' @details The data points provided in \code{x} are assumed to follow \deqn{X_t = f_t + \sigma\epsilon_t; t = 1,2,...,T,}
#'   where \eqn{T} is the total length of the data sequence, \eqn{X_t} are the observed
#'   data, \eqn{f_t} is a one-dimensional, deterministic signal with abrupt structural
#'   changes at certain points, and \eqn{\epsilon_t} is white noise. We denote by
#'   \eqn{r_1, r_2, ..., r_N} the elements in \code{cpt} and by \eqn{r_0 = 0} and
#'   \eqn{r_{N+1} = T}. Depending on the value that has been passed to \code{type}, the returned
#'   value is calculated as follows.
#'
#'   \itemize{
#'   \item {For \code{type = ``mean''}, in each segment \eqn{(r_j + 1, r_{j+1})}, \eqn{f_t} for
#'   \eqn{t \in (r_j + 1, r_{j+1})} is approximated by the mean of \eqn{X_t} calculated
#'   over \eqn{t \in (r_j + 1, r_{j+1})}.}
#'
#'   \item {For \code{type = ``slope''}, \eqn{f_t} is approximated by the linear spline fit with
#'   knots at \eqn{r_1, r_2, ..., r_N} minimising the \eqn{l_2} distance between the fit and the data.}
#'   }
#' @return
#'   A numeric vector with the estimated signal.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @examples
#' single.cpt.pcm <- c(rep(4,1000),rep(0,1000))
#' single.cpt.pcm.noise <- single.cpt.pcm + rnorm(2000)
#' cpt.single.pcm <- ID_pcm(single.cpt.pcm.noise)
#' fit.cpt.single.pcm <- est_signal(single.cpt.pcm.noise, cpt.single.pcm$cpt, type = "mean")
#'
#' three.cpt.pcm <- c(rep(4,500),rep(0,500),rep(-4,500),rep(1,500))
#' three.cpt.pcm.noise <- three.cpt.pcm + rnorm(2000)
#' cpt.three.pcm <- ID_pcm(three.cpt.pcm.noise)
#' fit.cpt.three.pcm <- est_signal(three.cpt.pcm.noise, cpt.three.pcm$pcm, type = "mean")
#'
#' single.cpt.plm <- c(seq(0,999,1),seq(998.5,499,-0.5))
#' single.cpt.plm.noise <- single.cpt.plm + rnorm(2000)
#' cpt.single.plm <- ID_cplm(single.cpt.plm.noise)
#' fit.cpt.single.plm <- est_signal(single.cpt.plm.noise, cpt.single.plm$cpt, type = "slope")
est_signal <- function(x, cpt, type = c("mean", "slope")) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which you want to estimate the underlying signal.")
  }
  x <- as.numeric(x)
  if (NA %in% x)
    stop("x vector cannot contain NA's")
  n <- length(x)
  if (type == "mean") {
    if (missing(cpt)){
      cpt <- ID_pcm(x)$cpt
    }
  if (!is.null(cpt)){
    if (any(is.na(cpt))){
      cpt <- cpt[!is.na(cpt)]
    }
  }
  cpt <- as.integer(cpt)
  len_cpt <- length(cpt)
  if (len_cpt) {
    if (min(cpt) < 0 || max(cpt) >= n)
      stop("change-points cannot be negative or greater than and n-1")
    cpt <- sort(cpt)
  }
   s <- e <- rep(0, len_cpt + 1)
   s[1] <- 1
   e[len_cpt + 1] <- n
   if (len_cpt) {
     s[2:(len_cpt + 1)] <- cpt + 1
     e[1:len_cpt] <- cpt
   }
   means <- rep(0, len_cpt + 1)
   for (i in 1:(len_cpt + 1)) {
     means[i] <- mean(x[s[i]:e[i]])
   }
   fit <- rep(means, e - s + 1)
  }
  if (type == "slope"){
    if (missing(cpt)){
      cpt <- ID_cplm(x)$cpt
    }
    if (!is.null(cpt)){
      if (any(is.na(cpt))){
        cpt <- cpt[!is.na(cpt)]
      }
    }
  cpt <- as.integer(cpt)
  len_cpt <- length(cpt)
  if (len_cpt) {
    if (min(cpt) < 0 || max(cpt) >= n)
      stop("change-points cannot be negative or greater than and n-1")
    cpt <- sort(cpt)
  }
   cpt <- sort(unique(c(cpt, 0, n)))
   fit <- rep(0, n)
   cpt <- setdiff(cpt, c(0, n))
   X <- splines::bs(1:n, knots = cpt, degree = 1, intercept = TRUE)
   fit <- stats::lm.fit(X, x)$fitted.values
   }
  return(fit)
}

#' Multiple change-point detection in the mean of a vector using the
#' Isolate-Detect methodology
#'
#' This function estimates the number and locations of multiple change-points in the mean
#' of the noisy piecewise-constant input vector \code{x}, using the Isolate-Detect methodology. The noise
#' is Gaussian. The estimated signal, as well as the solution path defined in \code{\link{sol_path_pcm}} are
#' also given. The function is a hybrid between the thresholding approach of \code{\link{win_pcm_th}} and the
#' information criterion approach of \code{\link{pcm_ic}} and estimates the change-points taking into
#' account both these approaches (see Details for more information and the relevant literature reference).
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param thr_id A positive real number with default value equal to 1. It is
#'   used to define the threshold, if the thresholding approach is to be followed; see \code{\link{pcm_th}}
#'   for more details.
#' @param th_ic_id  A positive real number with default value equal to 0.9. It is
#'   useful only if the model selection based Isolate-Detect method is to be followed.
#'   It is used to define the threshold value that will be used at the first step
#'   (change-point overestimation) of the model selection approach described in \code{\link{pcm_ic}}.
#' @param pointsth A positive integer with default value equal to 3. It is used only
#'   when the threshold based approach is to be followed and it defines the distance
#'   between two consecutive end- or start-points of the right- or left-expanding intervals,
#'   respectively.
#' @param pointsic A positive integer with default value equal to 10. It is used only
#'   when the information criterion based approach is to be followed and it defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @details Firstly, this function detects the change-points using \code{\link{win_pcm_th}}.
#'   If the estimated number of change-points is larger than 100, then the
#'   result is returned and we stop. Otherwise, \code{\link{ID_pcm}} proceeds to detect the
#'   change-points using \code{\link{pcm_ic}} and this is what is returned. To sum up,
#'   \code{\link{ID_pcm}} returns a result based on \code{\link{pcm_ic}} if the estimated number
#'   of change-points is less than 100. Otherwise, the result comes from thresholding.
#'   More details can be found in ``Detecting multiple generalized change-points by
#'   isolating single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'    \tabular{ll}{
#'    \cr \code{cpt} \tab A vector with the detected change-points.
#'    \cr \code{no_cpt}    \tab The number of change-points detected.
#'    \cr \code{fit} \tab A numeric vector with the estimated piecewise-constant signal.
#'    \cr \code{solution_path} \tab A vector containing the solution path.
#'   }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{win_pcm_th}} and \code{\link{pcm_ic}} which are the functions that \code{\link{ID_pcm}}
#' is based on. In addition, see \code{\link{ID_cplm}} for the case of detecting changes
#' in a continuous, piecewise-linear signal. The main function \code{\link{ID}}
#' of the package employs \code{\link{ID_pcm}}.
#' @examples
#' single.cpt <- c(rep(4,1000),rep(0,1000))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpts_detect <- ID_pcm(single.cpt.noise)
#'
#' three.cpt <- c(rep(4,500),rep(0,500),rep(-4,500),rep(1,500))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' cpts_detect_three <- ID_pcm(three.cpt.noise)
#'
#' multi.cpt <- rep(c(rep(0,50),rep(3,50)),20)
#' multi.cpt.noise <- multi.cpt + rnorm(2000)
#' cpts_detect_multi <- ID_pcm(multi.cpt.noise)
ID_pcm <- function(x, thr_id = 1, th_ic_id = 0.9, pointsth = 3, pointsic = 10) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
           which you would like to find change-points.")
  }
  if ( (thr_id <= 0) || (pointsth <= 0) || (th_ic_id <= 0) || (pointsic <= 0)){
    stop("The threshold constants as well as the `pointsth' and `pointsic' arguments
         should be positive numbers.")
  }
  if ( (abs(pointsth - round(pointsth)) > .Machine$double.eps ^ 0.5)
      || (abs(pointsic - round(pointsic)) > .Machine$double.eps ^ 0.5)){
    warning("The input values for `pointsth' and `pointsic' should be positive integers. If they are
            positive real numbers then the integer parts of the given numbers are used as the values
            of `pointsth' and `pointsic'.")
  }
  z <- win_pcm_th(x, thr_con  = thr_id, w_points = pointsth)
  if (length(z) >= 100) {
      cpt <- z
      nocpt <- length(z)
      mean <- est_signal(x, cpt, type = "mean")
      result <- list(cpt = sort(cpt), no_cpt = nocpt, fit = mean)
  } else {
      z <- pcm_ic(x, th_const = th_ic_id, points = pointsic)
      if (is.na(z$cpt_ic[[1]][1])) {
          cpt <- 0
          nocpt <- 0
      } else {
          cpt <- z$cpt_ic[[1]]
          nocpt <- z$no_cpt_ic[[1]]
      }
      mean <- est_signal(x, cpt, type = "mean")
      sol_path <- z$sol_path
      result <- list(cpt = sort(cpt), no_cpt = nocpt, fit = mean,
                     solution_path = sol_path)
  }
  return(result)
}

#' Multiple change-point detection in a continuous, piecewise-linear signal
#' via thresholding
#'
#' This function performs the Isolate-Detect methodology with the thresholding-based
#' stopping rule in order to detect multiple change-points in a continuous, piecewise-linear
#' noisy data sequence, with noise that is Gaussian. See Details for a brief explanation of the
#' Isolate-Detect methodology (with the relevant reference) and of the thresholding-based
#' stopping rule.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param sigma A positive real number. It is the estimate of the standard deviation
#'   of the noise in \code{x}. The default value is \code{mad(diff(diff(x)))/sqrt(6)}, where
#'   \code{mad(x)} denotes the median absolute deviation of \code{x} computed under the
#'   assumption that the noise is independent and identically distributed from the
#'   Gaussian distribution.
#' @param thr_const A positive real number with default value equal to 1.4. It is
#'   used to define the threshold; see \code{thr_fin}.
#' @param thr_fin With \code{T} the length of the data sequence, this is a positive real number
#'   with default value equal to \code{sigma * thr_const * sqrt(2 * log(T))}. It is the threshold,
#'   which is used in the detection process.
#' @param s,e Positive integers with \code{s} less than \code{e}, which indicate
#'   that you want to check for change-points in the data sequence with subscripts
#'   in \code{[s,e]}. The default values are \code{s} equal to 1 and
#'   \code{e} equal to \code{T}, with \code{T} the length of the data sequence.
#' @param points A positive integer with default value equal to 3. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively; see Details for more information.
#' @param k_l,k_r Positive integer numbers that get updated whenever the function
#'   calls itself during the detection process. They are not essential for the
#'   function to work, and we include them only to reduce the computational time.
#' @details The change-point detection algorithm that is used in \code{\link{cplm_th}} is the
#'   Isolate-Detect methodology described in ``Detecting multiple generalized
#'   change-points by isolating single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#'   The concept is simple and is split into two stages; firstly, isolation of each
#'   of the true change-points in subintervals of the data domain, and secondly their detection.
#'   ID first creates two ordered sets of \eqn{K = \lceil T/\code{points}\rceil} right- and left-expanding
#'   intervals as follows. The \eqn{j^{th}} right-expanding interval is \eqn{R_j = [1, j\times \code{points}]},
#'   while the \eqn{j^{th}} left-expanding interval is \eqn{L_j = [T - j\times \code{points} + 1, T]}.
#'   We collect these intervals in the ordered set \eqn{S_{RL} = \lbrace R_1, L_1, R_2, L_2, ... , R_K, L_K\rbrace}.
#'   For a suitably chosen contrast function, ID first identifies the point with the maximum contrast
#'   value in \eqn{R_1}. If its value exceeds a certain threshold, then it is taken as a change-point.
#'   If not, then the process tests the next interval in \eqn{S_{RL}} and repeats the above process.
#'   Upon detection, the algorithm makes a new start from estimated location.
#' @return
#'   A numeric vector with the detected change-points.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{win_cplm_th}}, \code{\link{ID_cplm}}, and \code{\link{ID}}, which employ
#' this function. In addition, see \code{\link{pcm_th}} for the case of detecting changes in
#' a piecewise-constant signal via thresholding.
#' @examples
#' single.cpt <- c(seq(0, 999, 1), seq(998.5, 499, -0.5))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single.th <- cplm_th(single.cpt.noise)
#'
#' three.cpt <- c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251,1249,2), seq(1248,749,-1))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' cpt.three.th <- cplm_th(three.cpt.noise)
#'
#' multi.cpt <- rep(c(seq(0,49,1), seq(48,0,-1)),20)
#' multi.cpt.noise <- multi.cpt + rnorm(1980)
#' cpt.multi.th <- cplm_th(multi.cpt.noise)
cplm_th <- function(x, sigma = stats::mad(diff(diff(x))) / sqrt(6),
                   thr_const = 1.4, thr_fin = sigma * thr_const * sqrt(2 * log(length(x))),
                   s = 1, e = length(x), points = 3, k_l = 1, k_r = 1) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  points <- as.integer(points)
  l <- length(x)
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, - points)
  chp <- 0
  if (e - s <= 2) {
    chp <- 0
    cpt <- chp
  } else {
      pos_r <- numeric()
      CUSUM_r <- numeric()
      pos_l <- numeric()
      CUSUM_l <- numeric()
      moving_points <- s_e_points(r_e_points, l_e_points, s, e)
      right_points <- moving_points[[1]]
      left_points <- moving_points[[2]]
      lur <- length(left_points)
      rur <- length(right_points)
      if (k_r < k_l) {
          while ( (chp == 0) & (k_r < min(k_l, rur))) {
              x_temp_r <- x[s:right_points[k_r]]
              ipcr <- cumsum_lin(x_temp_r)
              pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
              CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
              if (CUSUM_r[k_r] > thr_fin) {
                chp <- pos_r[k_r]
              } else {
                k_r <- k_r + 1
              }
          }
      }
      if (k_l < k_r) {
          while ( (chp == 0) & (k_l < min(k_r, lur))) {
              x_temp_l <- x[left_points[k_l]:e]
              ipcl <- cumsum_lin(x_temp_l)
              pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
              CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
              if (CUSUM_l[k_l] > thr_fin) {
                chp <- pos_l[k_l]
              } else {
                k_l <- k_l + 1
              }
          }
      }
      if (chp == 0) {
          while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
              x_temp_r <- x[s:right_points[k_r]]
              ipcr <- cumsum_lin(x_temp_r)
              pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
              CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
              if (CUSUM_r[k_r] > thr_fin) {
                 chp <- pos_r[k_r]
              } else {
                x_temp_l <- x[left_points[k_l]:e]
                ipcl <- cumsum_lin(x_temp_l)
                pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
                CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
                if (CUSUM_l[k_l] > thr_fin) {
                  chp <- pos_l[k_l]
                } else {
                  k_r <- k_r + 1
                  k_l <- k_l + 1
                }
              }
          }
      }
      if (chp != 0) {
          if (chp > ( (e + s) / 2)) {
              r <- cplm_th(x, s = s, e = chp, points = points,
                          thr_fin = thr_fin, k_r = k_r, k_l = 1)
          } else {
              r <- cplm_th(x, s = chp + 1, e = e, points = points,
                          thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
          }
          cpt <- c(chp, r)
      } else {
          cpt <- chp
      }
  }
  cpt <- cpt[cpt != 0]
  return(sort(cpt))
}

#' A windows-based approach for multiple change-point detection in a continuous,
#' piecewise-linear signal via thresholding
#'
#' This function performs the windows-based variant of the Isolate-Detect methodology
#' with the thresholding-based stopping rule in order to detect multiple change-points
#' in a continuous, piecewise-linear noisy data sequence, with the noise being Gaussian.
#' It is particularly helpful for very long data sequences, as due to applying Isolate-Detect
#' on moving windows, the computational time is reduced. See Details for a brief explanation of
#' this approach and for the relevant literature reference.
#'
#' @export
#' @param xd A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param sigma A positive real number. It is the estimate of the standard deviation
#'   of the noise in \code{xd}. The default value is \code{mad(diff(diff(xd)))/sqrt(6)}, where
#'   \code{mad(xd)} denotes the median absolute deviation of \code{xd} computed under the
#'   assumption that the noise is independent and identically distributed from the
#'   Gaussian distribution.
#' @param thr_con A positive real number with default value equal to 1.4. It is
#'   used to define the threshold. The change-points are estimated by thresholding
#'   with threshold equal to \code{sigma * thr_con * sqrt(2 * log(T))}, where
#'   \code{T} is the length of the data sequence \code{xd}.
#' @param c_win A positive integer with default value equal to 3000. It is the length
#'   of each window for the data sequence in hand. Isolate-Detect will be applied
#'   in segments of the form \code{[(i-1) * c_win + 1, i * c_win]}, for \eqn{i=1,2,...,K},
#'   where \eqn{K} depends on the length \code{T} of the data sequence.
#' @param w_points A positive integer with default value equal to 3. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @param l_win A positive integer with default value equal to 12000. If the length of
#'   the data sequence is less than or equal to \code{l_win}, then the windows-based approach
#'   will not be applied and the result will be obtained by the classical Isolate-Detect
#'   methodology based on thresholding.
#' @details The method that is implemented by this function is based on splitting the given
#'   data sequence uniformly into smaller parts (windows), to which Isolate-Detect, based on the
#'   thresholding stopping rule (see \code{\link{cplm_th}}), is then applied.
#' @return
#'   A numeric vector with the detected change-points.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{cplm_th}}, which is the function that \code{\link{win_cplm_th}} is based on. Also,
#' see \code{\link{ID_cplm}} and \code{\link{ID}}, which employ \code{\link{win_cplm_th}}. In addition,
#' see \code{\link{win_pcm_th}} for the case of detecting changes in a piecewise-constant signal via
#' thresholding.
#' @examples
#' single.cpt <- c(seq(0, 999, 1), seq(998.5, 499, -0.5))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single.th <- win_cplm_th(single.cpt.noise)
#'
#' three.cpt <- c(seq(0, 3999, 1), seq(3998.5, 1999, -0.5), seq(2001,9999,2), seq(9998,5999,-1))
#' three.cpt.noise <- three.cpt + rnorm(16000)
#' cpt.three.th <- win_cplm_th(three.cpt.noise)
win_cplm_th <- function(xd, sigma = stats::mad(diff(diff(xd))) / sqrt(6),
                        thr_con = 1.4, c_win = 3000, w_points = 3,
                        l_win = 12000) {
  if (!(is.numeric(xd))){
    stop("The input in `xd' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_con <= 0) || (w_points <= 0) || (c_win <= 0) || (l_win <= 0)){
    stop("The threshold constant as well as the `w_points', `c_win', `l_win' arguments should
         be positive numbers.")
  }
  if ( (abs(w_points - round(w_points)) > .Machine$double.eps ^ 0.5)
      || (abs(c_win - round(c_win)) > .Machine$double.eps ^ 0.5)
      || (abs(l_win - round(l_win)) > .Machine$double.eps ^ 0.5)){
    warning("The input values  for `w_points', `c_win', and  `l_win' should be positive integers.
            If either of them is a positive real number then the integer part of the given number
            is used to obtain the result.")
  }
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  l_win <- as.integer(l_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= l_win) {
      u <- cplm_th(x = xd, thr_const = thr_con, points = w_points)
      return(u)
  } else {
      K <- ceiling(lg / c_win)
      tsm <- list()
      u <- list()
      ufin <- numeric()
      uaddition <- list()
      tsm[[1]] <- xd[1:c_win]
      ufin <- cplm_th(tsm[[1]], thr_fin = t, points = w_points)
      uaddition[[1]] <- numeric()
      uaddition[[1]] <- cplm_th(x = xd[(max(1, c_win - (5 * w_points) + 1)):min(c_win + (5 * w_points), lg)], thr_fin = t, points = 2) + c_win - (5 * w_points)
      i <- 2
      while (i < K) {
          tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
          u[[i]] <- cplm_th(x = tsm[[i]], thr_fin = t, points = w_points) + (i - 1) * c_win
          uaddition[[i]] <- numeric()
          uaddition[[i]] <- cplm_th(x = xd[(max(1, i * c_win - (5 * w_points) + 1)):(min(i * c_win + (5 * w_points), lg))], thr_fin = t, points = 2) + i * c_win - (5 * w_points)
          ufin <- c(ufin, u[[i]], uaddition[[i]])
          i <- i + 1
      }
      tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
      u[[K]] <- cplm_th(tsm[[K]], thr_fin = t, points = w_points) + (K - 1) * c_win
      ufinl <- c(ufin, u[[K]])
      return(sort(unique(ufinl)))
  }
}

linear_contr_one <- function(x, s, e, b) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector.")
  }
  if ( (length(s) != length(b)) || (length(s) != length(e)) || (length(e) != length(b))){
    stop("The vectors s, b, e, should be of the same length")
  }
  if (any(s < 1) | any(b < 1) | any(e < 1)){
    stop("The entries of the vectors s, b, e should be positive integers.")
  }
  if (any(s > b) | any(b >= e)){
    stop("The value for b should be in the interval [s,e)")
  }
  if ( (any(abs( (s - round(s))) > .Machine$double.eps ^ 0.5))
      || (any(abs( (b - round(b))) > .Machine$double.eps ^ 0.5))
      || (any(abs( (e - round(e))) > .Machine$double.eps ^ 0.5))){
    stop("The input values  for s, b, and  e should be positive integers.")
  }
  r <- numeric()
  for (j in 1:length(b)) {
    x1 <- x[s[j]:e[j]]
    n <- length(x1)
    if ( (b[j] - s[j] + 1) == 1) {
      r[j] <- 0
    } else {
      y1 <- cumsum(x1 * (1:n))
      y <- cumsum(x1)
      a <- sqrt(6 / ( (n - 1) * n * (n + 1) * (2 - 2 * (b[j] - s[j] + 1) ^ 2 + 2 * (b[j] - s[j] + 1) * n - 1 + 2 * (b[j] - s[j] + 1) - n)))
      be <- sqrt( ( (n - (b[j] - s[j] + 1) + 1) * (n - (b[j] - s[j] + 1))) / ( (b[j] - s[j]) * (b[j] - s[j] + 1)))
      r[j] <- a * be * ( (2 * (b[j] - s[j] + 1) + n - 1) * y1[b[j] - s[j] + 1] - (n + 1) * (b[j] - s[j] + 1) * y[b[j] - s[j] + 1]) - (a / be) * ( (3 * n - 2 * (b[j] - s[j] + 1) + 1) * (y1[n] - y1[b[j] - s[j] + 1]) - (n + 1) * (2 * n - (b[j] - s[j] + 1)) * (y[n] - y[b[j] - s[j] + 1]))
    }
  }
  return(abs(r))
}

#' The solution path for the case of continuous piecewise-linear signals
#'
#' This function starts by over-estimating the number of true change-points.
#' After that, following an approach based on the values of a suitable contrast function,
#' it sorts the estimated change-points in a way that the estimation, which is
#' most-likely to be correct appears first, whereas the least likely to be correct,
#' appears last. The routine is typically not called directly by the user; it is
#' employed in \code{\link{cplm_ic}}. For more details, see References.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param thr_ic A positive real number with default value equal to 1.25. It is
#'   used to define the threshold. The change-points are estimated by thresholding
#'   with threshold equal to \code{sigma * thr_ic * sqrt(2 * log(T))}, where
#'   \code{T} is the length of the data sequence \code{x} and \code{sigma = mad(diff(diff(x)))/6}.
#'   Because, we would like to overestimate the number of the true change-points in \code{x}, it is
#'   suggested to keep \code{thr_ic} smaller than 1.4, which is the default value used as
#'   the threshold constant in the function \code{\link{win_cplm_th}}.
#' @param points A positive integer with default value equal to 3. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @return
#'   The solution path for the case of continuous piecewise-linear signals.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @references Anastasiou, A. and Fryzlewicz, P. (2018). Detecting multiple generalized change-points
#' by isolating single ones.
#' @examples
#' three.cpt <- c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(250.5,999,1.5), seq(998,499,-1))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' solution.path <- sol_path_cplm(three.cpt.noise)
sol_path_cplm <- function(x, thr_ic = 1.25, points = 3) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which the solution path will be calculated.")
  }
  if ( (thr_ic <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  lx_ic <- length(x)
  points <- as.integer(points)
  cpt_lower <- win_cplm_th(x, thr_con = thr_ic, w_points = points)
  lcpt_ic <- length(cpt_lower)
  if ( (lcpt_ic == 1) | (lcpt_ic == 0)) {
      return(cpt_lower)
  } else {
      seb_set <- c(unique(c(1, cpt_lower)), lx_ic)
      lseb_set <- length(seb_set)
      min_C <- numeric()
      while (lseb_set >= 3) {
          Rs <- linear_contr_one(x, seb_set[1:(lseb_set - 2)],
                                 seb_set[3:(lseb_set)],
                                 seb_set[2:(lseb_set - 1)])
          min_Rs <- seb_set[2:(lseb_set - 1)][which.min(Rs)]
          min_C <- c(min_C, min_Rs)
          seb_set <- seb_set[-which(seb_set == min_Rs)]
          lseb_set <- lseb_set - 1
      }
      return(min_C[length(min_C):1])
  }
}

#' Calculate the residuals related to the estimated signal
#'
#' This function returns the difference between \code{x} and the estimated signal
#' with change-points at \code{cpt}. The input in the argument \code{type_chg} will
#' indicate the type of changes in the signal.
#'
#' @export
#' @param x A numeric vector containing the data.
#' @param cpt A positive integer vector with the locations of the change-points.
#'   If missing, the \code{\link{ID}} function is called internally to detect any change-points
#'   that might be present in \code{x}.
#' @param type_chg A character string, which defines the type of the detected change-points.
#'   If \code{type_chg = ``mean''}, then the change-points represent the locations of changes
#'   in the mean of a piecewise-constant signal. If \code{type_chg = ``slope''}, then the
#'   change-points represent the locations of changes in the slope of a piecewise-linear
#'   and continuous signal.
#' @param type_res A choice of \code{``raw''} and \code{``standardised''} residuals.
#' @return
#'   If \code{type_res = ``raw''}, the function returns the difference between the data
#'   and the estimated signal. If \code{type_res = ``standardised''}, then the function
#'   returns the difference between the data and the estimated signal, divided by
#'   the estimated standard deviation.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @examples
#' single.cpt.pcm <- c(rep(4,1000),rep(0,1000))
#' single.cpt.pcm.noise <- single.cpt.pcm + rnorm(2000)
#' cpt_detect <- ID(single.cpt.pcm.noise, contrast = "mean")
#'
#' residuals_cpt_raw <- resid_ID(single.cpt.pcm.noise, cpt = cpt_detect$cpt, type_chg = "mean",
#' type_res = "raw")
#'
#' residuals_cpt_stand. <- resid_ID(single.cpt.pcm.noise, cpt = cpt_detect$cpt, type_chg = "mean",
#' type_res = "standardised")
#'
#' plot(residuals_cpt_raw)
#' plot(residuals_cpt_stand.)
resid_ID <- function(x, cpt, type_chg = c("mean", "slope"),
                  type_res = c("raw", "standardised")) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which the residuals will be calculated.")
  }
  x <- as.numeric(x)
  if (NA %in% x)
    stop("x vector cannot contain NA's")
  type_res <- match.arg(type_res, c("raw", "standardised"))
  type_chg <- match.arg(type_chg, c("mean", "slope"))
  if (!is.null(cpt))
    if (any(is.na(cpt)))
      cpt <- cpt[!is.na(cpt)]
  cpt <- as.integer(cpt)
  if (missing(cpt))
    cpt <- ID(x,  contrast = type_chg)$cpt
  res <- x - est_signal(x, cpt = cpt, type = type_chg)
  if (type_res == "raw")
      return(res) else return(res / stats::sd(res))
}

log_lik_slope <- function(x, cpt) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which the residuals will be calculated.")
  }
  x <- as.numeric(x)
  if (NA %in% x)
    stop("x vector cannot contain NA's")
  lx <- length(x)
  if (missing(cpt))
      cpt <- ID(x, contrast = "slope")$cpt
  if (!is.null(cpt))
      if (any(is.na(cpt)))
          cpt <- cpt[!is.na(cpt)]
  n_cpt <- length(cpt)
  res <- - lx / 2 * log(mean(resid_ID(x, cpt = cpt, type_chg = "slope") ^ 2))
  attr(res, "df") <- 2 * n_cpt + 3
  return(res)
}

#' Multiple change-point detection in a continuous piecewise-linear signal
#' via minimising an information criterion
#'
#' This function performs the Isolate-Detect methodology based on an information
#' criterion approach, in order to detect multiple change-points in a noisy, continuous,
#' piecewise-linear data sequence, with the noise being Gaussian. More information on
#' how this approach works as well as the relevant literature reference are given in Details.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param th_const A positive real number with default value equal to 1.25. It is
#'   used to define the threshold value that will be used at the first step of the
#'   model selection based Isolate-Detect method; see Details for more information.
#' @param Kmax A positive integer with default value equal to 200. It is the
#'   maximum allowed number of estimated change-points in the solution path; see
#'   \code{\link{sol_path_cplm}} for more details.
#' @param penalty A character vector with names of penalty functions used.
#' @param points A positive integer with default value equal to 10. It defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @details The approach followed in \code{\link{cplm_ic}} in order to detect the
#'   change-points is based on identifying the set of change-points that minimise an
#'   information criterion. At first, we employ \code{\link{sol_path_cplm}}, which
#'   overestimates the number of change-points using \code{th_const} in order to define the
#'   threshold and then sorts the obtained estimates in a way that the estimate,
#'   which is most likely to be correct appears first, whereas the least likely
#'   to be correct, appears last. Let \eqn{J} be the number of estimates
#'   that this overestimation approach returns. We will obtain a vector
#'   \eqn{b = (b_1, b_2, ..., b_J)}, with the estimates ordered as explained above. We
#'   define the collection \eqn{\left\{M_j\right\}_{j = 0,1,\ldots,J}}, where \eqn{M_0}
#'   is the empty set and \eqn{M_j = \left\{b_1,b_2,...,b_j\right\}}. Among the collection
#'   of models \eqn{M_j, j=0,1,...,J}, we select the one that minimises a predefined
#'   Information Criterion. The obtained set of change-points is apparently a subset of
#'   the solution path given in \code{\link{sol_path_cplm}}. More details can be found
#'   in ``Detecting multiple generalized change-points by isolating single ones'',
#'   Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'    \tabular{ll}{
#'    \cr \code{sol_path} \tab A vector containing the solution path.
#'    \cr \code{ic_curve} \tab A list with values of the chosen information criteria.
#'    \cr \code{cpt_ic} \tab A list with the change-points detected for each information
#'   criterion considered.
#'    \cr \code{no_cpt_ic} \tab The number of change-points detected for each information
#'   criterion considered.
#'    }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{ID_cplm}} and \code{\link{ID}}, which employ this function.
#' In addition, see \code{\link{pcm_ic}} for the case of detecting changes in
#' a piecewise-constant signal using the information criterion based approach.
#' @examples
#' single.cpt <- c(seq(0, 999, 1), seq(998.5, 499, -0.5))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single.ic <- cplm_ic(single.cpt.noise)
#'
#' three.cpt <- c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(250,1249,2), seq(1248,749,-1))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' cpt.three.ic <- cplm_ic(three.cpt.noise)
cplm_ic <- function(x, th_const = 1.25, Kmax = 200,
                       penalty = c("ssic_pen", "sic_pen"), points = 10) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you want to look for change-points.")
  }
  if ( (th_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  result <- list()
  penalty <- as.character(penalty)
  lx <- length(x)
  if (length(penalty)) {
      if (Kmax == 0 || lx <= 3) {
          result$ssic <- result$bic <- NA
          result$no_cpt_ic[c("ssic", "bic")] <- as.integer(c(0, 0))
          result$ic_curve <- list()
          result$ic_curve$bic <- result$ic_curve$ssic <- NA
          if (Kmax == 0) {
              stop("no change-points found, choose larger Kmax")
          } else {
              stop("sample size is too small")
          }
    } else {
          cpt_cand <- sol_path_cplm(x, thr_ic = th_const, points = points)
          if (length(cpt_cand) == 0){
          result$cpt_ic[c("ssic_pen", "sic_pen")] <- c(0, 0)
          result$no_cpt_ic[c("ssic_pen", "sic_pen")] <- as.integer(c(0, 0))
          result$ic_curve <- list()
          result$ic_curve$ssic_pen <- result$ic_curve$ssic_pen <- NA
          }
          if (length(cpt_cand) > min(Kmax, lx - 2))
              cpt_cand <- cpt_cand[1:min(Kmax, lx - 2)]
          len_cpt <- length(cpt_cand)
          result$sol_path <- cpt_cand
          result$ic_curve <- list()
          for (j in 1:length(penalty)) {
              result$ic_curve[[penalty[j]]] <- rep(0, len_cpt + 1)
          }
          if (len_cpt) {
              for (i in 1:len_cpt) {
                cur_cpt <- cpt_cand[1:i]
                lh <- log_lik_slope(x, cpt = cur_cpt)
                n_p <- attr(lh, "df")
                for (j in 1:length(penalty)) {
                  result$ic_curve[[penalty[j]]][i + 1] <- -2 * lh + do.call(penalty[j], list(n = lx, n_param = n_p))
                }
              }
          }
          for (j in 1:length(penalty)) {
              lh <- log_lik_slope(x, cpt = c())
              n_p <- attr(lh, "df")
              result$ic_curve[[penalty[j]]][1] <- -2 * lh + do.call(penalty[j], list(n = lx, n_param = n_p))
          }
          result$cpt_ic <- list()
          for (j in 1:length(penalty)) {
              tmp <- which.min(result$ic_curve[[penalty[j]]])
              if (tmp == 1) {
                result$cpt_ic[[penalty[j]]] <- NA
                result$no_cpt_ic[penalty[j]] <- as.integer(0)
              } else {
                result$cpt_ic[[penalty[j]]] <- sort(cpt_cand[1:(tmp - 1)])
                result$no_cpt_ic[penalty[j]] <- as.integer(tmp - 1)
              }
          }
      }
  }
  return(result)
}

#' Multiple change-point detection for a continuous, piecewise-linear signal
#' using the Isolate-Detect methodology
#'
#' This function estimates the number and locations of multiple change-points in the noisy,
#' continuous and piecewise-linear input vector \code{x}, using the Isolate-Detect methodology. The noise
#' follows the normal distribution. The estimated signal, as well as the solution path defined
#' in \code{\link{sol_path_cplm}} are also given. The function is a hybrid between the thresholding
#' approach of \code{\link{win_cplm_th}} and the information criterion approach of
#' \code{\link{cplm_ic}} and estimates the change-points taking into account both these
#' approaches (see Details for more information and the relevant literature reference).
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param thr_id A positive real number with default value equal to 1.4. It is
#'   used to define the threshold, if the thresholding approach is to be followed; see
#'   \code{\link{cplm_th}} for more details.
#' @param th_ic_id  A positive real number with default value equal to 1.25. It is
#'   useful only if the model selection based Isolate-Detect method is to be followed
#'   and it is used to define the threshold value that will be used at the first step
#'   (change-point overestimation) of the model selection approach described in
#'   \code{\link{cplm_ic}}.
#' @param pointsth A positive integer with default value equal to 3. It is used only
#'   when the threshold based approach is to be followed and it defines the distance
#'   between two consecutive end- or start-points of the right- or left-expanding intervals,
#'   respectively.
#' @param pointsic A positive integer with default value equal to 10. It is used only
#'   when the information criterion based approach is to be followed and it defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @details Firstly, this function detects the change-points using \code{\link{win_cplm_th}}.
#'   If the estimated number of change-points is larger than 100, then the
#'   result is returned and we stop. Otherwise, \code{\link{ID_cplm}} proceeds to detect the
#'   change-points using \code{\link{cplm_ic}} and this is what is returned. To sum up,
#'   \code{\link{ID_cplm}} returns a result based on \code{\link{cplm_ic}} if the estimated number
#'   of change-points is less than 100. Otherwise, the result comes from thresholding.
#'   More details can be found in ``Detecting multiple generalized change-points by
#'   isolating single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{cpt} \tab A vector with the detected change-points.
#'    \cr \code{no_cpt} \tab The number of change-points detected.
#'    \cr \code{fit} \tab A numeric vector with the estimated continuous piecewise-linear signal.
#'    \cr \code{solution_path} \tab A vector containing the solution path.
#'  }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{win_cplm_th}} and \code{\link{cplm_ic}} which are the functions that
#' \code{\link{ID_cplm}} is based on. In addition, see \code{\link{ID_pcm}} for the case of detecting changes
#' in the mean of a piecewise-constant signal. The main function \code{\link{ID}} of the package
#' employs \code{\link{ID_cplm}}.
#' @examples
#' single.cpt <- c(seq(0, 999, 1), seq(998.5, 499, -0.5))
#' single.cpt.noise <- single.cpt + rnorm(2000)
#' cpt.single <- ID_cplm(single.cpt.noise)
#'
#' three.cpt <- c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(250,1249,2), seq(1248,749,-1))
#' three.cpt.noise <- three.cpt + rnorm(2000)
#' cpt.three <- ID_cplm(three.cpt.noise)
#'
#' multi.cpt <- rep(c(seq(0,49,1), seq(48,0,-1)),20)
#' multi.cpt.noise <- multi.cpt + rnorm(1980)
#' cpt.multi <- ID_cplm(multi.cpt.noise)
ID_cplm <- function(x, thr_id = 1.4, th_ic_id = 1.25, pointsth = 3,
                   pointsic = 10) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_id <= 0) || (pointsth <= 0) || (th_ic_id <= 0) || (pointsic <= 0)){
    stop("The threshold constants as well as the `pointsth' and `pointsic' arguments
         should be positive numbers.")
  }
  if ( (abs(pointsth - round(pointsth)) > .Machine$double.eps ^ 0.5)
      || (abs(pointsic - round(pointsic)) > .Machine$double.eps ^ 0.5)){
    warning("The input values for `pointsth' and `pointsic' should be positive integers. If they are
            positive real numbers then the integer parts of the given numbers are used as the values
            of `pointsth' and `pointsic'.")
  }
  z <- win_cplm_th(x, thr_con = thr_id, w_points = pointsth)
  if (length(z) >= 100) {
      cpt <- z
      nocpt <- length(z)
      mean <- est_signal(x, cpt, type = "slope")
      result <- list(cpt = sort(cpt), no_cpt = nocpt, fit = mean)
  } else {
      z <- cplm_ic(x, th_const = th_ic_id, points = pointsic)
      if (is.na(z$cpt_ic[[1]][1])) {
          cpt <- 0
          nocpt <- 0
      } else {
          cpt <- z$cpt_ic[[1]]
          nocpt <- z$no_cpt_ic[[1]]
      }
      mean <- est_signal(x, cpt, type = "slope")
      sol_path <- z$sol_path
      result <- list(cpt = sort(cpt), no_cpt = nocpt, fit = mean,
                     solution_path = sol_path)
  }
  return(result)
}

#' Transform the noise to be closer to the Gaussian distribution
#'
#' This function pre-processes the given data in order to obtain a noise structure
#' that is closer to satisfying the Gaussianity assumption. See details for more information
#' and for the relevant literature reference.
#'
#' @export
#' @param x A numeric vector containing the data.
#' @param sc A positive integer number with default value equal to 3. It is
#'   used to define the way we pre-average the given data sequence.
#' @details For a given natural number \code{sc} and data \code{x} of length \eqn{T}, let us
#' denote by \eqn{Q = \lceil T/sc \rceil}. Then, \code{\link{normalise}} calculates
#' \deqn{\tilde{x}_q = 1/sc\sum_{t=(q-1) * sc + 1}^{q * sc}x_t,} for \eqn{q=1, 2, ..., Q-1}, while
#' \deqn{\tilde{x}_Q = (T - (Q-1) * sc)^{-1}\sum_{t = (Q-1) * sc + 1}^{T}x_t.}
#' More details can be found in the preprint ``Detecting multiple generalized
#' change-points by isolating single ones'', Anastasiou and Fryzlewicz (2018).
#' @return
#' The ``normalised'' vector \eqn{\tilde{x}} of length \eqn{Q}, as explained in Details.
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{ht_ID_pcm}} and \code{\link{ht_ID_cplm}}, which are
#'  functions that employ \code{\link{normalise}}.
#' @examples
#' t5 <- rt(n = 10000, df = 5)
#' n5 <- normalise(t5, sc = 3)
normalise <- function(x, sc = 3) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data.")
  }
  if ( (sc <= 0)){
    stop("The scale constant which is used to define the way we preaverage the data
         should be a positive number.")
  }
  if (abs(sc - round(sc)) > .Machine$double.eps ^ 0.5){
    warning("The input for `sc' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `sc'.")
  }
  sc <- as.integer(sc)
  res <- numeric()
  l <- length(x)
  l1 <- ceiling(l / sc)
  for (i in 1:(l1 - 1)) {
      res[i] <- (1 / sqrt(sc)) * sum(x[( (i - 1) * sc + 1):(i * sc)])
  }
  res[l1] <- (1 / sqrt(l - (l1 - 1) * sc)) * sum(x[( (l1 - 1) * sc + 1):l])
  return(res)
}

#' Apply the Isolate-Detect methodology for multiple change-point detection in the
#' mean of a vector with non Gaussian noise
#'
#' Using the Isolate-Detect methodology, this function estimates the number and locations
#' of multiple change-points in the mean of the noisy, piecewise-constant input vector \code{x},
#' with noise that is not normally distributed. It also gives the estimated signal, as well as
#' the solution path defined in \code{\link{sol_path_pcm}}. See Details for the relevant literature reference.
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param s.ht A positive integer number with default value equal to 3. It is
#'   used to define the way we pre-average the given data sequence (see Details).
#' @param q_ht A positive integer number with default value equal to 300. If the
#'   length of \code{x} is less than or equal to \code{q_ht}, then no pre-averaging
#'   will take place.
#' @param ht_thr_id A positive real number with default value equal to 1. It is
#'   used to define the threshold, if the thresholding approach is to be followed; see
#'   \code{\link{pcm_th}} for more details on the thresholding approach.
#' @param ht_th_ic_id  A positive real number with default value equal to 0.9. It is
#'   useful only if the model selection based Isolate-Detect method is to be followed
#'   and it is used to define the threshold value that will be used at the first step
#'   (change-point overestimation) of the model selection approach described in
#'   \code{\link{pcm_ic}}. It is applied to the new data, which are obtained after
#'   we pre-average \code{x}.
#' @param p_thr A positive integer with default value equal to 1. It is used only
#'   when the threshold based approach (as described in \code{\link{pcm_th}}) is to be followed
#'   and it defines the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @param p_ic A positive integer with default value equal to 3. It is used only
#'   when the information criterion based approach (described in \code{\link{pcm_ic}})
#'   is to be followed and it defines the distance between two consecutive end- or start-points
#'   of the right- or left-expanding intervals, respectively.
#' @details Firstly, in this function we call \code{\link{normalise}}, in order to
#'   create a new data sequence, \eqn{\tilde{x}}, by taking averages of observations in
#'   \code{x}. Then, we employ \code{\link{ID_pcm}} on \eqn{\tilde{x}_q} to obtain the
#'   change-points, namely \eqn{\tilde{r}_1, \tilde{r}_2, ..., \tilde{r}_{\hat{N}}} in
#'   increasing order. To obtain the original location of the change-points with,
#'   on average, the highest accuracy we define
#'   \eqn{\hat{r}_k = (\tilde{r}_{k}-1)*\code{s.ht} + \lfloor \code{s.ht}/2 + 0.5 \rfloor, k=1, 2,..., \hat{N}.}
#'   More details can be found in ``Detecting multiple generalized change-points by
#'   isolating single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{cpt} \tab A vector with the detected change-points.
#'    \cr \code{no_cpt} \tab The number of change-points detected.
#'    \cr \code{fit} \tab A numeric vector with the estimated piecewise-constant signal.
#'    \cr \code{solution_path} \tab A vector containing the solution path.
#'  }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{ID_pcm}} and \code{\link{normalise}}, which are functions that are
#' used in \code{\link{ht_ID_pcm}}. In addition, see \code{\link{ht_ID_cplm}} for the case
#' of continuous and piecewise-linear signals.
#' @examples
#' single.cpt <- c(rep(4,3000),rep(0,3000))
#' single.cpt.student <- single.cpt + rt(6000, df = 5)
#' cpts_detect <- ht_ID_pcm(single.cpt.student)
#'
#' three.cpt <- c(rep(4,2000),rep(0,2000),rep(-4,2000),rep(0,2000))
#' three.cpt.student <- three.cpt + rt(8000, df = 5)
#' cpts_detect_three <- ht_ID_pcm(three.cpt.student)
ht_ID_pcm <- function(x, s.ht = 3, q_ht = 300, ht_thr_id = 1,
                      ht_th_ic_id = 0.9, p_thr = 1, p_ic = 3) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (q_ht <= 0)){
    stop("The input for `q_ht' should be a positive number.")
  }
  if (abs(q_ht - round(q_ht)) > .Machine$double.eps ^ 0.5){
    warning("The input for `q_ht' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `q_ht'.")
  }
  if (length(x) <= q_ht) {
      result_new <- ID_pcm(x)
  } else {
      result_new <- list()
      nor <- normalise(x, sc = s.ht)
      result <- ID_pcm(x = nor, thr_id = ht_thr_id, th_ic_id = ht_th_ic_id,
                       pointsth = p_thr, pointsic = p_ic)
      r1 <- result$cpt
      r2 <- result$solution_path
      if (r1[1] == 0) {
          result_new$cpt <- 0
          result_new$no_cpt <- 0
        } else {
          result_new <- result
          result_new$cpt <- (r1 - 1) * (s.ht) + round(s.ht / 2)
          result_new$no_cpt <- length(result_new$cpt)
      }
      result_new$fit <- est_signal(x, result_new$cpt, type = "mean")
      result_new$solution_path <- (r2 - 1) * (s.ht) + round(s.ht / 2)
  }
  return(result_new)
}

#' Apply the Isolate-Detect methodology for multiple change-point detection in a
#' continuous, piecewise-linear vector with non Gaussian noise
#'
#' Using the Isolate-Detect methodology, this function estimates the number and locations
#' of multiple change-points in the noisy, continuous, piecewise-linear input vector \code{x},
#' with noise that is not normally distributed. It also gives the estimated signal, as well as
#' the solution path defined in \code{\link{sol_path_cplm}} (see Details for the relevant
#' literature reference).
#'
#' @export
#' @param x A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param s.ht A positive integer number with default value equal to 3. It is
#'   used to define the way we pre-average the given data sequence. For more information
#'   see Details.
#' @param q_ht A positive integer number with default value equal to 300. If the
#'   length of \code{x} is less than or equal to \code{q_ht}, then no pre-averaging
#'   will take place.
#' @param ht_thr_id A positive real number with default value equal to 1.4. It is
#'   used to define the threshold, if the thresholding approach (described in \code{\link{cplm_th}})
#'   is to be followed.
#' @param ht_th_ic_id  A positive real number with default value equal to 1.25. It is
#'   useful only if the model selection based Isolate-Detect method is to be followed
#'   and it is used to define the threshold value that will be used at the first step
#'   (change-point overestimation) of the model selection approach described in
#'   \code{\link{cplm_ic}}. It is applied to the new data, which are obtained
#'   after we take average values on \code{x}.
#' @param p_thr A positive integer with default value equal to 1. It is used only
#'   when the threshold based approach (described in \code{\link{cplm_th}}) is to
#'   be followed and it defines the distance between two consecutive end- or start-points
#'   of the right- or left-expanding intervals, respectively.
#' @param p_ic A positive integer with default value equal to 3. It is used only
#'   when the information criterion based approach (described in \code{\link{cplm_ic}}) is
#'   to be followed and it defines the distance between two consecutive end- or start-points
#'   of the right- or left-expanding intervals, respectively.
#' @details Firstly, in this function we call \code{\link{normalise}}, in order to
#'   create a new data sequence, \eqn{\tilde{x}}, by taking averages of observations in
#'   \code{x}. Then, we employ \code{\link{ID_cplm}} on \eqn{\tilde{x}_q} to obtain the
#'   change-points, namely \eqn{\tilde{r}_1, \tilde{r}_2, ..., \tilde{r}_{\hat{N}}} in
#'   increasing order. To obtain the original location of the change-points with,
#'   on average, the highest accuracy we define
#'   \deqn{\hat{r}_k = (\tilde{r}_{k}-1)*\code{s.ht} + \lfloor \code{s.ht}/2 + 0.5 \rfloor, k=1, 2,..., \hat{N}.}
#'   More details can be found in ``Detecting multiple generalized change-points by
#'   isolating single ones'', Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{cpt} \tab A vector with the detected change-points.
#'    \cr \code{no_cpt} \tab The number of change-points detected.
#'    \cr \code{fit} \tab A numeric vector with the estimated continuous piecewise-linear signal.
#'    \cr \code{solution_path} \tab A vector containing the solution path.
#'  }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{ID_cplm}} and \code{\link{normalise}}, which are functions that are
#' used in \code{\link{ht_ID_cplm}}. In addition, see \code{\link{ht_ID_pcm}} for the case
#' of piecewise-constant mean signals.
#' @examples
#' single.cpt <- c(seq(0, 1999, 1), seq(1998, -1, -1))
#' single.cpt.student <- single.cpt + rt(4000, df = 5)
#' cpt.single <- ht_ID_cplm(single.cpt.student)
#'
#' three.cpt <- c(seq(0, 3998, 2), seq(3996, -2, -2), seq(0,3998,2), seq(3996,-2,-2))
#' three.cpt.student <- three.cpt + rt(8000, df = 5)
#' cpt.three <- ht_ID_cplm(three.cpt.student)
ht_ID_cplm <- function(x, s.ht = 3, q_ht = 300, ht_thr_id = 1.4,
                      ht_th_ic_id = 1.25, p_thr = 1, p_ic = 3) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (q_ht <= 0)){
    stop("The input for `q_ht' should be a positive number.")
  }
  if (abs(q_ht - round(q_ht)) > .Machine$double.eps ^ 0.5){
    warning("The input for `q_ht' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `q_ht'.")
  }
  if (length(x) <= q_ht) {
      result_new <- ID_cplm(x)
  } else {
      result_new <- list()
      nor <- normalise(x, sc = s.ht)
      result <- ID_cplm(x = nor, thr_id = ht_thr_id, th_ic_id = ht_th_ic_id,
                       pointsth = p_thr, pointsic = p_ic)
      r1 <- sort(result$cpt)
      r2 <- result$solution_path
      if (r1[1] == 0) {
          result_new$cpt <- 0
          result_new$no_cpt <- 0
      } else {
          result_new <- result
          result_new$cpt <- (r1 - 1) * (s.ht) + floor(s.ht / 2)
          result_new$no_cpt <- length(result_new$cpt)
      }
      result_new$fit <- est_signal(x, result_new$cpt, type = "slope")
      result_new$solution_path <- (r2 - 1) * (s.ht) + floor(s.ht / 2)
  }
  return(result_new)
}

#' Multiple change-point detection in piecewise-constant or continuous, piecewise-linear
#' signals using the Isolate-Detect methodology
#'
#' This is the main, general function of the package. It employs more specialised functions in
#' order to estimate the number and locations of multiple change-points in the noisy, piecewise-constant
#' or continuous, piecewise-linear input vector \code{xd}. The noise can either follow the Gaussian
#' distribution or not. The approach that is followed is a hybrid between the thresholding approach
#' (explained in \code{\link{pcm_th}} and \code{\link{cplm_th}}) and the information criterion approach
#' (explained in \code{\link{pcm_ic}} and \code{\link{cplm_ic}}) and estimates the change-points
#' taking into account both these approaches. Further to the number and the location of the estimated
#' change-points, \code{\link{ID}}, returns the estimated signal, as well as the solution path.
#' For more information and the relevant literature reference, see Details.
#'
#' @export
#' @param xd A numeric vector containing the data in which you would like to find
#'   change-points.
#' @param th.cons A positive real number with default value equal to 1. It is
#'   used to define the threshold, if the thresholding approach (explained in \code{\link{pcm_th}})
#'   is to be followed to detect the change-points in the scenario of piecewise-constant signals.
#' @param th.cons_lin A positive real number with default value equal to 1.4. It is
#'   used to define the threshold, if the thresholding approach (explained in \code{\link{cplm_th}})
#'   is to be followed to detect the change-points in the scenario of continuous, piecewise-linear signals.
#' @param th.ic A positive real number with default value equal to 0.9. It is
#'   useful only if the model selection based Isolate-Detect method (described in
#'   \code{\link{pcm_ic}}) is to be followed for the scenario of piecewise-constant signals.
#'   It is used to define the threshold value that will be used at the first step (change-point
#'   overestimation) of the model selection approach.
#' @param th.ic.lin A positive real number with default value equal to 1.25. It is
#'   useful only if the model selection based Isolate-Detect method (described in
#'   \code{\link{cplm_ic}}) is to be followed for the scenario of continuous, piecewise-linear signals.
#'   It is used to define the threshold value that will be used at the first step (change-point
#'   overestimation) of the model selection approach.
#' @param lambda A positive integer with default value equal to 3. It is used only
#'   when the threshold based approach is to be followed and it defines the distance
#'   between two consecutive end- or start-points of the right- or left-expanding intervals,
#'   respectively.
#' @param lambda.ic A positive integer with default value equal to 10. It is used only
#'   when the information criterion based approach is to be followed and it defines
#'   the distance between two consecutive end- or start-points of the right- or
#'   left-expanding intervals, respectively.
#' @param contrast A character string, which defines the type of the contrast function to
#'   be used in the Isolate-Detect algorithm. If \code{contrast = ``mean''}, then the algorithm
#'   looks for changes in a piecewise-constant signal. If \code{contrast = ``slope''},
#'   then the algorithm looks for changes in a continuous, piecewise-linear signal.
#' @param ht A logical variable with default value equal to \code{FALSE}. If \code{FALSE},
#'   the noise is assumed to follow the Gaussian distribution. If \code{TRUE}, then the
#'   noise is assumed to follow a distribution that has tails heavier than those of the
#'   Gaussian distribution.
#' @param scale A positive integer number with default value equal to 3. It is
#'   used to define the way we pre-average the given data sequence only if
#'   \code{ht = TRUE}. See the Details in \code{\link{ht_ID_pcm}} for more information on
#'   how we pre-average.
#' @details The data points provided in \code{xd} are assumed to follow \deqn{X_t = f_t + \sigma\epsilon_t; t = 1,2,...,T,}
#'   where \eqn{T} is the total length of the data sequence, \eqn{X_t} are the observed
#'   data, \eqn{f_t} is a one-dimensional, deterministic signal with abrupt structural
#'   changes at certain points, and \eqn{\epsilon_t} are independent and identically
#'   distributed random variables with mean zero and variance one. In this function,
#'   the following scenarios for \eqn{f_t} are implemented.
#'   \itemize{
#'   \item {Piecewise-constant signal with Gaussian noise.
#'
#'   Use \code{contrast = ``mean''} and \code{ht = FALSE} here.}
#'
#'   \item {Piecewise-constant signal with heavy-tailed noise.
#'
#'   Use \code{contrast = ``mean''} and \code{ht = TRUE} here.}
#'
#'   \item {Continuous, piecewise-linear signal with Gaussian noise.
#'
#'   Use \code{contrast = ``slope''} and \code{ht = FALSE} here.}
#'
#'   \item {Continuous, piecewise-linear signal with heavy-tailed noise.
#'
#'   Use \code{contrast = ``slope''} and \code{ht = TRUE} here.}
#'   }
#'   In the case where \code{ht = FALSE}: the function firstly detects the change-points using
#'   \code{\link{win_pcm_th}} (for the case of piecewise-constant signal) or \code{\link{win_cplm_th}}
#'   (for the case of continuous, piecewise-linear signal). If the estimated number of change-points
#'   is greater than 100, then the result is returned and we stop. Otherwise, \code{\link{ID}} proceeds
#'   to detect the change-points using \code{\link{pcm_ic}} (for the case of piecewise-constant signal)
#'   or \code{\link{cplm_ic}} (for the case of continuous, piecewise-linear signal) and this is what is
#'   returned.\cr
#'   In the case where \code{ht = TRUE}: First we pre-average the given data sequence using \code{\link{normalise}}
#'   and then, on the obtained data sequence, we follow exactly the same procedure as the one when \code{ht = FALSE}
#'   above.\cr
#'   More details can be found in ``Detecting multiple generalized change-points by isolating single ones'',
#'   Anastasiou and Fryzlewicz (2018), preprint.
#' @return
#'   A list with the following components:
#'   \tabular{ll}{
#'    \cr \code{cpt} \tab A vector with the detected change-points.
#'    \cr \code{no_cpt} \tab The number of change-points detected.
#'    \cr \code{fit} \tab A numeric vector with the estimated signal.
#'    \cr \code{solution_path} \tab A vector containing the solution path.
#'  }
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}
#' @seealso \code{\link{ID_pcm}}, \code{\link{ID_cplm}}, \code{\link{ht_ID_pcm}}, and
#' \code{\link{ht_ID_cplm}}, which are the functions that are employed
#' in \code{\link{ID}}, depending on which scenario is imposed by the input arguments.
#' @examples
#' single.cpt.mean <- c(rep(4,3000),rep(0,3000))
#' single.cpt.mean.normal <- single.cpt.mean + rnorm(6000)
#' single.cpt.mean.student <- single.cpt.mean + rt(6000, df = 5)
#' cpt.single.mean.normal <- ID(single.cpt.mean.normal)
#' cpt.single.mean.student <- ID(single.cpt.mean.student, ht = TRUE)
#'
#' single.cpt.slope <- c(seq(0, 1999, 1), seq(1998, -1, -1))
#' single.cpt.slope.normal <- single.cpt.slope + rnorm(4000)
#' single.cpt.slope.student <- single.cpt.slope + rt(4000, df = 5)
#' cpt.single.slope.normal <- ID(single.cpt.slope.normal, contrast = "slope")
#' cpt.single.slope.student <- ID(single.cpt.slope.student, contrast = "slope", ht = TRUE)
ID <- function(xd, th.cons = 1, th.cons_lin = 1.4, th.ic = 0.9,
               th.ic.lin = 1.25, lambda = 3, lambda.ic = 10,
               contrast = c("mean", "slope"), ht = FALSE,
               scale = 3) {
    if (missing(contrast)){
      contrast <- "mean"
    }
    if (contrast == "mean") {
        if (ht == FALSE) {
            result <- ID_pcm(x = xd, thr_id = th.cons, th_ic_id = th.ic,
                             pointsth = lambda, pointsic = lambda.ic)
        } else {
            result <- ht_ID_pcm(x = xd, s.ht = scale,
                                p_thr = max(1, floor(lambda / scale)),
                                p_ic = max(1, floor(lambda.ic / scale)),
                                ht_thr_id = th.cons, ht_th_ic_id = th.ic)
        }
    }
    if (contrast == "slope") {
        if (ht == FALSE) {
            result <- ID_cplm(x = xd, thr_id = th.cons_lin,
                             th_ic_id = th.ic.lin, pointsth = lambda,
                             pointsic = lambda.ic)
        } else {
            result <- ht_ID_cplm(x = xd, s.ht = scale,
                                p_thr = max(1, floor(lambda / scale)),
                                p_ic = max(1, floor(lambda.ic / scale)),
                                ht_thr_id = th.cons_lin, ht_th_ic_id = th.ic.lin)
        }
    }
    return(result)
}
