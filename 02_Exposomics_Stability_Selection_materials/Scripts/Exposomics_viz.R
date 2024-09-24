#### LongITools Rome Showcase Utensils ####


#### plotting libraries ####

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(randomcoloR))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(aricode))

#### Custom plotting functions - Stability Selection ####

plot.selection_performance <- function(x,
                                       col = c("#6bcbb8","#ff0000", "#ffc500","#ff5c02"),
                                       col.axis = NULL,
                                       col.thr = "#A0AAA3",
                                       lty.thr = 2,
                                       n_predictors = NULL,
                                       theta = theta,
                                       ...) {
  
  
  # Storing extra arguments
  extra_args <- list(...)
  
  # Defining default settings
  if (!"type" %in% names(extra_args)) {
    extra_args$type <- "h"
  }
  if (!"xlab" %in% names(extra_args)) {
    extra_args$xlab <- ""
  }
  if (!"ylab" %in% names(extra_args)) {
    extra_args$ylab <- "Selection proportion"
  }
  if (!"las" %in% names(extra_args)) {
    extra_args$las <- 2
  }
  if (!"ylim" %in% names(extra_args)) {
    extra_args$ylim <- c(0, 1)
  }
  if (!"cex.lab" %in% names(extra_args)) {
    extra_args$cex.lab <- 1.2
  }
  if (!"cex.axis" %in% names(extra_args)) {
    extra_args$cex.axis <- 1
  }
  
  # Checking inputs
  if (length(col) != 4) {
    col <- rep(col[1], )
  }
  if (is.null(col.axis)) {
    col.axis <- col
  } else {
    if (length(col.axis) != 4) {
      col.axis <- rep(col.axis[1], 4)
    }
  }
  
  
  # Extracting selection proportions
  selprop <- SelectionProportions(x)
  
  # Defining colours
  mycolours <- ifelse(SelectedVariables(x) == 1, 
                      yes = ifelse(SelectedVariables(x) == theta, yes = col[1], no = col[2]),
                      no = ifelse(SelectedVariables(x) == theta, yes = col[3], no = col[4]))
  mycolours_axis <- mycolours
  
  # Re-ordering by decreasing selection proportion
  ids <- sort.list(selprop, decreasing = TRUE)
  if (is.null(n_predictors)) {
    n_predictors <- length(ids)
  }
  n_predictors <- min(n_predictors, length(ids))
  ids <- ids[1:n_predictors]
  selprop <- selprop[ids]
  mycolours <- mycolours[ids]
  mycolours_axis <- mycolours_axis[ids]
  
  # Extracting relevant extra arguments
  tmp_extra_args <- extra_args
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("xaxt", "col")]
  
  # Making plot
  do.call(base::plot, args = c(
    list(
      x = selprop,
      xaxt = "n",
      col = mycolours
    ),
    tmp_extra_args
  ))
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = graphics::axis)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("side", "at", "labels", "col.axis", "col.ticks", "type")]
  
  # Adding axis
  for (i in 1:length(selprop)) {
    do.call(graphics::axis, args = c(
      list(
        side = 1,
        at = i,
        labels = names(selprop)[i],
        col.axis = mycolours_axis[i],
        col.ticks = mycolours_axis[i]
      ),
      tmp_extra_args
    ))
  }
  graphics::abline(h = Argmax(x)[1, 2], lty = lty.thr, col = col.thr)
  
  # Extracting performance measures
  perf <- SelectionPerformance(x, theta)
  
  # Adding legend 
  
  legend(x = "right",
         lty = c(4,6),
         col= append(col,col.thr), 
         box.lty= 0,
         cex = 0.6,
         legend=c(paste0("True Positives N(",perf$TP,")"), 
                  paste0("False Positives N(",perf$FP,")"),
                  paste0("True Negatives N(",perf$TN,")"),
                  paste0("False Negatives N(",perf$FN,")"),
                  expression(hat(pi))))
}

plot.stab_paths <- function(stab_int, sim_int){
  
  stab_int <- stab
  
  sim_int <- simul
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # create dataframe of stability paths
  
  Pi <- stab_int$P
  
  Lambda <- stab_int$Lambda
  
  index_lambda <- 1:length(stab_int$Lambda)
  
  selprop <- t(sapply(index_lambda, FUN = function(x){SelectionProportions(stab_int, x)}))
  
  theta <- SelectedVariables(stab_int)
  
  theta_star <- sim_int$theta
  
  beta_true <- sim_int$beta
  
  stab_paths <- as.data.frame(cbind(Lambda, selprop))
  
  # define margins and layout
  
  par(oma = c(0.5,0.5,0.5,0.5))
  par(mar=c(4,4,4,8), mfrow = c(1,1), las = 1)
  
  # empty plot 
  
  plot(NA,
       xlim = (c(0, max(Lambda, na.rm = T))),
       ylim = c(0, 1),
       frame.plot =  FALSE, 
       axes = T,
       xlab = NA,
       ylab = NA,
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  # add axes
  
  axis(side = 1, at = append(seq(0,max(Lambda, na.rm = T), by = 1),ceiling(max(Lambda, na.rm = T))))
  axis(side = 2, at = seq(0,1, by = 0.1), las=1)
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = "Selection proportion as a function of lambda")
  
  title( 
    cex.lab=0.9, 
    ylab = "Selection proportion")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = "Lambda",)
  
  # add abline for calibrated parameters 
  abline(v = Argmax(stab_int)[1], lty = 2, col = addTrans("#A0AAA3", 255))
  abline(h = Argmax(stab_int)[2], lty = 2, col = addTrans("#A0AAA3", 255))
  lines(Lambda, Pi, col = addTrans("#A0AAA3", 150))
  
  # add stability paths 
  for (i in 2:ncol(stab_paths)){
    if (theta[colnames(stab_paths)[i]] == 1 & theta[colnames(stab_paths)[i]] == theta_star[colnames(stab_paths)[i],]){
      lines(stab_paths$V1, stab_paths[,i], col = addTrans("#6bcbb8", trans =  255), lwd = beta_true[colnames(stab_paths)[i],] + 1)
      text(x = Lambda[max(which(unique(stab_paths[,i]) != 1))], y = stab_paths[,i][max(which(unique(stab_paths[,i]) != 1))], labels =  colnames(stab_paths)[i],
           cex=0.65, col = addTrans("#6bcbb8", trans =  255))
    } else if(theta[colnames(stab_paths)[i]] == 0 & theta[colnames(stab_paths)[i]] == theta_star[colnames(stab_paths)[i],]){
      lines(stab_paths$V1, stab_paths[,i], col = addTrans("#ffc500", trans = 70), lty = 2, lwd = beta_true[colnames(stab_paths)[i],] + 1)
    } else if(theta[colnames(stab_paths)[i]] == 1 & theta[colnames(stab_paths)[i]] != theta_star[colnames(stab_paths)[i],]){
      lines(stab_paths$V1, stab_paths[,i], col = addTrans("#ff0000", trans = 255))
      text(x = Lambda[max(which(unique(stab_paths[,i]) != 1))], y = stab_paths[,i][max(which(unique(stab_paths[,i]) != 1))], labels =  colnames(stab_paths)[i],
           cex=0.65, col = addTrans("#ff0000", trans =  255))
    } else{ 
      lines(stab_paths$V1, stab_paths[,i], col = addTrans("#ff5c02", trans = 255), lty = 2, lwd = beta_true[colnames(stab_paths)[i],] + 1)
      text(x = Lambda[max(which(unique(stab_paths[,i]) != 1))], y = stab_paths[,i][max(which(unique(stab_paths[,i]) != 1))], labels =  colnames(stab_paths)[i],
           cex=0.65, col = addTrans("#ff5c02", trans =  255))
    }
  }
  
  # add gridlines 
  
  grid(ny = 10)
  
  # add stability as background?
  
  legend(x = "topright",
         legend = c("True Positives", "False Positives", "True Negatives", "False Negatives", substitute(hat(pi)), substitute(paste(hat(lambda), " , ", hat(pi)))),
         col = c("#6bcbb8","#ff0000", "#ffc500","#ff5c02", "#A0AAA3", "#A0AAA3"),
         lty=c(1,1,2,2,1,2),
         cex=0.5,
         inset = c(-0.1, 0),
         horiz= F,
         xpd = T,
         bty = "n",
         seg.len= 0.5,
         x.intersp = 0.5,
         y.intersp = 1.5)
  
  # line weight
  
}

#### Custom plotting functions - Clustering ####

rhc <- function (tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2, 
                 cluster = NULL) 
{
  if (length(h) > 1L | length(k) > 1L) 
    stop("'k' and 'h' must be a scalar")
  if (!is.null(h)) {
    if (!is.null(k)) 
      stop("specify exactly one of 'k' and 'h'")
    k <- min(which(rev(tree$height) < h))
    k <- max(k, 2)
  }
  else if (is.null(k)) 
    stop("specify exactly one of 'k' and 'h'")
  if (k < 2 | k > length(tree$height)) 
    stop(gettextf("k must be between 2 and %d", length(tree$height)), 
         domain = NA)
  if (is.null(cluster)) 
    cluster <- cutree(tree, k = k)
  clustab <- table(cluster)[unique(cluster[tree$order])]
  m <- c(0, cumsum(clustab))
  if (!is.null(x)) {
    if (!is.null(which)) 
      stop("specify exactly one of 'which' and 'x'")
    which <- x
    for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
  }
  else if (is.null(which)) 
    which <- 1L:k
  if (any(which > k)) 
    stop(gettextf("all elements of 'which' must be between 1 and %d", 
                  k), domain = NA)
  border <- rep_len(border, length(which))
  retval <- list()
  for (n in seq_along(which)) {
    rect(
      ybottom = m[which[n]] + 0.66,
      xright = par("usr")[3L],
      ytop = m[which[n] + 1] + 0.33,
      xleft = mean(rev(tree$height)[(k - 1):k]),
      border = border[n])
    retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
  }
  invisible(retval)
}

hclust_dist_plot <- function(myhclust){
  par(mar = c(4, 4, 4, 4))
  plot(as.dendrogram(myhclust),
       xaxt = "n", leaflab = "none", axes = FALSE,
       xlab = "", sub = "", main = "", ylab = "", horiz = T
  )
  axis(side = 3, at = seq(0, max(myhclust$height), by = 2), labels = NA)
  par(xpd = TRUE)
  text(seq(0, max(myhclust$height), by = 2),
       length(myhclust$order) * 1.11,
       seq(0, max(myhclust$height), by = 2), 
       cex = 0.7)
  text(mean(c(0, max(myhclust$height))),
       length(myhclust$order) * 1.15, 
       "Euclidean distance")
  rhc(myhclust, k = 3, border = "navy")
}

shclust_dist_plot <- function(shclust){
  par(mar = c(4, 4, 4, 2))
  plot_horiz.dendrogram(as.dendrogram(shclust), xaxt = "n", leaflab = "none", axes = FALSE,
                       sub = "", main = "")
  axis(side = 3, at = seq(0, 1, by = 0.2), labels = NA)
  par(xpd = TRUE)
  text(seq(0, 1, by = 0.2),
       111,
       seq(0, 1, by = 0.2), 
       cex = 0.7)
  text(0.5,
       115, 
       "1 - co-membership proportion")
  rhc(shclust, k = hat_N, border = "darkred")
}



# Comparison
hclust_comparison_plot <- function(myhclust, shclust){
  A <- matrix(0, nrow = length(myhclust$order), ncol = length(myhclust$order))
  rownames(A) <- myhclust$order
  colnames(A) <- rev(shclust$order)
  for (i in 1:nrow(A)) {
    hm <- cutree(myhclust, k = hat_N)
    hm <- names(hm[which(hm == hm[paste0("obs", i)])])
    
    shm <- cutree(shclust, k = hat_N)
    shm <- names(shm[which(shm == shm[paste0("obs", i)])])
    overlap <- intersect(hm, shm)
    if (length(overlap) / length(union(hm, shm)) > 0.5) {
      A[as.character(i), as.character(i)] <- 1
    } else {
      A[as.character(i), as.character(i)] <- 2
    }
  }
  rownames(A) <- paste0("h", rownames(A))
  colnames(A) <- paste0("sh", colnames(A))
  Abig <- Square(A)
  
  mycolours <- lighten(c("skyblue", "red", "forestgreen"), amount = 0.5)
  myedgecolours <- darken(c("lightgrey", "red"), amount = 0.1)
  g <- Graph(Abig, weighted = TRUE)
  V(g)$label <- gsub("h", "", gsub("sh", "", V(g)$label))
  V(g)$color <- mycolours[simul$theta[as.numeric(V(g)$label)]]
  V(g)$size <- 6
  V(g)$frame.color <- V(g)$color
  V(g)$label.color <- darken(V(g)$color, amount = 0.5)
  E(g)$color <- myedgecolours[E(g)$width]
  mylayout <- matrix(c(
    rep(0, length(myhclust$order)), rep(1, length(myhclust$order)),
    1:length(myhclust$order), rev(1:length(myhclust$order))
  ), ncol = 2)
  
  myasp <- 0.54
  {
    par(mar = rep(0, 4))
    plot(g, layout = mylayout, asp = 1 / myasp)
  }
}

plot_clustering_comparison <- function(myhclust, shclust){
  par(mfrow=c(1,3), xpd = T)
  hclust_dist_plot(myhclust)
  hclust_comparison_plot(myhclust, shclust)
  shclust_dist_plot(shclust)
}

### Cross-validation CalibrationPlot ----

plot.cv_calibration <- function(cv, model = "glinternet"){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # define variables
  
  if(model == "glinternet"){
    Lambda = cv$lambda
    cvErr = cv$cvErr
    cvErrStd = cv$cvErrStd
    cvErrmax = cvErr + cvErrStd
    cvErrmin = cvErr - cvErrStd
    lambdahat = cv$lambdaHat
    lambdahat1Std = cv$lambdaHat1Std
  } else if(model == "lasso"){
    Lambda = cv$lambda
    cvErr = cv$cvm
    cvErrStd = cv$cvsd
    cvErrmax = cv$cvup
    cvErrmin = cv$cvlo
    lambdahat = cv$lambda.min
    lambdahat1Std = cv$lambda.1se
  }
  
  
  # define margins and layout
  
  par(oma = c(1,1,1,1))
  par(mar=c(4,4,4,8), mfrow = c(1,1), las = 1)
  
  # empty plot 
  
  plot(NA,
       xlim = c(min(log(Lambda), na.rm = T), max(log(Lambda), na.rm = T)),
       ylim = c(0, max(cvErrmax)),
       frame.plot =  FALSE, 
       axes = T,
       xlab = NA,
       ylab = NA,
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  # add axes
  
  axis(side = 1, at = seq(floor(round(min(log(Lambda)), digits = 2)), ceiling(max(log(Lambda))), by = 1))
  axis(side = 2, at = seq(0,max(cvErrmax), by = 0.5))
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression(paste("Cross-validation Error as a Function of ", lambda)))
  
  title( 
    cex.lab=0.9, 
    ylab = "Cross-validation Error")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = expression(paste("log(",lambda,")")))
  
  # add abline for calibrated parameters
  abline(v = log(lambdahat), lty = 2, col = addTrans("#ff0000", 255))
  abline(v = log(lambdahat1Std), lty = 2, col = addTrans("#ffc500", 255))
  
  # plot cross-validation error
  lines(log(Lambda), cvErr)
  arrows(x0=log(Lambda), y0=cvErrmin, x1=log(Lambda), y1=cvErrmax, code=3, angle=90, col = addTrans("#A0AAA3", 150), length=0.05)
  
  # add legend
  
  legend(x = "topright",
         legend = c(substitute(paste(hat(lambda))), expression(paste(hat(lambda), " 1",sigma))),
         col = c("#ff0000", "#ffc500"),
         lty=c(2,2),
         cex=0.8,
         inset = c(-0.2, 0),
         horiz= F,
         xpd = T,
         bty = "n",
         seg.len= 1,
         x.intersp = 0.5,
         y.intersp = 1.5)
}

### Regularization Paths ----

plot.reg_paths <- function(reg_beta = NULL,
                           cv = cv,
                           calibration = "lambdahat_sd",
                           sim_int = sim_int){
  
  # custom functions
  addTrans <- function(color,trans)
  {
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
      hex <- unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
  }
  
  # retrieve plotting items
  
  if (is.null(reg_beta)){
    
    if(calibration == "lambdahat"){
      
      coef_lambdahat <- coef(cv, s = "lambda.min")
      
      theta <- ifelse(as.matrix(coef_lambdahat) == 0, 0, 1)[-1,]
      
      lambda_calibrate <- cv$lambda.min
      
    } else if (calibration == "lambdahat_sd"){
      
      coef_lambdahat_sd <- coef(cv, s = "lambda.1se")
      
      theta <- ifelse(as.matrix(coef_lambdahat_sd) == 0, 0, 1)[-1,]
      
      lambda_calibrate <- cv$lambda.1se
      
    }
    
    Lambda <- cv$lambda
    
    beta <- t(as.matrix(coef(cv, s = cv$lambda)))[,-1]
    
  } else {
    
    numLevels <- apply(sim_int$xdata,2, function(x) {
      if (is.factor(x))
        length(levels(x))
      else {1}})
    
    if(calibration == "lambdahat"){
      
      lambda_calibrate <- cv_int$lambdaHat
      
      glinternet_lambdahat <- glinternet::glinternet(sim_int$xdata, sim_int$ydata, lambda = cv_int$lambdaHat, numLevels = numLevels)
      
      theta <- glinternet.theta(mymodel = glinternet_lambdahat, xdata = sim_int$xdata)
      
    } else if (calibration == "lambdahat_sd"){
      
      lambda_calibrate <- cv_int$lambdaHat1Std
      
      glinternet_lambdahat_sd<- glinternet::glinternet(sim_int$xdata, sim_int$ydata, lambda = cv_int$lambdaHat1Std, numLevels = numLevels)
      
      theta <- glinternet.theta(mymodel = glinternet_lambdahat_sd, xdata = sim_int$xdata)
      
    }
    
    Lambda <- cv_int$lambda
    
    beta <- reg_beta
    
  }
  
  theta_star <- sim_int$theta
  
  beta_true <- sim_int$beta
  
  # define margins and layout
  
  par(oma = c(2,2,2,2))
  par(mar=c(4,4,4,10), mfrow = c(1,1), las = 1)
  
  # empty plot 
  
  plot(NA,
       xlim = c(min(log(Lambda)), ceiling(max(log(Lambda), na.rm = T))),
       ylim = c(floor((min(beta))),ceiling(max(beta))),
       frame.plot =  FALSE, 
       axes = T,
       xlab = NA,
       ylab = NA,
       xaxt = "n",
       yaxt = "n",
       xaxs = "i",
       yaxs = "i")
  
  # add axes
  
  axis(side = 1, at = seq(floor(round(min(log(Lambda)), digits = 2)), ceiling(max(log(Lambda))), by = 1))
  axis(side = 2, at = seq(floor((min(beta))),ceiling(max(beta)), by = 0.2), las=1)
  
  # Add titles
  
  title(
    adj = 0,
    line = 2,
    cex.axis=0.5, 
    cex.main=1, 
    main = expression(paste("Coefficients as a Function of ", lambda)))
  
  title( 
    cex.lab=0.9, 
    ylab = "Coefficients")
  
  title(
    cex.lab=0.9,
    adj = 1,
    xlab = expression(paste("log(",lambda,")")))
  
  # add abline for calibrated parameters 
  abline(v = log(lambda_calibrate), lty = 2, col = addTrans("#A0AAA3", 255))
  abline(h = 0, lty = 1, col = addTrans("black", 255))
  
  # add stability paths 
  for (i in 1:ncol(beta)){
    if (theta[colnames(beta)[i]] == 1 & theta[colnames(beta)[i]] == theta_star[colnames(beta)[i],]){
      lines(log(Lambda), beta[,i], col = addTrans("#6bcbb8", trans =  255), lwd = beta_true[colnames(beta)[i],] + 1)
      text(x = log(Lambda[max(which(unique(beta[,i]) != 1))]) + 0.5, y = beta[,i][max(which(unique(beta[,i]) != 1))], labels =  colnames(beta)[i],
           cex=0.65, col = addTrans("#6bcbb8", trans =  255))
    } else if(theta[colnames(beta)[i]] == 0 & theta[colnames(beta)[i]] == theta_star[colnames(beta)[i],]){
      lines(log(Lambda), beta[,i], col = addTrans("#ffc500", trans = 150), lty = 2, lwd = beta_true[colnames(beta)[i],] + 1)
    } else if(theta[colnames(beta)[i]] == 1 & theta[colnames(beta)[i]] != theta_star[colnames(beta)[i],]){
      lines(log(Lambda), beta[,i], col = addTrans("#ff0000", trans = 150))
    } else{ 
      lines(log(Lambda), beta[,i], col = addTrans("#ff5c02", trans = 255), lty = 2, lwd = beta_true[colnames(beta)[i],] + 1)
      text(x = log(Lambda[max(which(unique(beta[,i]) != 1))])  + 0.5, y = beta[,i][max(which(unique(beta[,i]) != 1))], labels =  colnames(beta)[i],
           cex=0.65, col = addTrans("#ff5c02", trans =  255))
    }
  }
  
  # add gridlines 
  
  grid(ny = 10)
  
  # add legend
  
  if (calibration == "lambdahat"){
    calibration_legend <- substitute(paste(hat(lambda)))
  } else if (calibration == "lambdahat_sd"){
    calibration_legend <- expression(paste(hat(lambda), " 1",sigma))
  }
  
  legend(x = "topright",
         legend = c("True Positives", "False Positives", "True Negatives", "False Negatives", calibration_legend),
         col = c("#6bcbb8","#ff0000", "#ffc500","#ff5c02", "#A0AAA3"),
         lty=c(1,1,2,2,2),
         cex=1,
         inset = c(-0.2, 0),
         horiz= F,
         xpd = T,
         bty = "n",
         seg.len= 1,
         x.intersp = 0.5,
         y.intersp = 1.5)
  
}




