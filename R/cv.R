## ================
## Cross Validation
## ================

#' Create folds for Crossfold Validation bu consecutive intervals
#'
#' @param folds Number of folds.
#' @param intervalSize The size of the consecutive intervals.
#' @param tileSize The size of one tile.
#' @return A list of as many elements as there are folds. Each element
#' contains the rows that are part of this fold.
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.leaveOutConsecutiveIntervals <- function(folds, intervalSize, tileSize){
    if(intervalSize >= 50) {
        warning("CV using large interval may only make sense if all experiments have replicates")
    }
 
    ## number of intervals
    n = round(tileSize/intervalSize)
    ## break points
    bps = round(seq(1, tileSize + 1, length.out = n + 1))
    leftout_intervals  = split(sample(n),rep(1:folds, length = n))
    res = lapply(
        leftout_intervals, 
        function(lo){
            ## we return the indices in increasing order to preserve
            ## as much as possible the Genomic ranges for subsequent subsetting
            sort(unlist(lapply(lo, function(i) bps[i]:(bps[i+1]-1))))
        }
        )
    return(res)
}

#' the function to be maximised for Cross Validation
#' Procedure:
#' 1. Fit model on all but the left out intervalls
#' 2. compute log-Likelihood on predicted values vs. real values
#' 3. compute the mean and return
#' @param pars The parameters to be optimized.
#' @param gtiles The datasets the likelihood is computed on.
#' @param CV_intervals A list of folds containing the rows to be left out for each fold.
#' @param formula A formula object.
#' @param experimentMatrix The experiment matrix with columnnames being the tracks and
#' rownames being the samples.
#' @param chunkCoords A GRanges object representing the chunk coordinates.
#' @param fixedpars A list of two elements - lambda and theta - being specified. If not
#' they are estimated.
#' @param verbose Print some more text?
#' @param logging Activate logging?
#' @param ... Additional parameters forwarded to mgcv::gam
#' @return A single numeric value - the likelihood for this combination of parameters.
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.loglik <- function(pars, gtiles, CV_intervals, formula,
                    fixedpars = list(lambda = NULL, theta = NULL), ...){
    coords <- attr(gtiles, "settings")$chunks
    if(is.null(fixedpars$lambda)) {
        fixedpars$lambda <- exp(pars[["lambda"]])
    }
    if(is.null(fixedpars$theta)) {
        fixedpars$theta <- exp(pars["theta"])
    }
    
    fullpred <- lapply(1:length(gtiles), function(y) {
        rep(NA, nrow(gtiles[[1]]))
    })
    names(fullpred) <- names(gtiles)
    ids <- expand.grid(folds = 1:length(CV_intervals), tiles = 1:length(gtiles))

    lambdaFun <- function(iter, ids, gtiles, CV_intervals, formula, 
                          fixedpars) {
        suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))
        id <- ids[iter,]
        
        testset <- gtiles[[id$tiles]][CV_intervals[[id$folds]],]
        trainset <- gtiles[[id$tiles]][-CV_intervals[[id$folds]],]
                    
        nvars <- length(.getVars(formula)) - 1
        ## the unnaming of theta is needed since it's a named vector and nb() transforms it soemhow into a list within when the name is not dropped
        mod <- mgcv::gam(formula, data = trainset, family = mgcv::nb(theta = unname(fixedpars$theta)), method = "REML",
                         sp = rep(fixedpars$lambda, nvars), ...)
            
        pred <- mgcv::predict.gam(mod, newdata = testset, type="response")
        return(pred)
    }
    
    if(fixedpars$theta < 1e-3) {
        fixedpars$theta <- 1e-3
    }
    
    cvs <- bplapply(1:nrow(ids), lambdaFun, ids = ids, gtiles = gtiles,
                    CV_intervals, formula = formula, fixedpars = fixedpars)

    for(ii in 1:length(cvs)) {
        id <- ids[ii,]
        fullpred[[id$tiles]][CV_intervals[[id$folds]]] <- cvs[[ii]]
    }

    res <- lapply(1:length(fullpred), function(y) {
        dens <- dnbinom(gtiles[[y]]$value, size = unname(fixedpars$theta), mu = fullpred[[y]], log = TRUE)
        gtiles[[y]]$yhat <- fullpred[[y]]
        gtiles[[y]]$ll <- dens
        return(gtiles[[y]])
    })
    names(res) <- names(fullpred)

    #############################################################
    ## out-of-sample log-likelihood of the fit on the chunk part of a tile
    
    ## ll: sum CV log-lik by regions and parameter combination
    ## we take the average i.e we assume all regions have the same length.
    ll <- mean(sapply(1:length(res), function(y) {
        id <- names(res)[y]
        chunks <- c(start(coords[mcols(coords)$id == as.integer(id),]),
                    end(coords[mcols(coords)$id == as.integer(id),]))
        sum(.untile(res[[id]], chunks)$ll, na.rm = TRUE)
    }))
    return(ll)
}

#' Main function to perform cross validation
#'
#' @param par The parameters to be optimized.
#' @param fn A function to be used for optimization.
#' @param data The datasets the likelihood is computed on.
#' @param formula A formula object.
#' @param experimentMatrix The experiment matrix with columnnames being the tracks and
#' rownames being the samples.
#' @param folds Number of k-folds.
#' @param intervalSize The size of the consecutive intervals.
#' @param fixedpars A list of two elements - lambda and theta - being specified. If not
#' they are estimated.
#' @param colData The sample specific values, that should be added to the data.
#' @param verbose Print some more text?
#' @param logging Activate logging?
.doCrossValidation <- function(par, fn, data, formula, folds, intervallSize,
                               fixedpars, experimentDesign, sf,
                               method = "Nelder-Mead",
                               control = list(maxit=100, fnscale=-1), ...) {
    
    settings <- metadata(data)
    meltedData <- .meltGTile(data, experimentDesign, sf, formula)
    attr(meltedData, "settings") <- settings
    cvint <- .leaveOutConsecutiveIntervals(folds, intervallSize, nrow(meltedData[[1]]))
    
    pars <- optim(par, fn, gtiles = meltedData, CV_intervals = cvint,
                  formula = formula, 
                  fixedpars = fixedpars,
                  method = method, control = control, ...)
    params <- exp(pars$par)
    
    if(length(params) == 1) {
        fixedpars[sapply(fixedpars, is.null)] <- params
        params <- unlist(fixedpars)
    }
    return(params)
}

## #' Shiny application to check fit of estimated lambda and theta
## #'
## #' A shiny application to check the fit of the estimated lambda and theta
## #'
## #' @param obj The return object of \code{\link{estimateLambda}}
## #' @return Starts a shiny application in webbrowser
## #' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
## #' @export
## plotCV <- function(obj){
##     if(!requireNamespace("shiny", quietly = TRUE)) {
##             stop("This is a shinyApp. Please install shiny.", call. = FALSE)
##         }
##     require(shiny)
##     lambdaTemp <- attr(obj,"lambdaTemp")
##     dataTemp <- attr(obj,"dataTemp")
##     bestLambdaTemp <- attr(obj,"bestLambdaTemp")

##     lambdaVector <- read.table(lambdaTemp, header = TRUE)[,1]
##     numLambda <- lambdaVector[1]
##     lambdaVector <- lambdaVector[-1]
##     bestLambda <- read.table(bestLambdaTemp, header = TRUE)[,1]
##     ll <- format(lambdaVector,scientific=TRUE,digits=3)
##     dr <- as.list(dataTemp)
##     data <- fread(as.character(dr[1]), header = TRUE)
##     data[,sample := as.factor(sample)]
##     tracks <- levels(data$sample)
##     ## regionNames <- saply(1:length(dr), function(y) paste("region", y, sep = ""))
##     ## names(dr) <- regionNames
##     if (numLambda > 1) {
##         chosenLambdas <- paste("Best Lambda chosen as", bestLambda[1])
##         for (ii in 2:numLambda) chosenLambdas <- paste(chosenLambdas, "and", bestLambda[ii])
##         app <- shinyApp(ui = shinyUI(fluidPage(
##                             titlePanel("Lambda debug plot"),
##                             sidebarLayout(position = "left",
##                                           sidebarPanel(
##                                               helpText(chosenLambdas),
##                                               sliderInput("lambda",
##                                                           label = "Lambda value:",
##                                                           min = 1, max = length(lambdaVector), value = which(lambdaVector == bestLambda[1]), step = 1,
##                                                           ##                              ticks = ll[seq(1,length(lambdaVector),length.out=10)],
##                                                           ticks = TRUE,
##                                                           width = '1000px'
##                                                           ),
##                                               sliderInput("lambda2",
##                                                           label = "Lambda value:",
##                                                           min = 1, max = length(lambdaVector), value = which(lambdaVector == bestLambda[2]), step = 1,
##                                                           ##                              ticks = ll[seq(1,length(lambdaVector),length.out=10)],
##                                                           ticks = TRUE,
##                                                           width = '1000px'
##                                                           ),
##                                               br(),
##                                               selectInput("region", 
##                                                           label = "Choose a region to display",
##                                                           choices = dr,
##                                                           selected = "region1"),
##                                               br(),
##                                               selectInput("inputTrack",
##                                                           label = "First track",
##                                                           choices = tracks,
##                                                           selected = tracks[1]),
##                                               br(),
##                                               selectInput("ipTrack",
##                                                           label = "Second track",
##                                                           choices = tracks,
##                                                           selected = tracks[2])),
##                                           mainPanel(
##                                               tableOutput("value1"),
##                                               tableOutput("value2"),
##                                               plotOutput("inputPlot"),
##                                               plotOutput("ipPlot")
##                                               )
##                                           )
##                             )),
                        
##                         server = shinyServer(function(input, output, session) {
##                             is.installed <- require(ggplot2)
##                             if (!is.installed) {
##                                 plotData <- function(data,title) {
##                                     plot(data$x,data$counts, col = "darkgray", xlab = "position", ylab = "counts and estimates", main = title)
##                                     lines(data$x,data$pred, col = "blue")
##                                 }
##                             }
##                             else {
##                                 plotData <- function(df,title) {
##                                     ggplot(data = df, aes(x=x,y=counts)) + geom_point(color = "darkgray") + 
##                                         geom_line(aes(x = x,y = pred), color = "blue") + ggtitle(title)
##                                 }
##                             }
##                             myLambda <- reactive({
##                                 as.character(lambdaVector[input$lambda])
##                             })
                            
##                             myLambda2 <- reactive({
##                                 as.character(lambdaVector[input$lambda2])
##                             })
                            
##                             dataset <- reactive({
##                                 df <- fread(as.character(input$region), header = TRUE)
##                                 df
##                             })

##                             track1 <- reactive({
##                                 as.character(input$inputTrack)
##                             })

##                             track2 <- reactive({
##                                 as.character(input$ipTrack)
##                             })
                            
##                             output$value1 <- renderTable({
##                                 data.frame(
##                                     Name = c("First Lambda Value", "Second Lambda Value"),
##                                     Value = c(as.character(c(myLambda())), as.character(c(myLambda2()))),
##                                     stringsAsFactors=FALSE
##                                     )
##                             })

##                             ## output$value2 <- renderTable({
##                             ##     data.frame(
##                             ##         Name = c("Lambda Value"),
##                             ##         Value = as.character(c(myLambda2())),
##                             ##         stringsAsFactors=FALSE
##                             ##         )
##                             ## })
                            
##                             ## only works with one lambda for now
##                             output$inputPlot <- renderPlot({
##                                 df <- dataset()
##                                 name <- track1()
##                                 dfInput <- df[sample == name & lambda1 == myLambda(),]
##                                 plotData(dfInput,name)
##                             })
                            
##                             output$ipPlot <- renderPlot({
##                                 df <- dataset()
##                                 name <- track2()
##                                 dfIP <- df[sample == name & lambda1 == myLambda2(),]
##                                 plotData(dfIP,name)
##                             })

##                             session$onSessionEnded(function() {
##                                 stopApp()
##                             })
##                         }))
##     }

##     else {

##         app <- shinyApp(ui = shinyUI(fluidPage(
##                             titlePanel("Lambda debug plot"),
##                             sidebarLayout(position = "left",
##                                           sidebarPanel(
##                                               helpText(paste("Best Lambda chosen as",bestLambda)),
##                                               sliderInput("lambda",
##                                                           label = "Lambda value:",
##                                                           min = 1, max = length(lambdaVector), value = which(lambdaVector == bestLambda[1]), step = 1,
##                                                           ##                              ticks = ll[seq(1,length(lambdaVector),length.out=10)],
##                                                           ticks = TRUE,
##                                                           width = '1000px'
##                                                           ),
##                                               br(),
##                                               selectInput("region", 
##                                                           label = "Choose a region to display",
##                                                           choices = dr,
##                                                           selected = "region1"),
##                                               br(),
##                                               selectInput("inputTrack",
##                                                           label = "First track",
##                                                           choices = tracks,
##                                                           selected = tracks[1]),
##                                               br(),
##                                               selectInput("ipTrack",
##                                                           label = "Second track",
##                                                           choices = tracks,
##                                                           selected = tracks[2])),
##                                           mainPanel(
##                                               tableOutput("value"),
##                                               plotOutput("inputPlot"),
##                                               plotOutput("ipPlot")
##                                               )
##                                           )
##                             )),
                        
##                         server = shinyServer(function(input, output, session) {
##                             is.installed <- require(ggplot2)
##                             if (!is.installed) {
##                                 plotData <- function(data,title) {
##                                     plot(data$x,data$counts, col = "darkgray", xlab = "position", ylab = "counts and estimates", main = title)
##                                     lines(data$x,data$pred, col = "blue")
##                                 }
##                             }
##                             else {
##                                 plotData <- function(df,title) {
##                                     ggplot(data = df, aes(x=x,y=counts)) + geom_point(color = "darkgray") + 
##                                         geom_line(aes(x = x,y = pred), color = "blue") + ggtitle(title)
##                                 }
##                             }
##                             myLambda <- reactive({
##                                 as.character(lambdaVector[input$lambda])
##                             })
                            
##                             dataset <- reactive({
##                                 df <- fread(as.character(input$region), header = TRUE)
##                                 df
##                             })

##                             track1 <- reactive({
##                                 as.character(input$inputTrack)
##                             })

##                             track2 <- reactive({
##                                 as.character(input$ipTrack)
##                             })
                            
##                             output$value <- renderTable({
##                                 data.frame(
##                                     Name = c("Lambda Value"),
##                                     Value = as.character(c(myLambda())),
##                                     stringsAsFactors=FALSE
##                                     )
##                             })
##                             ## only works with one lambda for now
##                             output$inputPlot <- renderPlot({
##                                 df <- dataset()
##                                 name <- track1()
##                                 dfInput <- df[sample == name & lambda1 == myLambda(),]
##                                 plotData(dfInput,name)
##                             })
                            
##                             output$ipPlot <- renderPlot({
##                                 df <- dataset()
##                                 name <- track2()
##                                 dfIP <- df[sample == name & lambda1 == myLambda(),]
##                                 plotData(dfIP,name)
##                             })

##                             session$onSessionEnded(function() {
##                                 stopApp()
##                             })
##                         }))
##     }
##     runApp(app)
## }

