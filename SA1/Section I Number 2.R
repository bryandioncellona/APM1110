library(shiny)
library(ggplot2)

# ==========================================
# USER INTERFACE
# ==========================================

ui <- fluidPage(
  
  titlePanel("Discrete Probability Distribution Calculator"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "type",
        "Select Distribution:",
        choices = c("Univariate", "Bivariate")
      ),
      
      conditionalPanel(
        condition = "input.type == 'Univariate'",
        
        textInput(
          "x", "Values of X:",
          value = "1,2,3,4"
        ),
        
        textInput(
          "p", "Probabilities:",
          value = "0.1,0.2,0.3,0.4"
        )
      ),
      
      conditionalPanel(
        condition = "input.type == 'Bivariate'",
        
        textInput(
          "bx", "Values of X:",
          value = "1,2"
        ),
        
        textInput(
          "by", "Values of Y:",
          value = "1,2"
        ),
        
        textAreaInput(
          "joint",
          "Joint Probabilities (rows separated by ;):",
          value = "0.1,0.2;0.3,0.4",
          rows = 4
        ),
        
        helpText(
          "Each row represents an X value.",
          "Each column represents a Y value."
        )
      ),
      
      actionButton(
        "calculate", "Calculate",
        class = "btn-primary"
      ),
      
      br(), br(),
      
      helpText(
        "All probabilities must be between 0 and 1",
        "and their total must equal 1."
      )
    ),
    
    mainPanel(
      
      h3("Results"),
      verbatimTextOutput("results"),
      
      conditionalPanel(
        condition = "input.type == 'Bivariate'",
        
        h3("Joint Probability Distribution"),
        tableOutput("jointTable"),
        
        h3("Marginal Distribution of X"),
        tableOutput("marginalX"),
        
        h3("Marginal Distribution of Y"),
        tableOutput("marginalY"),
        
        h3("Conditional Distribution of X given Y"),
        tableOutput("conditionalX"),
        
        h3("Conditional Distribution of Y given X"),
        tableOutput("conditionalY")
      ),
      
      h3("PDF Plot"),
      plotOutput("pdfPlot", height = "350px"),
      
      h3("CDF Plot"),
      plotOutput("cdfPlot", height = "350px")
    )
  )
)

# ==========================================
# SERVER
# ==========================================

server <- function(input, output, session) {
  
  # Parse comma-separated numbers
  parse_numbers <- function(s) {
    
    parts <- trimws(
      strsplit(s, ",", fixed = TRUE)[[1]]
    )
    
    if (length(parts) == 0 || any(parts == "")) {
      stop("Please enter complete numeric values.")
    }
    
    values <- suppressWarnings(as.numeric(parts))
    
    if (any(!is.finite(values))) {
      stop("All entries must be valid finite numbers.")
    }
    
    return(values)
  }
  
  # Validate probabilities
  check_probabilities <- function(p) {
    
    if (any(p < 0 | p > 1)) {
      stop("Probabilities must be between 0 and 1.")
    }
    
    if (abs(sum(p) - 1) > 1e-8) {
      stop(
        paste(
          "Probabilities must sum to 1.",
          "Current sum:", round(sum(p), 8)
        )
      )
    }
  }
  
  # ========================================
  # CALCULATIONS
  # ========================================
  
  calculated <- eventReactive(input$calculate, {
    
    tryCatch({
      
      if (input$type == "Univariate") {
        
        x <- parse_numbers(input$x)
        p <- parse_numbers(input$p)
        
        if (length(x) != length(p)) {
          stop(
            "The number of X values must match",
            "the number of probabilities."
          )
        }
        
        if (anyDuplicated(x)) {
          stop("X values must be unique.")
        }
        
        check_probabilities(p)
        
        # Sort values for correct CDF
        ord <- order(x)
        x <- x[ord]
        p <- p[ord]
        
        # Mean and variance
        meanX <- sum(x * p)
        varX <- sum((x - meanX)^2 * p)
        
        list(
          valid = TRUE,
          type = "Univariate",
          x = x,
          p = p,
          meanX = meanX,
          varX = varX
        )
        
      } else {
        
        x <- parse_numbers(input$bx)
        y <- parse_numbers(input$by)
        
        if (anyDuplicated(x) || anyDuplicated(y)) {
          stop("X and Y values must be unique.")
        }
        
        # Read joint probability matrix
        rows <- strsplit(
          input$joint, ";", fixed = TRUE
        )[[1]]
        
        parsed <- lapply(rows, parse_numbers)
        
        if (length(parsed) != length(x) ||
            any(lengths(parsed) != length(y))) {
          stop(
            "Joint probability matrix dimensions",
            "must match X and Y values."
          )
        }
        
        joint <- do.call(rbind, parsed)
        
        check_probabilities(as.vector(joint))
        
        # Sort X and Y values
        ox <- order(x)
        oy <- order(y)
        
        x <- x[ox]
        y <- y[oy]
        
        joint <- joint[ox, oy, drop = FALSE]
        
        # Marginal probabilities
        px <- rowSums(joint)
        py <- colSums(joint)
        
        # Means
        meanX <- sum(x * px)
        meanY <- sum(y * py)
        
        # Variances
        varX <- sum((x - meanX)^2 * px)
        varY <- sum((y - meanY)^2 * py)
        
        # Conditional P(X | Y)
        condX <- sweep(
          joint, 2,
          ifelse(py > 0, py, 1), "/"
        )
        
        condX[, py == 0] <- NA_real_
        
        # Conditional P(Y | X)
        condY <- sweep(
          joint, 1,
          ifelse(px > 0, px, 1), "/"
        )
        
        condY[px == 0, ] <- NA_real_
        
        list(
          valid = TRUE,
          type = "Bivariate",
          x = x,
          y = y,
          joint = joint,
          px = px,
          py = py,
          meanX = meanX,
          meanY = meanY,
          varX = varX,
          varY = varY,
          condX = condX,
          condY = condY
        )
      }
      
    }, error = function(e) {
      
      list(
        valid = FALSE,
        error = conditionMessage(e)
      )
      
    })
    
  })
  
  # Only use valid calculated data
  validData <- reactive({
    
    d <- calculated()
    
    validate(
      need(d$valid, d$error)
    )
    
    d
  })
  
  # ========================================
  # DISPLAY RESULTS
  # ========================================
  
  output$results <- renderPrint({
    
    req(input$calculate > 0)
    
    d <- calculated()
    
    if (!d$valid) {
      cat("ERROR:", d$error, "\n")
      return(invisible(NULL))
    }
    
    cat("========== RESULTS ==========\n\n")
    
    cat("Mean of X:", round(d$meanX, 6), "\n")
    cat("Variance of X:", round(d$varX, 6), "\n")
    
    if (d$type == "Bivariate") {
      
      cat("\nMean of Y:", round(d$meanY, 6), "\n")
      cat("Variance of Y:", round(d$varY, 6), "\n")
    }
    
    cat("\nProbability validation: PASSED\n")
  })
  
  # ========================================
  # JOINT PROBABILITY TABLE
  # ========================================
  
  output$jointTable <- renderTable({
    
    req(input$type == "Bivariate")
    d <- validData()
    req(d$type == "Bivariate")
    
    tab <- as.data.frame(d$joint)
    
    colnames(tab) <- paste0("Y = ", d$y)
    
    data.frame(
      X = d$x,
      tab,
      check.names = FALSE
    )
  })
  
  # ========================================
  # MARGINAL DISTRIBUTIONS
  # ========================================
  
  output$marginalX <- renderTable({
    
    req(input$type == "Bivariate")
    d <- validData()
    req(d$type == "Bivariate")
    
    data.frame(
      X = d$x,
      Probability = d$px,
      CDF = cumsum(d$px)
    )
  })
  
  output$marginalY <- renderTable({
    
    req(input$type == "Bivariate")
    d <- validData()
    req(d$type == "Bivariate")
    
    data.frame(
      Y = d$y,
      Probability = d$py,
      CDF = cumsum(d$py)
    )
  })
  
  # ========================================
  # CONDITIONAL DISTRIBUTIONS
  # ========================================
  
  output$conditionalX <- renderTable({
    
    req(input$type == "Bivariate")
    d <- validData()
    req(d$type == "Bivariate")
    
    tab <- as.data.frame(d$condX)
    
    colnames(tab) <- paste0("Y = ", d$y)
    
    data.frame(
      X = d$x,
      tab,
      check.names = FALSE
    )
  })
  
  output$conditionalY <- renderTable({
    
    req(input$type == "Bivariate")
    d <- validData()
    req(d$type == "Bivariate")
    
    tab <- as.data.frame(d$condY)
    
    colnames(tab) <- paste0("Y = ", d$y)
    
    data.frame(
      X = d$x,
      tab,
      check.names = FALSE
    )
  })
  
  # ========================================
  # PDF PLOT
  # ========================================
  
  output$pdfPlot <- renderPlot({
    
    d <- validData()
    
    if (d$type == "Univariate") {
      
      plotData <- data.frame(
        Value = d$x,
        Probability = d$p,
        Variable = "X"
      )
      
    } else {
      
      plotData <- rbind(
        
        data.frame(
          Value = d$x,
          Probability = d$px,
          Variable = "X"
        ),
        
        data.frame(
          Value = d$y,
          Probability = d$py,
          Variable = "Y"
        )
      )
    }
    
    ggplot(
      plotData,
      aes(
        x = factor(Value),
        y = Probability,
        fill = Variable
      )
    ) +
      
      geom_col(width = 0.6) +
      
      facet_wrap(
        ~Variable,
        scales = "free_x"
      ) +
      
      labs(
        title = "Probability Distribution (PDF)",
        x = "Random Variable Value",
        y = "Probability"
      ) +
      
      theme_minimal() +
      
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  })
  
  # ========================================
  # CDF PLOT
  # ========================================
  
  output$cdfPlot <- renderPlot({
    
    d <- validData()
    
    if (d$type == "Univariate") {
      
      plotData <- data.frame(
        Value = d$x,
        CDF = cumsum(d$p),
        Variable = "X"
      )
      
    } else {
      
      plotData <- rbind(
        
        data.frame(
          Value = d$x,
          CDF = cumsum(d$px),
          Variable = "X"
        ),
        
        data.frame(
          Value = d$y,
          CDF = cumsum(d$py),
          Variable = "Y"
        )
      )
    }
    
    ggplot(
      plotData,
      aes(x = Value, y = CDF)
    ) +
      
      geom_step(
        direction = "hv",
        color = "blue",
        linewidth = 1
      ) +
      
      geom_point(
        color = "blue",
        size = 2
      ) +
      
      facet_wrap(
        ~Variable,
        scales = "free_x"
      ) +
      
      coord_cartesian(ylim = c(0, 1)) +
      
      labs(
        title = "Cumulative Distribution Function (CDF)",
        x = "Random Variable Value",
        y = "Cumulative Probability"
      ) +
      
      theme_minimal() +
      
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  })
  
}

# ==========================================
# RUN APPLICATION
# ==========================================

shinyApp(ui = ui, server = server)