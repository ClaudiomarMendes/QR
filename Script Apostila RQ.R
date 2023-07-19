#-------------------------------------------------------------------------------------------------------------------------------------
# 1 Limpar a area de trabalho  ----
#-------------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())

#-------------------------------------------------------------------------------------------------------------------------------------
# 2 Importar os pacotes  ----
#-------------------------------------------------------------------------------------------------------------------------------------

library(quantreg)
library(Qtools)
library(dplyr)
library(ggplot2) 
library(gridExtra)
library(broom)
library(tcltk)

#-------------------------------------------------------------------------------------------------------------------------------------
# 3 Importar os dados  ----
#-------------------------------------------------------------------------------------------------------------------------------------

dados <- read.csv("dados.csv", header = T, sep = ";", dec = ",")

head(dados)

#-------------------------------------------------------------------------------------------------------------------------------------
# 4 Funcoes  ----
#-------------------------------------------------------------------------------------------------------------------------------------

# Plotagem dos coeficientes da Regressao Quantilica

coefRq <- function(modelo.rqs) {
  
  if (class(modelo.rqs)!= "rq")
    if (class(modelo.rqs)!= "rqs") {
      stop("Você deve usar essa função com objetos do tipo rq ou rqs")
    }
  
  graficos <- vector("list")
  
  for (i in 1:nrow(modelo.rqs$coefficients[, 0])) {
    graficos[[i]] <- broom::tidy(modelo.rqs) %>%
      filter(term == row.names(modelo.rqs$coefficients[, 0])[i]) %>%
      ggplot(aes(x = tau, y = estimate)) +
      labs(x = expression(R ^ {1} * (tau)),
           y = expression("Coeficiente"),
           title = paste("Coeficiente ", row.names(modelo.rqs$coefficients[, 0])[i])
      ) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
      geom_point(size = 2.5, shape = 18) +
      geom_hline(aes(yintercept = mean(estimate)), colour = "blue") +
      geom_hline(aes(yintercept = quantile(estimate, probs = 0.01)),
                 colour = "blue",
                 linetype = "dashed") +
      geom_hline(aes(yintercept = quantile(estimate, probs = 0.99)),
                 colour = "blue",
                 linetype = "dashed") +
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
  }
  
  gridGraficos <-
    do.call("grid.arrange", c(graficos, ncol = floor(sqrt(length(
      graficos
    )))))
  
  return(gridGraficos)
}

# Plotagem dos ajustes da Regressao Quantilica

ajusteRq <- function(dados, y, x, tau, method) {
  
  if (method != "rq")
    if (method != "rqss") {
      stop("use apenas 'rq' ou 'rqss' para o argumento method")
    }
  
  grafico <- ggplot(dados,
                    aes(x = eval(parse(text = x)), y = eval(parse(text = y)))) +
    labs(x = paste(x),
         y = paste(y),
         title = "Ajustes do Modelo de RQ") +
    geom_point(col = alpha("black", 0.7),
               size = 2.5,
               shape = 18) +
    geom_quantile(
      aes(colour = ..quantile..),
      quantiles = tau,
      method = method,
      lambda = 1
    ) +
    theme(
      plot.title = element_text(
        size = 16,
        vjust = 1,
        family = "serif"
      ),
      axis.title = element_text(
        size = 12,
        vjust = 1,
        family = "serif"
      )
    )
  return(grafico)
}

# Coeficiente de determinacao da Regressao Quantilica (R1)

R1 <- function(modelo.rqs) {
  
  if (class(modelo.rqs)!= "rq")
    if (class(modelo.rqs)!= "rqs") {
      stop("Você deve usar essa função com objetos do tipo rq ou rqs")
    }
  
  modelo.nulo =  suppressWarnings(rq(modelo.rqs$y ~ 1, tau = modelo.rqs$tau, method = modelo.rqs$method)) 
  
  res.abs.modelo.rqs = as.data.frame(abs(resid(modelo.rqs)))
  
  res.abs.modelo.nulo = as.data.frame(abs(resid(modelo.nulo)))
  
  r1 <- 1 - colSums(res.abs.modelo.rqs) / colSums(res.abs.modelo.nulo)
  
  return(r1)
  
}

# Plotagem do Coeficiente de determinacao da Regressao Quantilica

graphR1 = function(dataframe, modelo.rqs, trueScale = T) {
  
  if (class(modelo.rqs)!= "rq")
    if (class(modelo.rqs)!= "rqs") {
      stop("Você deve usar essa função com objetos do tipo rq ou rqs")
    }
  
  taus = modelo.rqs$tau
  R1OK = R1(modelo.rqs)
  
  data.graph = data.frame(taus, R1OK)
  graph = ggplot(data.graph , aes(x=taus , y=R1OK)) + ylim (c(0,1)) 
  graphOK = graph + geom_point() +geom_line(colour = "blue") + 
    labs(x = expression(tau), y = expression(R^{1}*(tau)), title = "Coeficiente de Determinacao") + 
    theme(
    plot.title = element_text(
      size = 16,
      vjust = 1,
      family = "serif"
    ),
    axis.title = element_text(
      size = 12,
      vjust = 1,
      family = "serif"
    )
  )
  
  return(graphOK)
}

# delta AIC

deltaAIC <- function(modelo.rqs) {
  
  if (class(modelo.rqs)!= "rq")
    if (class(modelo.rqs)!= "rqs") {
      stop("Você deve usar essa função com objetos do tipo rq ou rqs")
    }
  
  modelo.nulo =  suppressWarnings(rq(modelo.rqs$y ~ 1, tau = modelo.rqs$tau, method = modelo.rqs$method)) 
  
  delta = AIC(modelo.nulo) - AIC(modelo.rqs)
  
  return(delta)
  
}

# Funcao de distribuicao acumulada para a distribuicao de laplace assimetrica por SOUZA (2012)

palap <- function (q, mu = 0, sigma = 1, tau = 0.5) {
  
  saida <- vector(length = max(length(q), length(mu)))
  
  if (length(q) != length(mu))
  {
    for (k in 1:length(q))
    {
      if (!q[k] > mu)
        saida[k] <- tau * exp((1 / sigma) * (1 - tau) * (q[k] - mu))
      else
        saida[k] <- 1 - (1 - tau) * exp(-(tau / sigma) * (q[k] - mu))
    }
  }
  else
  {
    for (k in 1:length(q))
    {
      if (!q[k] > mu[k])
        saida[k] <- tau * exp((1 / sigma) * (1 - tau) * (q[k] - mu[k]))
      else
        saida[k] <- 1 - (1 - tau) * exp(-(tau / sigma) * (q[k] - mu[k]))
    }
  }
  return(saida)
}

# Grafico de residuos quantilicos em funcao dos valores ajustados por SOUZA (2012)

grafResiduosRQ <- function (model, scales = "fixed") {
  tau <- model$tau
  
  n <-
    ifelse(length(tau) == 1, length(residuals(model)), nrow(residuals(model)))
  
  preditos <- fitted(model)
  rho.hat <- model$rho / n
  
  if (length(tau) > 1)
  {
    rho.hat <- model$rho / n
    residuos <- list()
    
    for (k in 1:length(tau))
    {
      residuos[[k]] <-
        qnorm(palap(
          q = as.numeric(model$y),
          mu = preditos[, k],
          sigma = rho.hat[k],
          tau = tau[k]
        ))
    }
    
    residuos <- unlist(residuos)
    preditos <- as.vector(fitted(model))
    tau <- rep(tau, each = n)
    dados <- data.frame(preditos, residuos, tau)
    g <-
      ggplot(dados, aes(x = preditos, y = residuos, group = tau)) + geom_point() + facet_wrap( ~ tau, ncol = 3, scales = scales)
    g + geom_hline(aes(yintercept = qnorm(0.025)), colour = "blue") + geom_hline(aes(yintercept = qnorm(0.975)), colour = "blue") + 
      labs(x = "Valores preditos", y = "Resíduos Quantílicos", title = "Resíduos Quantílicos") + 
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
  }
  else
  {
    residuos <-
      qnorm(palap(as.numeric(model$y), preditos, rho.hat, tau))
    dados <- data.frame(preditos, residuos)
    g <-
      ggplot(dados, aes(x = preditos, y = residuos)) + geom_point()
    g + geom_hline(aes(yintercept = qnorm(0.025)), colour = "blue") + geom_hline(aes(yintercept = qnorm(0.975)), colour = "blue") + 
      labs(x = "Valores preditos", y = "Resíduos Quantílicos", title = "Resíduos Quantílicos") + 
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
  }
}

# Histogramas dos resíduos quantílicos por SOUZA (2012)

hist.ResiduosRQ <- function(model, scales, ncolunas) {
  
  tau <- model$tau
  rho <- model$rho
  
  if (length(tau) == 1)
  {
    n <- length(residuals(model))
    sigmahat <- model$rho / n
    predicted <- fitted(model)
    res.quant = qnorm(palap(
      q = as.numeric(model$y),
      mu = predicted,
      sigma =
        model$rho / n,
      tau = tau
    ))
    sample.quant <- sort(res.quant)
    db <- data.frame(sample.quant , Tau = paste("Tau = ", tau, sep = ""))
    g <-
      ggplot(db, aes(x = sample.quant, y = ..density..)) + geom_histogram()
    + labs(x = "Resíduos Quantílicos", y = "Densidade", title = "Histogramas dos Resíduos Quantílicos") 
    graph <-
      g + facet_wrap( ~ Tau) + geom_histogram(colour = "blue", fill = "white") + 
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
  }
  
  else
  {
    n <- nrow(residuals(model))
    residuos <- list()
    for (k in 1:length(tau))
    {
      residuos[[k]] <- qnorm(palap(
        q = as.numeric(model$y),
        mu = fitted(model)
        [, k],
        sigma = model$rho[k] / n ,
        tau = tau[k]
      ))
    }
    residuos <- lapply(residuos, sort)
    residuos <- unlist(residuos)
    predicted <- as.vector(fitted(model))
    tau.total <- rep(tau , each = n)
    dados <- data.frame(residuos, tau = tau.total)
    g <-
      ggplot(dados, aes (x = residuos, y = ..density.., group = tau)) + geom_histogram() + facet_wrap( ~
                                                                                                         tau , scales = scales, ncol = ncolunas)
    graph <-
      g + geom_histogram(colour = "blue", fill = "white") +
      labs(x = "Resíduos Quantílicos", y = "Densidade", title = "Histogramas dos Resíduos Quantílicos") + 
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
  }
  return (graph)
}

# Envelope dos resíduos quantílicos

envel.rq <- function(model,
                     data,
                     ncolunas = 1,
                     scales = "fixed") {
  ralap <-
    function(n,
             mu = 0,
             sigma = 1,
             tau = 0.5)
      # funcao geradora de numeros aleatorios da distribuicao lapace assimetrica
    {
      saida = vector(length = n)
      if (length(mu) == 1)
        mu = rep(mu, n)
      if (length(mu) != n)
        stop("Mu e n têm dimensões diferentes")
      for (k in 1:n)
      {
        u1 = rexp(1)
        u2 = rexp(1)
        saida.padrao <- u1 / tau - u2 / (1 - tau)
        saida[k] = mu[k] + sigma * saida.padrao
      }
      return (saida)
    }
  
  tau <- model$tau
  rho <- model$rho
  
  if (length(tau) == 1)
  {
    n <- length(residuals(model))
    sigmahat <- model$rho / n
    predicted <- fitted(model)
    res.quant = qnorm(palap(
      q = as.numeric(model$y),
      mu = predicted,
      sigma =
        model$rho / n ,
      tau = tau
    ))
    
    e <- matrix(0, n, 100)
    e1 <- numeric(n)
    e2 <- numeric(n)
    
    for (i in 1:100)
    {
      e[, i] <- ralap(n ,
                      mu = predicted,
                      sigma = sigmahat,
                      tau = tau)
      sim.model <-
        rq(as.formula(paste("e[,i] ~ " , model$formula[3])),
           data , tau = tau)
      e[, i] <-
        qnorm(palap(
          q = as.numeric(sim.model$y),
          mu = fitted(sim.model),
          sigma = sim.model$rho / n,
          tau = tau
        ))
      e[, i] <- sort(e[, i])
    } 
    
    for (i in 1:n)
    {
      eo <- sort(e[i,])
      e1[i] <- (eo[2] + eo[3]) / 2
      e2[i] <- (eo[97] + eo[98]) / 2
    }
    theoretical.quant <- qnorm(1:n / (n + 1))
    sample.quant <- sort(res.quant)
    db <-
      data.frame(theoretical.quant, sample.quant, Tau = paste ("Tau = ", tau , sep =
                                                                 ""))
    db1 <- data.frame(e1 = sort(e1), theoretical.quant)
    db2 <- data.frame(e2 = sort(e2), theoretical.quant)
    g <-
      ggplot(db, aes(x = theoretical.quant, y = sample.quant)) + geom_point() +
      labs(x = "Quantis teoricos", y = "Quantis amostrais", title = "Envelope dos Resíduos Quantílicos") + 
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
    graph <-
      g + geom_line(aes(y = e1), db1, colour = "blue") + geom_line(aes(y = e2), db2,colour = "blue") + facet_wrap( ~ Tau)
  }
  
  else
  {
    n <- nrow(residuals(model))
    residuos <- list()
    
    for (k in 1:length(tau))
    {
      residuos[[k]] <-
        qnorm(palap(
          q = as.numeric(model$y),
          mu = fitted(model)
          [, k],
          sigma = model$rho[k] / n,
          tau = tau[k]
        ))
    }
    
    residuos <- lapply(residuos, sort)
    residuos <- unlist(residuos)
    predicted <- as.vector(fitted(model))
    tau.total <- rep(tau, each = n)
    e <- list()
    e1 <- list(n)
    e2 <- list(n)
    
    for (k in 1:length(tau))
    {
      e[[k]] <- matrix(0, n, 100)
      e1[[k]] <- numeric(n)
      e2[[k]] <- numeric(n)
      
      for (i in 1:100)
      {
        e[[k]][, i] <- ralap(
          n ,
          mu = fitted(model)[, k],
          sigma = model$rho[k] / n ,
          tau = tau[k]
        )
        sim.model <-
          rq(as.formula(paste("e [[k]][,i] ~ ", model$formula[3])), data, tau = tau[k])
        e[[k]][, i] <-
          qnorm(palap(
            q = as.numeric(sim.model$y),
            mu = fitted(sim.model),
            sigma = sim.model$rho / n,
            tau = tau[k]
          ))
        e[[k]][, i] <- sort (e[[k]][, i])
      }
      
      for (i in 1:n)
      {
        eo <- sort(e[[k]][i,])
        e1[[k]][i] <- (eo[2] + eo[3]) / 2
        e2[[k]][i] <- (eo[97] + eo[98]) / 2
      }
    }
    
    e1 <- as.numeric(unlist(lapply(e1, sort)))
    e2 <- as.numeric(unlist(lapply(e2, sort)))
    quantis.teoricos <- vector(length = n * length(tau))
    quantis.teoricos <- qnorm(1:n / (n + 1))
    
    for (j in 2:length(tau))
    {
      quantis.teoricos <- c(quantis.teoricos, qnorm(1:n / (n + 1)))
    }
    
    dados <-
      data.frame(residuos, quantis.teoricos , tau = tau.total, e1, e2)
    g <-
      ggplot(dados, aes(x = quantis.teoricos, y = residuos , group = tau)) +
      geom_point() + facet_wrap(~ tau, scales = scales , ncol = ncolunas)
    graph <-
      g + geom_line(aes(y = e1, group = tau), dados, colour = "blue") + geom_line(aes(y = e2, group = tau), dados, colour = "blue") +
      labs(x = "Quantis teoricos", y = "Quantis amostrais", title = "Envelope dos Resíduos Quantílicos") + 
      theme(
        plot.title = element_text(
          size = 16,
          vjust = 1,
          family = "serif"
        ),
        axis.title = element_text(
          size = 12,
          vjust = 1,
          family = "serif"
        )
      )
  }
  return (graph)
}

#-------------------------------------------------------------------------------------------------------------------------------------
# 5 Resultados  ----
#-------------------------------------------------------------------------------------------------------------------------------------

modelo.rq = rq(Y ~ X, tau = 2:8/10, method = "br", data = dados)

suppressWarnings(summary(modelo.rq))

R1(modelo.rq)

Qtools::GOFTest(modelo.rq, alpha = 0.05, B = 1000) # Teste de falta de ajusde em modelos de regressao quantilica

deltaAIC(modelo.rq) # Normalmente, um deltaAIC > 10 sugere que o modelo tem algum poder explicativo

ajusteRq(dados, "Y", "X", 2:8/10, "rq")

ajusteRq(dados, "Y", "X", 2:8/10, "rqss")

coefs <- coefRq(modelo.rq)

graphR1(dados, modelo.rq)

grafResiduosRQ(modelo.rq)

hist.ResiduosRQ(modelo.rq, "fixed", 3)

envel.rq(modelo.rq, dados)
