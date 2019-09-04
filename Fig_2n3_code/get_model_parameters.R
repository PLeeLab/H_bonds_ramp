#------ modeling the pD behaviour ------#
get_model_parameters <- function(ELEMENT, START_AFTER_CODON_POSITION, END_AT_CODON_POSITION, RETURN="data"){ # ELEMENT is a string with the name of the element to subset (e.g. "Hydrogen bonds", "Oxygen", etc.)
  
  Sub_df <- Elements_df_boot %>% 
    filter(Position > START_AFTER_CODON_POSITION & Position <= END_AT_CODON_POSITION) %>% 
    filter(Element == ELEMENT)

  SS <- try(getInitial(Mean ~ SSlogis(Position, asymptote, xmid, scale), data=Sub_df), silent = TRUE)
  
  if(is(SS, "try-error")) {
    return(
      data.frame(Modeling = "Unsuccessful",
                 Element = ELEMENT,
                 Carrying_capacity = mean(Sub_df$Mean, trim=0.2), # Trim (filter out) 20% of the data at the extremes (100 pos * 20% = 20 codons are removed)
                 Rate = NA,
                 Initial_cost = min(Sub_df$Mean),
                 AIC_Uniform = NA, AIC_Linear = NA, AIC_Exponential = NA,
                 BIC_Uniform = NA, BIC_Linear = NA, BIC_Exponential = NA,
                 P_Carrying_capacity = NA, P_Rate = NA, P_Initial_cost=NA,
                 row.names = NULL)
    )
    
  } else {
    A_start <- SS["asymptote"]
    B_start <- 1/SS["scale"]
    C_start <- SS["asymptote"]/(exp(SS["xmid"]/SS["scale"])+1)
    names(A_start) <- NULL
    names(B_start) <- NULL
    names(C_start) <- NULL
    
    # uniform model 
    uniform_formula <- formula(Mean ~ A)
    uniform_model <- nls(uniform_formula, start=list(A = A_start), data = Sub_df)
    
    # linear model
    linear_model <- lm(Mean ~ Position, data = Sub_df) 
    
    # exponential model (https://en.wikipedia.org/wiki/Logistic_function) Section "Applications" In ecology: modeling population growth"
    #exponential_formula <- formula(Mean ~  A / (1+ ( (A-C) / C ) * exp(-B*Position) ) )
    exponential_formula <- formula(Mean ~  A*C*exp(B*Position)/(A+C*(exp(B*Position)-1)))
    
    exponential_model <- try(nls(exponential_formula, start=list(A = A_start, B = B_start, C = C_start), data = Sub_df), silent = TRUE)
    
    if(is(exponential_model, "try-error")) {
      return(
        data.frame(Modeling = "Unsuccessful",
                   Element = ELEMENT,
                   Carrying_capacity = mean(Sub_df$Mean, trim=0.2), # Trim (filter out) 20% of the data at the extremes (100 pos * 20% = 20 codons are removed)
                   Rate = NA,
                   Initial_cost = min(Sub_df$Mean),
                   AIC_Uniform = NA, AIC_Linear = NA, AIC_Exponential = NA,
                   BIC_Uniform = NA, BIC_Linear = NA, BIC_Exponential = NA,
                   P_Carrying_capacity = NA, P_Rate = NA, P_Initial_cost=NA,
                   row.names = NULL)
      )
    } else {
      Carrying_capacity <-  coefficients(exponential_model)["A"]
      Rate <- coefficients(exponential_model)["B"]
      Initial_cost <- coefficients(exponential_model)["C"]
    }
  }
  
  if (RETURN=="plot") { # Plot model fitness
    return(
      ggplot(data = cbind(Sub_df,
                          uniform_fit=predict(uniform_model),
                          linear_fit=predict(linear_model),
                          exponential_fit=predict(exponential_model)),
             aes(Position, Mean)) +
        geom_pointrange(mapping = aes(y=Mean, ymin=Lower, ymax=Upper), size=0.05) +
        #geom_point(shape=21, colour="black", fill="grey", size=0.25) +
        #geom_line(aes(Position, uniform_fit), colour="green", lwd=2, alpha=0.75) +
        #geom_line(aes(Position, linear_fit), colour="purple", lwd=2, alpha=0.75) +
        geom_line(aes(Position, exponential_fit), colour="red", lwd=0.5, alpha=1) +
        xlab(label = "Codon position (5'-end)") +
        ylab(label = "Hydrogen bonds per codon") +
        ggtitle(label = "Escherichia coli", subtitle = "str. K-12 substr. MG1655") +
        theme_new() +
        NULL
    )
  } else if (RETURN=="data") {
    return(
      data.frame(Modeling = "Successful",
                 Element = ELEMENT,
                 Carrying_capacity = Carrying_capacity,
                 Rate = Rate,
                 Initial_cost = Initial_cost,
                 AIC_Uniform = AIC(uniform_model), AIC_Linear = AIC(linear_model), AIC_Exponential = AIC(exponential_model),
                 BIC_Uniform = BIC(uniform_model), BIC_Linear = BIC(linear_model), BIC_Exponential = BIC(exponential_model),
                 P_Carrying_capacity = summary(exponential_model)$coefficients[1,4], P_Rate = summary(exponential_model)$coefficients[2,4], P_Initial_cost=summary(exponential_model)$coefficients[3,4],
                 row.names = NULL)
    )
  }
}