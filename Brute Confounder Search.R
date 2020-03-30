#Requires the initialization of the GenMatchRose function first.
genConfoundSearch <- function(Tr, X, Y, Z= X, V = rep(1, length(Y)), Gamma = 5, GammaInc = 0.25, BalanceMatrix = X, 
                              estimand = "ATT", M=1, weights = NULL, BiasAdjust = FALSE, Weight = 1, Var.calc = 0, sample = FALSE,
                              search.pop.size = 100, search.mutation.rate = 0.05, search.num.generations = 10, match.out = NULL,
                              version = "standard", match.pop.size = 100, match.max.generations = 100, match.wait.generations = 4, 
                              hard.generation.limit = FALSE, starting.values = rep(1,ncol(X)), MemoryMatrix = TRUE, exact = NULL,
                              caliper = NULL, replace = TRUE, ties = TRUE, CommonSupport = FALSE, nboots = 0, ks = TRUE,
                              verbose = FALSE, distance.tolerance = 1e-05, tolerance = sqrt(.Machine$double.eps), min.weight =0,
                              max.weight = 1000, Domains = NULL, print.level = 2, project.path = NULL, paired = TRUE, loss = 1,
                              data.type.integer = FALSE, restrict = NULL, cluster = FALSE, balance = TRUE){
  #Defines the main function to generate kappas. Has all inputs from GenMatchRose() function and new inputs to determine parameters
  #of the main genetic algorithm. Main function is a genetic algorithm built from scratch.
  
  base.genout <- GenMatchRose(Y=Y, Tr = Tr, X = X, BalanceMatrix = BalanceMatrix, estimand = estimand, M = M, weights = weights,
                         pop.size = match.pop.size, max.generation = match.max.generations, wait.generations = match.wait.generations,
                         hard.generation.limit = hard.generation.limit, starting.values = starting.values,
                         MemoryMatrix = MemoryMatrix, exact = exact, caliper = caliper, replace = replace, ties = ties,
                         CommonSupport = CommonSupport, nboots = nboots, ks = ks, verbose = verbose, distance.tolerance = distance.tolerance,
                         tolerance = tolerance, min.weight = min.weight, max.weight = max.weight, Domains= Domains,
                         print.level = print.level, project.path = project.path, paired = paired, loss = loss, 
                         data.type.integer = data.type.integer, restrict = restrict, cluster = cluster, balance = balance)
  base.mout <- Match(Y = Y, Tr = Tr, X= X, Z=Z, V=V, estimand = estimand, M = M, BiasAdjust = BiasAdjust, exact = exact, caliper = caliper,
              replace = replace, ties = ties, CommonSupport = CommonSupport, Weight = Weight, Weight.matrix = genout, weights = weights,
              Var.calc = Var.calc, sample = sample, restrict = restrict, match.out = match.out, version = version)
  base.sens <- psens(base.mout, Gamma = Gamma, GammaInc = GammaInc)
  #Runs the rosenbaum sensitivity analysis on the robustness optimized base variables to get a baseline of the robustness of the current data. 
  base.bounds <- base.sens$bounds[,3]
  base.len_sens <- length(base.bounds)
  base.fit.value <- base.len_sens + base.bounds[1]
  for(i in 1:base.len_sens){
    if(base.bounds[i] <= 0.05){
      base.fit.value <- base.len_sens - i + base.bounds[i]
    }
  }
  #Calculates a base fitness value using the same method as that described in the GenMatchRose function
  
  search.chromosome.length <- nrow(X)
  #Finds how many observations there are so each kappa can have the same amount of data points.
  search.num.survivors = trunc(search.pop.size/2)
  #Sets the number of survivors to half population size (truncated because it needs to be an integer)
  
  #Genetic Algorithms are composed of five functions: population initialization, fitness function, survivor selection, crossover, and mutation.
  
  init.pop <- function(){
    #defines the first function to initialize the population.
    gene_pool <- matrix(data = runif(search.pop.size*search.chromosome.length),nrow = search.chromosome.length, ncol = search.pop.size)
    #creates the gene pool, a matrix where each column is a kappa and so there are as many columns as initial population size, and rows are
    #the number of observations in the given data variables. All values are randomly selected and between 0 and 1 because scale should not matter.
    return(gene_pool)
    #returns created gene pool
  }
  
  fitness.value <- function(kap){
    #creates a function that returns the fitness value of the inputted kappa
    new_matrix <- cbind(X, kap)
    #combines the inputted kappa with the given observed variables
    new.genout <- GenMatchRose(Y=Y, Tr = Tr, X = new_matrix, BalanceMatrix = new_matrix, estimand = estimand, M = M, weights = weights,
                                pop.size = match.pop.size, max.generation = match.max.generations, wait.generations = match.wait.generations,
                                hard.generation.limit = hard.generation.limit, starting.values = starting.values,
                                MemoryMatrix = MemoryMatrix, exact = exact, caliper = caliper, replace = replace, ties = ties,
                                CommonSupport = CommonSupport, nboots = nboots, ks = ks, verbose = verbose, distance.tolerance = distance.tolerance,
                                tolerance = tolerance, min.weight = min.weight, max.weight = max.weight, Domains= Domains,
                                print.level = print.level, project.path = project.path, paired = paired, loss = loss, 
                                data.type.integer = data.type.integer, restrict = restrict, cluster = cluster, balance = balance)
    new.mout <- Match(Y = Y, Tr = Tr, X= new_matrix, Z=Z, V=V, estimand = estimand, M = M, BiasAdjust = BiasAdjust, exact = exact, caliper = caliper,
                       replace = replace, ties = ties, CommonSupport = CommonSupport, Weight = Weight, Weight.matrix = new.genout, weights = weights,
                       Var.calc = Var.calc, sample = sample, restrict = restrict, match.out = match.out, version = version)
    new.sens <- psens(new.mout, Gamma = Gamma, GammaInc = GammaInc)
    new.bounds <- new.sens$bounds[,3]
    new.len_sens <- length(new.bounds)
    new.fit.value <- new.len_sens + new.bounds[1]
    for(i in 1:new.len_sens){
      if(new.bounds[i] <= 0.05){
        new.fit.value <- new.len_sens - i + new.bounds[i]
      }
    }
    #Runs the same fitness value analysis done for the base and in the GenMatchRose function
    return(base.fit.value - new.fit.value)
    #returns the final fitness value which is the difference between the fitness of the this new dataset with the kappa and the base line fitness.
    #Difference between new and base is necessary because we want the kappas that improve robustness relative to the base.
  }
  
  fitness.function <- function(pop){
    #Defines the fitness function that assigns a fitness value to all generated chromosomes.
    list.fit <- c()
    #initializes a list for all the fitness values of the current population
    for(i in 1:length(pop[1,])){
      list.fit[i] <- fitness.value(kap = pop[,i])
    }
    #runs a for-loop over all the given population, finding the fitness value for all kappas in the population and assigning it to the fitness list
    return(list.fit)
    #returns the list with all the fitnesses
  }
  
  selection <- function(){
    #defines the survivor selection function
    ordered_fitness_vector <- order(fitness_vector, decreasing = TRUE)
    #Assigns the index of the fitness vector from high to low, with the index of the fitness vector being the sasme as that of the current population
    return(current.pop[,ordered_fitness_vector[1:search.num.survivors]])
    #Returns a surviving subset of the current population with the highest fitness values.
  }
  
  crossover <- function(){
    #Defines the crossover function, which is the reproduction part of the genetic metaphor.
    duplicate_size <- search.pop.size - search.num.survivors
    #finds the number of new chromosomes to produce
    duplicate_survivors <- survivors[,1:duplicate_size]
    #duplicates the existing survivors to fill up the population size
    for(i in 1:length(duplicate_size)){
      duplicate_survivors[,i] <- sample(duplicate_survivors[,i])
    }
    #reorganizes the order of the values in the duplicated chromosomes.
    return(duplicate_survivors)
    #returns the reorganized duplicates
  }
  
  mutation <- function(){
    #defines the mutation function
    for(i in 1:search.pop.size){
      #for-loop to go over all the population to see if they mutate
      if(runif(1) <= search.mutation.rate){
        #probability test to see if given chromosome gets to mutate based on given mutation rate
        for(z in 1:search.chromosome.length){
          #nested for-loop to go over all the genes within each chromosome
          if(runif(1) <= search.mutation.rate){
            #second probability test to see if given gene in chromoosome gets to mutate
            if(runif(1) <= 0.5){
              current.pop[z,i] <- current.pop[z,i]*1.25
            }
            else{
              current.pop[z,i] <- current.pop[z,i]*0.75
            }
            #test to see whether gene mutates in an increasing or decreasing direction (50% chance of either)
          }
        }
      }
    }}

  #start of te main genetic algorithm code
  
  current.pop <- init.pop()
  #creates the initial population
  
  for(i in 1:search.num.generations){
    #runs a for loop to go through the generations of the algorithm
    fitness_vector <- fitness.function(current.pop)
    #creates a fitness vector containing the fitness values of the current population
    survivors <- selection()
    #Finds the survivors of the current population
    crossed <- crossover()
    #Creates new crossed over chromosomes
    current.pop[,1:search.num.survivors] <- survivors
    current.pop[,(search.num.survivors+1):search.pop.size] <- crossed
    #Changes the current population to the survivors in the beginning of the current pop and the crossed over in the second part of the current pop
    mutation()
    #mutates the current population
  }
  
  fitness_vector <- fitness.function(current.pop)
  #finds the fitness vector on the last generation
  kap_index <- which(fitness_vector>0)
  #Gets the index of the kappas (chromsomes) who have a positive fitness value because we only want those that improve robustness
  
  if(length(kap_index) == 0){
    print("No Kappas found that reduce sensitivity to hidden")
    #if kap_index is zero (meaning no kappas were found that improved robustness) then returns a message saying none were found
  }
  if(length(kap_index)>0){
    #if statement for if there were kappas that improved robustness
    return.kappa <- current.pop[,kap_index]
    #assigns the kappas that improved robustness
    return.correlations <- matrix(data = NA, nrow = ncol(X), ncol = ncol(return.kappa))
    #creates a matrix to represent the correlations between the satisfactory kappas and the given, observed covariates
    for(i in 1:ncol(return.kappa)){
      #runs a for-loop over all the satisfactory kappas
      for(z in 1:ncol(X)){
        #runs a for-loop over all the given, observed covariates
        return.correlations[z,i] <- cor(return.kappa[,i], X[,z])
        #assigns the calculated correlation between the kappa and covariate to the appropriate position in the matrix
      }
    }
    return.fitness <- fitness_vector[kap_index]
    #gets the fitness values of the satisfactory kappas
    return.list <- list(kappa = return.kappa, correlations = return.correlations, fitness = return.fitness)
    #creates a list to be outputed by main function. List contains the kappas, the correlations, and the fitness values.
    return(return.list)
    #returns the return list
  }
}

gencon <- genConfoundSearch(X = X, Y = Y, Tr = treat, Gamma = 3, search.pop.size = 25, match.pop.size = 25)
gencon$correlations
gencon$fitness
gencon$kappa
