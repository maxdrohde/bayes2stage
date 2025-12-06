library(bayes2stage)
library(bench)

dataset <- bayes2stage::generate_data(N=5000)

f1 <- function(dataset){
  as.matrix(bayes2stage:::get_ods(dataset, "intercept"))
}

f2 <- function(dataset){
  as.matrix(bayes2stage:::get_ods_old(dataset, "intercept"))
}

benchmark_results <-
  bench::mark(`new` = f1(dataset),
              `old` = f2(dataset),
              iterations = 50,
              #time_unit = "m",
              memory = TRUE,
              check = TRUE,
              filter_gc = TRUE)

ggplot2::autoplot(benchmark_results)

