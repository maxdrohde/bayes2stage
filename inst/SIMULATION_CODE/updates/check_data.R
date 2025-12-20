raw_file <- "results/1/draws_1.parquet"
data <- arrow::read_parquet(raw_file)
print("Original names (subset):")
print(grep("rhat", names(data), value = TRUE))

cleaned <- janitor::clean_names(data)
print("Cleaned names (subset):")
print(grep("rhat", names(cleaned), value = TRUE))
