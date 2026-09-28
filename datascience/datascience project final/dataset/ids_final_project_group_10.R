data <- read.csv("biomedical_research_abstracts_2024_2026.csv", stringsAsFactors = FALSE)

names(data)
dim(data)

selected_data <- data[, c("pmid", "title", "abstract", "abstract_words", "journal", "pub_year", "country", "keywords", "major_topic")]

selected_data <- selected_data[!is.na(selected_data$abstract), ]
selected_data <- selected_data[selected_data$abstract != "", ]

selected_data$title <- trimws(selected_data$title)
selected_data$abstract <- trimws(selected_data$abstract)
selected_data$keywords <- trimws(selected_data$keywords)
selected_data$major_topic <- trimws(selected_data$major_topic)

selected_data <- selected_data[!duplicated(selected_data$pmid), ]
selected_data <- selected_data[!duplicated(selected_data$title), ]

selected_data <- selected_data[selected_data$abstract_words >= 50 & selected_data$abstract_words <= 1000, ]

dim(selected_data)
head(selected_data)




text_all <- paste(selected_data$title, selected_data$abstract, selected_data$keywords, selected_data$major_topic)

diabetes_data <- selected_data[grepl("diabetes|insulin|glucose", text_all, ignore.case = TRUE), ]
cancer_data <- selected_data[grepl("cancer|tumor|oncology|carcinoma", text_all, ignore.case = TRUE), ]
cardio_data <- selected_data[grepl("cardiovascular|heart|hypertension|stroke", text_all, ignore.case = TRUE), ]
neuro_data <- selected_data[grepl("neurological|brain|alzheimer|parkinson|dementia", text_all, ignore.case = TRUE), ]

diabetes_data$disease_group <- "Diabetes"
cancer_data$disease_group <- "Cancer"
cardio_data$disease_group <- "Cardiovascular"
neuro_data$disease_group <- "Neurological"

nrow(diabetes_data)
nrow(cancer_data)
nrow(cardio_data)
nrow(neuro_data)





set.seed(42)

diabetes_sample <- diabetes_data[sample(1:nrow(diabetes_data), 50), ]
cancer_sample <- cancer_data[sample(1:nrow(cancer_data), 50), ]
cardio_sample <- cardio_data[sample(1:nrow(cardio_data), 50), ]
neuro_sample <- neuro_data[sample(1:nrow(neuro_data), 50), ]

project_data <- rbind(diabetes_sample, cancer_sample, cardio_sample, neuro_sample)

project_data <- project_data[sample(1:nrow(project_data)), ]

table(project_data$disease_group)
dim(project_data)
head(project_data)



clean_text <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  return(x)
}

project_data$clean_abstract <- clean_text(project_data$abstract)
head(project_data[, c("disease_group", "abstract", "clean_abstract")])

write.csv(project_data, "ids_final_dataset_sample_group_10.csv", row.names = FALSE)


install.packages("tm")
install.packages("SnowballC")
install.packages("cluster")
install.packages("uwot")





library(tm)
library(SnowballC)
library(cluster)
library(uwot)

data <- read.csv("ids_final_dataset_sample_group_10.csv", stringsAsFactors = FALSE)

dim(data)
names(data)
table(data$disease_group)

data$text <- data$clean_abstract
data$text <- ifelse(is.na(data$text), "", data$text)

corpus <- VCorpus(VectorSource(data$text))

corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, stopwords("english"))
corpus <- tm_map(corpus, stripWhitespace)

dtm <- DocumentTermMatrix(corpus)

dtm <- removeSparseTerms(dtm, 0.97)

tfidf <- weightTfIdf(dtm)

tfidf_matrix <- as.matrix(tfidf)

dim(tfidf_matrix)

set.seed(42)

k_value <- 4

k_model <- kmeans(tfidf_matrix, centers = k_value, nstart = 25)

data$cluster <- as.factor(k_model$cluster)

table(data$cluster)
table(data$disease_group, data$cluster)

cluster_distribution <- as.data.frame(table(data$cluster))
names(cluster_distribution) <- c("Cluster", "Document_Count")
cluster_distribution

cross_table <- as.data.frame.matrix(table(data$disease_group, data$cluster))
cross_table

distance_matrix <- dist(tfidf_matrix)

silhouette_result <- silhouette(k_model$cluster, distance_matrix)

mean_silhouette <- mean(silhouette_result[, 3])
mean_silhouette

plot(silhouette_result, main = "Silhouette Plot for K-Means Clustering")

top_terms <- function(cluster_number, matrix_data, dataset, top_n = 10) {
  rows <- dataset$cluster == cluster_number
  cluster_words <- colMeans(matrix_data[rows, , drop = FALSE])
  cluster_words <- sort(cluster_words, decreasing = TRUE)
  names(cluster_words)[1:top_n]
}

cluster_1_words <- top_terms("1", tfidf_matrix, data)
cluster_2_words <- top_terms("2", tfidf_matrix, data)
cluster_3_words <- top_terms("3", tfidf_matrix, data)
cluster_4_words <- top_terms("4", tfidf_matrix, data)

cluster_keywords <- data.frame(
  Cluster_1 = cluster_1_words,
  Cluster_2 = cluster_2_words,
  Cluster_3 = cluster_3_words,
  Cluster_4 = cluster_4_words,
  stringsAsFactors = FALSE
)

cluster_keywords

set.seed(42)

umap_result <- umap(tfidf_matrix, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

umap_data <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Cluster = data$cluster,
  Disease = data$disease_group
)

plot(
  umap_data$UMAP1,
  umap_data$UMAP2,
  col = as.numeric(umap_data$Cluster),
  pch = 19,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "UMAP Visualization of Biomedical Abstract Clusters"
)

legend(
  "topright",
  legend = levels(umap_data$Cluster),
  col = 1:length(levels(umap_data$Cluster)),
  pch = 19,
  title = "Cluster"
)

plot(
  umap_data$UMAP1,
  umap_data$UMAP2,
  col = as.numeric(as.factor(umap_data$Disease)),
  pch = 19,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "UMAP Visualization by Disease Group"
)

legend(
  "topright",
  legend = levels(as.factor(umap_data$Disease)),
  col = 1:length(levels(as.factor(umap_data$Disease))),
  pch = 19,
  title = "Disease"
)

cluster_summary <- aggregate(
  pmid ~ cluster + disease_group,
  data = data,
  FUN = length
)

names(cluster_summary) <- c("Cluster", "Disease_Group", "Document_Count")

cluster_summary

write.csv(data, "ids_final_clustered_dataset_group_10.csv", row.names = FALSE)
write.csv(cluster_distribution, "cluster_distribution_group_10.csv", row.names = FALSE)
write.csv(cluster_keywords, "cluster_keywords_group_10.csv", row.names = FALSE)
write.csv(cluster_summary, "cluster_disease_summary_group_10.csv", row.names = FALSE)