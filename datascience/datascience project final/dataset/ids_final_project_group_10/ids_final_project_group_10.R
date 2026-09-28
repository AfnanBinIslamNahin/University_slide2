library(tm)
library(cluster)
library(uwot)
library(clusterSim)

raw_data <- read.csv("biomedical_research_abstracts_2024_2026.csv", stringsAsFactors = FALSE)
dim(raw_data )
head(raw_data )

selected_data <- raw_data[, c("pmid", "title", "abstract", "abstract_words", "journal", "pub_year", "country", "keywords", "major_topic")]

selected_data <- selected_data[!is.na(selected_data$abstract), ]
selected_data <- selected_data[selected_data$abstract != "", ]
dim(selected_data)
sum(is.na(selected_data$abstract))

selected_data$title <- trimws(selected_data$title)
selected_data$abstract <- trimws(selected_data$abstract)
selected_data$keywords <- ifelse(is.na(selected_data$keywords), "", selected_data$keywords)
selected_data$major_topic <- ifelse(is.na(selected_data$major_topic), "", selected_data$major_topic)
head(selected_data[, c("title", "keywords", "major_topic")])

selected_data <- selected_data[!duplicated(selected_data$pmid), ]
selected_data <- selected_data[!duplicated(selected_data$title), ]

selected_data <- selected_data[selected_data$abstract_words >= 80 & selected_data$abstract_words <= 350, ]

all_text <- paste(selected_data$title, selected_data$abstract, selected_data$keywords, selected_data$major_topic)
all_text <- tolower(all_text)

diabetes_pattern <- "diabetes|insulin|glucose|glycemic"
cancer_pattern <- "cancer|tumor|tumour|oncology|carcinoma|chemotherapy"
cardio_pattern <- "cardiovascular|heart disease|hypertension|myocardial|coronary"
neuro_pattern <- "neurological|brain|alzheimer|parkinson|dementia|stroke"

diabetes_data <- selected_data[
  grepl(diabetes_pattern, all_text) &
    !grepl(cancer_pattern, all_text) &
    !grepl(cardio_pattern, all_text) &
    !grepl(neuro_pattern, all_text),
]

cancer_data <- selected_data[
  grepl(cancer_pattern, all_text) &
    !grepl(diabetes_pattern, all_text) &
    !grepl(cardio_pattern, all_text) &
    !grepl(neuro_pattern, all_text),
]

cardio_data <- selected_data[
  grepl(cardio_pattern, all_text) &
    !grepl(diabetes_pattern, all_text) &
    !grepl(cancer_pattern, all_text) &
    !grepl(neuro_pattern, all_text),
]

neuro_data <- selected_data[
  grepl(neuro_pattern, all_text) &
    !grepl(diabetes_pattern, all_text) &
    !grepl(cancer_pattern, all_text) &
    !grepl(cardio_pattern, all_text),
]

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

clean_text <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  return(x)
}

project_data$clean_abstract <- clean_text(project_data$abstract)

table(project_data$disease_group)
dim(project_data)
head(project_data[, c("disease_group", "clean_abstract")])
write.csv(project_data, "ids_final_dataset_group_10.csv", row.names = FALSE)

data <- read.csv("ids_final_dataset_group_10.csv", stringsAsFactors = FALSE)

dim(data)
names(data)
table(data$disease_group)
head(data[, c("title", "disease_group", "clean_abstract")])

data$text <- data$clean_abstract
data$text <- ifelse(is.na(data$text), "", data$text)

corpus <- VCorpus(VectorSource(data$text))

corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, stopwords("english"))
corpus <- tm_map(corpus, removeWords, c("study","patient","patients","result","results","method","methods","conclusion","background","objective"))
corpus <- tm_map(corpus, stripWhitespace)

as.character(corpus[[1]])

dtm <- DocumentTermMatrix(corpus)

dtm <- removeSparseTerms(dtm, 0.95)

tfidf <- weightTfIdf(dtm)

tfidf_matrix <- as.matrix(tfidf)

dim(tfidf_matrix)
tfidf_matrix[1:5, 1:5]

set.seed(42)
pca_result <- prcomp(tfidf_matrix, center = TRUE, scale. = FALSE)
pca_data <- pca_result$x[, 1:20]
dim(pca_data)


set.seed(42)

k_values <- 2:8
sil_scores <- numeric(length(k_values))
db_scores <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  km <- kmeans(pca_data, centers = k_values[i], nstart = 50)
  sil <- silhouette(km$cluster, dist(pca_data))
  db <- index.DB(x = pca_data, cl = km$cluster, centrotypes = "centroids")
  sil_scores[i] <- mean(sil[, 3])
  db_scores[i] <- db$DB
}

evaluation_by_k <- data.frame(
  K = k_values,
  Silhouette_Score = round(sil_scores, 4),
  Davies_Bouldin_Index = round(db_scores, 4)
)

evaluation_by_k

best_k <- evaluation_by_k$K[which.max(evaluation_by_k$Silhouette_Score)]

best_k

plot(
  evaluation_by_k$K,
  evaluation_by_k$Silhouette_Score,
  type = "b",
  pch = 19,
  lwd = 2,
  xlab = "Number of Clusters (K)",
  ylab = "Average Silhouette Score",
  main = "Silhouette Score for Different K Values"
)
grid()

points(
  best_k,
  max(evaluation_by_k$Silhouette_Score),
  pch = 19,
  cex = 1.5,
  col="red"
)

text(
  best_k,
  max(evaluation_by_k$Silhouette_Score),
  labels = paste("Best K =", best_k),
  pos = 3,
  col="red"
)




best_k_db <- evaluation_by_k$K[which.min(evaluation_by_k$Davies_Bouldin_Index)]
plot(
  evaluation_by_k$K,
  evaluation_by_k$Davies_Bouldin_Index,
  type = "b",
  pch = 19,
  lwd = 2,
  col="darkgreen",
  xlab = "Number of Clusters (K)",
  ylab = "Davies-Bouldin Index",
  main = "Davies-Bouldin Index for Different K Values"
)
grid()
points(
  best_k_db,
  min(evaluation_by_k$Davies_Bouldin_Index),
  pch = 19,
  cex = 1.6,
  col = "red"
)

text(
  best_k_db,
  min(evaluation_by_k$Davies_Bouldin_Index),
  labels = paste("Lowest DBI K =", best_k_db),
  pos = 3,
  col = "red"
)

best_k <- evaluation_by_k$K[which.max(evaluation_by_k$Silhouette_Score)]
best_k



set.seed(42)

final_model <- kmeans(pca_data, centers = best_k, nstart = 50)

data$final_cluster <- as.factor(final_model$cluster)

table(data$final_cluster)

table(data$disease_group, data$final_cluster)






final_silhouette <- silhouette(final_model$cluster, dist(pca_data))

final_mean_silhouette <- mean(final_silhouette[, 3])

final_db <- index.DB(
  x = pca_data,
  cl = final_model$cluster,
  centrotypes = "centroids"
)

final_mean_silhouette

final_db$DB

evaluation_table <- data.frame(
  Metric = c("Silhouette Score", "Davies-Bouldin Index"),
  Value = c(round(final_mean_silhouette, 4), round(final_db$DB, 4)),
  Interpretation = c(
    "Higher value is better",
    "Lower value is better"
  ),
  stringsAsFactors = FALSE
)

evaluation_table




plot(
  final_silhouette,
  main = "Final Silhouette Plot",
  col = c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan"),
  border = NA,
  cex.names = 0.75,
  do.n.k = TRUE,
  do.clus.stat = TRUE
)

top_terms <- function(cluster_number, matrix_data, dataset, cluster_column, top_n = 10) {
  rows <- dataset[[cluster_column]] == cluster_number
  cluster_words <- colMeans(matrix_data[rows, , drop = FALSE])
  cluster_words <- sort(cluster_words, decreasing = TRUE)
  names(cluster_words)[1:top_n]
}

cluster_keywords <- data.frame(
  Cluster_1 = top_terms("1", tfidf_matrix, data, "final_cluster"),
  Cluster_2 = top_terms("2", tfidf_matrix, data, "final_cluster"),
  Cluster_3 = top_terms("3", tfidf_matrix, data, "final_cluster"),
  Cluster_4 = top_terms("4", tfidf_matrix, data, "final_cluster"),
  Cluster_5 = top_terms("5", tfidf_matrix, data, "final_cluster"),
  Cluster_6 = top_terms("6", tfidf_matrix, data, "final_cluster"),
  Cluster_7 = top_terms("7", tfidf_matrix, data, "final_cluster"),
  Cluster_8 = top_terms("8", tfidf_matrix, data, "final_cluster"),
  stringsAsFactors = FALSE
)

cluster_keywords

cluster_summary <- aggregate(
  pmid ~ final_cluster + disease_group,
  data = data,
  FUN = length
)

names(cluster_summary) <- c("Cluster", "Disease_Group", "Document_Count")

cluster_summary

set.seed(42)

umap_result <- umap(tfidf_matrix, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

umap_data <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Cluster = data$final_cluster,
  Disease = data$disease_group
)

plot(
  umap_data$UMAP1,
  umap_data$UMAP2,
  col = as.numeric(umap_data$Cluster),
  pch = 19,
  cex = 1.2,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "UMAP Visualization of Disease-Specific Clusters"
)

legend(
  "topright",
  legend = levels(umap_data$Cluster),
  col = 1:length(levels(umap_data$Cluster)),
  pch = 19,
  title = "Cluster",
  cex = 0.8
)

