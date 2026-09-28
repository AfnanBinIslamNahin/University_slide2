library(tm)
library(cluster)
library(uwot)
library(clusterSim)

raw_data <- read.csv("biomedical_research_abstracts_2024_2026.csv", stringsAsFactors = FALSE)

dim(raw_data)
head(raw_data)

selected_data <- raw_data[, c("pmid", "title", "abstract", "abstract_words", "journal", "pub_year", "country", "keywords", "major_topic")]

selected_data <- selected_data[!is.na(selected_data$abstract), ]
selected_data <- selected_data[selected_data$abstract != "", ]

selected_data$title <- trimws(selected_data$title)
selected_data$abstract <- trimws(selected_data$abstract)

selected_data$keywords <- ifelse(is.na(selected_data$keywords), "", selected_data$keywords)
selected_data$major_topic <- ifelse(is.na(selected_data$major_topic), "", selected_data$major_topic)

selected_data$keywords <- trimws(selected_data$keywords)
selected_data$major_topic <- trimws(selected_data$major_topic)

selected_data <- selected_data[selected_data$keywords != "" & selected_data$major_topic != "", ]

selected_data <- selected_data[!duplicated(selected_data$pmid), ]
selected_data <- selected_data[!duplicated(selected_data$title), ]

sum(selected_data$abstract_words >= 80 & selected_data$abstract_words <= 350)
sum(selected_data$abstract_words >= 100 & selected_data$abstract_words <= 400)
sum(selected_data$abstract_words >= 100 & selected_data$abstract_words <= 500)
sum(selected_data$abstract_words >= 150 & selected_data$abstract_words <= 400)

selected_data <- selected_data[selected_data$abstract_words >= 100 & selected_data$abstract_words <= 400, ]

dim(selected_data)
colSums(selected_data == "")
colSums(is.na(selected_data))

all_text <- paste(selected_data$title, selected_data$abstract, selected_data$keywords, selected_data$major_topic)
all_text <- tolower(all_text)

diabetes_pattern <- "diabetes|type 2 diabetes|type ii diabetes|insulin|glycemic|glucose|diabetic|hba1c|metabolic"
cancer_pattern <- "cancer|tumor|tumour|oncology|carcinoma|metastatic|chemotherapy|neoplasm"
cardio_pattern <- "cardiovascular|coronary|myocardial|heart failure|hypertension|cardiac|artery"
neuro_pattern <- "neurological|neurodegenerative|alzheimer|parkinson|dementia|stroke|brain"

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

sample_size <- 250

set.seed(42)

diabetes_sample <- diabetes_data[sample(1:nrow(diabetes_data), sample_size), ]
cancer_sample <- cancer_data[sample(1:nrow(cancer_data), sample_size), ]
cardio_sample <- cardio_data[sample(1:nrow(cardio_data), sample_size), ]
neuro_sample <- neuro_data[sample(1:nrow(neuro_data), sample_size), ]

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

colSums(project_data == "")
colSums(is.na(project_data))

head(project_data[, c("disease_group", "keywords", "major_topic", "clean_abstract")])

write.csv(project_data, "ids_final_dataset.csv", row.names = FALSE)


;;;;;;;;;;;;;;;;;
library(tm)
library(cluster)
library(uwot)
library(clusterSim)
library(dbscan)
library(Rtsne)




data <- read.csv("ids_final_dataset.csv", stringsAsFactors = FALSE)

dim(data)
names(data)
table(data$disease_group)
head(data[, c("title", "disease_group", "clean_abstract")])




data$text <- data$clean_abstract
data$text <- ifelse(is.na(data$text), "", data$text)

corpus <- VCorpus(VectorSource(data$text))




custom_stopwords <- c(
  stopwords("english"),
  "study", "studies", "patient", "patients",
  "result", "results", "method", "methods",
  "conclusion", "background", "objective",
  "clinical", "analysis", "data", "research",
  "significant", "associated", "using", "used",
  "based", "group", "groups", "effect", "effects",
  "treatment", "medical", "case", "cases",
  "review", "systematic", "among", "within",
  "included", "compared", "found", "showed"
)

corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, custom_stopwords)
corpus <- tm_map(corpus, stripWhitespace)

as.character(corpus[[1]])





dtm <- DocumentTermMatrix(corpus)

dtm <- removeSparseTerms(dtm, 0.98)

tfidf <- weightTfIdf(dtm)

tfidf_matrix <- as.matrix(tfidf)

dim(tfidf_matrix)
tfidf_matrix[1:5, 1:5]



set.seed(42)

pca_result <- prcomp(tfidf_matrix, center = TRUE, scale. = FALSE)

pca_data <- pca_result$x[, 1:40]

pca_data <- scale(pca_data)

dim(pca_data)




set.seed(42)

kmeans_model_4 <- kmeans(
  pca_data,
  centers = 4,
  nstart = 100
)

data$kmeans_cluster_4 <- as.factor(kmeans_model_4$cluster)

table(data$kmeans_cluster_4)

disease_cluster_table <- table(data$disease_group, data$kmeans_cluster_4)

disease_cluster_table
.........................................................

library(tm)
library(cluster)
library(uwot)
library(clusterSim)
library(dbscan)
library(Rtsne)

raw_data <- read.csv("biomedical_research_abstracts_2024_2026.csv", stringsAsFactors = FALSE)

dim(raw_data)
head(raw_data)

selected_data <- raw_data[, c("pmid", "title", "abstract", "abstract_words", "journal", "pub_year", "country", "keywords", "major_topic")]

selected_data <- selected_data[!is.na(selected_data$abstract), ]
selected_data <- selected_data[selected_data$abstract != "", ]

selected_data$title <- trimws(selected_data$title)
selected_data$abstract <- trimws(selected_data$abstract)

selected_data$keywords <- ifelse(is.na(selected_data$keywords), "", selected_data$keywords)
selected_data$major_topic <- ifelse(is.na(selected_data$major_topic), "", selected_data$major_topic)

selected_data$keywords <- trimws(selected_data$keywords)
selected_data$major_topic <- trimws(selected_data$major_topic)

selected_data <- selected_data[selected_data$keywords != "" & selected_data$major_topic != "", ]

selected_data <- selected_data[!duplicated(selected_data$pmid), ]
selected_data <- selected_data[!duplicated(selected_data$title), ]

range_check <- data.frame(
  Range = c("80-350", "100-400", "100-500", "150-400"),
  Count = c(
    sum(selected_data$abstract_words >= 80 & selected_data$abstract_words <= 350),
    sum(selected_data$abstract_words >= 100 & selected_data$abstract_words <= 400),
    sum(selected_data$abstract_words >= 100 & selected_data$abstract_words <= 500),
    sum(selected_data$abstract_words >= 150 & selected_data$abstract_words <= 400)
  )
)

range_check

selected_data <- selected_data[selected_data$abstract_words >= 100 & selected_data$abstract_words <= 400, ]

dim(selected_data)

blank_counts_selected <- sapply(selected_data, function(x) sum(trimws(as.character(x)) == "", na.rm = TRUE))
missing_counts_selected <- sapply(selected_data, function(x) sum(is.na(x)))

blank_counts_selected
missing_counts_selected

all_text <- paste(selected_data$title, selected_data$abstract, selected_data$keywords, selected_data$major_topic)
all_text <- tolower(all_text)

diabetes_pattern <- "diabetes|type 2 diabetes|type ii diabetes|insulin|glycemic|glucose|diabetic|hba1c"
cancer_pattern <- "cancer|tumor|tumour|oncology|carcinoma|metastatic|chemotherapy|neoplasm"
cardio_pattern <- "cardiovascular|coronary|myocardial|heart failure|hypertension|cardiac|artery"
neuro_pattern <- "neurological|neurodegenerative|alzheimer|parkinson|dementia|stroke"

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

available_documents <- data.frame(
  Disease_Group = c("Diabetes", "Cancer", "Cardiovascular", "Neurological"),
  Available_Documents = c(
    nrow(diabetes_data),
    nrow(cancer_data),
    nrow(cardio_data),
    nrow(neuro_data)
  )
)

available_documents

sample_size <- 250

if (nrow(diabetes_data) < sample_size) stop("Not enough Diabetes documents")
if (nrow(cancer_data) < sample_size) stop("Not enough Cancer documents")
if (nrow(cardio_data) < sample_size) stop("Not enough Cardiovascular documents")
if (nrow(neuro_data) < sample_size) stop("Not enough Neurological documents")

set.seed(42)

diabetes_sample <- diabetes_data[sample(1:nrow(diabetes_data), sample_size), ]
cancer_sample <- cancer_data[sample(1:nrow(cancer_data), sample_size), ]
cardio_sample <- cardio_data[sample(1:nrow(cardio_data), sample_size), ]
neuro_sample <- neuro_data[sample(1:nrow(neuro_data), sample_size), ]

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

blank_counts_project <- sapply(project_data, function(x) sum(trimws(as.character(x)) == "", na.rm = TRUE))
missing_counts_project <- sapply(project_data, function(x) sum(is.na(x)))

blank_counts_project
missing_counts_project

head(project_data[, c("disease_group", "keywords", "major_topic", "clean_abstract")])

write.csv(project_data, "final_dataset.csv", row.names = FALSE)

data <- read.csv("ids_final_dataset.csv", stringsAsFactors = FALSE)

dim(data)
names(data)
table(data$disease_group)
head(data[, c("title", "disease_group", "clean_abstract")])

data$text <- data$clean_abstract
data$text <- ifelse(is.na(data$text), "", data$text)

corpus <- VCorpus(VectorSource(data$text))

custom_stopwords <- c(
  stopwords("english"),
  "study", "studies", "patient", "patients",
  "result", "results", "method", "methods",
  "conclusion", "background", "objective",
  "clinical", "analysis", "data", "research",
  "significant", "associated", "using", "used",
  "based", "group", "groups", "effect", "effects",
  "treatment", "medical", "case", "cases",
  "review", "systematic", "among", "within",
  "included", "compared", "found", "showed"
)

corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, custom_stopwords)
corpus <- tm_map(corpus, stripWhitespace)

as.character(corpus[[1]])

dtm <- DocumentTermMatrix(corpus)

dtm <- removeSparseTerms(dtm, 0.98)

tfidf <- weightTfIdf(dtm)

tfidf_matrix <- as.matrix(tfidf)

dim(tfidf_matrix)
tfidf_matrix[1:5, 1:5]

set.seed(42)

pca_result <- prcomp(tfidf_matrix, center = TRUE, scale. = FALSE)

pca_data <- pca_result$x[, 1:40]

pca_data <- scale(pca_data)

dim(pca_data)

set.seed(42)

kmeans_model_4 <- kmeans(
  pca_data,
  centers = 4,
  nstart = 100
)

data$kmeans_cluster_4 <- as.factor(kmeans_model_4$cluster)

cluster_distribution <- as.data.frame(table(data$kmeans_cluster_4))
names(cluster_distribution) <- c("Cluster", "Document_Count")

cluster_distribution

disease_cluster_table <- table(data$disease_group, data$kmeans_cluster_4)

disease_cluster_table
........

.........

disease_cluster_percent <- round(
  prop.table(disease_cluster_table, margin = 1) * 100,
  2
)

disease_cluster_percent

cluster_disease_table <- table(data$kmeans_cluster_4, data$disease_group)

cluster_disease_table

majority_disease <- apply(cluster_disease_table, 1, function(x) {
  names(which.max(x))
})

majority_count <- apply(cluster_disease_table, 1, max)

cluster_size <- rowSums(cluster_disease_table)

cluster_purity <- round((majority_count / cluster_size) * 100, 2)

cluster_interpretation <- data.frame(
  Cluster = names(majority_disease),
  Majority_Disease = majority_disease,
  Majority_Count = majority_count,
  Cluster_Size = cluster_size,
  Purity_Percentage = cluster_purity,
  stringsAsFactors = FALSE
)

cluster_interpretation

overall_purity <- round((sum(majority_count) / nrow(data)) * 100, 2)

overall_purity

kmeans_4_silhouette <- silhouette(kmeans_model_4$cluster, dist(pca_data))

kmeans_4_mean_silhouette <- mean(kmeans_4_silhouette[, 3])

kmeans_4_db <- index.DB(
  x = pca_data,
  cl = kmeans_model_4$cluster,
  centrotypes = "centroids"
)

kmeans_4_mean_silhouette
kmeans_4_db$DB

kmeans_4_evaluation <- data.frame(
  Algorithm = "K-Means",
  K = 4,
  PCA_Components = 40,
  Nstart = 100,
  Silhouette_Score = round(kmeans_4_mean_silhouette, 4),
  Davies_Bouldin_Index = round(kmeans_4_db$DB, 4),
  Overall_Purity_Percentage = overall_purity,
  stringsAsFactors = FALSE
)

kmeans_4_evaluation

plot(
  kmeans_4_silhouette,
  main = "Silhouette Plot for K-Means Clustering with K = 4",
  col = c("red", "blue", "green", "purple"),
  border = NA,
  cex.names = 0.8,
  do.n.k = TRUE,
  do.clus.stat = TRUE
)

top_terms <- function(cluster_number, matrix_data, dataset, cluster_column, top_n = 10) {
  rows <- dataset[[cluster_column]] == cluster_number
  cluster_words <- colMeans(matrix_data[rows, , drop = FALSE])
  cluster_words <- sort(cluster_words, decreasing = TRUE)
  names(cluster_words)[1:top_n]
}

kmeans_4_keywords <- data.frame(
  Term_Rank = 1:10,
  Cluster_1 = top_terms("1", tfidf_matrix, data, "kmeans_cluster_4"),
  Cluster_2 = top_terms("2", tfidf_matrix, data, "kmeans_cluster_4"),
  Cluster_3 = top_terms("3", tfidf_matrix, data, "kmeans_cluster_4"),
  Cluster_4 = top_terms("4", tfidf_matrix, data, "kmeans_cluster_4"),
  stringsAsFactors = FALSE
)

kmeans_4_keywords

kmeans_cluster_summary <- aggregate(
  pmid ~ kmeans_cluster_4 + disease_group,
  data = data,
  FUN = length
)

names(kmeans_cluster_summary) <- c("Cluster", "Disease_Group", "Document_Count")

kmeans_cluster_summary

set.seed(42)

umap_result_4 <- umap(
  pca_data,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean"
)

umap_kmeans_4_data <- data.frame(
  UMAP1 = umap_result_4[, 1],
  UMAP2 = umap_result_4[, 2],
  Cluster = data$kmeans_cluster_4,
  Disease = data$disease_group
)

plot(
  umap_kmeans_4_data$UMAP1,
  umap_kmeans_4_data$UMAP2,
  col = as.numeric(umap_kmeans_4_data$Cluster),
  pch = 19,
  cex = 1.1,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "UMAP Visualization of K-Means Clusters with K = 4"
)

legend(
  "topright",
  legend = levels(umap_kmeans_4_data$Cluster),
  col = 1:length(levels(umap_kmeans_4_data$Cluster)),
  pch = 19,
  title = "Cluster",
  cex = 0.8
)

plot(
  umap_kmeans_4_data$UMAP1,
  umap_kmeans_4_data$UMAP2,
  col = as.numeric(as.factor(umap_kmeans_4_data$Disease)),
  pch = 19,
  cex = 1.1,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "UMAP Visualization by Actual Disease Group"
)

legend(
  "topright",
  legend = levels(as.factor(umap_kmeans_4_data$Disease)),
  col = 1:length(levels(as.factor(umap_kmeans_4_data$Disease))),
  pch = 19,
  title = "Disease Group",
  cex = 0.7
)

set.seed(42)

tsne_result_4 <- Rtsne(
  pca_data,
  dims = 2,
  perplexity = 30,
  verbose = TRUE,
  max_iter = 500
)

tsne_kmeans_4_data <- data.frame(
  TSNE1 = tsne_result_4$Y[, 1],
  TSNE2 = tsne_result_4$Y[, 2],
  Cluster = data$kmeans_cluster_4,
  Disease = data$disease_group
)

plot(
  tsne_kmeans_4_data$TSNE1,
  tsne_kmeans_4_data$TSNE2,
  col = as.numeric(tsne_kmeans_4_data$Cluster),
  pch = 19,
  cex = 1.1,
  xlab = "t-SNE 1",
  ylab = "t-SNE 2",
  main = "t-SNE Visualization of K-Means Clusters with K = 4"
)

legend(
  "topright",
  legend = levels(tsne_kmeans_4_data$Cluster),
  col = 1:length(levels(tsne_kmeans_4_data$Cluster)),
  pch = 19,
  title = "Cluster",
  cex = 0.8
)

set.seed(42)

hdbscan_data <- pca_result$x[, 1:10]
hdbscan_data <- scale(hdbscan_data)

minpts_values <- c(3, 5, 7, 10, 15)

hdbscan_results <- data.frame(
  MinPts = integer(),
  Cluster_Count = integer(),
  Noise_Points = integer(),
  Silhouette_Score = numeric(),
  Davies_Bouldin_Index = numeric(),
  stringsAsFactors = FALSE
)

for (mp in minpts_values) {
  hdb_model <- dbscan::hdbscan(hdbscan_data, minPts = mp)
  hdb_labels <- hdb_model$cluster
  
  non_noise_index <- hdb_labels != 0
  cluster_count <- length(setdiff(unique(hdb_labels), 0))
  noise_count <- sum(hdb_labels == 0)
  
  if (cluster_count >= 2 && sum(non_noise_index) > cluster_count) {
    hdb_sil <- silhouette(
      hdb_labels[non_noise_index],
      dist(hdbscan_data[non_noise_index, , drop = FALSE])
    )
    
    hdb_sil_score <- mean(hdb_sil[, 3])
    
    hdb_db <- index.DB(
      x = hdbscan_data[non_noise_index, , drop = FALSE],
      cl = hdb_labels[non_noise_index],
      centrotypes = "centroids"
    )
    
    hdb_db_score <- hdb_db$DB
  } else {
    hdb_sil_score <- NA
    hdb_db_score <- NA
  }
  
  hdbscan_results <- rbind(
    hdbscan_results,
    data.frame(
      MinPts = mp,
      Cluster_Count = cluster_count,
      Noise_Points = noise_count,
      Silhouette_Score = round(hdb_sil_score, 4),
      Davies_Bouldin_Index = round(hdb_db_score, 4)
    )
  )
}

hdbscan_results

write.csv(data, "ids_final_clustered_dataset_group_10.csv", row.names = FALSE)
write.csv(cluster_distribution, "cluster_distribution_group_10.csv", row.names = FALSE)
write.csv(as.data.frame.matrix(disease_cluster_table), "disease_cluster_table_group_10.csv")
write.csv(as.data.frame.matrix(disease_cluster_percent), "disease_cluster_percent_group_10.csv")
write.csv(cluster_interpretation, "cluster_interpretation_group_10.csv", row.names = FALSE)
write.csv(kmeans_4_evaluation, "kmeans_4_evaluation_group_10.csv", row.names = FALSE)
write.csv(kmeans_4_keywords, "kmeans_4_keywords_group_10.csv", row.names = FALSE)
write.csv(kmeans_cluster_summary, "kmeans_cluster_summary_group_10.csv", row.names = FALSE)
write.csv(hdbscan_results, "hdbscan_evaluation_group_10.csv", row.names = FALSE)



