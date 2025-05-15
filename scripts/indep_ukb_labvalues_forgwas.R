## find independent lab values for UKB GWAS

library(tidyverse)
library(data.table)
library(ggplot2)
library(caret) # for correlation analysis
library(SNPRelate) # helpful for genetic correlation
library(stringr)  # For string manipulation

# Load the data
biochem <- fread('/sc/arion/projects/paul_oreilly/data/shared/blood_biochem_trait.raw')
count <- fread('/sc/arion/projects/paul_oreilly/data/shared/blood_count_trait.raw')

# Examine the data structure
print("Biochem data structure:")
str(biochem)
print("Blood count data structure:")
str(count)

# Clean the data - remove missing values, output pheno file
biochem_clean <- biochem #%>% select_if(~sum(!is.na(.)) > 0.7 * nrow(biochem))
count_clean <- count #%>% select_if(~sum(!is.na(.)) > 0.7 * nrow(count))
master <- merge(biochem_clean, count_clean, by = "sample_id")
master$IID <- master$sample_id
master$FID <- master$sample_id
master$sample_id <- NULL
master <- master %>% select(FID, IID, everything()) 
write.table(master, file = "/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/master_labvalues.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

# Read in Neale lab h2 UKB file
h2 <- fread('/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/results/heritability/ukb31063_h2_topline.02Oct2019.tsv.gz')
# Fix typos in the file
h2$description <- gsub("Apoliprotein A", "Apolipoprotein A", h2$description, fixed=TRUE)
h2$description <- gsub("Apoliprotein B", "Apolipoprotein B", h2$description, fixed=TRUE)
# filter for nominally significant h2, h2>0.1
h2 <- h2[h2$h2_p < 0.05, ]
h2 <- h2[h2$h2_observed>0.1, ]


# Get list of all trait names in our master file (excluding FID, IID)
master_traits <- names(master)[!names(master) %in% c("FID", "IID")]

# Create a mapping between master trait names and Neale lab descriptions
# Convert underscore to space and remove units from Neale descriptions
  # Get all biochem and count trait names (excluding sample_id)
  biochem_traits <- names(biochem_clean)[names(biochem_clean) != "sample_id"]
  count_traits <- names(count_clean)[names(count_clean) != "sample_id"]
  
  # Create lookup table for matching between formats
  trait_lookup <- data.frame(
    master_name = c(biochem_traits, count_traits),
    master_name_spaces = gsub("_", " ", c(biochem_traits, count_traits)),
    stringsAsFactors = FALSE
  )
  
  # Add cleaned Neale lab descriptions (for matching)
  h2$clean_description <- gsub(" \\(.*\\)$", "", h2$description)  # Remove units in parentheses
  h2$clean_description <- tolower(h2$clean_description)  # Convert to lowercase
  
  # Convert master names to lowercase for matching
  trait_lookup$master_name_spaces_lower <- tolower(trait_lookup$master_name_spaces)
  
  # Create the final mapping
  mapping <- data.frame(
    master_name = character(),
    h2_phenotype = character(),
    h2_description = character(),
    h2_value = numeric(),
    h2_se = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each trait in our master file, find the best match in h2 data
  for(i in 1:nrow(trait_lookup)) {
    trait_name <- trait_lookup$master_name[i]
    trait_name_spaces_lower <- trait_lookup$master_name_spaces_lower[i]
    
    # Try exact matching first
    h2_match <- h2[tolower(h2$clean_description) == trait_name_spaces_lower, ]
    
    # If no exact match, try partial matching
    if(nrow(h2_match) == 0) {
      # Check if trait name is a substring of any description
      for(j in 1:nrow(h2)) {
        # Use fixed=TRUE to treat the strings as literal text, not regex patterns
        if(grepl(trait_name_spaces_lower, tolower(h2$clean_description[j]), fixed=TRUE) ||
           grepl(tolower(h2$clean_description[j]), trait_name_spaces_lower, fixed=TRUE)) {
          h2_match <- rbind(h2_match, h2[j,])
        }
      }
    }
    
    # If still no match, try even more fuzzy matching
    if(nrow(h2_match) == 0) {
      # Split trait name into words and check if all words appear in description
      words <- unlist(strsplit(trait_name_spaces_lower, " "))
      if(length(words) > 1) {
        for(j in 1:nrow(h2)) {
          # Use fixed=TRUE here as well
          if(all(sapply(words, function(w) grepl(w, tolower(h2$clean_description[j]), fixed=TRUE)))) {
            h2_match <- rbind(h2_match, h2[j,])
          }
        }
      }
    }
    
    # If we found matches, add to mapping
    if(nrow(h2_match) > 0) {
      # If multiple matches, take the one with highest h2
      if(nrow(h2_match) > 1) {
        h2_match <- h2_match[which.max(h2_match$h2_observed),]
      }
      
      mapping <- rbind(mapping, data.frame(
        master_name = trait_name,
        h2_phenotype = h2_match$phenotype,
        h2_description = h2_match$description,
        h2_value = h2_match$h2_observed,
        h2_se = h2_match$h2_liability_se,
        stringsAsFactors = FALSE
      ))
    }
  }

# Print the mapping results
cat("Mapped", nrow(mapping), "out of", length(master_traits), "traits to Neale lab heritability estimates\n")
print(head(mapping, 10))

# Save the mapping for reference
write.csv(mapping, "trait_heritability_mapping.csv", row.names = FALSE)

# Get numeric columns (excluding sample_id)
biochem_numeric <- biochem_clean %>% 
  select_if(is.numeric) %>% 
  select(-sample_id)

count_numeric <- count_clean %>% 
  select_if(is.numeric) %>% 
  select(-sample_id)

# Create h2 dataframes using the mapping
biochem_h2 <- mapping %>%
  filter(master_name %in% names(biochem_numeric)) %>%
  select(trait = master_name, h2 = h2_value, h2_se = h2_se) %>%
  # Remove traits with missing h2
  filter(!is.na(h2))

count_h2 <- mapping %>%
  filter(master_name %in% names(count_numeric)) %>%
  select(trait = master_name, h2 = h2_value, h2_se = h2_se) %>%
  # Remove traits with missing h2
  filter(!is.na(h2))

# Print summary of heritability data
cat("Biochemical traits with heritability estimates:", nrow(biochem_h2), "out of", ncol(biochem_numeric), "\n")
cat("Blood count traits with heritability estimates:", nrow(count_h2), "out of", ncol(count_numeric), "\n")

# Step 6: Clump traits based on correlation and heritability using combined approach
# First, combine the heritability datasets
combined_h2 <- bind_rows(biochem_h2, count_h2)

# Create a correlation matrix using the master dataframe
# Since master already has all variables combined with FID and IID columns
all_numeric_data <- master %>% 
  # Remove the ID columns
  select(-FID, -IID) %>%
  # Keep only columns that have h2 estimates
  select(any_of(combined_h2$trait))

# Calculate the correlation matrix
combined_cor <- cor(all_numeric_data, use = "pairwise.complete.obs")

# Convert to long format for clumping
combined_cor_long <- combined_cor %>%
  as.data.frame() %>%
  rownames_to_column("trait1") %>%
  pivot_longer(-trait1, names_to = "trait2", values_to = "correlation")

# Define the clumping function
clump_traits <- function(cor_data, h2_data, threshold) {
  # Start with all traits
  all_traits <- unique(h2_data$trait)
  selected_traits <- character()
  remaining_traits <- all_traits
  
  while(length(remaining_traits) > 0) {
    # Find trait with highest heritability among remaining traits
    remaining_h2 <- h2_data %>% filter(trait %in% remaining_traits)
    best_trait <- remaining_h2$trait[which.max(remaining_h2$h2)]
    
    # Add to selected traits
    selected_traits <- c(selected_traits, best_trait)
    
    # Find correlated traits
    correlated <- cor_data %>% 
      filter(trait1 == best_trait, trait2 %in% remaining_traits, abs(correlation) > threshold) %>%
      pull(trait2)
    
    # Remove best trait and correlated traits from remaining
    remaining_traits <- setdiff(remaining_traits, c(best_trait, correlated))
  }
  
  return(selected_traits)
}

# Apply clumping to the combined variables
combined_selected <- clump_traits(combined_cor_long, combined_h2, 0.1)

# Separate the results back to biochem and count for reporting
biochem_selected <- combined_selected[combined_selected %in% biochem_h2$trait]
count_selected <- combined_selected[combined_selected %in% count_h2$trait]

# Create final list of independent variables
final_traits <- combined_selected

# Print results
cat("Selected independent traits (combined approach) (n =", length(combined_selected), "):\n")
print(combined_h2 %>% filter(trait %in% combined_selected) %>% arrange(desc(h2)))

cat("\nOf which biochemistry traits: (n =", length(biochem_selected), "):\n")
print(biochem_h2 %>% filter(trait %in% biochem_selected) %>% arrange(desc(h2)))

cat("\nOf which blood count traits: (n =", length(count_selected), "):\n")
print(count_h2 %>% filter(trait %in% count_selected) %>% arrange(desc(h2)))

# Step 7: Create final list of independent variables for GWAS
final_traits <- c(biochem_selected, count_selected)

# Visualize selected traits and their heritability
all_selected_h2 <- bind_rows(
  biochem_h2 %>% filter(trait %in% biochem_selected) %>% mutate(category = "Biochemistry"),
  count_h2 %>% filter(trait %in% count_selected) %>% mutate(category = "Blood Count")
)

# Plot heritability of selected traits
pdf('test.pdf')
ggplot(all_selected_h2, aes(x = reorder(trait, h2), y = h2, fill = category)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = h2 - h2_se, ymax = h2 + h2_se), width = 0.2) +
  coord_flip() +
  labs(title = "Heritability of Selected Independent Traits",
       x = "Trait", y = "Heritability (h²)") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()
# Save the final list of traits for GWAS
write.csv(all_selected_h2, file = "selected_traits_for_gwas.csv", row.names = FALSE)