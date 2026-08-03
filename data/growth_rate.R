library(readxl)
library(dplyr)
library(here)

# --- Adjust these paths to your actual file locations ---
birdlife_file <- here("captain_test/captain2-main/generation_lengths.xlsx") 

# --- Read the Species sheet ---
# The main species table in BIRDBASE is usually named "Species"
sp <- read_excel(birdlife_file, sheet = "Generation lengths_v3.1_for_web")

# Inspect available columns first (run once)
# names(sp)

# Your 18 species (scientific names)
species <- c(
  "Apalis sharpii",
  "Aquila rapax",
  "Balearica pavonina",
  "Ceratogymna atrata",
  "Ceratogymna elata",
  "Chelictinia riocourii",
  "Circaetus beaudouini",
  "Criniger olivaceus",
  "Hylopsar cupreocauda",
  "Illadopsis rufescens",
  "Malaconotus lagdeni",
  "Merops mentalis",
  "Necrosyrtes monachus",
  "Parmoptila rubrifrons",
  "Phyllanthus atripennis",
  "Scotopelia ussheri",
  "Stephanoaetus coronatus",
  "Terathopius ecaudatus"
)

# --- Match species ---
# BIRDBASE uses a "Scientific name" column (adjust name if different)
# Common column names: "Scientific name", "Binomial", "sci_name"
result <- sp %>%
  filter(`Species name 2024` %in% species) %>%
  select(
    species = `Species name 2024`,
    generation_length = `Generation length`
  )

# --- Adjust these paths to your actual file locations ---
bird2020_file <- here("captain_test/captain2-main/bird2020.xlsx") 

# --- Read the Species sheet ---
# The main species table in BIRDBASE is usually named "Species"
bp <- read_excel(bird2020_file, sheet = "Table S2")

# Inspect available columns first (run once)
names(bp)

# Your 18 species (scientific names)
species <- c(
  "Apalis sharpii",
  "Aquila rapax",
  "Balearica pavonina",
  "Ceratogymna atrata",
  "Ceratogymna elata",
  "Chelictinia riocourii",
  "Circaetus beaudouini",
  "Criniger olivaceus",
  "Hylopsar cupreocauda",
  "Illadopsis rufescens",
  "Malaconotus lagdeni",
  "Merops mentalis",
  "Necrosyrtes monachus",
  "Parmoptila rubrifrons",
  "Phyllanthus atripennis",
  "Scotopelia ussheri",
  "Stephanoaetus coronatus",
  "Terathopius ecaudatus"
)

# --- Match species ---
# BIRDBASE uses a "Scientific name" column (adjust name if different)
# Common column names: "Scientific name", "Binomial", "sci_name"
resultbp <- bp %>%
  filter(`Scientific name` %in% species) %>%
  mutate(
    adult_survival = case_when(
      `Adult_survival_p` == 0 | is.na(`Adult_survival_p`) ~ `Adult_survival_m`,
      TRUE ~ `Adult_survival_p`
    )
  ) %>%
  select(
    species = `Scientific name`,
    adult_survival
  )
library(dplyr)

# Join the datasets
combined <- result %>%
  inner_join(resultbp, by = "species") %>%
  mutate(
    Gs = (1 / adult_survival)^(1 / generation_length)
  )

# Check results
print(combined)
summary(combined$Gs)

# Saves directly to your existing Documents folder
target_path <- file.path("~/captain_test/captain2-main/captaincode/data_ugf", "species_growth_rate.csv")
write.csv(combined, file = target_path, row.names = FALSE)
