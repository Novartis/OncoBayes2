args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("Usage: build-test-fixture.R <fixture-recipe> <output-rds>", call. = FALSE)
}

recipe <- args[[1L]]
output <- args[[2L]]

if (!file.exists(recipe)) {
  stop("Fixture recipe does not exist: ", recipe, call. = FALSE)
}


devtools::load_all(quiet = TRUE)

helper_files <- Sys.glob(file.path("tests", "testthat", "helper-*.R"))
for (helper_file in helper_files) {
  source(helper_file, local = FALSE)
}

recipe_env <- new.env(parent = globalenv())
source(recipe, local = recipe_env)

if (!exists("fixture", envir = recipe_env, inherits = FALSE)) {
  stop("Fixture recipe must create an object named `fixture`: ", recipe, call. = FALSE)
}

fixture <- get("fixture", envir = recipe_env, inherits = FALSE)
output_dir <- dirname(output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tmp <- tempfile(pattern = paste0(basename(output), "-"), tmpdir = output_dir)
on.exit(unlink(tmp), add = TRUE)

saveRDS(fixture, tmp)
if (!file.rename(tmp, output)) {
  stop("Could not move temporary fixture into place: ", output, call. = FALSE)
}

message("Wrote fixture: ", output)
