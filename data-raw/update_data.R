## code to prepare `words` dataset goes here
#
# Whenever words.txt or abbrev.txt is updated, run the code below.

update_word_data <- function() {
   message("Loading words.txt")
   words <- read.table(
      # "inst/extdata/words.txt",
      file=system.file(package="multienrichjam", "extdata", "words.txt"),
      sep="\t",
      col.names=c("from", "to"),
      quote="")
   message("Loading abbrev.txt")
   abbrev <- read.table(
      # "inst/extdata/abbrev.txt",
      file=system.file(package="multienrichjam", "extdata", "abbrev.txt"),
      sep="\t",
      col.names=c("from", "to"),
      quote="")
   
   suppressMessages(
      usethis::use_data(words, abbrev,
         # internal=TRUE,
         overwrite=TRUE))
}
update_word_data()
## For development work, refresh the package
# devtools::load_all()
