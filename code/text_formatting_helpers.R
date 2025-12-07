# Text formatting helper functions for result output
#
# These functions are used to format numbers and percentages for
# LaTeX tables and result snippet files.

# Format numbers with commas as thousand separators
format_num <- function(x) formatC(x, format = "d", big.mark = ",")

# Format percentages for LaTeX (with escaped percent sign)
format_pct <- function(x) sprintf("%.1f\\%%", x * 100)

# Format decimal numbers to specified digits
format_decimal <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), round(x, digits))

# Write a single value to a text file
write_result_snippet <- function(value, filename, path = RESULT_SNIPPETS_PATH) {
    writeLines(as.character(value), file.path(path, filename))
}