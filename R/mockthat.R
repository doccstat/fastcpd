utils_menu <- function() {
  utils::menu(
    c("yes", "no"),
    title = "ggplot2 is not installed. Do you want to install it?"
  )
}

require_namespace <- function(package_name) {
  requireNamespace(package_name, quietly = TRUE)
}

install_packages <- function(package_name) {
  utils::install.packages(
    package_name, repos = "https://cloud.r-project.org", quiet = TRUE
  )
}
