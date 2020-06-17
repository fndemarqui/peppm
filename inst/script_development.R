

usethis::use_build_ignore("inst/script_development.R")


usethis::use_travis()
devtools::install()
devtools::document()
devtools::check()
devtools::test()
devtools::check_win_devel()
rhub::check_for_cran()


devtools::build()
devtools::build_manual()
devtools::release_checks()
devtools::spell_check()
# devtools::test()
# devtools::reload()
# devtools::run_examples()


git add .
git commit

