




#usethis::use_travis()
devtools::install()
devtools::document()

usethis::use_build_ignore("inst/script_development.R")
usethis::use_build_ignore(".travis.yml")
# usethis::use_cran_comments()

devtools::check()
devtools::test()

# devtools::check_win_devel()
rhub::check_for_cran()


devtools::build()
devtools::build_manual()
devtools::release_checks()
devtools::spell_check()


git add .
git commit -m "cran submission"
git push --set-upstream origin master


git remote add origin https://github.com/fndemarqui/peppm.git