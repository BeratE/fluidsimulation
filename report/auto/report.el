(TeX-add-style-hook
 "report"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt" "letterpaper" "twocolumn")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "english") ("inputenc" "utf8") ("geometry" "margin=1.0in")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "babel"
    "inputenc"
    "geometry"
    "titlesec"
    "graphicx"
    "caption"
    "amsmath"
    "amssymb"
    "siunitx"
    "epigraph"
    "fancyhdr")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:methods"
    "eq:weight"
    "eq:pressure"
    "eg:pbfconstraint"
    "eg:deltax"
    "fig:mesh"
    "eq:localtriangle"
    "fig:trianglesampling"
    "fig:systemuml"
    "sec:results"
    "sec:conclusion"
    "sec:future")
   (LaTeX-add-bibliographies
    "references.bib"))
 :latex)

