(TeX-add-style-hook
 "GenoGAM"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "amsmath")
   (TeX-add-symbols
    "genogam"
    "chipseq")
   (LaTeX-add-labels
    "eq:mutantmodel")
   (LaTeX-add-bibliographies
    "bibliog"))
 :latex)

