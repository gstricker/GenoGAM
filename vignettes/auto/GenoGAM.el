(TeX-add-style-hook
 "GenoGAM"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "graphicx"
    "color"
    "framed"
    "alltt"
    "amsmath"
    "upquote")
   (TeX-add-symbols
    '("hlkwd" 1)
    '("hlkwc" 1)
    '("hlkwb" 1)
    '("hlkwa" 1)
    '("hlstd" 1)
    '("hlopt" 1)
    '("hlcom" 1)
    '("hlstr" 1)
    '("hlnum" 1)
    "genogam"
    "chipseq"
    "maxwidth"
    "hlipl"
    "FrameCommand")
   (LaTeX-add-labels
    "eq:mutantmodel")
   (LaTeX-add-environments
    "kframe"
    "knitrout")
   (LaTeX-add-bibliographies
    "bibliog"))
 :latex)

