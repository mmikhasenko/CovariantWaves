#!/bin/bash

#FORMULA="T_{\text{test}}"
#if (($#)) ; then
#     FORMULA=$1
#else
#     FORMULA=$(cat)
#fi

#It works nicely
FORMULA=${*:-`cat`}

TEXT="
\documentclass[preview]{standalone}
\usepackage{amsmath}
\begin{document}
\$
$FORMULA
\$
\end{document}
"
echo $TEXT | pdflatex &>/dev/null && rm standalone.log standalone.aux
#%Streamed expression:

