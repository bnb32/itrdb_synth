#!/bin/bash

if [ $# -ne 1 ]; then echo usage $0: file_name; exit 1; fi

FILE_NAME=$1

pdflatex ${FILE_NAME}
#latex ${FILE_NAME}
#latex ${FILE_NAME}
bibtex ${FILE_NAME}
#latex ${FILE_NAME}
#latex ${FILE_NAME}
#dvipdf ${FILE_NAME}.dvi
pdflatex ${FILE_NAME}
pdflatex ${FILE_NAME}
#${HOME}/gdrive upload ${FILE_NAME}.pdf
