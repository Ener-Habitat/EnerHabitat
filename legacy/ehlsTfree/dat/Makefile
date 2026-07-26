#
# Makefile
#

%.gnp:	phony

%.pdf:	%.gnp
	sed -e 's/^# set terminal.*$$/set terminal epsl   sol lw 4/' $(basename $<).gnp |\
	sed -e 's/^# set output.*$$/set output "'$(basename $<).tex'"/' |\
	gnuplot
	echo '\documentclass{article}' > fig.tex
	echo '\usepackage{graphicx}' >> fig.tex
	echo '\usepackage{amsmath}' >> fig.tex
	echo '\usepackage{amssymb}' >> fig.tex
	echo '\usepackage{latexsym}' >> fig.tex
	echo '' >> fig.tex
	echo '\pagestyle{empty}' >> fig.tex
	echo '' >> fig.tex
	echo '\\begin{document}' >> fig.tex
	echo '\input{'$(basename $<)'}' >> fig.tex
	echo '' >> fig.tex
	echo '\\end{document}' >> fig.tex
	latex fig
	dvips -E fig -o fig.eps
	epstopdf fig.eps -o $@

