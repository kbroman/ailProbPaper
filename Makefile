pdf: ailprob.pdf Revised/response.pdf

all: pdf

ailprob.pdf: ailprob.tex ailprob.bib genetics.bst Figs/fig1.eps
	pdflatex ailprob
	bibtex ailprob
	pdflatex ailprob
	pdflatex ailprob
	pdflatex ailprob

Figs/fig1.pdf: R/R_fig.R 
	cd R;R CMD BATCH '--args eps=FALSE' R_fig.R

Figs/fig1.eps: R/R_fig.R 
	cd R;R CMD BATCH '--args eps=TRUE' R_fig.R

clean: 
	\rm -f *.aux *.bbl *.blg *.log *.bak *~ *.Rout */*~ */*.Rout */*.aux */*.log

cleanall: 
	\rm -f *.aux *.bbl *.blg *.log *.pdf *.bak *~ *.Rout */*~ */*.Rout */*.pdf */*.aux */*.log
