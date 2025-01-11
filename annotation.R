##-----load a data------------------------------##


ANNOTATION= read.csv('/home/gabat/scratch/GENOM/annotations.tsv',sep='\t',header = TRUE,comment.char="*")  ##colonnes 3 : product_name represent les genes de GENOM_S




save.image(file = "/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/ANNOTATION.RData")

