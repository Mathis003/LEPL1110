import numpy as np
import os

path = "homework.c"

# Commenter les 4 lignes suivantes et decommenter la 5eme pour ne pas decompter les commentaires du fichier c
filename = "./homework_cleaned_comments.c"
if os.path.isfile(filename):
  	os.remove(filename)
os.system("gcc -fpreprocessed -dD -E " + path + " >> ./homework_cleaned_comments.c")
# filename = "homework.c"

lines = []
keywords = ["if", "else", "break" ,"int", "float", "char", "double", "long","for", "while", "do","void","goto","auto", "signed", "const", "extern", "register", "unsigned", "return","continue","enum","sizeof","struct", "typedef","union", "volatile","malloc","free"] 
with open(filename,"r") as file1:
	lines = file1.readlines()
lines_clean = [l.replace(" ","") for l in lines]
nbr_lines = len(lines) - lines_clean.count('\n')
a = ''.join(lines)
b = [a.count(c) for c in keywords]
nbr_pv = a.count(";")
nbr_acc = a.count("{")+a.count("}")
nbr_kw = np.sum(b)

print("Nombre de lignes :", nbr_lines, "\nNombre de points-virgules :", nbr_pv, "\nNombre d'accolades :", nbr_acc, "\nNombre de mots clés :", nbr_kw)

print("Index de longueur final : ", int(np.ceil((50 * nbr_lines + 750 * nbr_pv + 50 * nbr_acc + 150 * nbr_kw) / (nbr_lines + nbr_pv + nbr_acc + nbr_kw))))