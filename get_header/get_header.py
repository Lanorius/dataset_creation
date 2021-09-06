with open('/media/lanorius/Elements/data/BindingDB_All.tsv','r') as firstfile, open('BindingDB_ALL_header.txt','a') as secondfile:
	i = 0
	for line in firstfile:
             secondfile.write(line)
             if i == 5:
             	break
             i+=1
