with open('/media/daniel/Elements/Data/BindingDB_All_2021m7.tsv/BindingDB_All.tsv','r') as firstfile, open('BindingDB_ALL_header.txt','a') as secondfile:
     
    for line in firstfile:
             secondfile.write(line)
             break