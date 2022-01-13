import os
import sys



fileList= open("list.txt", "r")
files=fileList.readlines()
fileList.close()
wholeList = ""
subjectIDs = ""
lastID=""


for file in files:
    vcf = file.rstrip() + ".vcf"
    bed = file.rstrip() + ".bed"
    os.system("python filter.py %s %s" %(vcf,bed))
    

for file in files:
    subjectIDs = subjectIDs + file.rstrip() + "\t"
    lastID=file.rstrip()
    fullFile = file.rstrip() + ".bed"
    FileTempSorted = file.rstrip() + ".Hits_temp_sorted.bed"
    FileTempMerged = file.rstrip() + ".Hits_temp_sorted_merged.bed"
    os.system("bedtools sort -i %s   > %s" %(fullFile,FileTempSorted))
    os.system("bedtools merge -d 500 -i %s > %s" %(FileTempSorted,FileTempMerged))

    wholeList = wholeList + " " + FileTempMerged

#print (wholeList)
os.system("rm mergedST.bed")
os.system("rm mergedST.sorted.bed")

os.system("cat %s >> %s" %(wholeList,"mergedST.bed"))
os.system("bedtools sort -i mergedST.bed   > mergedST.sorted.bed")
os.system("bedtools merge -d 500 -i mergedST.sorted.bed > mergedAllST.sorted.bed")

os.system("bedtools window -w 0 -c -a mergedAllST.sorted.bed -b mergedST.sorted.bed > countsHERVs.bed")


####for every file do this:
#### bedtools window -w 0 -c -a mergedAll.sorted.bed -b LP6008462-DNA_F06.ProcessednovelHitsF.bed > del.txt
for file in files:
      fullFile = file.rstrip() + ".Hits_temp_sorted_merged.bed"
      file2 = file.rstrip() + ".Hits_counts.bed"
      os.system("bedtools window -w 0 -c -a mergedAllST.sorted.bed -b %s > %s" %(fullFile,file2))


##os.system("rm *novelHits_temp_sorted.bed")
##os.system("rm *novelHits_temp_sorted_merged.bed")

#NOW USE THE SAME CODE AS FOR KNOWN, BECAUSE NOW THERE IS A FILE WITH COUNTS FOR EACH SUBJECT
novel = file2
wholeList = ""
for file in files:
     fileCounts = file.rstrip() + ".Hits_counts.bed"
     wholeList = wholeList + " " + fileCounts

#print (wholeList)
#get the column 4 - which is 0 if that insertion has not been hit and 1 if it has - from each one of the files in
#the list
os.system("rm merged01nov.txt")
os.system("rm countsnov.txt")
matrix = open("matrix.txt","w")
matrix.write (subjectIDs + "\n")
matrix.close()
os.system("rm matrixIDs.txt")
os.system("python extractColumn.py 4 %s >> %s" %(wholeList,"matrix.txt"))
#print(novel)
os.system("python combineSidewise.py %s matrix.txt >> matrixIDs.txt" %(novel))

os.system("rm *Hits_temp_sorted*")
os.system("rm *Hits_counts.bed")


with open('matrixIDs.txt') as file:
    lis = [x.replace('\n', '').split('\t') for x in file]


for x in zip(*lis):
    for y in x:       
        print(y + '\t', end="")
    print('')
