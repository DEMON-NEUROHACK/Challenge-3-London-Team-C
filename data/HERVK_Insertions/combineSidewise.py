import sys

file1 = open (sys.argv[1], "r")#novel
#print ("novel file is " + sys.argv[1])
file2 = open (sys.argv[2], "r")#counts
#print ("counts file is "	+ sys.argv[2])
lines1=file1.readlines()
file1.close()

lines2=file2.readlines()
file2.close()

#if not (len(lines1) == len(lines2)) or not (len(lines1)+1 == len(lines2)) or not (len(lines1)+1 == len(lines2)) :
if abs(len(lines2)-len(lines1))>1:
    print ("different niumber of lines, possible error\nsys.argv[1]\tsys.argv[2]")
    print (len(lines1))
    print (len(lines2))
if (len(lines1)+1 == len(lines2)):
    lines1.insert(0,"Chr\tStart\tEnd\t")
    
for i in range(0,len(lines1)):
    line1=lines1[i]
    #line1 = line1.strip()
    bits=line1.split()
    part1=bits[0] + "_" + bits[1] + "_" + bits[2]
    line2=lines2[i]
    line2 = line2.strip()
    line3=part1 + "\t" + line2
    print (line3)
