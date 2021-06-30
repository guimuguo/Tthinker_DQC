import sys

fname = sys.argv[1]
fwn = fname + "_null"
fw = open(fwn, "w")
file1 = open(fname, 'r') 
Lines = file1.readlines()
count = 0
for line in Lines:
    cur = int(line.split()[0])
    diff = cur - count
    for i in range(diff):
        fw.write("\n")
        count += 1
    fw.write(line)
    count += 1
file1.close()
fw.close()
