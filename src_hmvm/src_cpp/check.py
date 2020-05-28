import sys

# arg1 = dir1
# arg2 = dir2
# arg3 = filename
args = sys.argv

fname1 = args[1] + args[3]
fname2 = args[2] + args[3]
print("diff: " + fname1 + " " + fname2)

f1 = open(fname1)
l1 = f1.readlines()
f2 = open(fname2)
l2 = f2.readlines()
f1.close()
f2.close()

#print(l1)

len1 = len(l1)
len2 = len(l2)

#print(len1)
#print(len2)

nf = 0
for i in range(len(l1)):
    str1 = l1[i]
    num1 = float(str1[0:5])
    str2 = l2[i]
    num2 = float(str1[0:5])
    # print(num1,num2)
    ret = (num1 == num2)
    if(ret==False):
        print("False: " + i + " " + num1 + " " + num2)
        nf = nf + 1
if(nf==0):
    print("no error")
else:
    print(str(nf) + " error(s) found")


