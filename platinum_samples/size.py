import sys
#import numpy


hej=[]
for line in open(sys.argv[1]):	
	hej.append(int(line.strip()))
	

print(sum(hej)/len(hej))
print(min(hej))
print(max(hej))

