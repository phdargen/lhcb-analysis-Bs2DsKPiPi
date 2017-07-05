import sys
jobNumber = int(sys.argv[1])
length = len(jobs(jobNumber).subjobs)
file = open('LFNs.txt','w')
for i in range(0,length):
    output = jobs(jobNumber).subjobs(i).backend.getOutputDataLFNs()
    if( output.hasLFNs() == True):
        s = str(jobs(jobNumber).subjobs(i).backend.getOutputDataLFNs().getReplicas().keys())
        print s
        file.write(s)
        file.write('\n')
file.close()


