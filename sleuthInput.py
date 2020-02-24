
def SleuthInput():
    SRRs = 'SRR5660030','SRR5660033','SRR5660044','SRR5660045'
    covFile = open('cov.txt','w')
    condition1 = "2dpi"
    condition2 = "6dpi"
    covFile.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    for i in SRRs:
        path = '/data/kdelany/compBio_miniProject/'+i
        if int(i[3:])%2==0:
            covFile.write(str(i)+ '\t' + condition1 + '\t'+ str(path)+ '\n')
        else:
            covFile.write(str(i)+ '\t' + condition2 + '\t'+ str(path)+ '\n')
    covFile.close()
