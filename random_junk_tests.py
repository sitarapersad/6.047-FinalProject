f = open('scz_bip.log')
for line in f:
    print line
    if line =='Summary of Genetic Correlation Results\n':
        print 'HOLLA'

f.close()