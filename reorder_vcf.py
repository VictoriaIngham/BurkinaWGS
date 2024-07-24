import os, sys
Argument = []
Argument = sys.argv[1:]

if (len(Argument)) < 1:
        print("Usage:Input_vcf Outputfile")
        sys.exit()

output = open(Argument[1],"w")
input = open(Argument[0])

def numeric_compare(x, y):
        x1 = int(x)
        y1 = int(y)
        return x1 - y1
Chromosome = ["2R","3R","2L","UNKN","3L","X","Y_unplaced","Mt"]
VCF = {}
for line in input:
        if line.startswith("#"):
                output.write(str(line))
                continue
        v = []
        v = line.strip("\n").split("\t")

        if v[0] not in VCF:
                VCF[v[0]] = {}
                VCF[v[0]][v[1]] = line
        else:
                VCF[v[0]][v[1]] = line
for chr in Chromosome:
        for pos in sorted(VCF[chr].keys(),cmp=numeric_compare):
                output.write(str(VCF[chr][pos]))
                output.flush()
output.close()