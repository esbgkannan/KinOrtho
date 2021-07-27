import sys

#I/O files
file_input = sys.argv[1]
file_output = file_input[0:file_input.rindex('.')] + '.rels.raw'

#Define column
col_p1 = 1
col_p2 = 4
col_relation = 10

#Unique protein ID pairs
unique_p1_p2 = {}

#Read KinOrtho output file
file = open(file_input, 'r')
lines = file.readlines()
file.close()
for i in range(1, len(lines)):
    lines[i] = lines[i].replace('\n', '')
    toks = lines[i].split('\t')
    if toks[col_relation] == 'Ortholog':
        #UniProt format
        token_p1 = toks[col_p1].split('|')
        p1 = token_p1[1]
        token_p2 = toks[col_p2].split('|')
        p2 = token_p2[1]
        p1_p2 = p1 + '\t' + p2
        p2_p1 = p2 + '\t' + p1
        if p2_p1 not in unique_p1_p2:
            unique_p1_p2[p1_p2] = ''

out = open(file_output, 'w')
for p1_p2 in unique_p1_p2:
    out.write(p1_p2 + '\n')
out.close()
