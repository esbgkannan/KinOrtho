import math
import os
import sys

def avg_evalue(ev_1, ev_2):
    if ev_1 == 0:
        ev_1 = float(min_e)
    if ev_2 == 0:
        ev_2 = float(min_e)
    avg = pow(10, ((math.log10(ev_1)+math.log10(ev_2))/2))
    return avg

def build_blast_db(seq_file, db_name):
    command = 'makeblastdb -in ' + seq_file
    command += ' -dbtype prot -out ' + db_name
    os.system(command)

def build_graph(ortholog_file, out_file, is_full_pipeline):
    sys.stdout.write('\tReading ortholog file')
    species_pair_weight_sum = {}
    species_pair_weight_count = {}
    species_weight_sum = {}
    species_weight_count = {}
    file = open(ortholog_file, 'r')
    lines = file.readlines()
    file.close()
    for i in range(0, len(lines)):
        toks = lines[i].split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        weight = -math.log10(float(toks[2]))
        species_1 = ''
        species_2 = ''
        if is_full_pipeline:
            species_1 = protein_to_species[protein_1]
            species_2 = protein_to_species[protein_2]
        else:
            species_1 = protein_to_species[protein_1[:protein_1.rindex('|')]]
            species_2 = protein_to_species[protein_2[:protein_2.rindex('|')]]
        species_pair = species_1 + '\t' + species_2
        if species_1 > species_2:
            species_pair = species_2 + '\t' + species_1
        
        # Ortholog and Co-ortholog
        if species_1 != species_2:
            if species_pair not in species_pair_weight_sum:
                species_pair_weight_sum[species_pair] = weight
                species_pair_weight_count[species_pair] = 1
            else:
                species_pair_weight_sum[species_pair] += weight
                species_pair_weight_count[species_pair] += 1
        # In-paralog
        else:
            if species_1 not in species_weight_sum:
                species_weight_sum[species_1] = weight
                species_weight_count[species_1] = 1
            else:
                species_weight_sum[species_1] += weight
                species_weight_count[species_1] += 1

        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    
    # Normalization
    sys.stdout.write('\tNormalizing weight')
    out = open(out_file, 'w')
    for i in range(0, len(lines)):
        toks = lines[i].split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        weight = -math.log10(float(toks[2]))
        species_1 = ''
        species_2 = ''
        if is_full_pipeline:
            species_1 = protein_to_species[protein_1]
            species_2 = protein_to_species[protein_2]
        else:
            species_1 = protein_to_species[protein_1[:protein_1.rindex('|')]]
            species_2 = protein_to_species[protein_2[:protein_2.rindex('|')]]
        species_pair = species_1 + '\t' + species_2
        if species_1 > species_2:
            species_pair = species_2 + '\t' + species_1
        
        # Ortholog and Co-ortholog
        norm_weight = weight
        if species_1 != species_2:
            norm_weight /= (species_pair_weight_sum[species_pair]/species_pair_weight_count[species_pair])
        # In-paralog
        else:
            norm_weight /= (species_weight_sum[species_1]/species_weight_count[species_1])
        out.write(protein_1 + '\t' + protein_2 + '\t' + str(norm_weight) + '\n')
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    out.close()

def clustering(graph_file, out_name):
    command = 'mcxload -abc ' + graph_file
    command += ' -o ' + out_name + '.mci'
    command += ' -write-tab ' + out_name + '.tab'
    os.system(command)
    command = 'mcl ' + out_name + '.mci'
    command += ' -I ' + inflat + ' -te ' + thread
    command += ' -o ' + out_name + '.out'
    os.system(command)

def combin_result(full_result, domain_result, out_file):
    
    # Full-length result
    sys.stdout.write('\tReading full-length result')
    protein_pair_to_info = {}
    file = open(full_result, 'r')
    lines = file.readlines()
    file.close()
    for i in range(1, len(lines)):
        toks = lines[i].rstrip().split('\t')
        protein_1 = toks[1]
        protein_2 = toks[3]
        ev = toks[4]
        weight = toks[5]
        relationship = toks[6]
        info = ev + '\t' + weight + '\t' + relationship
        protein_pair_to_info[protein_1 + '\t' + protein_2] = info
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    
    # Domain-based result
    sys.stdout.write('\tReading domain-based result')
    out = open(out_file, 'w')
    header = 'Species_1' + '\t' + 'Protein_1' + '\t' + 'Domain_1' + '\t'
    header += 'Species_2' + '\t' + 'Protein_2' + '\t' + 'Domain_2' + '\t'
    header += 'E-value_Full' + '\t' + 'E-value_Domain' + '\t'
    header += 'Weight_Full' + '\t' + 'Weight_Domain' + '\t' + 'Relationship'
    out.write(header + '\n')
    file = open(domain_result, 'r')
    lines = file.readlines()
    file.close()
    for i in range(1, len(lines)):
        toks = lines[i].rstrip().split('\t')
        species_1 = toks[0]
        protein_1 = toks[1]
        domain_1 = toks[2]
        species_2 = toks[3]
        protein_2 = toks[4]
        domain_2 = toks[5]
        ev = toks[6]
        weight = toks[7]
        relationship = toks[8]
        protein1_protein2 = protein_1 + '\t' + protein_2
        protein2_protein1 = protein_2 + '\t' + protein_1
        
        #Overlapping
        if protein1_protein2 in protein_pair_to_info:
            info = protein_pair_to_info[protein1_protein2].split('\t')
            ev_full = info[0]
            weight_full = info[1]
            relationship_full = info[2]
            if relationship_full == relationship:
                print_line = species_1 + '\t' + protein_1 + '\t' + domain_1 + '\t'
                print_line += species_2 + '\t' + protein_2 + '\t' + domain_2 + '\t'
                print_line += ev_full + '\t' + ev + '\t'
                print_line += weight_full + '\t' + weight + '\t' + relationship
                out.write(print_line + '\n')
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    out.close()

def filtering(homolog_file, ortholog_file, graph_file, cluster_name, out_file, is_full_pipeline):
    
    # Query sequences (closest)
    query_seq = {}
    query_seq_best_hit = {}
    if homolog_file != '':
        sys.stdout.write('\tReading homology search result')
        file = open(homolog_file, 'r')
        lines = file.readlines()
        file.close()
        for i in range(0, len(lines)):
            toks = lines[i].split('\t')
            if toks[0] not in query_seq:
                query_seq[toks[0]] = ''
                query_seq_best_hit[toks[1]] = ''
            
            if i%(int(len(lines)/10)) == 0:
                sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
                sys.stdout.flush()
        sys.stdout.write('\n')
    
    # Index-protein mapping
    index_to_protein = {}
    protein_to_index = {}
    file = open(cluster_name + '.tab', 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        toks = line.rstrip().split('\t')
        index_to_protein[toks[0]] = toks[1]
        protein_to_index[toks[1]] = toks[0]
    
    # Index-group mapping
    index_to_group = {}
    group_contain_query_seq = {}
    file = open(cluster_name + '.out', 'r')
    lines = file.readlines()
    file.close()
    begin = 0
    now_group = ''
    for line in lines:
        toks = line.rstrip().split()
        if toks[0] == 'begin':
            begin = 1
        elif begin:
            first_index = 0
            if line[0:1] != ' ':
                now_group = toks[0]
                first_index = 1
            for i in range(first_index, len(toks)):
                if toks[i] != '$':
                    index_to_group[toks[i]] = now_group
                    protein = index_to_protein[toks[i]]
                    if protein in query_seq_best_hit:
                        group_contain_query_seq[now_group] = ''
                    elif protein[:protein.rindex('|')] in query_seq_best_hit:
                        group_contain_query_seq[now_group] = ''
    
    # Weight between protein 1 and protein 2
    sys.stdout.write('\tReading graph file')
    protein_pair_weight = {}
    file = open(graph_file, 'r')
    lines = file.readlines()
    file.close()
    for i in range(0, len(lines)):
        toks = lines[i].rstrip().split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        weight = toks[2]
        protein_pair_weight[protein_1 + '\t' + protein_2] = weight
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    
    # Print results
    sys.stdout.write('\tPrinting result')
    out = open(out_file, 'w')
    if is_full_pipeline:
        header = 'Species_1' + '\t' + 'Protein_1' + '\t' + 'Species_2' + '\t' + 'Protein_2'
        header += '\t' + 'E-value' + '\t' + 'Weight' + '\t' + 'Relationship'
        out.write(header + '\n')
    else:
        header = 'Species_1' + '\t' + 'Protein_1' + '\t' + 'Domain_1' + '\t'
        header += 'Species_2' + '\t' + 'Protein_2' + '\t' + 'Domain_2' + '\t'
        header += 'E-value' + '\t' + 'Weight' + '\t' + 'Relationship'
        out.write(header + '\n')
    file = open(ortholog_file, 'r')
    lines = file.readlines()
    file.close()
    for i in range(0, len(lines)):
        toks = lines[i].rstrip().split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        domain_1 = protein_1[protein_1.rindex('|')+1:]
        domain_2 = protein_2[protein_2.rindex('|')+1:]
        protein1_protein2 = protein_1 + '\t' + protein_2
        protein2_protein1 = protein_2 + '\t' + protein_1
        weight = ''
        if protein1_protein2 in protein_pair_weight:
            weight = protein_pair_weight[protein1_protein2]
        else:
            weight = protein_pair_weight[protein2_protein1]
        ev = toks[2]
        relationship = toks[3]
        index_1 = protein_to_index[protein_1]
        index_2 = protein_to_index[protein_2]
        group_1 = index_to_group[index_1]
        group_2 = index_to_group[index_2]
        if not is_full_pipeline:
            protein_1 = protein_1[:protein_1.rindex('|')]
            protein_2 = protein_2[:protein_2.rindex('|')]
        species_1 = protein_to_species[protein_1]
        species_2 = protein_to_species[protein_2]
        
        # Filtering
        if group_1 == group_2 and (group_1 in group_contain_query_seq or homolog_file == ''):
            if is_full_pipeline:
                print_line = species_1 + '\t' + protein_1 + '\t' + species_2 + '\t' + protein_2
                print_line += '\t' + ev + '\t' + weight + '\t' + relationship
                out.write(print_line + '\n')
            else:
                print_line = species_1 + '\t' + protein_1 + '\t' + domain_1 + '\t'
                print_line += species_2 + '\t' + protein_2 + '\t' + domain_2 + '\t'
                print_line += ev + '\t' + weight + '\t' + relationship
                out.write(print_line + '\n')
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    out.close()

def homology_search(seq_file, db_name, out_file):
    command = 'blastp -outfmt "6 qseqid sseqid sstart send evalue bitscore positive" -max_hsps 1'
    command += ' -query ' + seq_file
    command += ' -db ' + db_name
    command += ' -num_threads ' + thread
    command += ' -evalue ' + e_value
    command += ' > ' + out_file
    os.system(command)

def orthology_inference(allvsall_file, out_file, is_full_pipeline):
    
    # Ortholog
    sys.stdout.write('\tInferring ortholog')
    best_protein_species = {}
    best_protein_species_ev = {}
    is_ortholog = {}
    has_orthologs = {}
    file = open(allvsall_file, 'r')
    lines = file.readlines()
    file.close()
    out = open(out_file, 'w')
    for i in range(0, len(lines)):
        toks = lines[i].split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        ev = float(toks[4])
        species_1 = ''
        species_2 = ''
        if is_full_pipeline:
            species_1 = protein_to_species[protein_1]
            species_2 = protein_to_species[protein_2]
        else:
            species_1 = protein_to_species[protein_1[:protein_1.rindex('|')]]
            species_2 = protein_to_species[protein_2[:protein_2.rindex('|')]]
        protein1_species2 = protein_1 + '\t' + species_2
        protein2_species1 = protein_2 + '\t' + species_1
        
        # Identify bidirectional best hit
        if species_1 != species_2:
            if protein1_species2 not in best_protein_species:
                best_protein_species[protein1_species2] = protein_2
                best_protein_species_ev[protein1_species2] = ev
                if protein2_species1 in best_protein_species:
                    if best_protein_species[protein2_species1] == protein_1:
                        is_ortholog[protein_1 + '\t' + protein_2] = ''
                        has_orthologs[protein_1] = ''
                        has_orthologs[protein_2] = ''
                        avg_ev = avg_evalue(ev, best_protein_species_ev[protein2_species1])
                        out.write(protein_1 + '\t' + protein_2 + '\t' + str(avg_ev) + '\t' + 'Ortholog' + '\n')
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    
    best_protein_species = {}
    best_protein_species_ev = {}
    
    # In-paralog
    sys.stdout.write('\tInferring in-paralog')
    in_other_species = {}
    candidate_inparalog = {}
    protein2inparalogs = {}
    for i in range(0, len(lines)):
        toks = lines[i].split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        ev = float(toks[4])
        species_1 = ''
        species_2 = ''
        if is_full_pipeline:
            species_1 = protein_to_species[protein_1]
            species_2 = protein_to_species[protein_2]
        else:
            species_1 = protein_to_species[protein_1[:protein_1.rindex('|')]]
            species_2 = protein_to_species[protein_2[:protein_2.rindex('|')]]
        protein1_protein2 = protein_1 + '\t' + protein_2
        protein2_protein1 = protein_2 + '\t' + protein_1
        
        # Identify bidirectional better hit
        if species_1 != species_2:
            in_other_species[protein_1] = ''
        elif protein_1 not in in_other_species:
            candidate_inparalog[protein1_protein2] = ev
            if protein2_protein1 in candidate_inparalog:
                if protein_1 in has_orthologs or protein_2 in has_orthologs:
                    is_inparalog = 0
                    if is_full_pipeline:
                        if protein_1 != protein_2:
                            is_inparalog = 1
                    else:
                        if protein_1[:protein_1.rindex('|')] != protein_2[:protein_2.rindex('|')]:
                            is_inparalog = 1
                    if is_inparalog:
                        avg_ev = avg_evalue(ev, candidate_inparalog[protein2_protein1])
                        out.write(protein_1 + '\t' + protein_2 + '\t' + str(avg_ev) + '\t' + 'In-paralog' + '\n')
                        if protein_1 not in protein2inparalogs:
                            protein2inparalogs[protein_1] = protein_2
                        else:
                            protein2inparalogs[protein_1] += "\t" + protein_2
                        if protein_2 not in protein2inparalogs:
                            protein2inparalogs[protein_2] = protein_1
                        else:
                            protein2inparalogs[protein_2] += "\t" + protein_1
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    
    in_other_species = {}
    candidate_inparalog = {}
    
    # Co-ortholog
    sys.stdout.write('\tInferring co-ortholog')
    candidate_coortholog = {}
    for i in range(0, len(lines)):
        toks = lines[i].split('\t')
        protein_1 = toks[0]
        protein_2 = toks[1]
        ev = float(toks[4])
        species_1 = ''
        species_2 = ''
        if is_full_pipeline:
            species_1 = protein_to_species[protein_1]
            species_2 = protein_to_species[protein_2]
        else:
            species_1 = protein_to_species[protein_1[:protein_1.rindex('|')]]
            species_2 = protein_to_species[protein_2[:protein_2.rindex('|')]]
        protein1_protein2 = protein_1 + '\t' + protein_2
        protein2_protein1 = protein_2 + '\t' + protein_1
    
        if species_1 != species_2:
            if protein2_protein1 not in candidate_coortholog:
                # Case 1: A1 -> In-paralog -> A2 -> Ortholog -> B2
                if protein_1 in protein2inparalogs:
                    inparalogs = protein2inparalogs[protein_1].split('\t')
                    for j in inparalogs:
                        inparalog_protein2 = j + '\t' + protein_2
                        protein2_inparalog = protein_2 + '\t' + j
                        if inparalog_protein2 in is_ortholog or protein2_inparalog in is_ortholog:
                            candidate_coortholog[protein1_protein2] = ev
                            break
                if protein1_protein2 not in candidate_coortholog and protein_2 in protein2inparalogs:
                    inparalogs = protein2inparalogs[protein_2].split('\t')
                    for j in inparalogs:
                        inparalog_protein1 = j + '\t' + protein_1
                        protein1_inparalog = protein_1 + '\t' + j
                        if inparalog_protein1 in is_ortholog or protein1_inparalog in is_ortholog:
                            candidate_coortholog[protein1_protein2] = ev
                            break
                # Case 2: A1 -> In-paralog -> A2 -> Ortholog -> B2 -> In-paralog -> B1
                if protein_1 in protein2inparalogs and protein_2 in protein2inparalogs:
                    inparalogs_1 = protein2inparalogs[protein_1].split('\t')
                    inparalogs_2 = protein2inparalogs[protein_2].split('\t')
                    for j in inparalogs_1:
                        for k in inparalogs_2:
                            inparalog1_inparalog2 = j + '\t' + k
                            inparalog2_inparalog1 = k + '\t' + j
                            if inparalog1_inparalog2 in is_ortholog or inparalog2_inparalog1 in is_ortholog:
                                candidate_coortholog[protein1_protein2] = ev
                                break
            else:
                if protein1_protein2 not in is_ortholog and protein2_protein1 not in is_ortholog:
                    avg_ev = avg_evalue(ev, candidate_coortholog[protein2_protein1])
                    out.write(protein_1 + '\t' + protein_2 + '\t' + str(avg_ev) + '\t' + 'Co-ortholog' + '\n')
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    out.close()

def print_manual():
    print('')
    print('Usage:')
    print('python kinortho.py [-i input] [-f full_seq] [-d domain_seq] [-o output] [-E e_value]')
    print('\t\t[-t thread] [-I inflat] [-e min_e] [-s start_step] [-S stop_step]')
    print('')
    print('Options:')
    print('-i: the folder of input proteomes (required)')
    print('-f: the fasta file of full-length query sequences')
    print('-d: the fasta file of domain-based query sequences')
    print('-o: output folder (defalut: ./kinortho_out/)')
    print('-E: E-value threshold (default: 1e-5)')
    print('-t: number of thread (default: 1)')
    print('-I: inflation value (default: 1.5)')
    print('-e: minimum E-value (default: 1e-200)')
    print('-s: start step (default: 0; value: {0,1,2,3,4,5,6})')
    print('-S: stop step (default: 6; value: {0,1,2,3,4,5,6})')
    print('')
    print('Steps:')
    print('0: Initiation')
    print('1: Homology Search')
    print('2: Building BLAST DB')
    print('3: All-vs-all Homology Search')
    print('4: Orthology Inference')
    print('5: Cluster Analysis')
    print('6: Combining Results')
    print('')

def rebuild_blast_db(proteomes_seq, blast_result, out_file, is_full_pipeline):
    
    # Get BLAST result
    sys.stdout.write('\tGetting BLAST result')
    homologs_id = {}
    homologs_id_domain_pos = {}
    file = open(blast_result, 'r')
    lines = file.readlines()
    file.close()
    for i in range(0, len(lines)):
        toks = lines[i].split('\t')
        if is_full_pipeline:
            homologs_id[toks[1]] = ''
        else:
            if toks[1] not in homologs_id_domain_pos:
                homologs_id_domain_pos[toks[1]] = toks[2] + '-' + toks[3]
            else:
                homologs_id_domain_pos[toks[1]] += ';' + toks[2] + '-' + toks[3]
        
        if i%(int(len(lines)/10)) == 0:
            sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
            sys.stdout.flush()
    sys.stdout.write('\n')
    
    # Generate homolog sequences
    sys.stdout.write('\tGenerating homolog sequences')
    out = open(out_file, 'w')
    file = open(proteomes_seq, 'r')
    lines = file.readlines()
    file.close()
    is_homolog = 0
    now_title = ''
    now_seq = ''
    if is_full_pipeline:
        for i in range(0, len(lines)):
            if lines[i][0:1] == '>':
                toks = lines[i][1:].split(' ')
                seq_id = toks[0]
                if seq_id in homologs_id:
                    is_homolog = 1
                    out.write(lines[i])
                else:
                    is_homolog = 0
            elif is_homolog:
                out.write(lines[i])
            
            if i%(int(len(lines)/10)) == 0:
                sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
                sys.stdout.flush()
        sys.stdout.write('\n')
    else:
        for i in range(0, len(lines)):
            if lines[i][0:1] == '>':
                toks = lines[i][1:].split(' ')
                seq_id = toks[0]
                if seq_id in homologs_id_domain_pos:
                    is_homolog = 1
                    now_title = lines[i][1:]
                    now_seq = ''
                else:
                    is_homolog = 0
                    now_title = ''
                    now_seq = ''
            elif is_homolog:
                now_seq += lines[i].rstrip()
                if i == len(lines)-1 or lines[i+1][0:1] == '>':
                    # Determine domain boundary
                    in_domain_pos = {}
                    toks = now_title.split(' ', 1)
                    seq_id = toks[0]
                    domain_pos = homologs_id_domain_pos[seq_id].split(';')
                    for j in domain_pos:
                        start, end = j.split('-')
                        for k in range(int(start), int(end)+1):
                            in_domain_pos[k] = ''
                    # Print domain sequences
                    domain_index = 1
                    pos = sorted(in_domain_pos)
                    now_start = pos[0]
                    now_end = pos[0]
                    for j in range(0, len(pos)):
                        now_end = pos[j]
                        title = '>' + seq_id
                        title += '|dom' + str(domain_index) + '_' + str(now_start) + '_' + str(now_end)
                        title += ' ' + toks[1]
                        if j != len(pos)-1:
                            if pos[j+1]-pos[j] > 1:
                                out.write(title)
                                out.write(now_seq[now_start-1:now_end] + '\n')
                                domain_index += 1
                                now_start = pos[j+1]
                                now_end = pos[j+1]
                        else:
                            out.write(title)
                            out.write(now_seq[now_start-1:now_end] + '\n')
            
            if i%(int(len(lines)/10)) == 0:
                sys.stdout.write('...' + str((i/(int(len(lines)/10))*10)) + '%')
                sys.stdout.flush()
        sys.stdout.write('\n')
    out.close()
    
    # Build BLAST DB
    build_blast_db(out_file, out_file[:out_file.rindex('.')])

##############################
# Arguments
##############################

argv = sys.argv
args = {'-i': '', '-f': '', '-d': '', '-o': './kinortho_out/', '-E': '1e-5',
        '-t': '1', '-I': '1.5', '-e': '1e-200', '-s': 0, '-S': 6}

if argv[1] == '-h' or argv[1] == '-help':
    print_manual()
    sys.exit()

for i in range(1, len(argv), 2):
    key = argv[i]
    if key in args:
        if i+1 < len(argv):
            args[key] = argv[i+1]
        else:
            print('Error: no value for the option "' + key + '"')
            print_manual()
            sys.exit()
    else:
        print('Error: there is no option "' + key + '"')
        print_manual()
        sys.exit()

if args['-i'] == '':
    print('Error: input is required')
    print_manual()
    sys.exit()

if int(args['-s']) not in [0, 1, 2, 3, 4, 5, 6] or int(args['-S']) not in [0, 1, 2, 3, 4, 5, 6]:
    print('Error: start/Stop step must be {0,1,2,3,4,5,6}')
    print_manual()
    sys.exit()

try:
    proteomes, full_seq, domain_seq, output_folder, e_value, thread, inflat, min_e, start_step, stop_step = list(map(
        args.get, ['-i', '-f', '-d', '-o', '-E', '-t', '-I', '-e', '-s', '-S']))
except:
    print_manual()
    sys.exit()

proteomes = os.path.abspath(proteomes) + '/'
output_folder = os.path.abspath(output_folder) + '/'
if full_seq != '':
    full_seq = os.path.abspath(full_seq)
if domain_seq != '':
    domain_seq = os.path.abspath(domain_seq)    

##############################
# Step 0: Initiation
##############################

# Create a protein-species mapping (cannot skip it)
protein_to_species = {}
proteome_names = []
input_listdir = os.listdir(proteomes)
input_listdir.sort()
for file_name in input_listdir:
    if file_name[-6:] == '.fasta':
        proteome_names.append(file_name)

print('Step 0: Initiation')
for i in range(0, len(proteome_names)):
    print('\tReading proteomes (' + str((i+1)) + '/' + str(len(proteome_names)) + '): ' + proteome_names[i])
    sys.stdout.flush()
    file = open(proteomes + proteome_names[i], 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        if line[0:1] == '>':
            toks = line[1:].split(' ')
            seq_id = toks[0]
            protein_to_species[seq_id] = proteome_names[i]

step = int(start_step)
if step == 0:
    
    # Create output folder
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Combine proteomes
    out = open(output_folder + 'proteomes.fasta', 'w')
    for i in proteome_names:
        file = open(proteomes + i, 'r')
        out.write(file.read())
        file.close()
    out.close()

    # Build a searchable database
    build_blast_db(output_folder + 'proteomes.fasta', output_folder + 'proteomes')
    
    if step == int(stop_step):
        sys.exit()
    step += 1

##############################
# Step 1: Homology Search
##############################
if step == 1:
    
    # Full-length pipeline
    if full_seq != '':
        seq_file = full_seq
        db_name = output_folder + 'proteomes'
        out_file = output_folder + 'homologs_full.blast'
        print('Step 1: Homology search\n\tDB: ' + db_name + '\n\tSequence: ' + seq_file)
        sys.stdout.flush()
        homology_search(seq_file, db_name, out_file)

    # Domain-based pipeline
    if domain_seq != '':
        seq_file = domain_seq
        db_name = output_folder + 'proteomes'
        out_file = output_folder + 'homologs_domain.blast'
        print('Step 1: Homology search\n\tDB: ' + db_name + '\n\tSequence: ' + seq_file)
        sys.stdout.flush()
        homology_search(seq_file, db_name, out_file)

    if step == int(stop_step):
        sys.exit()
    step += 1

##############################
# Step 2: Build BLAST DB
##############################
if step == 2:
    
    # Full-length pipeline
    if full_seq != '':
        proteomes_seq = output_folder + 'proteomes.fasta'
        blast_result = output_folder + 'homologs_full.blast'
        out_file = output_folder + 'homologs_full.fasta'
        is_full_pipeline = 1
        print('Step 2: Building BLAST DB')
        rebuild_blast_db(proteomes_seq, blast_result, out_file, is_full_pipeline)

    # Domain-based pipeline
    if domain_seq != '':
        proteomes_seq = output_folder + 'proteomes.fasta'
        blast_result = output_folder + 'homologs_domain.blast'
        out_file = output_folder + 'homologs_domain.fasta'
        is_full_pipeline = 0
        print('Step 2: Building BLAST DB')
        rebuild_blast_db(proteomes_seq, blast_result, out_file, is_full_pipeline)

    if step == int(stop_step):
        sys.exit()
    step += 1

##############################
# Step 3: All-vs-all
##############################
if step == 3:

    # No-query
    if full_seq == '' and domain_seq == '':
        seq_file = output_folder + 'proteomes.fasta'
        db_name = output_folder + 'proteomes'
        out_file = output_folder + 'all-vs-all.blast'
        print('Step 3: All-vs-all Homology Search\n\tDB: ' + db_name + '\n\tSequence: ' + seq_file)
        sys.stdout.flush()
        homology_search(seq_file, db_name, out_file)

    # Full-length pipeline
    if full_seq != '':
        seq_file = output_folder + 'homologs_full.fasta'
        db_name = output_folder + 'homologs_full'
        out_file = output_folder + 'all-vs-all_full.blast'
        print('Step 3: All-vs-all Homology Search\n\tDB: ' + db_name + '\n\tSequence: ' + seq_file)
        sys.stdout.flush()
        homology_search(seq_file, db_name, out_file)

    # Domain-based pipeline
    if domain_seq != '':
        seq_file = output_folder + 'homologs_domain.fasta'
        db_name = output_folder + 'homologs_domain'
        out_file = output_folder + 'all-vs-all_domain.blast'
        print('Step 3: All-vs-all Homology Search\n\tDB: ' + db_name + '\n\tSequence: ' + seq_file)
        sys.stdout.flush()
        homology_search(seq_file, db_name, out_file)
    
    if step == int(stop_step):
        sys.exit()
    step += 1

##############################
# Step 4: Orthology Inference
##############################
if step == 4:
    
    # No-query
    if full_seq == '' and domain_seq == '':
        allvsall_file = output_folder + 'all-vs-all.blast'
        out_file = output_folder + 'orthology_inference.txt'
        is_full_pipeline = 1
        print('Step 4: Orthology Inference')
        orthology_inference(allvsall_file, out_file, is_full_pipeline)

    # Full-length pipeline
    if full_seq != '':
        allvsall_file = output_folder + 'all-vs-all_full.blast'
        out_file = output_folder + 'orthology_inference_full.txt'
        is_full_pipeline = 1
        print('Step 4: Orthology Inference')
        orthology_inference(allvsall_file, out_file, is_full_pipeline)

    # Domain-based pipeline
    if domain_seq != '':
        allvsall_file = output_folder + 'all-vs-all_domain.blast'
        out_file = output_folder + 'orthology_inference_domain.txt'
        is_full_pipeline = 0
        print('Step 4: Orthology Inference')
        orthology_inference(allvsall_file, out_file, is_full_pipeline)

    if step == int(stop_step):
        sys.exit()
    step += 1

##############################
# Step 5: Cluster Analysis
##############################
if step == 5:
    
    # No-query
    if full_seq == '' and domain_seq == '':
        ortholog_file = output_folder + 'orthology_inference.txt'
        out_file = output_folder + 'graph.abc'
        is_full_pipeline = 1
        print('Step 5: Cluster Analysis')
        build_graph(ortholog_file, out_file, is_full_pipeline)

        graph_file = out_file
        out_name = output_folder + 'clustering'
        clustering(graph_file, out_name)

        homolog_file = ''
        cluster_name = out_name
        out_file = output_folder + 'Results_No-query.txt'
        filtering(homolog_file, ortholog_file, graph_file, cluster_name, out_file, is_full_pipeline)

    # Full-length pipeline
    if full_seq != '':
        ortholog_file = output_folder + 'orthology_inference_full.txt'
        out_file = output_folder + 'graph_full.abc'
        is_full_pipeline = 1
        print('Step 5: Cluster Analysis')
        build_graph(ortholog_file, out_file, is_full_pipeline)

        graph_file = out_file
        out_name = output_folder + 'clustering_full'
        clustering(graph_file, out_name)

        homolog_file = output_folder + 'homologs_full.blast'
        cluster_name = out_name
        out_file = output_folder + 'Results_Full-length.txt'
        filtering(homolog_file, ortholog_file, graph_file, cluster_name, out_file, is_full_pipeline)

    # Domain-based pipeline
    if domain_seq != '':
        ortholog_file = output_folder + 'orthology_inference_domain.txt'
        out_file = output_folder + 'graph_domain.abc'
        is_full_pipeline = 0
        print('Step 5: Cluster Analysis')
        build_graph(ortholog_file, out_file, is_full_pipeline)

        graph_file = out_file
        out_name = output_folder + 'clustering_domain'
        clustering(graph_file, out_name)

        homolog_file = output_folder + 'homologs_domain.blast'
        cluster_name = out_name
        out_file = output_folder + 'Results_Domain-based.txt'
        filtering(homolog_file, ortholog_file, graph_file, cluster_name, out_file, is_full_pipeline)

    if step == int(stop_step):
        sys.exit()
    step += 1

##############################
# Step 6: Combine Results
##############################
if step == 6:
    if full_seq != '' and domain_seq != '':
        full_result = output_folder + 'Results_Full-length.txt'
        domain_result = output_folder + 'Results_Domain-based.txt'
        out_file = output_folder + 'Results_Overlapping.txt'
        print('Step 6: Combining Results')
        combin_result(full_result, domain_result, out_file)

