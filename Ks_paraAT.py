import sys, os, glob
import argparse



def cds_maker (htm_gb, afg_out, gb_out, dicSeq, protein_gb, cds_gb, seqments) :
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    fwcds = open(afg_out + '.cds', 'w')
    for seq in SeqIO.parse(open(afg_out), 'fasta') :
        tem_seq = seq.seq
        cdsPos = []
        for block in seqments :
            tem_cds_pos = []
            for bl_ion in block.split() :
                bl_st = int(bl_ion)
                gap_count = seq.seq[:int(bl_st)].count("-")
                if tem_cds_pos == [] :
                    gap_count += 1
                cd_st = (int(bl_st) - gap_count) * 3
                tem_cds_pos.append(cd_st)
            cdsPos.append(tem_cds_pos)
        new_seq_id = '---'.join(seq.id.split('---')[:-1])
        full_cds = dicSeq[new_seq_id]['cds']
        new_cds = '>%s\n'%seq.id
        for cds_ion in cdsPos :
            new_cds += full_cds[cds_ion[0]:cds_ion[1]]
        new_cds += '\n'
        new_cds_line = str(new_cds)
        fwcds.write(new_cds_line)
    fwcds.close()
    

    fgb_n = open(gb_out + '.new', 'w')
    for fgb_seq in SeqIO.parse(open(gb_out), 'fasta') :
        new_seq = str(fgb_seq.seq).replace(' ', '')
        new_line = '>%s\n%s\n'%(fgb_seq.id, new_seq)
        fgb_n.write(new_line)
    fgb_n.close()
    
    cat_cmd = 'cat %s.new >> %s\n'%(gb_out, protein_gb)
    cat_cmd += 'cat %s.cds >> %s\n'%(afg_out, cds_gb)
    cat_cmd += 'rm %s %s %s %s %s.new %s.cds'%(afg_out[:-4], gb_out, htm_gb, afg_out, gb_out, afg_out)
    for cmd_line in cat_cmd.split('\n') :
        os.system(cmd_line)
    return (protein_gb, cds_gb)



def seq_renewer (htm_gb, afg_out, gb_out, dicSeq, protein_gb, cds_gb) :
    fhgb = open(htm_gb, 'r')
    line = fhgb.readline()
    while line != '' :  
        if line.find("Flanks:") != -1 :
            seqments = []
            for se_unit in line.split('[')[1:]:
                seqments.append(se_unit.split(']')[0])
        line = fhgb.readline()
    fhgb.close()

    if seqments != [] :
        protein_gb, cds_gb = cds_maker(htm_gb, afg_out, gb_out, dicSeq, protein_gb, cds_gb, seqments)
    return (protein_gb, cds_gb)



def make_input(in_file, sp_list, dicSeq, cpu_num) :
    h_file_list = []
    fr = open(in_file, 'r')
    line = fr.readline()
    dic = {}
    while line != '' :
        unit = line.strip().split('\t')
        #if unit[-1] == 'Ref' :
        sp, clu, gene_id = unit[1], unit[2], unit[3]
        dic.setdefault(clu, {}).setdefault(sp, {}).setdefault(gene_id, 1)
        line = fr.readline()
    fr.close()

    h_pa_out = '%s_%s.homolog'%(sp_list[0], sp_list[1])
    info_h_pa_out = '%s_%s.homolog.info'%(sp_list[0], sp_list[1])
    protein_gb = '%s_%s.protein'%(sp_list[0], sp_list[1])
    cds_gb = '%s_%s.cds'%(sp_list[0], sp_list[1])
    dic_clu = {}
    h_file_list.append(h_pa_out)
    fw = open(h_pa_out, 'w')
    fw2 = open(info_h_pa_out, 'w')
    for clu_ion in dic :
        pare_list = []
        try :
            one_id_list = dic[clu_ion][sp_list[0]].keys()
            try :
                two_id_list = dic[clu_ion][sp_list[1]].keys()
            except KeyError :
                two_id_list = []
        except KeyError :
            one_id_list = []
        if len(one_id_list) > 0 and  len(two_id_list) > 0 :
            print one_id_list, two_id_list
            fa = open(clu_ion, 'w')
            for one_gene_ion in one_id_list :
                for two_gene_ion in two_id_list :
                    fa = open (clu_ion, 'w')
                    one_seq = '>%s---%s\n%s\n'%(one_gene_ion, two_gene_ion, dicSeq[one_gene_ion]['protein'])
                    fa.write(one_seq)
                    two_seq = '>%s---%s\n%s\n'%(two_gene_ion, one_gene_ion,  dicSeq[two_gene_ion]['protein'])
                    fa.write(two_seq)
                    info_new_line = '%s\t%s\t%s\n'%(clu_ion, one_gene_ion, two_gene_ion)
                    match_name = '%s---%s-%s---%s'%(one_gene_ion, two_gene_ion, two_gene_ion, one_gene_ion)
                    dic_clu[match_name] = clu_ion
                    fw2.write(info_new_line)
                    homolo_line = '%s---%s\t%s---%s\n'%(one_gene_ion, two_gene_ion, two_gene_ion, one_gene_ion)
                    fa.close()
        
                    afg_out = '%s.afg'%clu_ion
                    gb_out = '%s-gb'%(afg_out)
                    htm_gb = '%s.htm'%(gb_out)

                    mafft_cmd = 'mafft --ep 0.123 --thread %s --auto %s > %s'%(cpu_num, clu_ion, afg_out)
                    os.system(mafft_cmd)

                    gblock_cmd = 'Gblocks %s -t=p -b3=8'%(afg_out)
                    os.system(gblock_cmd)
                    
                    if dicSeq[one_gene_ion].has_key('cds') and dicSeq[two_gene_ion].has_key('cds'):
                        protein_gb, cds_gb = seq_renewer(htm_gb, afg_out, gb_out, dicSeq, protein_gb, cds_gb)
                        fw.write(homolo_line)
    fw.close()
    fw2.close()
    return (h_file_list, dic_clu, protein_gb, cds_gb)



def paraAT_executer(h_file_list, dic_clu, protein_gb, cds_gb):
    for homo_file in h_file_list :
        out_put = homo_file+ '.kaks.out'
        cmd = "perl ParaAT.pl -h %s -n %s -a %s -p proc -kaks -f axt -o %s -g -m mafft"%(homo_file, cds_gb, protein_gb, out_put)
        os.system(cmd)
        out_file = './%s/*_aln.axt.kaks'%out_put
        final_out = (homo_file.split('.'))[0] + '.Ks.txt'

        title_line = "Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites(0:2:4)\tSubstitutions\tS-Substitutions\tN-Substitutions\tFold-S-Substitutions(0:2:4)\tFold-N-Substitutions(0:2:4)\tDivergence-Time\tSubstitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\tGC(1:2:3)\tML-Score\tAICc\tAkaike-Weight Model\n"
        fwo = open(final_out, 'w')
        fwo.write(title_line)
        dicKs = {}
        for out_f in glob.glob(out_file):
            #print out_f
            match_acc = out_f.split('/')[-1].split('.cds_aln.axt')[0]
            match_clu_id = dic_clu[match_acc]
            f = open(out_f, 'r')
            line = f.readline()
            line = f.readline()
            while line != '' :
                unit = line.split('\t')
                if int(line.split('\t')[6]) > 99 :
                    try :
                        pre_ks = float(dicKs[match_clu_id].split('\t')[3]) 
                    except KeyError :
                        pre_ks = 0 
                    if unit[3] != 'NA' and pre_ks < float(unit[3]) :
                        dicKs[match_clu_id] = line
                            
                    #fwo.write(line)
                line = f.readline()
            f.close()
        for match_line in dicKs :
            fwo.write(dicKs[match_line])
    fwo.close()
    return (1)

def make_dicSeq(cds, protein):
    dicSeq = {}
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    for fcp in [cds, protein] :
        if fcp.lower().find('protein') != -1 :
            file_index = 'protein'
        else :
            file_index = 'cds'
        for seq in SeqIO.parse(open(fcp), 'fasta') :
            seq_acc = seq.id.split()[0].split('|')[-1]
            dicSeq.setdefault(seq_acc, {}).setdefault(file_index, seq.seq)
    return (dicSeq)
        
        

if __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'clustered file with gene list', required = True)
    parser.add_argument('-s', help = 'sp_list with , deliminated format', required = True)
    parser.add_argument('-p', help = 'protein seq', required = True)
    parser.add_argument('-n', help = 'cds seq', required = True)
    parser.add_argument('-a', help = 'cpu number', required = True)

    args = parser.parse_args()
    in_file = args.i
    sp_name = args.s
    sp_list = sp_name.split(',')
    cds = args.n
    protein = args.p
    cpu_num = args.a

    dicSeq = make_dicSeq (cds, protein)
    h_file_list, dic_clu, protein_gb, cds_gb = make_input(in_file, sp_list, dicSeq, cpu_num)
    ans = paraAT_executer(h_file_list, dic_clu, protein_gb, cds_gb)

