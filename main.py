#!/usr/bin/python3
"""
1: verificar se é single ou pair-end
2: calcular GC e gerar csv.
3: passar csv no WEKA pra gerar o arff com os clusteres
4: scprit python para dividir os clusteres lembrando se é single ou pair end.
5: calculo do Kmer com kmerStream
6: ordenar a coluna F0 e pergar os 5 primeiros da coluna k.
"""
import sys, os, subprocess, _thread as thread, time, multiprocessing
from shutil import copyfile
from pandas import DataFrame
from sklearn.cluster import KMeans
import time
global path
import json

path = os.path.dirname(os.path.realpath(__file__))

def w_log(log_file , log, start_time):
    elapsed_time = time.time() - start_time
    final = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    Logs = [
        '0\t    start',
        '1\t    upload complete', 
        '2\t    load gc content', 
        '3\t    Cluster',
        '4\t    std assembly',
        '5\t    cluster assembly',
        '6\t    map contigs',
        '7\t    close gaps',
        '8\t    quast'
        ]
        
    with open(log_file, 'a+') as l:
        l.write(Logs[log] + '\t' + str(final) + '\n')


def scaffolds_to_contigs(scaffolds, out):
    file = open(scaffolds)
    contig = ''
    id_fasta = file.readline()
    sequences = file.read()
    sequences = sequences.split('N')
    for seq in sequences:
        if seq.strip() != '':
            contig += id_fasta + seq + '\n'
    

    with open(out, 'w') as o:
        o.write(contig)

    return out


def calc_gaps(contig):
    gaps = 0 
    
    seq = ""
    with open(contig) as fasta:
        for line in fasta:
            if line[0] != ">":
                seq += line.strip()

    Ns = seq.count("N")
    start = 0
    for base in seq:
        if base == "N":
            if start == 0:
                gaps +=1
                start = 1
        else:
            start = 0
    
    return [gaps, Ns]

def dedupe(contig, out, identity = '100'):
    cmd = path + "/bbmap/dedupe.sh in=" + contig + \
    " out=" + out + " s=1 e=2 minidentity=" + identity + " overwrite=true"

    os.system(cmd)
    return out


def fastq_rename(fastq, outfile):
    cmd = "awk \'{if( (NR-1)%4 ) print; else printf(\"@id_%d\\n\",cnt++)}\' " + fastq + " > " + outfile
    os.system(cmd)

    return os.path.abspath(outfile)

def fasta_rename(fasta, outfile, name = 'id_'):
    cmd = "awk \'/^>/{print \">" + name + "\" ++i; next}{print}\' < " + fasta + " > " + outfile
    os.system(cmd)
    
    return os.path.abspath(outfile)

def bioawk(std_contigs, project, t_contigs, u_contigs):
    t_contigs = project + "/" + t_contigs
    u_contigs = project + "/" + u_contigs
    
    cmd = path + "/bioawk -c fastx '(length($seq) > 500){ print \">\"$name\"\\n\"$seq }' " + std_contigs + " > " + t_contigs
    os.system(cmd)
    
    cmd = path + "/bioawk -c fastx '(length($seq) < 500){ print \">\"$name\"\\n\"$seq }' " + std_contigs + " > " + u_contigs
    os.system(cmd)
    
    return t_contigs,u_contigs

def cat_files(files, out):
    cmd = "cat "
    for f in files:
        cmd += f + " " 
    cmd +=  " > " + out

    os.system(cmd)
    
    return out

def gaa(contig_list, project, out,):
    tmp = project + 'merge/'
    ref = contig_list.pop(0)

    while contig_list:
        print(ref, contig_list[0])
        cmd = path + "/gaa/gaa.pl -t " + ref + " -q  " + contig_list.pop(0) + " -o " + tmp + " > " + project + "gaa.tmp "
        os.system(cmd)
        os.system("cat " + tmp + "*.fasta.fa > " + out)
        ref = out
        os.system("rm -r " + tmp)

    return out




def quast(contig_list, ref_fasta, project, gff=None):
    out = project + 'quast/'
    if not os.path.exists(out):
        os.mkdir(out)
    cmd = "quast.py -R " + ref_fasta + \
        " -o " + out + " "
    
    for contig in contig_list:
        if contig:
            cmd += contig + ' '

    if gff:
        cmd += " -G " + gff

    cmd += " > " + out + "quast.log"

    os.system(cmd)

    return out


def calc_gc(seq):
    """
    calcula a porcentagem de gc da sequencia  (GC/ACGT) * 100
    retorna a porcentagem.
    """
    counts = {'C' : 0, 'G' : 0}
    for bases in seq:
        if(bases == 'G' or bases == 'C'):
            counts[bases] += 1
    gc = ((counts['G'] + counts['C']) / len(seq)) * 100
    
    return  int(gc)

def gc_content(read):
    x = open(read)
    gc = []
    while True:
        line = x.readline()
        if line == '':
            break
        gc.append(int(calc_gc(x.readline().strip())))
        x.readline()
        x.readline()

    return gc

def read_cluster(cluster):
    """
    Lê a saida do WEKA com as informações dos clustererers e salva em uma lista;
    [cluster] = saida do weka, out.arff
    """
    alist = []
    x = open(cluster).readlines()
    while True:
        x.remove(x[0])
        if x[0].find("@data") > -1 :
            x.remove(x[0])
            x.remove(x[0])
            break
    for line in x:
        line = line.rstrip().split(',')
        alist.append(line[2][-1])

    return alist

def read_kmer(filename):
    """
    Lê o arquivo de saida do KmerStream orderna decrescente pela coluna F0 e pega os 5 primeiros K
    """
    kmer = open(filename).readlines()
    kmer.remove(kmer[0])
    x = []
    k = ''
    for line in kmer:
        l = line.split("\t")
        x.append([int(l[2]),int(l[1])])
        i = 0

    for line in (sorted(x,reverse=True)[0:5]):
        k+= str(line[1]) + ","

    return k[:-1]
    
def get_argx(argx):
    x = {'read1'    : None , 
        'read2'     : None, 
        'reference' : None, 
        'scaffolds' : None, 
        'project'   : None, 
        'email'     : None, 
        'gff'       : None,
        'gft'       : 'contiguator',
        'threads'   : '24', 
        'outcontig' : None, 
        'identity'  : '100',
        'clusters'  : 4,
        'assembler' : 'spades'
        }
    
    for a in argx:
        if a in x:
            x[a] = argx[argx.index(a) + 1]

    return x

def contiguator(contigs, reference, project, threads=16):
    out = project + "contiguator/"
    

    if not os.path.isdir(out):
        os.mkdir(out)
    
    run = 'python ' + path + '/CONTIGuator/CONTIGuator.py'
    cmd = 'cd ' + out + ' && '
    cmd += run + ' -c ' + contigs + ' -r ' +\
    reference + ' -t ' + str(threads) + ' > contiguat.log'
    
    os.system(cmd)
    
    cmd = "find " + out + " -name 'PseudoContig.fsa' -exec cat {} + > " + project + "/PseudoContig.fsa"
    os.system(cmd)
    
    
    return project + "/PseudoContig.fsa"
    
def merge_fastq(fastq_a,fastq_b, out_name):
    """
    converte os dois fastq de entrada para interlivered
    """
    a = open(fastq_a)
    b = open(fastq_b)
    fastq_out = ''
    i = 1
    while True:
        line = a.readline()
        if line == "":
            break
        fastq_out += '@id-' + str(i) + '\n'
        i+=1
        b.readline() #id b
        fastq_out += a.readline().strip()+ b.readline() 
        fastq_out += "+\n"
        a.readline()#sep a
        b.readline()#sep b

        fastq_out += a.readline().strip() + b.readline()

        if line == "":
            break
    
    write_file(out_name, fastq_out)
      
    
    return out_name


def edena(read1, project, prefix,  read2=None, threads='24'):
    print('edena...')

    out = project + prefix
    if not os.path.isdir(out):
        os.mkdir(out)


    cmd = path + "/edena -nThreads " + str(threads)

    if not read2:
        cmd += " -singleEnd " + read1

    else:
        cmd += " -paired " + read1 + " " + read2

    cmd = "cd " + out + " && " + cmd + " > edena.log"
    
    print(cmd)
    os.system(cmd)
    cmd = "cd " + out + " && " +  path + "/edena -e out.ovl > edena.log"
    print(cmd)
    os.system(cmd)

    return out + '/out_contigs.fasta'


def spades(read1, project, prefix,  read2=None, trusted_contigs=None, S=None, threads='16', untrusted_contigs = None, kmer = None, interlaced = False):
    out = project + prefix + 'spades/'
    if not os.path.exists(out):
        os.mkdir(out)

    cmd =  path + '/spades/bin/spades.py -o ' + out + ' -t ' + str(threads) + ' '
    if read2:
        cmd += '-1 ' + read1 + ' -2 ' + read2 + ' '
        if S:
            cmd += '-s ' + S + ' '
    else:
        if interlaced:
            cmd += '--12 ' + read1
        else:
            cmd += '-s ' + read1

    if trusted_contigs:
        cmd += ' --trusted-contigs ' + trusted_contigs + ' '
        
    if untrusted_contigs:
        cmd += ' --untrusted-contigs ' + untrusted_contigs + ' '

    if kmer:
        cmd+= " -k " + kmer

    cmd += ' > ' + out + 'execution.log'

    write_file(out + "comandline.txt", cmd)
    os.system(cmd)

    if os.path.isfile(out + 'contigs.fasta'):
        return os.path.abspath(out + 'contigs.fasta')
    else:
        return False
       
def fgap(scaffolds, contigs, out_prefix, threads = '16'):
    cmd = path + "/fgap/run_fgap.sh " + path + "/mcr/v717 " +\
    "-d " + scaffolds + " -a " + contigs + \
    " -b " + path + "/fgap/blast/ -o " + out_prefix + " -t " +\
    str(threads) + " > " + out_prefix + ".cmd"
        
    os.system(cmd)

    out_file = out_prefix + '.final.fasta'

    if os.path.exists(out_file):
        return out_file
    else:
        return False

def gapblaster(scaffolds, contigs, project, threads = '24'):
    out = project  + 'gapblaster/'
    if not os.path.isdir(out):
        os.mkdir(out)
    
    run = "./jre/bin/java -Xms6000m -Xmx12000m -jar "

    cmd =  run + path + "/gapblaster-cli.jar " + scaffolds + " " +\
    contigs + " " +  out + " 2>&1 | tee " + project + "gapblasterlog_terminal.txt "
    os.system(cmd)
    
    out_file = out + 'scaffolds.fasta'

    if os.path.exists(out_file):
        return  out_file
    else:
        return False


def gmcolser(scaffolds, contigs, project, threads = '24'):
    cmd =  path + "/GMcloser/gmcloser -t " + path + "/" + scaffolds + " -q " + path + "/" \
    + contigs + " -c -n " + str(threads) + " -p " + "gmcolser > gmcloser.log"
    print(cmd)
    os.system("cd " + project + " && " + cmd)

    outfile = project + "gmcolser.gapclosed.fa"

    if os.path.exists(outfile):
        return outfile
    else:
        return False

def send_mail(email, project):
    os.system("curl 'http://200.239.92.140/ClustGClosed/PHP/send_mail.php?email=" + email + "&project=" + project + "'")
    
def medusa(contigs, reference, project, threads=16):
    ref_dir = project + 'reference/'
    if not os.path.isdir(ref_dir):
        os.mkdir(ref_dir)
    
    copyfile(reference, ref_dir + 'reference.fna')

    cmd = 'java -jar ' + path + '/medusa/medusa.jar -d -f ' + ref_dir +\
    ' -i ' + contigs + ' -o ' + project + '/PseudoContig.fsa ' +\
    ' -random 5 -scriptPath ' + path + '/medusa/medusa_scripts/' +\
    ' > ' + project + 'medusa.log'
    
    os.system(cmd)

    return project + 'PseudoContig.fsa'
        
def write_file(filename, data, mode="w"):
    with open(filename, mode) as out:
        out.write(data)

def kmer_stream(read1, project,read2=None, threads=16):
    cmd = path + "/KmerStream/KmerStream -k \
    11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,55,59,61,65,67,69,71,75,79,81,85,89,91,93,95 \
    --tsv -t " + str(threads) + " -o " + project + "result.kmer " + read1
    if read2:
        cmd += " " + read2
    os.system(cmd)
    kmer = read_kmer(project + 'result.kmer')

    return kmer

def weka_kmeans(input_csv, output, n_cluster):
        cmd = "java -classpath " + path + "/weka/weka.jar:. \-Xms16000m -Xmx18000m \
        weka.filters.unsupervised.attribute.AddCluster -W 'weka.clusterers.SimpleKMeans \
        -N " + str(n_cluster) +"'  -i " + input_csv + " -o " + output +\
        " -I 1-1"
        
        write_file(output + '.cmd', cmd)
        os.system(cmd)

        return output

def split_clusters(c_labels, read1, project, n_cluster, read2 = None):  
    clusters = {}  
    if read2:
        fastq1 = open(read1).read()[1:].split('\n@id')   
        fastq2 = open(read2).read()[1:].split('\n@id')   
        for i in range(len(c_labels)):
            write_file(project + 'cluster_A_' + str(c_labels[i]) + '.fastq', '@id' + fastq1[i] + '\n', 'a+')
            write_file(project + 'cluster_B_' + str(c_labels[i]) + '.fastq', '@id' + fastq2[i] + '\n', 'a+')

        clusters['A'] = [project + 'cluster_A_' + str(i) + '.fastq' for i in range(n_cluster)]
        clusters['B'] = [project + 'cluster_B_' + str(i) + '.fastq' for i in range(n_cluster)]

    else:
        fastq1 = open(read1).read()[1:].split('\n@id')   
        for i in range(len(c_labels)):
            write_file(project + 'cluster_A_' + str(c_labels[i]) + '.fastq', '@id' + fastq1[i] + '\n', 'a+')
        
        clusters['A'] = [project + 'cluster_A_' + str(i) + '.fastq' for i in range(n_cluster)]
        clusters['B'] = [None for i in range(n_cluster)]

    return clusters

def save_project(project, pname):
    tmp = project + pname
    if not os.path.isdir(tmp):
        os.mkdir(tmp)

    file_list = ['std_assembly.fasta',
                 'std_mod.fasta',
                 'quast',
                 'contigs_final.fasta',
                 'scaffolds_final.fasta',
                 'log.txt',
                 'gaps.json'
                 ]
    for f in file_list:
        os.system('cp -r ' + project + f + ' ' + tmp)
    

    os.system('cd ' + project + ' && tar -czvf ../../download/' + pname + '.tar.gz ' + pname)
    #os.system('rm -r ' + project) #apaga projeto


def main(argx):
    read1       = argx['read1']
    read2       = argx['read2']
    gff         = argx['gff']
    gft         = argx['gft']
    project_name = argx['project']
    reference   = argx['reference']
    threads     = argx['threads']
    scaffolds    = argx['scaffolds']
    n_cluster   = int(argx['clusters'])
    project     = path + '/projects/' + argx['project'] + '/'
    
    log_file = project + 'log.txt'
    s_time = time.time()
    

    if not os.path.isdir(project):
        os.mkdir(project)

    w_log(log_file, 1, s_time)

    print('rename fastq')
    read1 = fastq_rename(read1, project + 'read1.fastq')
    if read2:
        read2 = fastq_rename(read2, project + 'read2.fastq')


    print("loading gc content..")
    if read2:
        merged_fastq = merge_fastq(read1, read2, project + 'merged.fastq')
        gc_list = gc_content(merged_fastq)
    else:
        gc_list = gc_content(read1)
    
    w_log(log_file, 2, s_time)

    print("clustering Kmeans...")        
    kmeans = KMeans(n_clusters = n_cluster).fit(DataFrame(gc_list))
    print("separating files by groups...")        
    clusters = split_clusters(kmeans.labels_ , read1, project, n_cluster, read2)
    w_log(log_file, 3, s_time)
    cluster_contig = []

    if not scaffolds:
        print("Std assembly..")
        if argx['assembler'] == 'spades':
            kmer = kmer_stream(read1, project, read2, threads = threads)
            std_contigs = spades(read1, project, 'std_', threads=threads, kmer = kmer, read2 = read2)
        
        else:
            std_contigs = edena(read1, project, 'std_', threads=threads, read2 = read2)
        
        std_contigs = fasta_rename(std_contigs, project + 'std_mod.fasta', name = 'id_')
        
        t_contigs, u_contigs = bioawk(std_contigs, project, 'trusted_contigs.fasta', 'untrusted_contigs.fasta')

    else:
        scaffolds = fasta_rename(scaffolds, project + 'scaffolds.fasta')
        std_contigs = scaffolds
        t_contigs = scaffolds

    w_log(log_file, 3, s_time)

    for i in range(len(clusters['A'])):
        print("assembling cluster " + str(i + 1))

        if argx['assembler'] == 'spades':
            kmer = kmer_stream(clusters['A'][i], project, clusters['B'][i], threads = threads)
            contig = spades(clusters['A'][i], project, 'ass_cluster_' + str(i), threads=threads, kmer = kmer, read2= clusters['B'][i], trusted_contigs = t_contigs, untrusted_contigs = u_contigs)
        
        else:
            contig = edena(clusters['A'][i], project, 'ass_cluster_' + str(i), threads=threads, read2= clusters['B'][i])
        
        if contig:
            contig = fasta_rename(contig, project + 'contig_' + str(i) + '.fasta', name =  'Contig0.')         
            cluster_contig.append(contig)
        i+=1
    
    w_log(log_file, 5, s_time)
    cluster_contigs = cat_files(cluster_contig, project + 'clusters_contigs.fasta')
    #cluster_contigs = gaa(cluster_contig, project, project + 'clusters_contigs.fasta')
    print('removing redundance...')
    #cluster_contigs = dedupe(cluster_contigs, project + 'clusters_nr_contigs.fasta')

    

    print("mapping contigs...")
    if reference:
        reference =  fasta_rename(reference, project + 'reference.fasta', name = 'id_')
        if 'contiguator' in gft:
            print("contiguator")
            std_contigs = contiguator(std_contigs, reference, project)
        else:
            print("medusa")
            std_contigs = medusa(std_contigs, reference, project)

    else:
        reference = std_contigs            

    w_log(log_file, 6, s_time)

    
    cluster_contigs = fasta_rename(cluster_contigs, project + 'clusters_mod.fasta')

    print("fgap")
    fgap_contigs = fgap(std_contigs, cluster_contigs, project + 'fgap', threads = threads)
    if not fgap_contigs:
        fgap_contigs = std_contigs

    print("gapblaster")
    gp_contigs = gapblaster(fgap_contigs, cluster_contigs, project, threads = threads)
    if not gp_contigs:
        gp_contigs = fgap_contigs
    

    print("gmcloser")
    gmcloser_contigs = fgap(gp_contigs, cluster_contigs, project + 'combined', threads = threads)
    if not gmcloser_contigs:
        gmcloser_contigs = gp_contigs

    os.system("cp " + gmcloser_contigs + " " + project + "scaffolds_final.fasta")

    gaps = {}
    gaps['before'] = calc_gaps(std_contigs)
    gaps['after'] = calc_gaps(gmcloser_contigs)
    gaps = json.dumps(gaps)
    write_file(project + 'gaps.json', gaps)

    gmcloser_contigs = scaffolds_to_contigs(gmcloser_contigs, project + 'contigs_final.fasta')    
    std_contigs = scaffolds_to_contigs(std_contigs, project + 'std_assembly.fasta')

    print("quast")
    w_log(log_file, 7, s_time)

    quast([std_contigs, gmcloser_contigs], reference, project, gff = gff)
    save_project(project, project_name)
    w_log(log_file, 8, s_time)
    
    if 'email' in argx:
        send_mail(argx['email'], argx['project'])

    
if 'read1' in sys.argv:
    argx = get_argx(sys.argv)
    main(argx)
    