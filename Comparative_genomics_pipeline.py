#!/usr/bin/python3
"""
@author: Kirti Deepak Chhatlani, Katherine Duchesneau, Xu Qiu, Ashika Ramesh, Huy Tran, Joseph Luper Tsenum, Cheng Zhang
"""
import subprocess
import os
import argparse
import csv



def fastANI(input_dir,output_folder):
    compGenFolder = os.path.dirname(__file__)
    for filename in os.listdir(input_dir):
        f = os.path.join(input_dir, filename)
        k = f.split("/")[-1]
        print(f,file= open(output_folder+'/newfile.txt','a'))

    print('**********Running fastANI**********')
    os.system("fastANI --ql "+output_folder+"/newfile.txt --rl "+output_folder+"/newfile.txt -o "+output_folder+"/fastani.out")

    print('**********Making a distance matrix in phylip format**********')
    os.system(compGenFolder+"/pairwise_identities_to_distance_matrix.py "+output_folder+"/fastani.out > "+output_folder+"/fastani.phylip")

    print('**********Converting the phylip output to newick format**********')
    os.system("Rscript "+compGenFolder+"/bionj_tree.R "+output_folder+"/fastani.phylip "+output_folder+"/fastani.newick")
    os.system('rm '+output_folder+'/newfile.txt')

def ksnp3(input_dir,output_folder):
    col1 = []
    col2 = []
    for filename in os.listdir(input_dir):
        path = os.path.abspath(input_dir)
        f = os.path.join(path, filename)
        col1.append(f)
        name=os.path.basename(f)
        name2=name[ : 7]
        col2.append(name2)

    with open('in_list.txt', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(col1, col2))
    print('**********kSNP running**********')
    subprocess.run(['MakeFasta', 'in_list.txt', 'kchooser.fasta'])
    subprocess.run(['Kchooser','kchooser.fasta'])
    subprocess.run(['rm','kchooser.fasta'])
    with open("Kchooser.report") as f:
        lines = f.read()
        first = lines.split('\n', 1)[0]
    K=first[22: 24]
    subprocess.run(['kSNP3', '-in', 'in_list.txt', '-k', K,'-outdir',output_folder+"_ksnp", '-ML', '-NJ'])
    krm=["rm Kchooser.report"]
    subprocess.run(krm,shell = True)
    print('**********The output files are places in the output directory**********')

def chewbecca(input_dir,output_folder):
    print('**********Running chewbecca**********')
    subprocess.call(["mkdir", output_folder+"/chewBBACA"], universal_newlines=True)
    subprocess.call(["mkdir", output_folder+"/chewBBACA/Listeria"], universal_newlines=True)


    subprocess.check_output(["chewBBACA.py", "DownloadSchema","-sp","6","-sc","1","-o", \
        output_folder+"/chewBBACA/Listeria"], universal_newlines=True)

    subprocess.check_output(["mv", output_folder+"/chewBBACA/Listeria/lmonocytogenes_Pasteur_cgMLST", \
        output_folder+"/chewBBACA/Listeria_db"], universal_newlines=True)

    subprocess.check_output(["rm","-r", output_folder+"/chewBBACA/Listeria"], universal_newlines=True)

    subprocess.check_output(["chewBBACA.py", "AlleleCall", "-i", input_dir, "-g", output_folder+"/chewBBACA/Listeria_db", "-o", \
        output_folder+"/chewBBACA/Allele", "--cpu", "4", "--ptf", output_folder+"/chewBBACA/Listeria_db/Listeria_monocytogenes.trn"], universal_newlines=True)


    directory = output_folder+"/chewBBACA/Allele"
    for filename in os.listdir(directory):
        fname = os.path.join(directory, filename)
        # print(fname)
        subprocess.check_output(["mv", fname, output_folder+"/chewBBACA/Alleleresults"], universal_newlines=True)


    subprocess.check_output(["rm","-r", output_folder+"/chewBBACA/Allele"], universal_newlines=True)

    subprocess.check_output(["chewBBACA.py", "TestGenomeQuality", "-i", output_folder+"/chewBBACA/Alleleresults/results_alleles.tsv", \
        "-n", "1", "-t", "100", "-s", "5", "-o", output_folder+"/chewBBACA/Test"], universal_newlines=True)

    subprocess.check_output(["chewBBACA.py", "ExtractCgMLST", "-i", output_folder+"/chewBBACA/Alleleresults/results_alleles.tsv", "--r", \
        output_folder+"/chewBBACA/Alleleresults/RepeatedLoci.txt", "--g", output_folder+"/chewBBACA/Test/removedGenomes.txt", \
        "-o", output_folder+"/chewBBACA/Extract", "--t","0.95"], universal_newlines=True)

    f = open(output_folder+"/chewBBACA/chewbbaca.nwk","w")
    subprocess.call(["grapetree", "-p", output_folder+"/chewBBACA/Extract/cgMLST.tsv", "-m", "NJ"], stdout=f)
    print('**********The output files are placed in the output folder**********')
  
def abricate(input_dir, output_folder):
    print('**********Searching for AMR genes**********')
    os.system("abricate --db card " + input_dir + "/*.fasta > "+output_folder+"/results.tab")
    os.system("abricate --summary "+output_folder+"/results.tab > "+output_folder+"/antibiotic_summary.tab")
    print('**********Antibiotic profile written**********')
    print('**********Searching for virulence genes**********')
    os.system("abricate --db vfdb " + input_dir + "/*.fasta > "+output_folder+"/results.tab")
    os.system("abricate --summary "+output_folder+"/results.tab > "+output_folder+"/virulence_summary.tab")
    print('**********Virulence profile written**********')
    print('ABRicate DONE')

def main():
    parser = argparse.ArgumentParser(description='A generalized pipeline for comparative genomics.')
    parser.add_argument('-i', '--input', type=str, required=True, help="Path for directory containing contig FASTA files.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output folder name.")
    parser.add_argument('-f',action = 'store_true', help="Option for using fastani as a tool")
    parser.add_argument('-ks',action = 'store_true',help="Option for using ksnp3 as a tool")
    parser.add_argument('-c', action = 'store_true', help ="Option for using chewbecca as a tool")
    parser.add_argument('-a',action = 'store_true', help="Option for using the tool abricate")
    args = parser.parse_args()

    input_dir = args.input
    output_folder = args.output

    if args.f:
        fastANI(input_dir,output_folder)

    if args.ks:
        ksnp3(input_dir,output_folder)

    if args.c:
        chewbecca(input_dir,output_folder)
    
    if args.a:
        abricate(input_dir,output_folder)

if __name__ == "__main__":
    main()
