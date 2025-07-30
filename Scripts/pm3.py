
# -*- coding: UTF-8 -*-
import sys
import subprocess
import argparse
import re
import os
from multiprocessing import Pool

def parse_args(dir_cache_default, transcript_file_default):
    dir_cache_path_default = os.path.expanduser(dir_cache_default)
    print('dir_cache_default===',dir_cache_default)
    print('dir_cache_path_default====',dir_cache_path_default)
    transcript_path_default = os.path.expanduser(transcript_file_default)
    print('transcript_file_default===',transcript_file_default)
    print('transcript_path_default====',transcript_path_default)

    parse = argparse.ArgumentParser(description="Extract AR genes and case/family samples, VEP annotation , gnomAD_joint_grpmax_AF annotation, variants filtration, RGnet construction, PM3 tagging.")
    parse.add_argument('--vcf',type=str,required=True,help='The vcf file.')
    parse.add_argument('--bed',type=str,required=True,help='The bed file of AR genes.')
    parse.add_argument('--sam',type=str,required=True,help='The case or family sample file.')
    parse.add_argument('--cache',type=str,default=dir_cache_path_default,help='The path of VEP cache files.')
    parse.add_argument('--gv',type=str,required=True,help='The version of genome.')
    parse.add_argument('--fa',type=str,required=True,help='The reference fasta file.')
    parse.add_argument('--conf',type=str,required=True,help='The configure file of vcfanno.')
    parse.add_argument('--af',type=float,default=0.005,help='The allele frequency of PM2.')
    parse.add_argument('--plp',type=str,required=True,help='The P/LP file.')
    parse.add_argument('--blb',type=str,default=None,help='The B/LB file.')
    parse.add_argument('--incis', type=str, default=None, help='The in cis file.')
    parse.add_argument('--transcript',type=str,default=transcript_path_default,help='ENST IDs file of the transcript of the target gene.')
    parse.add_argument('--ncd',type=str,default='No',choices=['Yes','No'],help='Whether to include non-coding region variants.')
    parse.add_argument('--intrans',type=str,default=None,help='The in trans file.')
    parse.add_argument('--scorePy',type=str,required=True,help='The scoring script.')
    parse.add_argument('--output',type=str,required=True,help='The name of output file.')

    args = parse.parse_args()
    return args

def filt_sample_AR(input_vcf, bed_file, sample_file, paral_num, output_vcf):
    #print(input_vcf)
    print('filt_sample_AR+++++++++++++++++++++++')
    print(bed_file)
    print(sample_file)
    print(output_vcf)
    # 根据ARgene的bed文件和case(+family)样本文件，提取指定样本和AR gene的位点
    # 将多个等位基因变异拆分为单条记录，bcftools norm -m-any
    bcftools_command = f"bcftools view -S {sample_file} --force-samples -R {bed_file} {input_vcf} --threads {paral_num} | bcftools norm -m-any -o {output_vcf}"
    subprocess.run(bcftools_command, shell=True, check=True)
    print('filt_sample_AR over==============')



def vep(input_vcf, output_vcf, cache, paral_num, genome_version, fasta_file):
    print('vep+++++++++++++++++++++++')
    print(input_vcf)
    print(output_vcf)
    print(cache)
    print(paral_num)
    print(genome_version)
    print(fasta_file)
    vep_command = f"vep -i {input_vcf} -o {output_vcf} --dir_cache {cache} --fork {paral_num} --cache --cache_version 113 --offline --assembly {genome_version} --fasta {fasta_file} --hgvs --symbol --biotype  --force_overwrite --vcf --plugin SpliceRegion,Extended"
    subprocess.run(vep_command, shell=True, check=True)
    print('vep over==============')

def split_csq(input_vcf, output_vcf):
    print('split_csq++++++++++++++++++++++++')
    print(input_vcf)
    print(output_vcf)
    with open(output_vcf,'w') as out_f:
        with open(input_vcf) as in_f:
            for line in in_f:
                line = line.strip()
                if (line.startswith("##") and 'INFO=<ID=CSQ,' in line):
                    ####���按照冒号分割INFO=<ID=CSQ行，在用|分割，最后得到csq的所有字段名
                    csq_key = line.split(":")[1].strip().strip(">\"").split('|')
                    for csq_key_each in csq_key:
                        print('csq_key_each===',csq_key_each)
                        # print("##INFO=<ID=VEP_" + csq_key_each + ",Number=.,Type=String,Description=\"NA\">")
                        out_f.write("##INFO=<ID=VEP_" + csq_key_each + ",Number=.,Type=String,Description=\"NA\">" + "\n")
                    #####为新生成的VEP CSQ字段添加header
                    # print(line)
                    out_f.write(line + '\n')
                elif (line.startswith("##")):
                    # print(line)
                    out_f.write(line + '\n')
                elif (line.startswith("#")):
                    #######找到一个#好开头的那一行，得到INFO所在的字段坐标
                    info_start = line.strip().strip("#").split("\t").index("INFO")
                    # print(line)
                    out_f.write(line + '\n')
                else:
                    ######开始处理variant行
                    ######用制表符分割variant行，并取出info信息，用分号分割，保存到variant_value中
                    variant_value = line.strip().split("\t") #每个变异的所有信息
                    # print('variant_value====',variant_value)
                    info_value = line.strip().split("\t")[info_start].split(";")  ##INFO列的信息
                    # print(info_value)
                    for each in range(len(info_value)):
                        ######找到CSQ所在字段
                        # print('info_value[each]====',info_value[each])
                        if (info_value[each].startswith("CSQ=")):

                            # print('info_value[each]=======', info_value[each])
                            csq_value_all = []
                            #######构造一个空的csq_value数组，用于存储csq value
                            for tmp_index in range(len(csq_key)):
                                csq_value_all.append("")
                            ############用逗号分割CSQ的各个value组
                            csq_value = info_value[each].strip("CSQ=").split(',')
                            for index in range(len(csq_value)):
                                ########用|分割每一组CSQ信息，并保存到对应的CSQ INFO字段下面，并用|分割
                                csq_value_each = csq_value[index].split('|')
                                for index_each in range(len(csq_value_each)):
                                    csq_value_all[index_each] += '|' + csq_value_each[index_each]
                            out_vep_info = ""
                            #########按照key=value的形式输出VEP的所有CSQ字段
                            for index in range(len(csq_value_all)):
                                out_vep_info += "VEP_" + csq_key[index] + "=" + csq_value_all[index][1:] + ";"
                            info_value[each] = out_vep_info
                    # out_vep_info = ""
                    # #########按照key=value的形式输出VEP的所有CSQ字段
                    # for index in range(len(csq_value_all)):
                    #     out_vep_info += "VEP_" + csq_key[index] + "=" + csq_value_all[index][1:] + ";"
                    # info_value[each] = out_vep_info
                    #########连接新的csq info value到variant行中
                    variant_value[info_start] = ';'.join(info_value).strip(";")
                    #print("\t".join(variant_value))
                    out_f.write("\t".join(variant_value) + '\n')
    print('split over================')
    print('split over2222================')


def vcfanno(input_vcf, output_vcf, paral_num, conf_file):
    print('vcfanno++++++++++++++++++')
    print(input_vcf)
    print(output_vcf)
    print(paral_num)
    print(conf_file)
    #conf_file需要根据实际情况，（提取gnomAD v4.1 特定基因的vcf文件/指定对应chr的vcf文件，并建索引），作为注释文件
    vcfanno_command = f"vcfanno -p {paral_num} {conf_file} {input_vcf} > {output_vcf}"
    subprocess.run(vcfanno_command, shell=True, check=True)
    print('vcfanno over======================')

def extract(input_vcf,out_tsv):
    print('extract+++++++++++++++++++++++')
    print(input_vcf)
    print(out_tsv)
    extract_command = f'bcftools query -H -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/VEP_SYMBOL\t%INFO/VEP_Feature\t%INFO/VEP_IMPACT\t%INFO/VEP_SpliceRegion\t%INFO/VEP_Consequence\t%INFO/gnomAD_AF_grpmax_joint\t%INFO/VEP_HGVSc\t%INFO/VEP_HGVSp[\t%GT]\n" {input_vcf} > {out_tsv}'
    subprocess.run(extract_command, shell=True, check=True)
    print('extract over=====================')


def filt_var(input_tsv, output_tsv, bed, plp, linkage, af, ncd, transcript):
    #PLP.txt中位点直接保留，未在PLP中的位点，根据AF值（可选参数）过滤。
    #是否考虑non-coding，默认不考虑-会根据VEPIMPCAT VEPconsequence过滤；如果考虑，则不过滤
    print('filtAF++++++++++++++++++++++++')
    print(input_tsv)
    print(output_tsv)
    print(plp)
    print(af)
    print(ncd)
    tmp_check_tsv = 'tmp_filtcheck.tsv'
    tmp_total_tsv = 'tmp_filttotal.tsv'
    # read gene_symbol 对应到VEP_Symbol
    with open(bed, 'r') as f:
        gene_symbol = []
        for line in f:
            head = line.strip().split('\t')[3]
            # print('symbol======',head)
            gene_symbol.append(head)
            # print(gene_symbol)

    with open(plp, 'r') as f:
        linenum = 0
        snp_p = []
        for line in f:
            linenum += 1
            if (linenum == 1):
                continue
            else:
                head = line.strip().split('\t')[0:4]
                snp_p.append('_'.join(head))

    with open(transcript, 'r') as f:
        transcripts = []
        for line in f:
            transcripts.append(line.strip().split('\t')[0])

    # linkage位点记录
    with open(linkage, 'r') as f:
        lineNum = 0
        link = {}
        for line in f:
            lineNum += 1
            if (lineNum == 1):
                continue
            else:
                head = line.strip().split('\t')
                snp = '_'.join(head[1:5])
                if snp not in link:
                    link[snp] = [head[0]]
                else:
                    link[snp].append(head[0])

    with (open(input_tsv,'r') as f):
        linenum = 0
        header = []
        out = {}
        for line in f:
            linenum += 1
            if (linenum == 1):
                temp = line.strip().split('\t')
                header = [i.split("]")[1].split(":")[0] for i in temp]
                with open(output_tsv, 'w') as f_out:
                    # f_out.write('\t'.join(temp[0:7]) + '\t' + temp[9] + '\t' + '\t'.join(header[10:]) + '\n')
                    f_out.write(line)
                with open(tmp_check_tsv,'w') as f:
                    f.write(line)
                with open(tmp_total_tsv,'w') as f:
                    f.write(line)

            else:
                sampleIDs = header[12:]
                tmp = line.strip().split('\t')
                # genotypes=[i.split(";")[0] for i in tmp[17:]]
                genotypes = tmp[12:]
                snp = "_".join(tmp[0:4])
                ## 标记连锁位点
                if snp in link.keys():
                    # print(snp)
                    for i in range(len(sampleIDs)):
                        sid = sampleIDs[i]
                        if sid in link[snp]:
                            # print(genotypes[i])
                            # print(sid)
                            genotypes[i] = 'rm'
                ## 拆 Feature,判断基因
                VEP_SYMBOL = tmp[4].split('|')
                VEP_Feature = tmp[5].split('|')
                VEP_IMPACT = tmp[6].split('|')
                # VEP_MANE_SELECT = tmp[7].split('|')
                # VEP_MANE_PLUS_CLINICAL = tmp[8].split('|')
                # VEP_Gene=tmp[7].split('|')
                VEP_SpliceRegion = tmp[7].split('|')
                VEP_Consequence = tmp[8].split('|')
                VEP_AF = tmp[9]
                VEP_HGVSc = tmp[10].split('|')
                VEP_HGVSp = tmp[11].split('|')

                if snp in snp_p:
                    for i in range(len(VEP_Feature)):
                        key = '\t'.join(tmp[0:4]) + '\t' + VEP_SYMBOL[i] + '\t' + VEP_Feature[i] + '\t' + \
                              VEP_IMPACT[i] + '\t' + VEP_SpliceRegion[i] + '\t' + VEP_Consequence[i] +  '\t' + VEP_AF + '\t' + VEP_HGVSc[i] + '\t' + VEP_HGVSp[i] + '\t' + '\t'.join(genotypes)
                        with open(tmp_total_tsv,'a') as f: #check，写入总文件
                            f.write(key + '\n')
                        if VEP_SYMBOL[i] in gene_symbol and  VEP_Feature[i].split('.')[0] in transcripts: #MANE转录本过滤, gene symbol过滤
                                # print(snp)
                                # print('\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' +'\t'.join(genotypes))
                                # with open(output, "a") as f:
                            with open(output_tsv, 'a') as f_out:
                                f_out.write(key + '\n')
                        else:  #check,因为symbol或transcript被过滤的位点
                            with open(tmp_check_tsv,'a') as f:
                                f.write(key + '\n')

                        # out[key]=VEP_IMPACT[i] + '\t'+ VEP_SYMBOL[i] + '\t'+VEP_Gene[i] + '\t'+'\t'.join(genotypes)

                        # f.write('\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' +VEP_IMPACT[i] + '\t' +VEP_SYMBOL[i] + '\t'+VEP_Gene[i] + '\t'+'\t'.join(genotypes)+'\n')
                else:
                    for i in range(len(VEP_Feature)):
                        key = '\t'.join(tmp[0:4]) + '\t' + VEP_SYMBOL[i] + '\t' + VEP_Feature[i] + '\t' + \
                              VEP_IMPACT[i] + '\t' + VEP_SpliceRegion[i] + '\t' + VEP_Consequence[i] + '\t' + VEP_AF  + '\t' + VEP_HGVSc[i] + '\t' + VEP_HGVSp[i] + '\t' + '\t'.join(genotypes)
                        with open(tmp_total_tsv,'a') as f: #check，写入总文件
                            f.write(key + '\n')
                    gnomad_AF = tmp[9]  ### 判断频率
                    filter_flag = 1
                    if (gnomad_AF == '.'):
                        filter_flag = 1
                    elif float(gnomad_AF) > af:  # 频率＞0.005过滤掉
                        filter_flag = 0
                    if filter_flag == 1:
                        for i in range(len(VEP_Feature)):
                            key = '\t'.join(tmp[0:4]) + '\t' + VEP_SYMBOL[i] + '\t' + VEP_Feature[i] + '\t' + \
                                  VEP_IMPACT[i] + '\t' + VEP_SpliceRegion[i] + '\t' + VEP_Consequence[i] + '\t' + VEP_AF + '\t' + VEP_HGVSc[i] + '\t' + VEP_HGVSp[i] + '\t' + '\t'.join(genotypes)
                            if (VEP_Feature[i].split('.')[0] in transcripts and VEP_SYMBOL[i] in gene_symbol): #MANE转录本过滤, gene symbol过滤
                                if ncd == 'Yes': #如果考虑non-coding,则不过滤位点
                                    with open(output_tsv, 'a') as f_out:
                                        f_out.write(key + '\n')
                                else:   ##纳入所有CDS区并考虑UTR
                                    if VEP_IMPACT[i] == 'MODERATE' or VEP_IMPACT[i] == 'HIGH' \
                                            or re.findall(r"extended_intronic_splice_region_variant", VEP_SpliceRegion[i]) \
                                            or re.findall(r"splice_region_variant", VEP_Consequence[i]) \
                                            or re.findall(r"synonymous_variant",VEP_Consequence[i]) \
                                            or re.findall(r"start_retained_variant",VEP_Consequence[i]) \
                                            or re.findall(r"stop_retained_variant",VEP_Consequence[i]) \
                                            or re.findall(r"5_prime_UTR_variant", VEP_Consequence[i]) \
                                            or re.findall(r"3_prime_UTR_variant", VEP_Consequence[i]) :
                                        # and VEP_IMPACT[i]=='MODERATE' or VEP_IMPACT[i]=='HIGH':
                                        with open(output_tsv, 'a') as f_out:
                                            f_out.write(key + '\n')
                                    else: #check，因为non-conding被过滤的位点
                                        with open(tmp_check_tsv, 'a') as f: #记录被过滤的位点
                                            f.write(key + '\n')
                            else: #check，因为转录本或symbol被过滤的位点
                                with open(tmp_check_tsv,'a') as f:
                                    f.write(key + '\n')
                    else: #check，因为AF被过滤的位点
                        for i in range(len(VEP_Feature)):
                            key = '\t'.join(tmp[0:4]) + '\t' + VEP_SYMBOL[i] + '\t' + VEP_Feature[i] + '\t' + \
                                  VEP_IMPACT[i] + '\t' + VEP_SpliceRegion[i] + '\t' + VEP_Consequence[i] + '\t' + VEP_AF + '\t' + VEP_HGVSc[i] + '\t' + VEP_HGVSp[i] + '\t' + '\t'.join(genotypes)
                            with open(tmp_check_tsv,'a') as f:
                                f.write(key + '\n')

    print('filtAF over========================')


def uniq_transcripts(input_tsv, output_txt):
    uniq_command = f"cut -f6 {input_tsv} |sort -u > {output_txt}"
    subprocess.run(uniq_command, shell=True, check=True)

def call_scoring(args):
    uniq_transc, score_py, network, plp, blb, output, intrans = args
    command = ['python',score_py, uniq_transc, network, plp, blb, output]
    print('intrans===',intrans)
    if intrans:
        command.append(intrans)
    subprocess.run(command)
    # call_command = f"cat {uniq_transc} | xargs -P 20 -n 1 python {score_py} {network_tsv} {plp} {output} {intrans}"
    # subprocess.run(call_command,shell=True,check=True)

if __name__ == '__main__':
    paral_num = 10
    dir_cache_default = '~/.vep'
    transcript_file_default = '~/PM3/database/MANE_ENST_NM.txt'
    tmp_filt_out = 'tmp_filt.vcf'
    tmp_vep_out = 'tmp_vep.vcf'
    tmp_splitcsq_out = 'tmp_splitcsq.vcf'
    tmp_vcfanno_out = 'tmp_vcfanno.vcf'
    tmp_extract_out = 'tmp_extract.tsv'
    tmp_filtvar_out = 'tmp_filtvar.tsv'
    tmp_uniqtranscripts_out = 'tmp_uniqtranscripts.txt'

    args = parse_args(dir_cache_default,transcript_file_default)

    filt_sample_AR(args.vcf, args.bed, args.sam, paral_num, tmp_filt_out)
    vep(tmp_filt_out, tmp_vep_out, args.cache, paral_num, args.gv, args.fa)
    split_csq(tmp_vep_out, tmp_splitcsq_out)
    vcfanno(tmp_splitcsq_out, tmp_vcfanno_out, paral_num, args.conf)
    extract(tmp_vcfanno_out,tmp_extract_out)
    filt_var(tmp_extract_out, tmp_filtvar_out, args.bed, args.plp, args.incis, args.af, args.ncd, args.transcript)
    uniq_transcripts(tmp_filtvar_out,tmp_uniqtranscripts_out)

    with open(tmp_uniqtranscripts_out,'r') as f:
        uniq_transcs = [line.strip() for line in f.readlines()]
    print('uniq_transc',uniq_transcs)

    # args_list = ['scoring.py','tmp_filtvar.tsv', 'PLP_chr1_9gene.txt','output.tsv','intrans.txt']
    # args_list = [(uniq_transc,'scoring.py','tmp_filtvar.tsv', 'PLP_chr1_9gene.txt','output.tsv') for uniq_transc in uniq_transcs]
    args_list = [(uniq_transc, args.scorePy, tmp_filtvar_out, args.plp, args.blb, args.output, args.intrans) for uniq_transc in uniq_transcs]
    with Pool(processes = paral_num) as pool:
        pool.map(call_scoring,args_list)

    # 合并结果文件
    with open(args.output, 'w') as outfile:
        outfile.write('CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'Point_Hom' + '\t' + \
                'Symbol' + '\t' + 'Transcript' + '\t' + 'VEP_IMPACT' + '\t' + 'VEP_Consequence' + '\t' + \
                'AF_Max' + '\t' + 'VEP_HGVSc' + '\t' + 'VEP_HGVSp' + '\t' + 'HomozygousNum' + '\t' + 'Score' + '\t' + 'PM3_Level' + '\n')

        for sub_file in [f for f in os.listdir() if f.endswith('_'+ args.output)]:
            with open(sub_file) as infile:
                outfile.write(infile.read())
            os.remove(sub_file)

    with open('allFreqSmps.txt', 'w') as outfile:
        outfile.write(
            'PV' + '\t' + 'CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'Point_Hom' + '\t' + 'Transcript' + '\t' + 'Var2' + '\t' + 'Sample' + '\t' + 'Weight' + '\n')
        for sub_file in [f for f in os.listdir() if f.endswith('_allFreqSmps.txt')]:
            with open(sub_file) as infile:
                outfile.write(infile.read())
            os.remove(sub_file)

    with open('multiVars.tsv', 'w') as outfile:
        outfile.write(
            'SampleID' + '\t' + 'CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'Point_Hom' + '\t' + 'Transcript' + '\n')
        for sub_file in [f for f in os.listdir() if f.endswith('_multiVars.tsv')]:
            with open(sub_file) as infile:
                outfile.write(infile.read())
            os.remove(sub_file)