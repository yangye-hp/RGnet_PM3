#!/usr/bin/python
#coding:utf-8
import sys
import argparse
import gzip
out_str=""
input_file= sys.argv[1]
with open(input_file) as in_f:
  for line in in_f:
    line=line.strip()
######INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE"
    if(line.startswith("##") and 'INFO=<ID=CSQ,' in line ):
####���按照冒号分割INFO=<ID=CSQ行，在用|分割，最后得到csq的所有字段名
      csq_key=line.split(":")[1].strip().strip(">\"").split('|')
      for csq_key_each in csq_key:
        print("##INFO=<ID=VEP_"+csq_key_each+",Number=.,Type=String,Description=\"NA\">")
#####为新生成的VEP CSQ字段添加header
      print (line)
    elif(line.startswith("##")):
      print (line)
    elif(line.startswith("#")):
#######找到一个#好开头的那一行，得到INFO所在的字段坐标
      info_start=line.strip().strip("#").split("\t").index("INFO")
      print (line)
    else:
######开始处理variant行
######用制表符分割variant行，并取出info信息，用分号分割，保存到variant_value中
      variant_value=line.strip().split("\t")
      info_value=line.strip().split("\t")[info_start].split(";")
      for each in range(len(info_value)):
######找到CSQ所在字段
        if(info_value[each].startswith("CSQ=")):
          csq_value_all=[]
#######构造一个空的csq_value数组，用于存储csq value
          for tmp_index in range(len(csq_key)):
            csq_value_all.append("")
############用逗号分割CSQ的各个value组
          csq_value=info_value[each].strip("CSQ=").split(',')
          for index in range(len(csq_value)):
########用|分割每一组CSQ信息，并保存到对应的CSQ INFO字段下面，并用|分割
            csq_value_each=csq_value[index].split('|')
            for index_each in range(len(csq_value_each)):
                csq_value_all[index_each]+='|'+csq_value_each[index_each]
      out_vep_info=""
#########按照key=value的形式输出VEP的所有CSQ字段
      for index in range(len(csq_value_all)):
        out_vep_info+="VEP_"+csq_key[index]+"="+csq_value_all[index][1:]+";"
      info_value[each]=out_vep_info
#########连接新的csq info value到variant行中
      variant_value[info_start]=';'.join(info_value).strip(";")
      print ("\t".join(variant_value))

