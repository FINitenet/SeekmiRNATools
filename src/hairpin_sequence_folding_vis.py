#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :hairpin_sequence_folding_vis.py
# @Time      :2024/12/01 16:50:00
# @Author    :Yuchen@rlab

import pandas as pd
import subprocess
import os

# GFF 文件路径
gff_file = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_merge/converted.gff"
gff_data = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                       names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

# 选择感兴趣的行，例如 gene 类型
gene_data = gff_data[gff_data['type'] == 'siRNA']

# BAM 和基因组文件路径
bam_file = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/R_0/merged.bam"
genome_file = "/bios-store1/chenyc/Reference_Source/Artemisia_Reference/Aannua.fa"

# 输出文件夹路径
output_folder = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_folding/wildtype"
os.makedirs(output_folder, exist_ok=True)

# 遍历每个基因，生成并运行 `strucVis` 命令
for _, row in gene_data.iterrows():
    seqid = row['seqid']
    start = row['start']
    end = row['end']
    strand = row['strand']
    
    # 从 attributes 中提取 Locus 名称
    attributes = row['attributes']
    locus_name = None
    for attr in attributes.split(";"):
        if attr.startswith("ID="):
            locus_name = attr.split("=")[1]
            break
    
    if locus_name:
        # 构建输出 PS 文件路径
        output_ps = os.path.join(output_folder, f"{locus_name}_{strand}.ps")
        
        # 构建 strucVis 命令，生成 PS 文件
        command = [
            "strucVis",
            "-b", bam_file,
            "-g", genome_file,
            "-c", f"{seqid}:{start}-{end}",
            "-s", strand,
            "-p", output_ps,
            "-n", locus_name
        ]
        subprocess.run(command)
        
        # 转换 PS 文件为 PDF 文件
        output_pdf = os.path.join(output_folder, f"{locus_name}_{strand}.pdf")
        convert_command = ["ps2pdf", output_ps, output_pdf]
        subprocess.run(convert_command)
        
        # 可选择删除原始 PS 文件
        os.remove(output_ps)

print("strucVis commands executed and PS files converted to PDF successfully.")