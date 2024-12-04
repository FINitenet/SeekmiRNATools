#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :utils.py
# @Time      :2024/12/03 20:03:12
# @Author    :Yuchen@rlab

import gffutils
import datetime
import subprocess
from Bio import SeqIO

def print_current_time(message=""):
    """
    Print the current time with an optional message
    """
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if message:
        print(f"{message} - Current time: {current_time}")
    else:
        print(f"Current time: {current_time}")

def rnafold_predict(sequence):
    # 使用 RNAfold 通过 subprocess 调用
    process = subprocess.Popen(
        ['RNAfold'], 
        stdin=subprocess.PIPE, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        text=True
    )
    
    # 传递序列给 RNAfold 并获取输出
    stdout, stderr = process.communicate(input=sequence)

    # 分析 RNAfold 输出，获取 ΔG (自由能)
    # RNAfold 的输出通常是两行，第二行为结构和 ΔG
    lines = stdout.strip().split('\n')
    if len(lines) > 1:
        structure_info = lines[1].split()[-1]  # 结构末尾有 ΔG
        structure = lines[1].split()[0]
        delta_g = float(structure_info.strip('()'))  # 去掉括号并转为浮点数
    else:
        structure = None
        delta_g = None

    return delta_g, structure

def get_parent_sequences_from_fasta(fasta_file, parent_names):
    """
    从FASTA文件中获取父类序列

    参数:
    fasta_file: str, FASTA文件的路径
    parent_names: list, 包含父类的名称或ID列表

    返回:
    parent_sequences: dict, 父类名称或ID为键，对应序列为值的字典
    """
    # 初始化一个空字典存储父类的序列
    parent_sequences = {}

    # 使用 SeqIO 解析 FASTA 文件
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # 如果该序列的ID或Name在父类列表中
            if record.id in parent_names:
                # 将该序列添加到字典中，父类ID作为键，序列作为值
                parent_sequences = str(record.seq)

    return parent_sequences


def parse_gff_for_mirna(db, flag):
    # 初始化空字典存储结果
    result_dict = {}

    for feature in db.features_of_type("miRNA_primary_transcript"):
        parent_id = feature.attributes['ID'][0]

        # 获取父特征对象
        parent_feature = db[parent_id]
        
        # 获取父特征的 Name 属性、位置信息和链信息
        parent_name = parent_feature.attributes.get('Name', [parent_feature.id])[0]
        parent_start = parent_feature.start
        parent_end = parent_feature.end
        parent_strand = parent_feature.strand
        # 初始化字典条目，如果父类的 Name 还未添加到字典中，则创建新的条目
        if parent_name not in result_dict:
            result_dict[parent_name] = {
                "Parent_contig": parent_feature.seqid,
                "Parent_start": parent_start,
                "Parent_end": parent_end,
                "Parent_strand": parent_strand,  # 添加父类的链信息
                "Children": []  # 创建一个空列表来存储子类信息
            }
    
    if flag:
        print("The miRNA_primary_transcript annotation information has been successfully parsed.")
    else:
        for feature in db.features_of_type("miRNA"):
            parent_name = feature.attributes['Derives_from'][0]
            # 获取子特征的 Name 属性、位置信息和链信息
            child_name = feature.attributes.get('Name', [feature.id])[0]
            child_contig = feature.seqid
            child_start = feature.start
            child_end = feature.end
            result_dict[parent_name]["Children"].append({
                            "Child_name": child_name,
                            "Child_contig": child_contig,
                            "Child_start": child_start,
                            "Child_end": child_end
                        })

    return result_dict