#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :mirna_annotation_update.py
# @Time      :2024/12/02 00:51:19
# @Author    :Yuchen@rlab
# @Version   :20241202
# @Description: Update miRNA annotation information based on NGS data
# @Update    :20241202

from Bio import SeqIO
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def Long_to_Short(long_sequence_file, short_sequences_file):
    # 读取长序列文件（多个序列）
    long_sequences = {}
    for record in SeqIO.parse(long_sequence_file, "fasta"):
        long_sequences[record.id] = str(record.seq)
    print(f"Read {len(long_sequences)} long sequences from {long_sequence_file}")
    # 打印序列数量，确保读取成功
    print(f"Read {len(long_sequences)} long sequences")

    # 读取短序列文件
    short_sequences = pd.read_csv(short_sequences_file)

    # 确保短序列文件包含必要的列
    required_columns = {
        "MIR", "Peak1_most_common_read", "Peak1_RPM",
        "Peak2_most_common_read", "Peak2_RPM"
    }
    if not required_columns.issubset(short_sequences.columns):
        missing_columns = required_columns - set(short_sequences.columns)
        raise ValueError(f"The short sequences file is missing the following columns: {', '.join(missing_columns)}")

    # 查找短序列在长序列中的位置，并标记mature和star
    positions = []
    for _, row in short_sequences.iterrows():
        mir_id = row['MIR']
        peak1_seq = row['Peak1_most_common_read']
        peak2_seq = row['Peak2_most_common_read']
        peak1_rpm = row['Peak1_RPM']
        peak2_rpm = row['Peak2_RPM']

        # 验证短序列是否为有效字符串
        peak1_valid = isinstance(peak1_seq, str) and len(peak1_seq) > 0
        peak2_valid = isinstance(peak2_seq, str) and len(peak2_seq) > 0
        if not peak1_valid and not peak2_valid:
            print(f"Warning: Both short sequences for {mir_id} in file {short_sequences_file} are empty or invalid, skipping processing. Context: MIR={mir_id}, Peak1={peak1_seq}, Peak2={peak2_seq}")
            continue

        # 获取当前长序列
        long_sequence = long_sequences[mir_id]

        # 只有 Peak1 有效
        if peak1_valid and not peak2_valid:
            mature_seq, mature_label = peak1_seq, f"{mir_id}_mature"
            mature_start = long_sequence.find(mature_seq)
            if mature_start != -1:
                mature_end = mature_start + len(mature_seq) - 1
                positions.append((mir_id, mature_seq, mature_start, mature_end, mature_label))

        # 只有 Peak2 有效
        elif peak2_valid and not peak1_valid:
            mature_seq, mature_label = peak2_seq, f"{mir_id}_mature"
            mature_start = long_sequence.find(mature_seq)
            if mature_start != -1:
                mature_end = mature_start + len(mature_seq) - 1
                positions.append((mir_id, mature_seq, mature_start, mature_end, mature_label))

        # 如果两个序列都有效
        elif peak1_valid and peak2_valid:
            if peak1_rpm >= peak2_rpm:
                mature_seq, star_seq = peak1_seq, peak2_seq
                mature_label, star_label = f"{mir_id}_mature", f"{mir_id}_star"
            else:
                mature_seq, star_seq = peak2_seq, peak1_seq
                mature_label, star_label = f"{mir_id}_mature", f"{mir_id}_star"

            # 定位mature序列
            mature_start = long_sequence.find(mature_seq)
            if mature_start != -1:
                mature_end = mature_start + len(mature_seq) - 1
                positions.append((mir_id, mature_seq, mature_start, mature_end, mature_label))

            # 定位star序列
            star_start = long_sequence.find(star_seq)
            if star_start != -1:
                star_end = star_start + len(star_seq) - 1
                positions.append((mir_id, star_seq, star_start, star_end, star_label))

    # 转为DataFrame
    positions_df = pd.DataFrame(positions, columns=['ID', 'Sequence', 'Start', 'End', 'Label'])
    return positions_df

def Update_GFF(gff_file, positions_df, output_gff):
    # 读取GFF文件
    gff_columns = ['Seqname', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute']
    gff_data = pd.read_csv(gff_file, sep="\t", header=None, names=gff_columns, comment='#')

    # 读取positions_df
    positions_df = positions_df

    # 创建新的子类注释
    new_annotations = []

    for _, pos_row in positions_df.iterrows():
        # 获取父类信息
        logging.info(f"Processing row: {pos_row}")
        parent_id = pos_row['ID']
        parent_rows = gff_data[gff_data['Attribute'].str.contains(f"ID={parent_id}(?:;|$)", regex=True)]
        
        if parent_rows.empty:
            print(f"Warning: No annotation found for parent {parent_id}")
            continue

        for _, parent_row in parent_rows.iterrows():
            # 计算子类的起始和终止位点
            parent_start = parent_row['Start']
            parent_end = parent_row['End']
            child_start = parent_start + pos_row['Start']
            child_end = parent_start + pos_row['End']

            # 构造子类的Attribute字段
            child_id = pos_row['Label']
            child_attribute = f"ID={child_id};Name={pos_row['Label']};Derives_from={parent_id}"

            # 添加子类注释
            new_annotations.append([
                parent_row['Seqname'],
                parent_row['Source'],
                "miRNA",  # 子类特征名
                child_start,
                child_end,
                ".",  # 默认Score为空
                parent_row['Strand'],
                ".",  # 默认Frame为空
                child_attribute
            ])

    # 转换新注释为DataFrame
    new_annotations_df = pd.DataFrame(new_annotations, columns=gff_columns)

    # 合并原始GFF数据和新注释
    updated_gff = pd.concat([gff_data, new_annotations_df], ignore_index=True)
    gff_data_sorted = updated_gff.sort_values(by=['Seqname', 'Start', 'End'], ascending=[True, True, True])


    # 保存更新后的GFF文件
    # Write the DataFrame to CSV in chunks to avoid memory issues
    with open(output_gff, 'w') as f:
        gff_data_sorted.to_csv(f, sep="\t", index=False, header=False, chunksize=10000)

    print(f"The updated GFF file has been saved to {output_gff}")

if __name__ == "__main__":

    long_sequence_file = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_extract/loci.u.fa"  # 长序列FASTA文件
    short_sequences_file = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/R_0/miRNA_processing_results.csv"  # 短序列CSV文件
    annotation_file = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_merge/converted.gff"  # 原注释文件
    output_file = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_merge/updated_genome_annotation.gff"  # 输出文件

    gff_file = "/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_merge/converted.gff"  # 已有的GFF文件
    output_gff = "/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_merge/updated_genome_annotation.gff"  # 输出的GFF文件

    positions_df = Long_to_Short(long_sequence_file, short_sequences_file)
    Update_GFF(gff_file, positions_df, output_gff)  # 更新GFF文件