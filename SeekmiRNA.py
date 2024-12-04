#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :SeekmiRNA.py
# @Time      :2024/12/04 13:50:06
# @Author    :Yuchen@rlab
# @Description: This script is used to predict miRNA from GFF file and calculate the peak and strand specificity of miRNA from BAM file.
# @Software  :PyCharm
# @Environment:python3


import os
import pysam
import subprocess
import gffutils
import logging
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict, Counter
from src.utils import get_parent_sequences_from_fasta, rnafold_predict, print_current_time, parse_gff_for_mirna

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("SeekmiRNA.log"),
        logging.StreamHandler()
    ]
)

def extract_fasta_from_gff(hairpin_sequences, GENOME_FILE, GFF_FILE):
        command = [
            "gffread",
            "-g", GENOME_FILE,
            "-w", hairpin_sequences,
            GFF_FILE
        ]

        subprocess.run(command, check=True)

        logging.info(f"FASTA file has been successfully generated: {hairpin_sequences}")

def calculate_peaks_and_rnafold(hairpin_sequences, hairpin_dict, bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # 获取 BAM 文件中的总reads数，用于RPM计算
        bam_total_reads = bam.count()
        
        # 初始化用于保存结果的列表
        peak_data = []

        # 遍历字典中的父类信息
        for parent_name, parent_info in hairpin_dict.items():
            parent_contig = parent_info['Parent_contig']
            parent_start = parent_info['Parent_start']
            parent_end = parent_info['Parent_end']

            # 获取父类序列并计算 RNAfold 预测结果
            parent_sequence = get_parent_sequences_from_fasta(hairpin_sequences, parent_name)
            if parent_sequence:
                delta_g, structure = rnafold_predict(parent_sequence)
            else:
                structure = None
                delta_g = None

            # 使用 defaultdict 统计父类区域中每个位置的reads数
            position_counts = defaultdict(int)
            read_sequences_by_position = defaultdict(list)

            # 遍历父类区域内的所有 reads，并统计每个位置的reads数和序列
            for read in bam.fetch(region=None, start=parent_start, end=parent_end, contig=parent_contig):
                read_start = read.reference_start
                seq = read.get_forward_sequence()

                read_sequences_by_position[read_start].append(seq)
                position_counts[read_start] += 1

            peaks = []
            current_peak_count = 0
            current_peak_start = None
            current_peak_reads = []

            sorted_positions = sorted(position_counts.keys())
            for i, pos in enumerate(sorted_positions):
                if current_peak_start is None:
                    current_peak_start = pos

                if i > 0 and pos > sorted_positions[i - 1] + 24:
                    peaks.append((current_peak_start, sorted_positions[i - 1], current_peak_count, current_peak_reads))
                    current_peak_start = pos
                    current_peak_count = 0
                    current_peak_reads = []

                current_peak_count += position_counts[pos]
                current_peak_reads.extend(read_sequences_by_position[pos])

            # 记录最后一个peak
            if current_peak_start is not None:
                peaks.append((current_peak_start, sorted_positions[-1], current_peak_count, current_peak_reads))

            # 如果没有找到 peaks，则跳过输出
            if not peaks:
                peak_data.append({
                    "MIR": parent_name,
                    "Peak1_most_common_read": None,
                    "Peak1_most_common_length": None,
                    "Peak1_RPM": None,
                    "Peak2_most_common_read": None,
                    "Peak2_most_common_length": None,
                    "Peak2_RPM": None,
                    "distance": None,
                    "deltaG": delta_g,
                    "structure": structure
                })
                continue
            
            # 按照peak reads数进行排序，找出丰度最高的两个peaks
            peaks = sorted(peaks, key=lambda x: x[2], reverse=True)

            if len(peaks) > 1:
                # 选择丰度最高的两个peaks
                peak1 = peaks[0]
                peak2 = peaks[1]

                # 计算RPM
                peak1_rpm = (peak1[2] / bam_total_reads) * 1e6
                peak2_rpm = (peak2[2] / bam_total_reads) * 1e6

                # 计算两个peak的物理距离（即第一个peak的结束位置和第二个peak的开始位置的差值）
                physical_distance = peak2[0] - peak1[1]

                # 找出每个peak中最主要的reads序列及其长度
                peak1_reads = peak1[3]
                peak2_reads = peak2[3]

                # 计算最常见的reads序列
                peak1_most_common_read, peak1_most_common_count = Counter(peak1_reads).most_common(1)[0]
                peak2_most_common_read, peak2_most_common_count = Counter(peak2_reads).most_common(1)[0]

                # 获取最常见reads序列的长度
                peak1_most_common_length = len(peak1_most_common_read)
                peak2_most_common_length = len(peak2_most_common_read)

            elif len(peaks) == 1:
                peak1 = peaks[0]
                peak2 = None

                # 计算RPM
                peak1_rpm = (peak1[2] / bam_total_reads) * 1e6
                peak2_rpm = None

                # 计算两个peak的物理距离（即第一个peak的结束位置和第二个peak的开始位置的差值）
                physical_distance = None

                # 找出每个peak中最主要的reads序列及其长度
                peak1_reads = peak1[3]
                peak2_reads = None

                # 计算最常见的reads序列
                peak1_most_common_read, peak1_most_common_count = Counter(peak1_reads).most_common(1)[0]
                peak2_most_common_read, peak2_most_common_count = None, None

                # 获取最常见reads序列的长度
                peak1_most_common_length = len(peak1_most_common_read)
                peak2_most_common_length = None
            
            # 保存结果
            peak_data.append({
                "MIR": parent_name,
                "Peak1_most_common_read": peak1_most_common_read,
                "Peak1_most_common_length": peak1_most_common_length,
                "Peak1_RPM": peak1_rpm,
                "Peak2_most_common_read": peak2_most_common_read,
                "Peak2_most_common_length": peak2_most_common_length,
                "Peak2_RPM": peak2_rpm,
                "distance": physical_distance,
                "deltaG": delta_g,
                "structure": structure
            })

        df_peaks = pd.DataFrame(peak_data)

        return df_peaks
    

def calculate_strand_specificity(mirna_dict, bam_file):
    with pysam.AlignmentFile(bam_file, "rb", threads=12) as bam:
        # 初始化用于保存链特异性结果的列表
        strand_specificity_data = []

        # 遍历字典中的父类和子类信息
        for parent_name, parent_info in tqdm(mirna_dict.items()):
            parent_contig = parent_info['Parent_contig']
            parent_start = parent_info['Parent_start']
            parent_end = parent_info['Parent_end']

            # 计数变量
            forward_reads = 0  # 正链 reads 计数
            reverse_reads = 0  # 反链 reads 计数
            total_reads_count = 0
    
            # 遍历父类区域内的所有 reads
            for read in bam.fetch(region=None, start=parent_start, end=parent_end, contig=parent_contig):
                total_reads_count += 1
                if read.is_reverse:
                    reverse_reads += 1
                else:
                    forward_reads += 1

            if reverse_reads > forward_reads:
                parent_strand = '-'
                same_strand_count = reverse_reads
            else:
                parent_strand = '+'
                same_strand_count = forward_reads

            # 避免除以 0 的情况
            if total_reads_count > 0:
                strand_specificity = same_strand_count / total_reads_count
            else:
                strand_specificity = 0  # 如果没有 reads，链特异性为 0

            # 将结果保存到 strand_specificity_data 列表
            strand_specificity_data.append({
                "MIR": parent_name,
                "Strand": parent_strand,
                "Strand_specificity": strand_specificity
            })

    # 将结果转换为 DataFrame
    df_strand_specificity = pd.DataFrame(strand_specificity_data)

    return df_strand_specificity
    


def Long_to_Short(long_sequence_file, short_sequences_file, df_peaks):
    # 读取长序列文件（多个序列）
    long_sequences = {}
    for record in SeqIO.parse(long_sequence_file, "fasta"):
        long_sequences[record.id] = str(record.seq)
    logging.info(f"Read {len(long_sequences)} long sequences from {long_sequence_file}")

    # 读取短序列文件
    short_sequences = df_peaks

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
    for _, row in tqdm(short_sequences.iterrows()):
        mir_id = row['MIR']
        peak1_seq = row['Peak1_most_common_read']
        peak2_seq = row['Peak2_most_common_read']
        peak1_rpm = row['Peak1_RPM']
        peak2_rpm = row['Peak2_RPM']

        # 验证短序列是否为有效字符串
        peak1_valid = isinstance(peak1_seq, str) and len(peak1_seq) > 0
        peak2_valid = isinstance(peak2_seq, str) and len(peak2_seq) > 0
        if not peak1_valid and not peak2_valid:
            logging.info(f"Warning: Both short sequences for {mir_id} in file {short_sequences_file} are empty or invalid, skipping processing. Context: MIR={mir_id}, Peak1={peak1_seq}, Peak2={peak2_seq}")
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
    df_positions = pd.DataFrame(positions, columns=['ID', 'Sequence', 'Start', 'End', 'Label'])
    return df_positions

def Update_GFF(gff_file, df_positions, output_gff):
    # 读取GFF文件
    # 读取 GFF 文件
    gff_columns = ['Seqname', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute']
    gff_data = pd.read_csv(gff_file, sep="\t", header=None, names=gff_columns, comment='#')

    for _, pos_row in tqdm(df_positions.iterrows()):
        # 获取父类 ID
        parent_id = pos_row['ID']
        new_strand = pos_row['Strand']  # 新的 Strand 信息

        # 查找父类记录
        parent_rows = gff_data[gff_data['Attribute'].str.contains(f"ID={parent_id}(?:;|$)", regex=True)]
        if parent_rows.empty:
            logging.info(f"警告：未找到父类 {parent_id} 的注释")
            continue

        # 修改父类 Strand 信息
        gff_data.loc[parent_rows.index, 'Strand'] = new_strand

        # 添加子类注释
        for _, parent_row in parent_rows.iterrows():
            parent_start = parent_row['Start']
            parent_end = parent_row['End']
            parent_strand = new_strand

            # 计算子类的起始和终止位置
            if parent_strand == "+":
                child_start = parent_start + pos_row['Start']
                child_end = parent_start + pos_row['End']
            elif parent_strand == "-":
                child_start = parent_end - pos_row['End']
                child_end = parent_end - pos_row['Start']
            else:
                logging.info(f"Unknown strand information {parent_strand}, skipping {parent_id}")
                continue
    
            # 构造子类的 Attribute 字段
            child_id = pos_row['Label']
            child_attribute = f"ID={child_id};Name={pos_row['Label']};Derives_from={parent_id}"

            # 构造新的一行
            new_row = {
                'Seqname': parent_row['Seqname'],
                'Source': parent_row['Source'],
                'Feature': "miRNA",
                'Start': child_start,
                'End': child_end,
                'Score': ".",
                'Strand': new_strand,  # 子类继承更新后的 Strand
                'Frame': ".",
                'Attribute': child_attribute
            }
            # 插入新行
            gff_data = pd.concat([gff_data, pd.DataFrame([new_row])], ignore_index=True)

    # 按照 Seqname、Start 和 End 排序
    gff_data_sorted = gff_data.sort_values(by=['Seqname', 'Start', 'End'], ascending=[True, True, True])
    # 保存更新后的 GFF 文件
    gff_data_sorted.to_csv(output_gff, sep="\t", index=False, header=False)
    logging.info(f"Updated GFF file has been saved to {output_gff}")

def compute_final_precision_for_mirna(mirna_dict, bam_file):
    """
    •	父类和子类位置信息：从之前生成的字典中提取父类的 start 和 end 位置，以及子类的 start 和 end 位置。
    •	使用 pysam 读取 BAM 文件：通过 bam.count(region=None, start=start, end=end, contig=contig) 方法，计算指定区域内的 reads 数。
    •	子类区域±2：子类的 start 和 end 各向外扩展 2 个碱基（child_start - 2 和 child_end + 2）。
    •	计算精度：累加所有子类±2区域内的 reads 数，除以父类区域内的 reads 总数，计算精度。
    •	处理父类 reads 为 0 的情况：如果父类区域内没有 reads，直接将精度设置为 0，避免除以 0 的错误。
    """
    # 初始化用于保存精度结果的列表
    precision_data = []

    with pysam.AlignmentFile(bam_file, "rb", threads=12) as bam:
        # 遍历字典中的父类和子类信息
        for parent_name, parent_info in tqdm(mirna_dict.items()):
            parent_contig = parent_info['Parent_contig']
            parent_start = parent_info['Parent_start']
            parent_end = parent_info['Parent_end']
            
            # 计算父类区域内的 reads 总数
            parent_reads_count = bam.count(region=None, start=parent_start, end=parent_end, contig=parent_contig)

            # 遍历父类的每个子类
            child_total_reads = 0
            for child_info in parent_info['Children']:
                child_contig = child_info['Child_contig']
                child_start = child_info['Child_start']
                child_end = child_info['Child_end']

                # 在子类的区域加上 ±2 范围
                child_start_range = max(parent_start, child_start - 2)
                child_end_range = min(parent_end, child_end + 2)

                # 计算子类 ±2 范围内的 reads 数
                child_reads_count = bam.count(region=None, start=child_start_range, end=child_end_range, contig=parent_contig)
                
                # 累加子类范围内的 reads 数
                child_total_reads += child_reads_count

            # 避免父类 reads 总数为 0 的情况
            if parent_reads_count > 0:
                precision = child_total_reads / parent_reads_count
            else:
                precision = 0  # 如果父类区域没有 reads，则精度为 0
            
            # 检查精度是否大于 1
            if precision > 1:
                logging.info(f"异常: Precision > 1\n"
                    f"Parent: {parent_name}, Reads: {parent_reads_count}\n"
                    f"Children Total Reads: {child_total_reads}")
                
            # 将结果保存到 precision_data 列表
            precision_data.append({
                "MIR": parent_name,
                "Precision": precision
            })

    # 将结果转换为 DataFrame
    df_precision = pd.DataFrame(precision_data)
    
    return df_precision

def main():
    OUTDIR = "/bios-store1/chenyc/scripts/SeekmiRNATools/test"
    GENOME_FILE = "/bios-store1/chenyc/Reference_Source/Artemisia_Reference/Aannua.fa"
    GFF_FILE = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/res_loci_merge/converted.gff"
    BAM_FILE = "/home/chenyc/bioinfo/bios-store2/Cooperative_Project/Project_chensusu_240807N/R_0/merged.bam"

    logging.info(f"BAM file used: {BAM_FILE}")
    logging.info(f"GFF file used: {GFF_FILE}")

    hairpin_sequences = os.path.join(OUTDIR, "miRNA_hairpin_sequences.fa")
    mature_sequences = os.path.join(OUTDIR, "miRNA_mature_sequences.fa")
    gff_update_file = os.path.join(OUTDIR, "updated.gff")

    db = gffutils.create_db(GFF_FILE, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

    # 从GFF文件中提取miRNA的hairpin序列
    extract_fasta_from_gff(hairpin_sequences, GENOME_FILE, GFF_FILE)

    # peak search and calculate 
    logging.info(f"Starting miRNA Peak Search Processing")
    hairpin_dict = parse_gff_for_mirna(db, flag=True)
    df_peaks = calculate_peaks_and_rnafold(hairpin_sequences, hairpin_dict, BAM_FILE)
    logging.info(f"Caculate strand specificity")
    df_strand_specificity = calculate_strand_specificity(hairpin_dict, BAM_FILE)
    df = pd.merge(df_peaks, df_strand_specificity, on="MIR")

    # update gff file
    logging.info(f"Processing sequence file and update GFF file")
    df_positions = Long_to_Short(hairpin_sequences, mature_sequences, df)
    df_positions = pd.merge(df_positions, df_strand_specificity, left_on="ID", right_on="MIR", how="left")
    Update_GFF(GFF_FILE, df_positions, gff_update_file)

    # compute final precision for miRNA
    logging.info(f"Caculate final precision for miRNA")
    db = gffutils.create_db(gff_update_file, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
    mature_dict = parse_gff_for_mirna(db, flag=False)
    df_precision = compute_final_precision_for_mirna(mature_dict, BAM_FILE)

    # save results
    logging.info(f"Save results to csv file")
    df = pd.merge(df, df_precision, on="MIR")
    df.to_csv(os.path.join(OUTDIR, "miRNA_results.csv"), index=False)
    extract_fasta_from_gff(mature_sequences, GENOME_FILE, gff_update_file)

if __name__ == "__main__":
    logging.info(print_current_time("Program started"))
    main()
    logging.info(print_current_time("Program ended"))

