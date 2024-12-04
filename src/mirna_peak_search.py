#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :mirna_peak_search.py
# @Time      :2024/10/23 16:38:48
# @Author    :Yuchen@rlab
# @Description: This script is used to call miRNA peaks from the input bam file. And calculate the miRNA precision and recall.

import gffutils
import pysam
import argparse
import datetime
import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, Counter

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

def parse_gff_for_mirna(gff_file):
    # 创建数据库或加载已有数据库
    db = gffutils.create_db(gff_file, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

    # 初始化空字典存储结果
    result_dict = {}

    # 遍历所有特征
    for feature in db.all_features():
        # 如果特征包含 'Derives_from' 属性
        if 'Derives_from' in feature.attributes:
        # if 'Derives_from' in feature.attributes or 'ID' in feature.attributes:
            # 获取 Derives_from 的父特征 ID
            parent_id = feature.attributes['Derives_from'][0]

            # 获取父特征对象
            parent_feature = db[parent_id]
            
            # 获取父特征的 Name 属性、位置信息和链信息
            parent_name = parent_feature.attributes.get('Name', [parent_feature.id])[0]
            parent_start = parent_feature.start
            parent_end = parent_feature.end
            parent_strand = parent_feature.strand
            
            # 获取子特征的 Name 属性、位置信息和链信息
            child_name = feature.attributes.get('Name', [feature.id])[0]
            child_contig = feature.seqid
            child_start = feature.start
            child_end = feature.end
            
            # 初始化字典条目，如果父类的 Name 还未添加到字典中，则创建新的条目
            if parent_name not in result_dict:
                result_dict[parent_name] = {
                    "Parent_contig": parent_feature.seqid,
                    "Parent_start": parent_start,
                    "Parent_end": parent_end,
                    "Parent_strand": parent_strand,  # 添加父类的链信息
                    "Children": []  # 创建一个空列表来存储子类信息
                }

            # 将子特征的信息（Name、位置信息和链信息）添加到父类的 Children 列表中
            result_dict[parent_name]["Children"].append({
                "Child_name": child_name,
                "Child_contig": child_contig,
                "Child_start": child_start,
                "Child_end": child_end
            })
    return result_dict

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
        for parent_name, parent_info in mirna_dict.items():
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
                print(f"异常: Precision > 1\n"
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

def calculate_strand_specificity(mirna_dict, bam_file):
    """
    •	链信息匹配：
        •	read.is_reverse 属性用于确定比对是否是反向链的。通过检查 read.is_reverse，我们可以确定某个 read 是正链还是反链。
        •	如果父类链信息是 '+'（正链），且 read 不是反向链（not read.is_reverse），则认为这个 read 与父类链相同。
        •	如果父类链信息是 '-'（负链），且 read 是反向链（read.is_reverse），则认为这个 read 与父类链相同。
        •	总 reads 数量：使用 total_reads_count 来统计父类区域内的所有 reads 数。
	•	链特异性计算：
        •	如果父类区域内没有 reads，链特异性为 0。
        •	否则，链特异性等于与父类链相同的 reads 数量除以父类区域内的 reads 总数。
    """

    with pysam.AlignmentFile(bam_file, "rb", threads=12) as bam:
        # 初始化用于保存链特异性结果的列表
        strand_specificity_data = []

        # 遍历字典中的父类和子类信息
        for parent_name, parent_info in mirna_dict.items():
            parent_contig = parent_info['Parent_contig']
            parent_start = parent_info['Parent_start']
            parent_end = parent_info['Parent_end']
            parent_strand = parent_info['Parent_strand']  # 父类的链信息

            # 计数变量
            same_strand_count = 0  # 与父类链相同的reads数
            total_reads_count = 0  # 父类区域内的总reads数

            # 遍历父类区域内的所有 reads
            for read in bam.fetch(region=None, start=parent_start, end=parent_end, contig=parent_contig):
                total_reads_count += 1  # 总reads数 +1
                
                # 检查read的链信息是否与父类相同
                if (parent_strand == '+' and not read.is_reverse) or (parent_strand == '-' and read.is_reverse):
                    same_strand_count += 1  # 与父类链相同的reads数 +1

            # 避免除以 0 的情况
            if total_reads_count > 0:
                strand_specificity = same_strand_count / total_reads_count
            else:
                strand_specificity = 0  # 如果没有 reads，链特异性为 0

            # 将结果保存到 strand_specificity_data 列表
            strand_specificity_data.append({
                "MIR": parent_name,
                "Strand_specificity": strand_specificity
            })

    # 将结果转换为 DataFrame
    df_strand_specificity = pd.DataFrame(strand_specificity_data)

    return df_strand_specificity

def calculate_peaks_and_rnafold(mirna_dict, bam_file):
    # 打开 BAM 文件
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # 获取 BAM 文件中的总reads数，用于RPM计算
        bam_total_reads = bam.count()
        
        # 初始化用于保存结果的列表
        peak_data = []

        # 遍历字典中的父类信息
        for parent_name, parent_info in mirna_dict.items():
            parent_contig = parent_info['Parent_contig']
            parent_start = parent_info['Parent_start']
            parent_end = parent_info['Parent_end']

            # 获取父类序列并计算 RNAfold 预测结果
            parent_sequence = get_parent_sequences_from_fasta(FASTA_FILE, parent_name)
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

                 # 检查 read 是否是反向互补链
                if read.is_reverse:
                    # 获取反向互补序列
                    seq = str(Seq(read.query_sequence).reverse_complement())
                else:
                    # 正向链序列
                    seq = read.query_sequence

                # 记录序列
                read_sequences_by_position[read_start].append(seq)
                position_counts[read_start] += 1


            # 寻找相邻位置±2nt的peak
            peaks = []
            current_peak_count = 0
            current_peak_start = None
            current_peak_reads = []

            sorted_positions = sorted(position_counts.keys())
            for i, pos in enumerate(sorted_positions):
                if current_peak_start is None:
                    current_peak_start = pos

                if i > 0 and pos > sorted_positions[i - 1] + 25:
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

def process_mirna_peak_search(FASTA_FILE, BAM_FILE, GFF_FILE, OUTPUT_FILE):
    print(f"Starting miRNA Peak Search Processing")
    print(f"FASTA file used: {FASTA_FILE}")
    print(f"BAM file used: {BAM_FILE}")
    print(f"GFF file used: {GFF_FILE}")
    print(f"Output file: {OUTPUT_FILE}")

    mirna_dict = parse_gff_for_mirna(GFF_FILE)
    df_precision = compute_final_precision_for_mirna(mirna_dict, BAM_FILE)
    df_strand_specificity = calculate_strand_specificity(mirna_dict, BAM_FILE)
    df_peaks = calculate_peaks_and_rnafold(mirna_dict, BAM_FILE)
    merged_df = pd.merge(df_peaks, df_strand_specificity, on='MIR', how='outer')
    merged_df = pd.merge(merged_df, df_precision, on='MIR', how='outer')
    merged_df.to_csv(OUTPUT_FILE, index=False)
    
    print(f"miRNA Peak Search Processing completed, results saved to {OUTPUT_FILE}")

# Argument parsing / help message / version
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument("-v", "--version", action="version",
                    version='%(prog)s v1.0.20241023')
parser.add_argument("-i", "--inputfile", type= str, default= os.getcwd(),
                    help="inputfile is a sorted bam file, mapped reads to the genome")
parser.add_argument("-o", "--outdir", type= str, default= os.getcwd(),
                    help="outputfile is a csv file, which contains the miRNA processing results")
parser.add_argument("-f", "--fastafile", type= str, default= "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_hairpin_bowtie_index/ath_hairpin_miRBase.fa",
                    help="fastafile is a fasta file, which contains the miRNA sequences")
parser.add_argument("-g", "--gfffile", type= str, default= "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_miRNA_miRBase.gff3",
                    help="gfffile is a gff file, which contains the miRNA annotations")
args = parser.parse_args()

if __name__ == "__main__":
    # FASTA_FILE = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_hairpin_bowtie_index/ath_hairpin_miRBase.fa"
    # BAM_FILE = "/bios-store1/chenyc/Project_miRNA_reannotation/PBOX_test/PBOX_2024_analysis3/merged_alignments.bam"
    # GFF_FILE = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_miRNA_miRBase.gff3"
    # OUTPUT_FILE = "/bios-store1/chenyc/test/test4precision/miRNA_processing_results.csv"
    
    # Program start
    print_current_time("Program started")

    # Assign file paths
    FASTA_FILE = args.fastafile
    BAM_FILE = args.inputfile
    GFF_FILE = args.gfffile
    OUTPUT_FILE = os.path.join(args.outdir, "miRNA_processing_results.csv")

    # Call the miRNA peak search processing function
    process_mirna_peak_search(FASTA_FILE, BAM_FILE, GFF_FILE, OUTPUT_FILE)

    # Program end
    print_current_time("Program ended")