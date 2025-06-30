#!/usr/bin/env python
# coding: utf-8

import os
import re
import argparse
import pandas as pd

def removeChr(contig):
	return contig.replace("chr", "")

def parseGtfAttribute(attribute, gtf):
    if attribute == "exon_number":
        attribute_pattern = rf'{attribute} (\d+)'
    else:
        attribute_pattern = rf'{attribute} "([^"]+)"'
    parsed = gtf["attributes"].apply(lambda x: re.search(attribute_pattern, x).group(1) if re.search(attribute_pattern, x) else "")
    failedToParse = parsed == ""
    if any(failedToParse):
        print(f"Failed to parse '{attribute}' attribute of {sum(failedToParse)} record(s).")
    return parsed

def read_annotations(annotation_file, ref_genome):
    print("Loading annotation for genome: ", ref_genome)
    exonsFile = annotation_file
    exons = pd.read_csv(exonsFile, sep="\t", comment="#", quotechar='"', names=["contig", "src", "type", "start", "end", "score", "strand", "frame", "attributes"], low_memory=False)
    exons = exons[exons["type"].isin(["CDS"])]
    exons = exons[["contig", "type", "start", "end", "strand", "attributes"]]
    exons["contig"] = exons["contig"]
    exons["geneID"] = parseGtfAttribute("gene_id", exons)
    exons["geneName"] = parseGtfAttribute("gene_name", exons)
    exons["transcript"] = parseGtfAttribute("transcript_id", exons)
    exons["exonNumber"] = parseGtfAttribute("exon_number", exons)
    exons.drop(["attributes"], axis=1, inplace=True)
    exons.drop_duplicates(inplace=True, ignore_index=True)
    return exons

def query_annotations(annotations, contig, position):
    filtered_annotations = annotations[(annotations['contig'] == str(contig)) & (annotations['start'] <= position) & (annotations['end'] >= position)].copy()
    if not filtered_annotations.empty:
        filtered_annotations['distance'] = 0
        return filtered_annotations[['geneID', 'geneName', 'transcript', 'exonNumber', 'distance']].reset_index(drop=True)
    else:
        annotations_contig = annotations[annotations['contig'] == str(contig)].copy()
        if annotations_contig.empty:
            return None
        annotations_contig['distance'] = annotations_contig.apply(lambda x: min(abs(x['start'] - position), abs(x['end'] - position)), axis=1)
        closest_records = annotations_contig.nsmallest(1, 'distance')
        return closest_records[['geneID', 'geneName', 'transcript', 'exonNumber', 'distance']].reset_index(drop=True)

def parse_svrvor(sample, input_file, ref_genome, annotation_file):
    annotations = read_annotations(annotation_file, ref_genome)
    print("Processing sample: ", sample)
    in1 = input_file
    sv0 = pd.read_csv(in1, sep="\t")
    l_svbc = []
    l_gen1 = []
    l_gen2 = []
    l_str1 = []
    l_str2 = []
    l_brk1 = []
    l_brk2 = []
    l_sit1 = []
    l_sit2 = []
    l_type = []
    l_spl1 = []
    l_spl2 = []
    l_dism = []
    l_cov1 = []
    l_cov2 = []
    l_conf = []
    l_rfra = []
    l_tags = []
    l_rpds = []
    l_clo1 = []
    l_clo2 = []
    l_gid1 = []
    l_gid2 = []
    l_tra1 = []
    l_tra2 = []
    l_dir1 = []
    l_dir2 = []
    l_filt = []
    l_ftra = []
    l_pseq = []
    l_rdid = []
    l_chr1 = []
    l_sta1 = []
    l_end1 = []
    l_chr2 = []
    l_sta2 = []
    l_end2 = []
    for index, row in sv0.iterrows():
        l_svbc.append(row['SV_BC'])
        l_gen1.append(row['gene1'])
        l_gen2.append(row['gene2'])
        if row['STRANDS'] == '++':
            l_str1.append('+/+')
            l_str2.append('+/+')
        elif row['STRANDS'] == '--':
            l_str1.append('-/-')
            l_str2.append('-/-')
        elif row['STRANDS'] == '+-':
            l_str1.append('+/+')
            l_str2.append('-/-')    
        elif row['STRANDS'] == '-+':
            l_str1.append('-/-')
            l_str2.append('+/+')
        else:
            l_str1.append('.')
            l_str2.append('.') 
        l_brk1.append(str(row['CHROM']) + ':' + str(row['POS']))
        l_brk2.append(str(row['CHR2']) + ':' + str(row['END']))
        l_sit1.append(row['site1'])
        l_sit2.append(row['site2'])
        if row['SVTYPE'] == 'TRA':
            l_type.append('Translocation')
        elif row['SVTYPE'] == 'DUP':
            l_type.append('Duplication')
        elif row['SVTYPE'] == 'INV':
            l_type.append('Inversion')
        elif row['SVTYPE'] == 'DEL':
            l_type.append('Deletion')
        else:
            l_type.append('.')
        l_spl1.append('.')
        l_spl2.append('.')
        l_dism.append('.')
        l_cov1.append('.')
        l_cov2.append('.')
        l_conf.append('High')
        l_tags.append('.')
        l_rpds.append('.')
        atr1 = query_annotations(annotations, row['CHROM'], row['POS'])
        if atr1 is not None:
            l_clo1.append('exon ' + str(atr1['exonNumber'][0]) + ' ' + str(atr1['distance'][0]) + 'bp away')
            l_gid1.append(atr1['geneID'][0])
            l_tra1.append(atr1['transcript'][0])                    
        else:
            l_clo1.append('.')
            l_gid1.append('.')
            l_tra1.append('.')
        atr2 = query_annotations(annotations, row['CHR2'], row['END'])
        if atr2 is not None:
            l_clo2.append('exon ' + str(atr2['exonNumber'][0]) + ' ' + str(atr2['distance'][0]) + 'bp away')
            l_gid2.append(atr2['geneID'][0])
            l_tra2.append(atr2['transcript'][0])
        else:
            l_clo2.append('.')
            l_gid2.append('.')
            l_tra2.append('.')
        if (atr1['distance'][0] != 0) or (atr2['distance'][0] != 0):
            l_rfra.append('out-of-frame')
        else:
            l_rfra.append('in-frame')
        if row['STRANDS'] == '++':
            l_dir1.append('downstream')
            l_dir2.append('downstream')
        elif row['STRANDS'] == '--':
            l_dir1.append('upstream')
            l_dir2.append('upstream')
        elif row['STRANDS'] == '+-':
            l_dir1.append('downstream')
            l_dir2.append('upstream')
        elif row['STRANDS'] == '-+':
            l_dir1.append('upstream')
            l_dir2.append('downstream')
        l_filt.append('PASS')
        l_ftra.append('.')
        l_pseq.append('.')
        l_rdid.append('.')
        l_chr1.append(row['CHROM'])
        l_sta1.append(row['POS'])
        l_end1.append(row['POS'] + 1 + abs(row['SVLEN']))
        l_chr2.append(row['CHR2'])
        l_sta2.append(row['END'])
        l_end2.append(row['END'] + 1 + abs(row['SVLEN']))
    df1 = pd.DataFrame()
    df1['SV_BC'] = pd.Series(l_svbc)
    df1['gene1'] = pd.Series(l_gen1)
    df1['gene2'] = pd.Series(l_gen2)
    df1['strand1(gene/SV)'] = pd.Series(l_str1)
    df1['strand2(gene/SV)'] = pd.Series(l_str2)
    df1['Breakpoint1'] = pd.Series(l_brk1)
    df1['Breakpoint2'] = pd.Series(l_brk2)
    df1['site1'] = pd.Series(l_sit1)
    df1['site2'] = pd.Series(l_sit2)
    df1['type'] = pd.Series(l_type)
    df1['split_reads1'] = pd.Series(l_spl1) 
    df1['split_reads2'] = pd.Series(l_spl2) 
    df1['discordant_mates'] = pd.Series(l_dism)
    df1['coverage1'] = pd.Series(l_cov1)
    df1['coverage2'] = pd.Series(l_cov2)
    df1['confidence'] = pd.Series(l_conf)
    df1['reading_frame'] = pd.Series(l_rfra)
    df1['tags'] = pd.Series(l_tags)
    df1['retained_protein_domains'] = pd.Series(l_rpds)
    df1['closest_genomic_breakpoint1'] = pd.Series(l_clo1)
    df1['closest_genomic_breakpoint2'] = pd.Series(l_clo2)
    df1['gene_id1'] = pd.Series(l_gid1)
    df1['gene_id2'] = pd.Series(l_gid2)
    df1['transcript_id1'] = pd.Series(l_tra1)
    df1['transcript_id2'] = pd.Series(l_tra2)
    df1['direction1'] = pd.Series(l_dir1)
    df1['direction2'] = pd.Series(l_dir2)
    df1['filters'] = pd.Series(l_filt)
    df1['SV_transcript'] = pd.Series(l_ftra)
    df1['peptide_sequence'] = pd.Series(l_pseq)
    df1['read_identifiers'] = pd.Series(l_rdid)
    ot1 = sample + '_DrawSVs.tsv'
    df1.to_csv(ot1, sep="\t", index=False)
    df2 = pd.DataFrame()
    df2['chr1'] = pd.Series(l_chr1)
    df2['start1'] = pd.Series(l_sta1)
    df2['end1'] = pd.Series(l_end1)
    df2['chr2'] =  pd.Series(l_chr2)
    df2['start2'] = pd.Series(l_sta2)
    df2['end2'] = pd.Series(l_end2)
    df2['type'] = pd.Series(l_type)
    ot2 = sample + '_links.csv'
    df2.to_csv(ot2, index=False)
    
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process structural variant files and generate output.")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--input", help="Input file")
    parser.add_argument("--genome", default="HG19VS", help="Reference genome")
    parser.add_argument("--annotations", help="Genome corresponding annotation GTF file")
    args = parser.parse_args()

    parse_svrvor(args.sample, args.input, args.genome, args.annotations)