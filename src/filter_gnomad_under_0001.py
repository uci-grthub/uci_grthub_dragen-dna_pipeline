#!/usr/bin/env python

import gzip
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Filter variants from Nirvana JSON output based on gnomAD allele frequency.")
parser.add_argument("json_path", type=str, help="Path to the Nirvana JSON output file (.json.gz)")
args = parser.parse_args()

json_path = Path(args.json_path)

# Read Nirvana JSON output by lines

# Nirvana / Illumina Connected Annotations writes a line-delimited JSON:
#   line 1               {"header":{...},"positions":[
#   position objects     {"chromosome":...},   (one per line, trailing comma)
#   ]                     closes the positions array
#   ,"genes":[            optional genes section (may be its own line or
#                         concatenated as  ],"genes":[  )
#   gene objects          {"name":...},
#   ,"samples":[ / ]}     other sections we ignore
# Section markers can appear standalone or fused onto a data line, so detect
# them by substring rather than exact match.
header = ''
positions = []
genes = []
section = None  # 'positions' | 'genes' | None
with gzip.open(json_path, 'rt') as f:
    for line in f:
        t = line.strip()
        if not t:
            continue
        if t.startswith('{"header":'):
            # header object + opening of the positions array
            head_part = t[:-len(',"positions":[')] if t.endswith(',"positions":[') else t
            header = head_part[len('{"header":'):]
            section = 'positions'
            continue
        # section-boundary lines (close an array and/or open the next one)
        if t.startswith(']') or t.startswith(',"'):
            if '"genes":[' in t:
                section = 'genes'
            elif '"positions":[' in t:
                section = 'positions'
            else:
                section = None  # samples or final }
            continue
        # data line belonging to the current section
        obj = t.rstrip(',')
        if section == 'positions':
            positions.append(obj)
        elif section == 'genes':
            genes.append(obj)

position_count = len(positions)
gene_count = len(genes)
print('header object:', header)
print('number of positions:', position_count)
print('number of genes:', gene_count)

# Retrieve variants under a gnomAD allele frequency threshold

import pandas as pd
import json

variants_field = 'variants'
gnomad_field = 'gnomad'
freq_threshold = 0.0001
freq_data = {'variant_id': [], 'gnomAD_allele_freq': []}
for position in positions:
    position_dict = json.loads(position)
    if variants_field in position_dict:
        for variant_dict in position_dict[variants_field]:
            if gnomad_field in variant_dict:
                freq = variant_dict[gnomad_field]['allAf']
                if freq < freq_threshold:
                    freq_data['variant_id'].append(variant_dict['vid'])
                    freq_data['gnomAD_allele_freq'].append(freq)

freq_df = pd.DataFrame(data=freq_data)
freq_df

# Retrieve all relevant genes and their OMIM gene names

gene_data = {'gene': [], 'OMIM_gene_name': []}
for gene in genes:
    gene_dict = json.loads(gene)
    gene_data['gene'].append(gene_dict['name'])
    omim_gene_name = ''
    if 'omim' in gene_dict:
        omim_dict = gene_dict['omim'][0]
        if 'geneName' in omim_dict:
            omim_gene_name = omim_dict['geneName']
    gene_data['OMIM_gene_name'].append(omim_gene_name)

gene_df = pd.DataFrame(data=gene_data)
gene_df

# Retrieve variants under a gnomAD allele frequency threshold

import pandas as pd
import json

variants_field = 'variants'
gnomad_field = 'gnomad'
transcripts_field = 'transcripts'
freq_threshold = 0.0001
freq_data = {'variant_id': [], 'gnomAD_allele_freq': [], 'transcripts': []}
for position in positions:
    position_dict = json.loads(position)
    if variants_field in position_dict:
        for variant_dict in position_dict[variants_field]:
            if gnomad_field in variant_dict and transcripts_field in variant_dict:
                freq = variant_dict[gnomad_field]['allAf']
                if freq < freq_threshold:
                    freq_data['variant_id'].append(variant_dict['vid'])
                    freq_data['gnomAD_allele_freq'].append(freq)
                    freq_data['transcripts'].append(variant_dict['transcripts'])

freq_df2 = pd.DataFrame(data=freq_data)

df_exploded = freq_df2.explode('transcripts')
df_unnested = pd.json_normalize(df_exploded['transcripts']).set_index(df_exploded['variant_id']).reset_index()

# Derive output CSV filename from the input JSON path so the output matches the input base name
# Handle .json.gz specially to avoid producing filenames like sample.json.csv when using Path.with_suffix on .json.gz
if json_path.name.endswith('.json.gz'):
    output_name = json_path.name[:-len('.json.gz')] + '.csv'
else:
    # Path.with_suffix removes only the last suffix; this works for .json, .gz, or no suffix
    output_name = json_path.with_suffix('.csv').name

output_path = json_path.with_name(output_name)

df_unnested.to_csv(output_path)
