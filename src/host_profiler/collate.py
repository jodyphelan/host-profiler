import argparse
import os
import sys
import json
import subprocess as sp
import pandas as pd
from tqdm import tqdm
import pathogenprofiler as pp
from collections import defaultdict


def collate(args: argparse.Namespace):
    """Main function for the collate cli subcommand."""

    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    if len(samples) == 0:
        raise ValueError("No samples found in directory: %s" % args.dir)


    fields = ["sample_id","chrom","pos","db-snp-id","gene_name","gene_id","change","nucleotide_change","protein_change","depth","freq","forward_reads","reverse_reads","filter","type","genotype","description"]
    variant_rows = []
    variant_info = {}
    coverage_info = []

    for s in tqdm(samples, desc="Loading reports"):
        # Data has the same structure as the .result.json files
        if not os.path.isfile("%s/%s%s" % (args.dir,s,args.suffix)):
            sys.stderr.write(f"Can't find sample {s}\n")
            continue
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for region in data['qc']['target_qc']:
            region['sample_id'] = s
            coverage_info.append(region)
        for var in data["variants"] + data['fail_variants']:
            var['sample_id'] = s
            if 'drugs' in var:
                var['drugs'] = ",".join([a["drug"] for a in var["drugs"]])
            else:
                var['drugs'] = ''
            var['genotype'] = None if var['filter'] == "soft_fail" else 1
            for ann in var['annotation']:
                var.update(ann)
            variant_rows.append({f:var.get(f,None) for f in fields})
            variant_info[var['gene_name'],var['change']] = var

    df = pd.DataFrame(variant_rows)

    variants = set([(v.chrom,v.pos) for v in list(df[['chrom','pos']].value_counts().reset_index(name='count').itertuples(index=False))])
    print(variants)

    mutation2genomepos = {}
    samples2mutations = defaultdict(set)
    for v in df.itertuples(index=False):
        mutation2genomepos[(v.gene_name,v.change)] = (v.chrom,v.pos)
        samples2mutations[v.sample_id].add((v.gene_name,v.change))

    conf = pp.get_db(db_dir=args.db_dir, db_name = args.db)
    for g in conf['json_db']:
        for m in conf['json_db'][g]:
            mut = conf['json_db'][g][m]
            if (g,m) not in mutation2genomepos:
                always_report = False

                for ann in mut['annotations']:
                    if 'always_report' in ann:
                        always_report = True
                if always_report:
                    mutation2genomepos[(g,m)] = (mut['chromosome'],mut['genome_positions'][0])
                    variants.add((mut['chromosome'],mut['genome_positions'][0]))
                    variant_info[(g,m)] = {
                        "chrom":mut['chromosome'],
                        "pos": mut['genome_positions'][0],
                        "gene_name":g,
                        "gene_id": g,
                        "change": '',
                        "nucleotide_change": '',
                        "protein_change": '',
                        "type": '',
                    }
                    variant_info[(g,m)].update(ann)

    depth = defaultdict(lambda: defaultdict(int))
    for s in tqdm(samples, desc="Calculating depth"):
        bam = pp.filecheck(f"{args.dir}/{s}.bam")
        for l in sp.Popen(f'samtools depth {bam}',shell=True,stdout=sp.PIPE).stdout:
            chrom,pos,dp = l.decode().strip().split()
            if (chrom,int(pos)) in variants:
                depth[(chrom,int(pos))][s] = dp




    for s in samples:
        for m in mutation2genomepos:
            gpos = mutation2genomepos[m]
            if m in samples2mutations[s]:
                pass
            else:
                v = dict(variant_info[m])
                    
                variant_rows.append({
                    "sample_id": s,
                    "chrom":gpos[0],
                    "pos": gpos[1],
                    "db-snp-id": v['db-snp-id'],
                    "gene_name":v['gene_name'],
                    "gene_id": v['gene_id'],
                    "change": v['change'],
                    "nucleotide_change": v['nucleotide_change'],
                    "protein_change": v['protein_change'],
                    "depth": depth[gpos][s],
                    "freq": 0,
                    "forward_reads": None,
                    "reverse_reads": None,
                    "filter":'reference' if int(depth[gpos][s]) >= args.min_depth else 'depth_fail',
                    "type": v['type'],
                    "genotype": 0 if int(depth[gpos][s]) >= args.min_depth else 'NA',
                    "description": v['description'],
                })

    df = pd.DataFrame(variant_rows)

    df.to_csv(args.out+'.variants.csv',index=False,na_rep="NA")

    covdf = pd.DataFrame(coverage_info)[['sample_id','target','percent_depth_pass','median_depth']]
    covdf.to_csv(args.out+'.coverage.csv',index=False,na_rep="NA")