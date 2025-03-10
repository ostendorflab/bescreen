import os
import sys
import argparse
import sqlite3
import polars as pl
import gzip
import csv

def arguments():

    parser = argparse.ArgumentParser(description='rsid_db',
                                     prog='python rsid_db.py')

    # required
    required = parser.add_argument_group("required arguments")

    required.add_argument('-v', '--dbsnp-vcf',
                          help='Path to the input VCF from dbSNP',
                          metavar='<PATH>',
                          required=True)
    required.add_argument('-d', '--rsid-db',
                          help='Path to the output rsID database',
                          metavar='<PATH>',
                          required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if args.dbsnp_vcf and not os.path.isfile(args.dbsnp_vcf):
        sys.exit('Input VCF from dbSNP not found!')

    if args.rsid_db and os.path.isfile(args.rsid_db):
        sys.exit('Output rsID database already exists!')


    snp_vcf_path = args.dbsnp_vcf
    rsid_db_path = args.rsid_db

    return (snp_vcf_path,
            rsid_db_path)


def generate_rsid_db(snp_vcf_path,
                     rsid_db_path):

    with sqlite3.connect(rsid_db_path) as rsid_db:

        rsid_cursor = rsid_db.cursor()

        rsid_cursor.execute(
        """
        CREATE TABLE rsids (
        rsid TEXT PRIMARY KEY,
        chrom TEXT,
        pos INTEGER,
        ref TEXT,
        alt TEXT
        ) WITHOUT ROWID
        """
        )

        with gzip.open(snp_vcf_path, 'rt') as vcf:
            reader = csv.reader(vcf, delimiter="\t")
            i = 0
            for row in reader:
                if not ''.join(row).startswith('#'):
                    rsid_cursor.execute(
                    f"""
                    INSERT INTO rsids (rsid, chrom, pos, ref, alt)
                    VALUES ('{row[2]}', '{row[0]}', '{row[1]}', '{row[3]}', '{row[4]}');
                    """
                    )
                    i += 1
                    if i % 1000000 == 0:
                        print(f'{i} entries written')
                        rsid_db.commit()

        rsid_db.commit()


def query_rsids(rsid_list,
                rsid_db_path):

    with sqlite3.connect(rsid_db_path) as rsid_db:

        rsid_cursor = rsid_db.cursor()
        rsid_cursor.execute(
        f"""
        SELECT * FROM rsids
        WHERE rsid IN ('{"', '".join(rsid_list)}');
        """
        )
        rsid_locations = rsid_cursor.fetchall()

    return rsid_locations


def query_rsid(rsid,
               rsid_db_path): # is this function of any use?

    with sqlite3.connect(rsid_db_path) as rsid_db:

        rsid_cursor = rsid_db.cursor()
        rsid_cursor.execute(
        f"""
        SELECT * FROM rsids
        WHERE rsid == '{rsid}';
        """
        )
        rsid_locations = rsid_cursor.fetchall()

    return rsid_locations


def transform_locations(rsid_locations):

    rsids = []
    chroms = []
    poss = []
    refs = []
    alts = []

    for rsid_location in rsid_locations:
        rsids.append(rsid_location[0])
        chroms.append(rsid_location[1])
        poss.append(rsid_location[2])
        refs.append(rsid_location[3].split(','))
        alts.append(rsid_location[4].split(','))

    rsidvars = pl.DataFrame({'rsid': rsids,
                            'chrom': chroms,
                            'pos': poss,
                            'ref': refs,
                            'alt': alts},
                            schema={'rsid': pl.String,
                                    'chrom': pl.String,
                                    'pos': pl.Int64,
                                    'ref': pl.List(pl.String),
                                    'alt': pl.List(pl.String)})

    rsidvars = rsidvars.explode('ref').explode('alt')
    rsidvars = rsidvars.with_columns(variant = pl.concat_str([pl.col('chrom'), pl.col('pos'), pl.col('ref'), pl.col('alt')], separator='_'))

    return rsidvars


if __name__ == '__main__':

    snp_vcf_path_arg, \
    rsid_db_path_arg = arguments()

    generate_rsid_db(snp_vcf_path_arg,
                     rsid_db_path_arg)
