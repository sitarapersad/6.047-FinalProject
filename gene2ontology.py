import MySQLdb
import pandas as pd
import operator

def get_conn():
    db = MySQLdb.connect(host="mysql-amigo.ebi.ac.uk",    # your host, usually localhost
                         user="go_select",         # your username
                         passwd="amigo",  # your password
                         db="go_latest",
                         port=4085)        # name of the data base

    # you must create a Cursor object. It will let
    #  you execute all the queries you need
    cur = db.cursor()
    return db, cur


# Use all the SQL you like
def get_annotations(genes):
    db, cur = get_conn()
    format_strings = ','.join(['%s'] * len(genes))

    sql="""
        SELECT DISTINCT
            gene_product.symbol,
            term.acc,
            term.name,
            term.term_type,
            gene_product.full_name

        FROM
            gene_product
            INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)
            INNER JOIN species ON (gene_product.species_id=species.id)
            INNER JOIN association ON (gene_product.id=association.gene_product_id)
            INNER JOIN evidence ON (association.id=evidence.association_id)
            INNER JOIN term ON (association.term_id=term.id)
        WHERE
            gene_product.symbol IN (%s)
            AND species.genus = 'homo'
            AND species.species = 'sapiens'
        ORDER BY gene_product.symbol, term.term_type, term.acc
        ;
        """ % format_strings

    df = pd.read_sql(sql, params=tuple(genes), con=db)
    db.close()
    return df


def get_genes(disease1, disease2):
    f = open('genes_of_interest/unique_genes_'+disease1+'_'+disease2+'.txt', 'r')
    genes = []
    for line in f:
        genes.append(line[:-1])  # remove \n

    return genes

def get_annotations_from_disease(disease1, disease2):
    df = get_annotations(get_genes(disease1, disease2))
    df.to_csv('annotations/annot_'+disease1+'_'+disease2+'.csv')
    return df

def analysis():
    diseases = ['aut', 'add', 'bip', 'mdd', 'scz']
    go_dict = {}
    for i in xrange(len(diseases) - 1):
        for j in xrange(i + 1, len(diseases)):
            df = pd.read_csv('annotations/annot_'+diseases[i]+'_'+diseases[j]+'.csv')
            df = df[df.term_type == 'biological_process']
            for index, row in df.iterrows():
                go = row['acc']
                go_dict[go] = go_dict.get(go, {})
                go_dict[go]['diseases'] = go_dict.get('diseases', set()) | {diseases[i], diseases[j]}
                go_dict[go]['description'] = go_dict.get('description', set()) | {row['name']}
    df = pd.DataFrame()
    df['go'] = pd.Series([go for go in go_dict])
    df['diseases'] = pd.Series([go_dict[go]['diseases'] for go in go_dict])
    df['counts'] = pd.Series(len(go_dict[go]['diseases']) for go in go_dict)
    df['descriptions'] = pd.Series([go_dict[go]['description'] for go in go_dict])

    df.to_csv('test', sep='\t')
    return df.sort_values(['counts'], ascending=False)

if __name__ == '__main__':
    diseases = ['aut', 'add', 'bip', 'mdd', 'scz']
    for i in xrange(len(diseases) - 1):
        for j in xrange(i + 1, len(diseases)):
            print diseases[i], diseases[j]
            get_annotations_from_disease(diseases[i], diseases[j])