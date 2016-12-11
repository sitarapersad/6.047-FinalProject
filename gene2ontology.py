import MySQLdb

db = MySQLdb.connect(host="mysql-amigo.ebi.ac.uk",    # your host, usually localhost
                     user="go_select",         # your username
                     passwd="amigo",  # your password
                     db="go_latest",
                     port=4085)        # name of the data base

# you must create a Cursor object. It will let
#  you execute all the queries you need
cur = db.cursor()

# Use all the SQL you like
sql="""
SELECT
 term.name AS superterm_name,
 term.acc AS superterm_acc,
 term.term_type AS superterm_type,
 association.*,
 gene_product.symbol AS gp_symbol,
 gene_product.symbol AS gp_full_name,
 dbxref.xref_dbname AS gp_dbname,
 dbxref.xref_key AS gp_acc,
 species.genus,
 species.species,
 species.ncbi_taxa_id,
 species.common_name
FROM term
 INNER JOIN graph_path ON (term.id=graph_path.term1_id)
 INNER JOIN association ON (graph_path.term2_id=association.term_id)
 INNER JOIN gene_product ON (association.gene_product_id=gene_product.id)
 INNER JOIN species ON (gene_product.species_id=species.id)
 INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)
WHERE
 species.genus = 'Homo'
 AND
 species.species = 'sapiens'
 AND
 gene_product.symbol = 'SLC35F3';
"""

cur.execute(sql)

for row in cur.fetchall():
    print row

db.close()