#### clin-workflow

ETL of the workflow to create a mutation-centric index in Elasticsearch 6.7.x

see clin-mutation-centric.json

0..1 => mutation 0..n => donor => 0..n phenotypes

donor and phenotypes is of "nested type"

## For description of the annotations
See https://github.com/cr-ste-justine/clin-project/wiki/SNPEFF
See https://github.com/cr-ste-justine/clin-project/wiki/SNPSIFT-et-annotations
And
https://usf.app.box.com/s/cdws8yx5occ603ccbknwyamz5reapdug

##### ElasticSearch Set-up
To create the index 'mutations', run the following command line where ElasticSearch is available (ssh -L or ssh thru environment)

```curl -XPUT "http://localhost:9200/mutations_genes" -H 'Content-Type: application/json' -d @clin-mutation-centric.json```

