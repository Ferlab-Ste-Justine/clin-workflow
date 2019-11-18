#### clin-workflow

ETL of the workflow to create a mutation-centric index in Elasticsearch 6.7.x

see clin-mutation-centric.json

= 0..1 > mutation  = 0..n > donor

donors is of "nested type"

## For description of the annotations
See https://github.com/cr-ste-justine/clin-project/wiki/SNPEFF
See https://github.com/cr-ste-justine/clin-project/wiki/SNPSIFT-et-annotations
And
https://usf.app.box.com/s/cdws8yx5occ603ccbknwyamz5reapdug

##### ElasticSearch Set-up
To create the index 'mutations', run the following command line where ElasticSearch is available (ssh -L or ssh thru environment)

```curl -XPUT "http://localhost:9200/mutations" -H 'Content-Type: application/json' -d @clin-mutation-centric.json```

#### ETL Genomic Algorithm 
     
    Initiate Spark
    Open VCF and split into partitions
    For each partition:
     
        Read a line / Prepare ES payload from line
            clinical data fetching (per specimen/donor)
            Lire les ligne du VCF
            Lire les annotations de VEP (comma; pipe; &)
            external db analysis
            
            transcripts grouping (gene, aaChange, consequence, CDNAChange, strand)
            VEP imapct scoring 0-4 
            Gene analysis (Redis)
            donor/specimen analysis (FHIR)
            internal & laboratory frequencies
         Group by a number of ES transactions (parameters)
         Bulk operation to ES
         
### Predictions
reference: https://onlinelibrary.wiley.com/doi/full/10.1002/humu.21517

#### FATHMM
- D : DAMAGING
- T : TOLERATED

#### POLYPHEN2_HVAR
“benign”, “possibly damaging”, “probably damaging”
- B : benign
- P : possibly damaging
- D : probably damaging

#### LRT

- N :  predicted N(eutral) 
- D :  predicted D(eleterious) 
- U :  U(nknown)

#### SIFT

- T : T(olerated)
- D : D(amaging)

