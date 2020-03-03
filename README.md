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
            family & transmission analysis (AD, AR, DeNovo, XD, XR)
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

### Frequencies

#### internal cohort & by Lab (LDx) & eventually by studies
- any result with a dot (.) i.e ./. or 1/. is discraded
- PN est le nombre de patient ayant une mutation (1/0, 0/1, 1/1)  sur l'allèle en question
- PN is the number of patient that have a mutation on the allele (1/0, 0/1, 1/1)
- AC est le nombre d'allèle muté (1 seulement)
- AC is the count of mutated allele found (1)
- AN est le nombre total d'allèle rencontré (0 ou 1)
- AN is the total number of allele found (0 or 1)
- HC est le nombre de Homozygote rencontré (1/1; 1|1)
- HC is the total number of Homozygous individual (1/1; 1|1)
- AF = AC / AN
#####NOTA BENE
- Les vcfs on été normalisés...  donc, on pas de 2+
- We treat normalized VCFs, we no longer have value greater than 1 on the genotype of a patient

### Donor annotations
- adAlt : Allelic depths for the alt alleles
- adTotal : Total Allelic depths for the ref and alt alleles
- gq : Genotype Quality (integer)
- gt : Genotype
- adFreq : Ratio between adAlt and adTotal 
- qd : Variant Confidence/Quality by Depth (float)

### Others
- Ensembl_transcriptid=(from dbNSFP) Ensembl transcript ids (Multiple entries separated by ";")
- FeatureId (from VEP) Feature - Ensembl stable ID of feature 


### Notes
- Fev 21, 2020
on ne garde pas le ensemblTranscriptId  provenant de vep/dbNSFP; décision prise par Vincent sur Slack.
Enlever ensemblTranscriptID -- next version de l’index - check
Puisqu’on utilise pas le champ Picked, on devrait l’enlever. On va essayer d’épurer un peu l’index
