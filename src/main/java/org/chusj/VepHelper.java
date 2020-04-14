package org.chusj;


import org.apache.http.HttpHost;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.time.LocalDate;
import java.util.*;

import static org.chusj.PatientHelper.loadPedigree;
import static org.chusj.RedisGeneSetHelper.getMembersForEnsId;
import static org.chusj.VEPSparkDriverProgram.getMD5Hash;


public class VepHelper {

    private static float avgFuncAnnoPerMutation = 0.0f;
    private static int countFuncAnnoPerMutation = 0;
    private static int countMutation =0;
    private static int countRedisCalls =0;
    private static NumberFormat NF = NumberFormat.getInstance();
    private static JSONObject lastOne;
    private static boolean toPrint = false;
    private static final int CLINVAR = 0;
    private static final int OMIM = 1;
    private static final int ENSEMBL = 2;
    private static final int ORPHANET = 3;
    private static final int DBSNP = 4;
    private static final int PUBMED = 5;
    private static final int CLINVAR_SIG = 6;
    private static final int CLINVAR_TRAIT = 7;
    private static final int TYPES = 8;
    private static final int GENES = 9;
    private static final int PHENO = 10;
    private static final int EXPECTED_VEP_ANNOTATION_QTY = 66;
    private static final Map<String,String> FATHMM_PRED_VALUE
            = new HashMap<String, String>(2){{
        put("D", "DAMAGING");
        put("T", "TOLERATED");
    }};
    private static final Map<String,String> POLYPHEN2_HVAR_PRED_VALUE
            = new HashMap<String, String>(3){{
        put("B", "BENIGN");
        put("P", "POSSIBLY DAMAGING");
        put("D", "PROBABLY DAMAGING");
    }};
    private static final Map<String,String> LRT_PRED_VALUE
            = new HashMap<String, String>(3){{
        put("N", "NEUTRAL");
        put("D", "DELETERIOUS");
        put("U", "UNKNOWN");
    }};
    private static final Map<String,String> SIFT_PRED_VALUE
            = new HashMap<String, String>(2){{
        put("D", "DAMAGING");
        put("T", "TOLERATED");
    }};

    /*
    This main method is only used to test locally the JSON object creation
     */
    public static void main(String[] args) throws Exception {

        if (args.length != 3) {
            args = new String[]{"FAM_C3_92_new5k.txt", "pedigree.properties", "pedigree.ped"};
//            args = new String[]{"batch5000.txt",  "pedigree.properties", "pedigreeTest1.ped"};
        }

        String extractFile = args[0];
        String pedigrePropsFile = args[1];
        String pedFile = args[2];

        List<String> specimenList = getSpecimenList(extractFile);

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {


            PatientHelper.client = clientTry;
            Properties pedigreeProps = getPropertiesFromFile(pedigrePropsFile);
            List<Pedigree> pedigrees = loadPedigree(pedFile);
//            Map<String, Patient> patientMap = PatientHelper.preparePedigreeFromPedAndFHIR(pedigrees);
            Map<String, Patient> patientsMap = PatientHelper.preparePedigree(specimenList);
//        Map<String, Patient> patientMap = PatientHelper.preparePedigreeFromProps(pedigreeProps);

            specimenList.forEach(System.out::println);

            pedigrees.forEach(System.out::println);
//            patientMap.forEach((k, v) -> System.out.println(k + "\n\t" + v));
//            patientsMap.forEach((k, v) -> System.out.println(k + "\n\t" + v));


            Map<String,Family> familyMap = getFamilyMap(specimenList, pedigrees);

            familyMap.forEach((k, v) -> System.out.println(k + "\n\t" + v));

            //System.exit(0);

            try (BufferedReader buf = new BufferedReader(new FileReader(extractFile))) {

                String fetchedLine;

                buf.readLine();

                //System.out.println(pedigreeProps.toString());
//                List<String> hpoTerms = new ArrayList<>(Arrays.asList("HP:0005280","HP:0001773"));

                while (true) {
                    fetchedLine = buf.readLine();
                    if (fetchedLine == null) {
                        break;
                    } else {
                        Variant oneVariant = processVcfDataLine(fetchedLine, pedigreeProps, patientsMap,
                                pedigrees, specimenList, familyMap);

                        if (toPrint && oneVariant != null) {
                            JSONObject propertiesOneMutation = new JSONObject(oneVariant.getJsonObjInString());
                            System.out.println(propertiesOneMutation.toString(0));
                            //extractGenesFromMutation(propertiesOneMutation, "0", true).stream().forEach((gene) -> System.out.println("\t"+gene.toString(0)));
                            toPrint = false;
                        }

                        //System.err.println(propertiesOneMutation.toString(0));
                    }
                }
            }
            if (countMutation > 0) {
                avgFuncAnnoPerMutation = (float) countFuncAnnoPerMutation / countMutation;
            }

            if (lastOne != null) {
//            extractGenesFromMutation(lastOne, "0", false).stream().forEach((gene) -> System.out.println("gene="+gene.toString(0)));
                System.out.println("lastOne=" + lastOne.toString(2));
//            extractGenesFromMutation(lastOne).stream().forEach((gene) -> System.out.println("gene="+gene.toString(0)));
            }

            System.out.println("\navgFuncAnnoPerMutation=" + avgFuncAnnoPerMutation + "  mutationCount=" + countMutation);
            System.out.println("redisCallsCount=" + countRedisCalls);
        }
    }



    static Variant processVcfDataLine(String extractedLine,  Properties pedigreeProps,
                                         Map<String, Patient> patientMap, List<Pedigree> pedigrees,
                                         List<String> specimenList, Map<String,Family> familyMapRef) {

        // CHROM

        String[] lineValueArray = extractedLine.split("\t");
        // dynamic positioning system -- pos counter++
//        toPrint = false;
        int pos = 0;
        int impactScore;
        String chrom = lineValueArray[pos++];
        // CHROM <- line from snpsift extract fields execution
        // #CHROM <- line from extraction utilities for knowning list of specimen id
        if ("CHROM".equalsIgnoreCase(chrom) || "#CHROM".equalsIgnoreCase(chrom)) {
            return null; // Meta data line
        }

        boolean[] dbExt = {false, false, false, false, false, false, false, false, false, false, false};
        List<Set<String>> dbExtId = new ArrayList<>();

        for (int i=0; i<dbExt.length; i++) {
            dbExtId.add(i, new HashSet<>());
        }

        Set<Gene> geneSet = new HashSet<>();
        Set<Clinvar> clinvarSet = new HashSet<>();

        Map<Integer, FunctionalAnnotation> faMap = new HashMap<>();
        Map<String, Transcript> deStructuredTranscript = new HashMap<>();

        JSONObject propertiesOneMutation = new JSONObject();
        JSONArray donorArray = new JSONArray();
        JSONArray specimenArray = new JSONArray();
        JSONObject bdExtObj = new JSONObject();
        JSONObject availBdExtObj = new JSONObject();
        JSONObject clinvarObj = new JSONObject();
        JSONObject geneObj;
        JSONArray geneArray = new JSONArray();

        String build = pedigreeProps.getProperty("assemblyVersion");
        propertiesOneMutation.put("assemblyVersion", build);
        propertiesOneMutation.put("annotationTool", pedigreeProps.getProperty("annotationTool"));

        LocalDate localDate = LocalDate.now();
        propertiesOneMutation.put("lastAnnotationUpdate", localDate);


        Map<String,Family> familyMap = new HashMap<>();
        familyMapRef.forEach((k, v) -> familyMap.put(k,v.clone()));

        int nbDonor = specimenList.size();
        JSONObject[] arrayDonor = new JSONObject[nbDonor];
        for (int i=0; i<nbDonor; i++) {
            arrayDonor[i] = new JSONObject();
        }

        // POS	REF	ALT	FILTER	DP	MQ	QD	RPA	RU	STR
        countMutation++;
        String position = lineValueArray[pos++];

//        String dbSNP_ID = lineValueArray[pos++];
//        if (addStrToJsonObject("dbSNP_ID", dbSNP_ID, propertiesOneMutation, false)) {
//            // bdExtObj.put( new JSONObject().put("dbSNP",
//            dbExt[DBSNP] =  true;
//            addStrToListSet(dbExtId, DBSNP, dbSNP_ID);
//        }
        String reference = lineValueArray[pos++];

        String alt = lineValueArray[pos++].replace(",<NON_REF>", ""); // CT,<NON_REF> or G,TGG,<NON_REF>
//        String qual = lineValueArray[pos++];
        String filter = lineValueArray[pos++];

        String dpS = lineValueArray[pos++];
        String mqS = lineValueArray[pos++];
        String qdS = lineValueArray[pos++];
//        String rpaS = lineValueArray[pos++];
//        String ruS = lineValueArray[pos++];
//        String strS = lineValueArray[pos++];
//        String mqRankSum = lineValueArray[pos++];
//        String readPosRankSum = lineValueArray[pos++];
//        String lod = lineValueArray[pos++];
//        String fractionInformativeReads = lineValueArray[pos++];

        propertiesOneMutation.put("altAllele", alt);
        //GT:AD:AF:DP:F1R2:F2R1:FT:GP:GQ:PL:PP
        // GEN[*].GT	GEN[*].AD	GEN[*].AF	GEN[*].DP	GEN[*].FT	GEN[*].GQ	GEN[*].SB	CSQ

        String[] gt = lineValueArray[pos++].split(",");
        String adStr = lineValueArray[pos++];
        String[] ad = adStr.split(",", -1);
        String af = lineValueArray[pos++];//.split(",");
        //String[] dn = lineValueArray[pos++].split(",");
        String genDP = lineValueArray[pos++];//.split(",");
//        String[] f1r2 = lineValueArray[pos++].split(",");
//        String[] f2r1 = lineValueArray[pos++].split(",");
        String ft = lineValueArray[pos++];//.split(",");
//        String[] gp = lineValueArray[pos++].split(",");
        String[] gq = lineValueArray[pos++].split(",");
//        String[] dq = lineValueArray[pos++].split(",");
        String sb = lineValueArray[pos++];//.split(",");
//        String[] pl = lineValueArray[pos++].split(",");
//        String[] pp = lineValueArray[pos++].split(",");

//        String hiConfDeNovo = lineValueArray[pos++];

        String csq = lineValueArray[pos];

        String[] csqArray = csq.split(",");  // return vep functionalAnnotation array
        countFuncAnnoPerMutation += csqArray.length;

        String chrPos = chrom.substring(3); // remove 'chr'
        /*
        Comment est-ce que l'on dealer avec des chromosomes comme ceux-ci :
        14_GL000225v1_random; Un_KI270442v1; etc...  ???
        Alex DL 10:44 AM:
        En clinique on peut les jeter a la poubelle. PAr exemple exomiser ne retient que les variants dans les chromosomes chr1-chr22 et chrX chrY chrM
         */
        if (chrPos.length() > 2) {
            return null;
        }
        Variant variant = new Variant();
        String dnaChange = reference + ">" + alt.split(",")[0];
        String mutationId = "chr" + chrPos + ":g." + position + dnaChange;
        String uid = getMD5Hash(mutationId +"@"+build);

        propertiesOneMutation.put("id", uid);
        variant.setId(uid);
        propertiesOneMutation.put("mutationId", mutationId);
        variant.setMutationId(mutationId);
        propertiesOneMutation.put("dnaChange", dnaChange);
        propertiesOneMutation.put("chrom", chrPos);
        propertiesOneMutation.put("refAllele", reference);
        propertiesOneMutation.put("start", Long.valueOf(position));

        JSONArray functionalAnnotations = new JSONArray();
        JSONObject frequencies = null;
        JSONObject funcAnnotation;


        // vep annotation analysis (column CSQ)
        impactScore = 0;
        for (String s : csqArray) {

            if (s.trim().isEmpty()) {
                System.err.println("#############\n############\n"+
                        "ID="+mutationId+
                        "#############\n############\n");
                break;
            }

            funcAnnotation = processVepAnnotations(s, dbExtId, dbExt, geneSet, deStructuredTranscript,
                    reference, alt, clinvarSet);

            if (funcAnnotation == null) {
                System.err.println("#############\n############\n"+
                        "ID="+mutationId+
                        "#############\n############\n");
                break;
            }

            if (!funcAnnotation.isNull("frequencies")) {
                frequencies = (JSONObject) funcAnnotation.remove("frequencies");

            }
            functionalAnnotations.put(funcAnnotation);
            String funcAnnotationImpact = (String) funcAnnotation.remove("impact");
            int funcAnnotationScore = getImpactScore(funcAnnotationImpact);
            if (funcAnnotationScore > impactScore ) impactScore = funcAnnotationScore;

            String gene = "";
            String geneId = "";
            if (!funcAnnotation.isNull("geneAffectedSymbol")) {
                gene = (String) funcAnnotation.remove("geneAffectedSymbol");
                geneId = (String) funcAnnotation.remove("geneAffectedId");
            }
            String consequence = (String) funcAnnotation.remove("consequence");
            boolean canonical = (!funcAnnotation.isNull("canonical") && (boolean) funcAnnotation.get("canonical"));
            boolean isPick = (funcAnnotation.isNull("pick") ? false : true);
            Long strand = null;
            if (!funcAnnotation.isNull("strand")) {
                strand = (long) funcAnnotation.remove("strand");
            }
            String aaChange = "";
            if (!funcAnnotation.isNull("aaChange")) {
                aaChange = (String) funcAnnotation.remove("aaChange");
            }
            String cdnaChange = "";
            if (!funcAnnotation.isNull("cdnaChange")) {
                cdnaChange = (String) funcAnnotation.remove("cdnaChange");
            }
            //FunctionalAnnotation(String gene, String aaChange, String consequence, String codingDNAChange, long strand)
            JSONObject scores = (JSONObject) funcAnnotation.remove("conservationsScores");
            JSONObject predictions = null;
            if (!funcAnnotation.isNull("predictions")) {
                predictions = (JSONObject) funcAnnotation.remove("predictions");
            }
            String biotype = (String) funcAnnotation.remove("biotype");
            // Vincent from Slack on Fev 21, 2020 -> on va garder juste le featureID (donc pas le ensemblTranscriptId)
            funcAnnotation.remove("ensemblTranscriptId");

            FunctionalAnnotation functionalAnnotation = new FunctionalAnnotation(gene, aaChange, consequence, cdnaChange, strand);
//            Set<String> consequenceSet = new HashSet(Arrays.asList(consequence.split("&")));
//            if (consequenceSet.size()>1) {
//                toPrint=true;
//            }
            functionalAnnotation.setGeneId(geneId);
            functionalAnnotation.setImpact(funcAnnotationImpact);
            functionalAnnotation.setScores(scores);
            if (predictions !=null) functionalAnnotation.setPredictions(predictions);
            functionalAnnotation.setBiotype(biotype);

            Integer hash = functionalAnnotation.hashCode();
            if (faMap.containsKey(hash)) {
                FunctionalAnnotation prevFA = faMap.get(hash);
                if (isPick) {
                    prevFA.setPick(isPick);
                }
                if (canonical) {
                    prevFA.setCanonical(canonical);
                }
                prevFA.getTheRest().put(funcAnnotation);
            } else {
                JSONArray ja = new JSONArray();
                functionalAnnotation.setTheRest(ja);
                if (isPick) {
                    functionalAnnotation.setPick(isPick);
                }
                if (canonical) {
                    functionalAnnotation.setCanonical(canonical);
                }
                ja.put(funcAnnotation);
                faMap.put(hash, functionalAnnotation);
            }
        }

        JSONArray consequences = new JSONArray();


        propertiesOneMutation.put("impactScore", impactScore);
        //propertiesOneMutation.put("type",  variant_class.get("type"));
        String types = toStringList(dbExtId.get(TYPES));
        if (!types.trim().isEmpty()) {
            propertiesOneMutation.put("type", types);
        } else {
            System.err.println("empty type for :" + mutationId);
        }


        // Genes analysis
        for (Gene gene: geneSet) {
            geneObj = new JSONObject();
            if (gene.getGeneSymbol().trim().isEmpty()) continue;

            geneObj.put("geneSymbol", gene.getGeneSymbol());
            geneObj.put("ensemblId", gene.getEnsemblId());
            geneObj.put("biotype", gene.getBiotype());
            addGeneSetsToObjs(gene.getEnsemblId(), geneObj, availBdExtObj, patientMap, pedigrees);

//            String tableRow = "gene:" + gene.getGeneSymbol() + "@" +
//                    "type:" + gene.getBiotype() + "@" +
//                    "location:" + (geneObj.isNull("location") ? "" : geneObj.get("location")) + "@" +
//                    "ensemblId:" + gene.getEnsemblId();
//            geneObj.put("tableRow", tableRow);
            if (!geneObj.isNull("oldGeneSymbol")) {
                String newAlias = (String) geneObj.remove("oldGeneSymbol");
                String newSymbol = (String) geneObj.get("geneSymbol");
                gene.setNewSymbol(newSymbol);
                gene.setNewAlias(newAlias);

                faMap.forEach((k,v) -> {
                    if (v != null && v.getGeneId().equalsIgnoreCase(gene.getEnsemblId())) {
                        v.setGene(newSymbol);
                    }
                });
            }
            geneArray.put(geneObj);
            variant.addGene(gene);
        }
        faMap.forEach((k,v) -> {
            if (v != null) consequences.put(v.getJsonObj());
        });
        propertiesOneMutation.put("consequences", consequences);

        // Patient and donor analysis
        Map<String, Frequencies> frequenciesPerLabos = new HashMap<>();
        patientMap.forEach((id, patient) -> {
            String labName = patient.getLabAlias();
            if (!frequenciesPerLabos.containsKey(labName)) {
                Frequencies freqenceLabo = new Frequencies();
                frequenciesPerLabos.put(labName, freqenceLabo);
            }
        });

        int patientNb = 0;
        float alleleCount = 0, alleleNumber = 0, homozygoteCount = 0;
        float alleleFrequencies = 0f;
        // snpSift Extract does not keep correctly the qty of unknown value . (only give 1 number (0) in GEN[*].AD
        int adPointer = 0;

        for (int i=0; i< nbDonor; i++) {

            //Pedigree currentDonor = pedigrees.get(i);
            Patient currentPatient = patientMap.get(specimenList.get(i));
            String labo = currentPatient.getLabAlias();
            Frequencies freqenceLabo = frequenciesPerLabos.get(labo);
            String zygosity = zygosity(gt[i]);
            if (!"UNK".equalsIgnoreCase(zygosity)) {
                int countAltAllele = countAlternativeAllele(gt[i]);
                alleleCount += countAltAllele;
                freqenceLabo.setAc(freqenceLabo.getAc() + countAltAllele);
                int alleleNumberCount = countNumberOfAllele(gt[i]);
                alleleNumber += alleleNumberCount;
                freqenceLabo.setAn(freqenceLabo.getAn() + alleleNumberCount);
                if (!"HOM REF".equalsIgnoreCase(zygosity)) {
                    freqenceLabo.setPn(freqenceLabo.getPn() + 1);
                    patientNb++;
                }
                if ("HOM".equalsIgnoreCase(zygosity)) {
                    homozygoteCount++;
                    freqenceLabo.setHc(freqenceLabo.getHc()+1);
                }
                freqenceLabo.updateAf();
            }

            //if ("HOM REF".equalsIgnoreCase(zygosity)) continue;
            arrayDonor[i].put("lastUpdate", localDate);
            //arrayDonor[i].put("phenotypes", phenotypesArray);
//            addNumberToJsonObject("quality", qual, arrayDonor[i], false, 'f');
            arrayDonor[i].put("filter", filter);
            arrayDonor[i].put("gt", gt[i]);
            arrayDonor[i].put("zygosity", zygosity);
            addNumberToJsonObject("gq", gq[i], arrayDonor[i], false, 'l');
            String adS;
            // snpSift Extract does not keep correctly the qty of unknown value . (only give 1 number (0) in GEN[*].AD

            if (ad.length < nbDonor * 2) {
                //toPrint = true;
                if ("./.".equalsIgnoreCase(gt[i])) {
                    adS = ad[adPointer] + "," + ad[adPointer++]; // increment only by 1
                } else {
                    adS = ad[adPointer++] + "," + ad[adPointer++];
                }
            } else {
                adS = ad[i * 2] + "," + ad[(i * 2) + 1];
            }

            String[] adSArray = adS.split(",");
//            arrayDonor[i].put("ad", adS);
            int adRef = Integer.parseInt(adSArray[0]);
            int adAlt = Integer.parseInt(adSArray[1]);
            int adTotal = adRef+adAlt ;
            float adFreq = 0f;
            if (adTotal > 0) {
                adFreq = (float) adAlt / adTotal;
            }
//            arrayDonor[i].put("adRef", adRef);
            arrayDonor[i].put("adAlt", adAlt);
            arrayDonor[i].put("adTotal", adTotal);
            arrayDonor[i].put("adFreq", adFreq);
//            arrayDonor[i].put("af", );
//            addNumberToJsonObject("af", af[i], arrayDonor[i], false, 'f');
//            arrayDonor[i].put("f1r2", f1r2[i*2] + "," + f1r2[i*2+1]);
//            arrayDonor[i].put("f2r1", f2r1[i*2] + "," + f2r1[i*2+1]);
//            arrayDonor[i].put("dp", );
//            addNumberToJsonObject("dp", genDP[i], arrayDonor[i], false, 'l');
//            if (Integer.parseInt(genDP[i]) != adTotal) {
//                toPrint = true;
//            }
//            arrayDonor[i].put("sb", sb[i]);
//            arrayDonor[i].put("mb", mb[i]);
            //addNumberToJsonObject("mq", mqS, arrayDonor[i], false, 'l'); //.E0
            //addNumberToJsonObject("mqRankSum", mqRankSum, arrayDonor[i], false, 'f'); //".8116E1"
            //addNumberToJsonObject("depth", dpS, arrayDonor[i], false, 'l');
//            addNumberToJsonObject("readPosRankSum", readPosRankSum, arrayDonor[i], false, 'f');
            addNumberToJsonObject("qd", qdS, arrayDonor[i], false, 'f');

            arrayDonor[i].put("specimenId", currentPatient.getSpecimenId());

            //addStrToJsonObject("dn", dn[i], arrayDonor[i], false);
            //addNumberToJsonObject("dq", dq[i], arrayDonor[i], false, 'f');

            arrayDonor[i].put("patientId", currentPatient.getPatientId());
            arrayDonor[i].put("familyId", currentPatient.getFamilyId());
            arrayDonor[i].put("relation", currentPatient.getRelation());
            arrayDonor[i].put("studyId", currentPatient.getStudyId());
            arrayDonor[i].put("practitionerId", currentPatient.getRequesterId());
            arrayDonor[i].put("organizationId", currentPatient.getReqOrgId());
            arrayDonor[i].put("sequencingStrategy", currentPatient.getSequencingStrategy());
            //arrayDonor[i].put("exomiserScore", 0f);
            familyMap.get(currentPatient.getFamilyId()).addFamily(i);

        }
        if (alleleNumber > 0) {
            alleleFrequencies = alleleCount /alleleNumber;
        }

        // transmission and familial analysis
        // Per Family - need a grouping from PED structure
        // Ped need to match donor array <---

        //familyMap.forEach((k, v) -> System.out.println(k + "\n\t" + v));
        for (Map.Entry<String, Family> entry : familyMap.entrySet()) {
            String familyId = entry.getKey();
            Family family = entry.getValue();
            int familyCompositionIndex = -1;
            List<Integer> familyCompositionArray = family.getFamilyComposition();
            int familySize = familyCompositionArray.size();

            for (Integer index : familyCompositionArray) {
                familyCompositionIndex++;
                //Pedigree currentDonor = pedigrees.get(index);
                //Pedigree currentDonor = pedigrees.get(index);
                Pedigree currentDonor  =  getPedigreeBySpecimenId(specimenList.get(index), pedigrees);
                Patient currentPatient = patientMap.get(currentDonor.getId());
                String relation = currentPatient.getRelation();
                String zygosity = (String) arrayDonor[index].get("zygosity");
                // discard 0/0 and ./. but we keep a reference for the proband of the family...
                if ("HOM REF".equalsIgnoreCase(zygosity) || "UNK".equalsIgnoreCase(zygosity)) continue;
                if ("Proband".equalsIgnoreCase(relation)) {
                    if (familySize > 1) {

                        List<String> genotypeOfAllMembers = new ArrayList<>(familySize);
                        List<Pedigree> familyPed = new ArrayList<>(familySize);
                        //genotypeOfAllMembers.add((String) arrayDonor[index].get("zygosity"));
                        genotypeOfAllMembers.add((String) arrayDonor[index].get("gt"));
                        familyPed.add(currentDonor);
                        StringBuilder genotypeFamily = new StringBuilder();
                        // put genotype of all other family member inside a field named genotypeFamily...
                        // since we discard 0/0 and ./. of all donor for clinical
                        for (int j = familyCompositionIndex + 1; j < familySize; j++) {
                            Integer otherRelationIndex = familyCompositionArray.get(j);
                            genotypeFamily.append(arrayDonor[otherRelationIndex].get("relation")).append(":").append(gt[otherRelationIndex]).append((j == familySize - 1) ? "" : ",");
                            genotypeOfAllMembers.add((String) arrayDonor[otherRelationIndex].get("gt"));
                            familyPed.add(getPedigreeBySpecimenId(specimenList.get(otherRelationIndex), pedigrees));

                            if (genotypeFamily.length() > 0) {
                                arrayDonor[index].put("genotypeFamily", genotypeFamily.toString());
                            }
                        }
                        // analysis of transmission - X & autosomal (dominant & recessif)
                        boolean strictMode = true;
                        JSONArray transmission = new JSONArray();
                        if ( !( chrPos.equalsIgnoreCase("X") || chrPos.equalsIgnoreCase("Y") ) ) {
                            if (isAutosomalDominant(genotypeOfAllMembers, familyPed, strictMode)) {
                                transmission.put("AD");
                            }
                            if (isAutosomalRecessif(genotypeOfAllMembers, familyPed, strictMode)) {
                                transmission.put("AR");
                            }
                        } else if (chrPos.equalsIgnoreCase("X")){
                            if (isXDominant(genotypeOfAllMembers, familyPed, strictMode)) {
                                transmission.put("XD");
                            }
                            if (isXRecessif(genotypeOfAllMembers, familyPed, strictMode)) {
                                transmission.put("XR");
                            }
                        }
                        String denovo = isDeNovo(genotypeOfAllMembers, familyPed, strictMode);
                        if (!denovo.equalsIgnoreCase("NO")) {
                            transmission.put(denovo);
                            //System.out.println(genotypeOfAllMembers);
                            //toPrint=true;
                        }
                        if (transmission.length() > 0) {
                            arrayDonor[index].put("transmission", transmission);
                        }
                    }
                }
                if (currentPatient.isAffected()) {
                    if (currentPatient.getQtyHposTermsFound() > 0) {
                        List<String> hposTerms = new ArrayList<>(currentPatient.getHposTerms());
                        if (currentDonor.getPhenotype().length() > 2) {
                            hposTerms.add(currentDonor.getPhenotype());
                        }
                        //arrayDonor[index].put("HpoTerms", hposTerms);
                        arrayDonor[index].put("nbHpoTerms", currentPatient.getQtyHposTermsFound());
                        //toPrint = true;
                        // setback to 0 for new variant
                        currentPatient.setQtyHposTermsFound(0);
                    }
                }
                specimenArray.put(currentDonor.getId());
                donorArray.put(arrayDonor[index]);
                //String labo = currentPatient.getLabAlias();
                //Frequencies freqenceLabo = frequenciesPerLabos.get(labo);
                //freqenceLabo.setPn(freqenceLabo.getPn() + 1f);
                //patientNb++;
            }
        }

        if (frequencies == null ) frequencies = new JSONObject();
        JSONObject freqenceInterne = new JSONObject();
        JSONArray labNames = new JSONArray();

        addNumberToJsonObject("AC", new BigDecimal(alleleCount), freqenceInterne, false, 'l');
        addNumberToJsonObject("PN", new BigDecimal(patientNb), freqenceInterne, false, 'l');
        addNumberToJsonObject("AN", new BigDecimal(alleleNumber), freqenceInterne, false, 'l');
        addNumberToJsonObject("HC", new BigDecimal(homozygoteCount), freqenceInterne, false, 'l');
        addNumberToJsonObject("AF", new BigDecimal(alleleFrequencies), freqenceInterne, false, 'f');

        frequencies.put("interne", freqenceInterne);

        JSONObject finalFrequencies = frequencies;
        frequenciesPerLabos.forEach((lab, frequence) -> {
            JSONObject freqenceLabo = new JSONObject();
            addNumberToJsonObject("AC", new BigDecimal(frequence.getAc()), freqenceLabo, false, 'l');
            addNumberToJsonObject("PN", new BigDecimal(frequence.getPn()), freqenceLabo, false, 'l');
            addNumberToJsonObject("AN", new BigDecimal(frequence.getAn()), freqenceLabo, false, 'l');
            addNumberToJsonObject("HC", new BigDecimal(frequence.getHc()), freqenceLabo, false, 'l');
            addNumberToJsonObject("AF", new BigDecimal(frequence.getAf()), freqenceLabo, false, 'f');
            finalFrequencies.put(lab, freqenceLabo);
            labNames.put(lab);
        });

        propertiesOneMutation.put("labos", labNames);

//        JSONObject frequencies = (!propertiesOneMutation.isNull("frequencies")) ? (JSONObject) propertiesOneMutation.get("frequencies") : new JSONObject();
        //
        propertiesOneMutation.put("frequencies", frequencies);

        addSetsToArrayToObjs(dbExtId, DBSNP, bdExtObj, "dbSNP", availBdExtObj);
        String clinvarids = toStringList(dbExtId.get(CLINVAR));
        if (!clinvarids.isEmpty()) {
            clinvarObj.put("clinvar_id", clinvarids);
            propertiesOneMutation.put("clinvar", clinvarObj);
        }
        addSetsToArrayToObjs(dbExtId, CLINVAR, bdExtObj, "clinvar", availBdExtObj);
//        addSetsToArrayToObj(dbExtId, ENSEMBL, bdExtObj, "ensembl", "id");
        addSetsToArrayToObjs(dbExtId, OMIM, bdExtObj, "omim", availBdExtObj);
        addSetsToArrayToObjs(dbExtId, ORPHANET, bdExtObj, "orphanet", availBdExtObj);
        addSetsToArrayToObjs(dbExtId, PUBMED, bdExtObj, "pubmed", availBdExtObj);

        String clsig = toStringList(dbExtId.get(CLINVAR_SIG));
        if (!clsig.isEmpty()) {
            clinvarObj.put( "clinvar_clinsig", clsig.split(",") ); }
        String cltraits = toStringList(dbExtId.get(CLINVAR_TRAIT));

        addSetsToArrayToObjs(dbExtId, CLINVAR_TRAIT, clinvarObj, "clinvar_trait", null);

        if (bdExtObj.length() >0) {
            propertiesOneMutation.put("bdExt", bdExtObj);
            if (!bdExtObj.isNull("omim")) {
//                JSONArray omimArray = (JSONArray) bdExtObj.get("omim");
//                if (omimArray.length() > 1) {
//                    toPrint = true;
//                }
            }
        }
        if (availBdExtObj.length()>0) {
            propertiesOneMutation.put("availableDbExt", availBdExtObj);
        }

        propertiesOneMutation.put("donors", donorArray);
        propertiesOneMutation.put("specimenList", specimenArray);
        if (geneArray.length()>0) {
            propertiesOneMutation.put("genes", geneArray);
        }

        if ( dbExt[CLINVAR]  || dbExt[OMIM] || dbExt[ORPHANET]  )//|| dbExt[DBSNP] )
            lastOne = propertiesOneMutation;


        variant.setJsonObjInString(propertiesOneMutation.toString(0));
        //variant.setGenes();
        //return propertiesOneMutation;
        return variant;
    }

    private static List<Pedigree> getFamilyPed(String familyId, List<Pedigree> pedigrees) {
        return null;
    }

    static Pedigree getPedigreeBySpecimenId(String specimen, List<Pedigree> pedigrees) {
        //System.out.println("specimen"+specimen);
        for (Pedigree ped : pedigrees) {
//            System.out.println("ped="+ped);
            if (ped.getId().equalsIgnoreCase(specimen)) {
                return ped;
            }
        }
        return null;

    }
    /*
        VEP specific annotation processing (CSQ)
     */
    private static JSONObject processVepAnnotations(String csqLine, List<Set<String>> dbExtId, boolean[] dbExt,
                                                    Set<Gene> geneSet, Map<String, Transcript> deStructuredTranscript,
                                                    String ref, String alt, Set<Clinvar> clinvarSet) {

        String[] functionalAnnotationArray = csqLine.split("[|]", -1);
        // dynamic positioning system -- pos counter++

        //System.out.print(id+"("+functionalAnnotationArray.length+")");
        //CHROM	POS	ID	REF	ALT	QUAL	FILTER	DP	MQ	MQRankSum	ReadPosRankSum	LOD	FractionInformativeReads	SNP	MNP	INS	DEL	MIXED	HOM	HET	GEN[*].GT	GEN[*].GQ	GEN[*].AD	GEN[*].AF	GEN[*].F1R2	GEN[*].F2R1	GEN[*].DP	GEN[*].SB	GEN[*].MB	GEN[*].DN	GEN[*].DQ	hiConfDeNovo	CSQ
        //CHROM	POS	ID	REF	ALT	QUAL	FILTER	DP	MQ	MQRankSum	ReadPosRankSum	LOD	FractionInformativeReads	SNP	MNP	INS	DEL	MIXED	HOM	HET	GEN[*].GT	            GEN[*].AD	GEN[*].AF	GEN[*].F1R2	GEN[*].F2R1	GEN[*].DP	GEN[*].SB	GEN[*].MB	CSQ

        int pos = 0;

        if (functionalAnnotationArray.length != EXPECTED_VEP_ANNOTATION_QTY) {
            System.out.println(">" + csqLine + "\n unexpected size, qty=" + functionalAnnotationArray.length);
            return null;
        }
        JSONObject funcAnnotation = new JSONObject();
        JSONObject frequencies = new JSONObject();
        JSONObject frequencyExAc = new JSONObject();
        JSONObject frequency1000Gp3 = new JSONObject();
        JSONObject frequencyUk10k = new JSONObject();
        JSONObject frequencyGnomadEx = new JSONObject();
        JSONObject frequencyGnomadGen = new JSONObject();
        JSONObject prediction;
        JSONObject conservation = new JSONObject();

        String cdnaChange = "";

        //    Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|

        String Allele = functionalAnnotationArray[pos++];
        String Consequence = functionalAnnotationArray[pos++];
        addStrToJsonObject("consequence", Consequence, funcAnnotation, false);
        String impact = functionalAnnotationArray[pos++];
        addStrToJsonObject("impact", impact, funcAnnotation, false);
        String symbol = functionalAnnotationArray[pos++];
        String Gene = functionalAnnotationArray[pos++];
        addStrToJsonObject("ensemblTranscriptId", Gene, funcAnnotation, false);
        String Feature_type = functionalAnnotationArray[pos++];
        addStrToJsonObject("featureType", Feature_type, funcAnnotation, false);
        String featureId = functionalAnnotationArray[pos++];
        addStrToJsonObject("featureId", featureId , funcAnnotation, false);
        String biotype = functionalAnnotationArray[pos++];
        addStrToJsonObject("biotype", biotype, funcAnnotation, false);
        String EXON = functionalAnnotationArray[pos++];
//        addStrToJsonObject("exon", functionalAnnotationArray[pos++], funcAnnotation, false);
        String INTRON = functionalAnnotationArray[pos++];
//        addStrToJsonObject("intron", functionalAnnotationArray[pos++], funcAnnotation, false);

        // N10- HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|

        String HGVSc = functionalAnnotationArray[pos++];
        String HGVSp = functionalAnnotationArray[pos++];
//        addStrToJsonObject("hgvsC", functionalAnnotationArray[pos++], funcAnnotation, false);
//        addStrToJsonObject("hgvsP", HGVSp, funcAnnotation, false);
        String cdnaPos = functionalAnnotationArray[pos++];
//        addStrToJsonObject("cdnaPos", cdnaPos , funcAnnotation, false);
        String CDS_position = functionalAnnotationArray[pos++];
//        addStrToJsonObject("cdsPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
        String Protein_position = functionalAnnotationArray[pos++];
//        addStrToJsonObject("ProteinPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
        String Amino_acids = functionalAnnotationArray[pos++];
//        addStrToJsonObject("aminoAcids", Amino_acids , funcAnnotation, false);
        String Codons = functionalAnnotationArray[pos++];
//        addStrToJsonObject("codons", functionalAnnotationArray[pos++], funcAnnotation, false);
        String Existing_variation = functionalAnnotationArray[pos++];
        String DISTANCE = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("distance", functionalAnnotationArray[pos++] , funcAnnotation, false, 'l');
        addNumberToJsonObject("strand", functionalAnnotationArray[pos++] , funcAnnotation, false, 'l'); // can be empty

        //N20       FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|
        String FLAGS = functionalAnnotationArray[pos++];
        String PICK = functionalAnnotationArray[pos++];
        addStrToJsonObject("pick", PICK, funcAnnotation, false);
        String VARIANT_CLASS = functionalAnnotationArray[pos++];
        //JSONObject variant_class = new JSONObject();
        //addStrToJsonObject("type", VARIANT_CLASS, variant_class, false);
        if (!VARIANT_CLASS.isEmpty()) {
            dbExt[TYPES] = true;
            addStrToListSet(dbExtId, TYPES, VARIANT_CLASS);
        }
        String SYMBOL_SOURCE = functionalAnnotationArray[pos++];

        // N23 HGNC_ID|CANONICAL|GIVEN_REF|USED_REF|BAM_EDIT|  DOMAINS

        String HGNC_ID = functionalAnnotationArray[pos++];
        addBooleanToJsonObject("canonical", functionalAnnotationArray[pos++], funcAnnotation, false);

        String RefSeq = functionalAnnotationArray[pos++];
        addStrToJsonObject("refSeqId", RefSeq, funcAnnotation, false);
//        String GIVEN_REF = functionalAnnotationArray[pos++];
//        String USED_REF = functionalAnnotationArray[pos++];
//        String BAM_EDIT = functionalAnnotationArray[pos++];
//        String domains = functionalAnnotationArray[pos++];
        //addStrToJsonObject("domains", domains, funcAnnotation, false);

        // N29 HGVS_OFFSET|      CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|
        // N26 HGVS_OFFSET|HGVSg|CLIN_SIG|SOMATIC|PHENO|PUBMED|

        String HGVS_OFFSET = functionalAnnotationArray[pos++];
        String HGVSg = functionalAnnotationArray[pos++];
//        addStrToJsonObject("HGVSg", HGVSg, funcAnnotation, false);
        String CLIN_SIGStr = functionalAnnotationArray[pos++];

        String SOMATIC = functionalAnnotationArray[pos++];
        //addStrToJsonObject("SOMATIC", SOMATIC, funcAnnotation, false);

        String PHENOStr = functionalAnnotationArray[pos++];
        if (!PHENOStr.isEmpty()) {
            dbExt[PHENO] = true;
            addStrToListSet(dbExtId, PHENO, PHENOStr);
        }
        String pubmed = functionalAnnotationArray[pos++];

        if (!pubmed.isEmpty()) {
            dbExt[PUBMED] = true;
            addStrToListSet(dbExtId, PUBMED, pubmed);
        }
//        String MOTIF_NAME = functionalAnnotationArray[pos++];
//        String MOTIF_POS = functionalAnnotationArray[pos++];
//        String HIGH_INF_POS = functionalAnnotationArray[pos++];
//        String MOTIF_SCORE_CHANGE = functionalAnnotationArray[pos++];

        // N38+12
        // 1000Gp3_AC|1000Gp3_AF|
        // // 1000Gp3_AFR_AC|1000Gp3_AFR_AF|1000Gp3_AMR_AC|1000Gp3_AMR_AF|1000Gp3_EAS_AC|1000Gp3_EAS_AF|1000Gp3_EUR_AC|1000Gp3_EUR_AF|1000Gp3_SAS_AC|1000Gp3_SAS_AF|
        // 1000Gp3_AC|1000Gp3_AF|

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'l');
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'f');
//        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'l');
//        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'f');
//        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'l');
//        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'f');
//        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'l');
//        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'f');
//        addNumberToJsonObject("EUR_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'l');
//        addNumberToJsonObject("EUR_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'f');
//        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'l');
//        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequency1000Gp3, true, 'f');

//        String ALSPAC_AC = functionalAnnotationArray[pos++];
//        String ALSPAC_AF = functionalAnnotationArray[pos++];
//        String APPRIS = functionalAnnotationArray[pos++];
//        String Aloft_Confidence = functionalAnnotationArray[pos++];
//        String Aloft_Fraction_transcripts_affected = functionalAnnotationArray[pos++];
//        String Aloft_pred = functionalAnnotationArray[pos++];
//        String Aloft_prob_Dominant = functionalAnnotationArray[pos++];
//        String Aloft_prob_Recessive = functionalAnnotationArray[pos++];
//        String Aloft_prob_Tolerant = functionalAnnotationArray[pos++];
//        String AltaiNeandertal = functionalAnnotationArray[pos++];
//        String Ancestral_allele = functionalAnnotationArray[pos++];

        // CADD_phred|CADD_raw_rankscore|DANN_rankscore|DEOGEN2_pred|DEOGEN2_rankscore
        // N50 - CADD_phred|CADD_raw_rankscore|DANN_rankscore|DEOGEN2_pred|DEOGEN2_rankscore|
        //                  CADD_raw_rankscore|DANN_rankscore|
//        String CADD_phred = functionalAnnotationArray[pos++];
//        String CADD_raw = functionalAnnotationArray[pos++];
        String CADD_raw_rankscore = functionalAnnotationArray[pos++];
        String DANN_rankscore = functionalAnnotationArray[pos++];
//        addStrToJsonObject("CADD_raw_rankscore", CADD_raw_rankscore, funcAnnotation, false);
//        addStrToJsonObject("DANN_rankscore", DANN_rankscore, funcAnnotation, false);
//        String DANN_score = functionalAnnotationArray[pos++];
//        String DEOGEN2_pred = functionalAnnotationArray[pos++];
//        String DEOGEN2_rankscore = functionalAnnotationArray[pos++];
//        String DEOGEN2_score = functionalAnnotationArray[pos++];
//        String Denisova = functionalAnnotationArray[pos++];
        // N55 - // ESP6500_AA_AC|ESP6500_AA_AF|
        // N55 - ESP6500_EA_AC|ESP6500_EA_AF|
//        addNumberToJsonObject("AA_AC", functionalAnnotationArray[pos++] , frequencyEsp6500, true, 'l');
//        addNumberToJsonObject("AA_AF", functionalAnnotationArray[pos++] , frequencyEsp6500, true, 'f');
//        addNumberToJsonObject("EA_AC", functionalAnnotationArray[pos++] , frequencyEsp6500, true, 'l');
//        addNumberToJsonObject("EA_AF", functionalAnnotationArray[pos++] , frequencyEsp6500, true, 'f');

        // Eigen-PC-phred_coding|Eigen-PC-raw_coding_rankscore|Eigen-pred_coding|Eigen-raw_coding_rankscore|
        // N59 - Eigen-PC-phred_coding|Eigen-PC-raw_coding_rankscore|Eigen-pred_coding|Eigen-raw_coding_rankscore|

//        String Eigen_PC_phred_coding = functionalAnnotationArray[pos++];
//        String Eigen_PC_raw_coding = functionalAnnotationArray[pos++];
//        String Eigen_PC_raw_coding_rankscore = functionalAnnotationArray[pos++];
//        String Eigen_pred_coding = functionalAnnotationArray[pos++];
//        String Eigen_raw_coding = functionalAnnotationArray[pos++];
//        String Eigen_raw_coding_rankscore = functionalAnnotationArray[pos++];

        // Ensembl_geneid | Ensembl_proteinid | Ensembl_transcriptid
        // N63 Ensembl_geneid|Ensembl_transcriptid|

        String Ensembl_geneid = functionalAnnotationArray[pos++];
//        if (Ensembl_geneid.trim().isEmpty() && Gene.startsWith("ENSG")) {
//            Ensembl_geneid = Gene;
//            biotype = "specialPseudo";
//            //System.err.println("Ensembl_geneid empty but found from Gene "+Gene);
//            //toPrint = true;
//        }

        if (!Gene.isEmpty()) {
            dbExt[ENSEMBL] = true;
            addStrToListSet(dbExtId, ENSEMBL, Gene);
            addStrToJsonObject("geneAffectedId", Gene, funcAnnotation, false);
        }
//        String Ensembl_proteinid = functionalAnnotationArray[pos++];
        //addStrToJsonObject("Ensembl_proteinid_full", Ensembl_proteinid, funcAnnotation, false);
        String Ensembl_transcriptid = functionalAnnotationArray[pos++];
        //addStrToJsonObject("Ensembl_transcriptid_full", Ensembl_transcriptid, funcAnnotation, false);

        // ExAC_AC|ExAC_AF|ExAC_AFR_AC|ExAC_AFR_AF|ExAC_AMR_AC|ExAC_AMR_AF|ExAC_Adj_AC|ExAC_Adj_AF|ExAC_EAS_AC|ExAC_EAS_AF|ExAC_FIN_AC
        //|ExAC_FIN_AF|ExAC_NFE_AC|ExAC_NFE_AF|ExAC_SAS_AC|ExAC_SAS_AF|ExAC_nonTCGA_AC|ExAC_nonTCGA_AF|ExAC_nonTCGA_AFR_AC|ExAC_nonTCGA_AFR_AF|ExAC_nonTCGA_AMR_AC|ExAC_nonTCGA_AMR_AF
        //|ExAC_nonTCGA_Adj_AC|ExAC_nonTCGA_Adj_AF|ExAC_nonTCGA_EAS_AC|ExAC_nonTCGA_EAS_AF|ExAC_nonTCGA_FIN_AC|ExAC_nonTCGA_FIN_AF|ExAC_nonTCGA_NFE_AC|ExAC_nonTCGA_NFE_AF|ExAC_nonTCGA_SAS_AC
        //|ExAC_nonTCGA_SAS_AF|ExAC_nonpsych_AC|ExAC_nonpsych_AF|ExAC_nonpsych_AFR_AC|ExAC_nonpsych_AFR_AF|ExAC_nonpsych_AMR_AC|ExAC_nonpsych_AMR_AF|ExAC_nonpsych_Adj_AC|ExAC_nonpsych_Adj_AF
        //|ExAC_nonpsych_EAS_AC|ExAC_nonpsych_EAS_AF|ExAC_nonpsych_FIN_AC|ExAC_nonpsych_FIN_AF|ExAC_nonpsych_NFE_AC|ExAC_nonpsych_NFE_AF|ExAC_nonpsych_SAS_AC|ExAC_nonpsych_SAS_AF
        //
        // N65 ExAC_AC|ExAC_AF|ExAC_nonTCGA_AC|ExAC_nonTCGA_AF|ExAC_nonpsych_AC|ExAC_nonpsych_AF|
        // ExAC_AC|ExAC_AF|

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("Adj_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("Adj_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("FIN_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("FIN_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequencyExAc, true, 'l');
//        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequencyExAc, true, 'f');
//        String ExAC_nonTCGA_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_AF = functionalAnnotationArray[pos++];

//        String ExAC_nonTCGA_AFR_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_AFR_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_AMR_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_AMR_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_Adj_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_Adj_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_EAS_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_EAS_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_FIN_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_FIN_AF = functionalAnnotationArray[pos++];

//        String ExAC_nonTCGA_NFE_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_NFE_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_SAS_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonTCGA_SAS_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_AFR_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_AFR_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_AMR_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_AMR_AF = functionalAnnotationArray[pos++];
        //
//        String ExAC_nonpsych_Adj_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_Adj_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_EAS_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_EAS_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_FIN_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_FIN_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_NFE_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_NFE_AF = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_SAS_AC = functionalAnnotationArray[pos++];
//        String ExAC_nonpsych_SAS_AF = functionalAnnotationArray[pos++];

        //N71  FATHMM_converted_rankscore|FATHMM_pred|GERP++_NR|GERP++_RS_rankscore|GTEx_V7_tissue|
        // FATHMM_converted_rankscore|FATHMM_pred|GTEx_V7_tissue|

        String FATHMM_converted_rankscore = functionalAnnotationArray[pos++];
        String FATHMM = functionalAnnotationArray[pos++];
//        String FATHMM_score = functionalAnnotationArray[pos++];
//        String GENCODE_basic = functionalAnnotationArray[pos++];
//        String GERPpp_NR = functionalAnnotationArray[pos++];
//        String GERPpp_RS = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("GERP", GERPpp_RS, conservation, true, 'f');
//        String GERPpp_RS_rankscore = functionalAnnotationArray[pos++];
//        String GM12878_confidence_value = functionalAnnotationArray[pos++];
//        String GM12878_fitCons_rankscore = functionalAnnotationArray[pos++];
//        String GM12878_fitCons_score = functionalAnnotationArray[pos++];
//        String GTEx_V7_gene = functionalAnnotationArray[pos++];
        String GTEx_V7_tissue = functionalAnnotationArray[pos++];
        //addStrToJsonObject("GTEx_V7_tissue", GTEx_V7_tissue, funcAnnotation, false);
        //N76 GenoCanyon_rankscore|Interpro_domain|LRT_converted_rankscore|LRT_pred|M-CAP_pred|M-CAP_rankscore|
        //    Interpro_domain|LRT_converted_rankscore|LRT_pred|

//        String GenoCanyon_rankscore = functionalAnnotationArray[pos++];
//        String GenoCanyon_score = functionalAnnotationArray[pos++];

        String Interpro_domain = functionalAnnotationArray[pos++];
//        addStrToJsonObject("Interpro_domain", Interpro_domain, funcAnnotation, false);
//        String LRT_Omega = functionalAnnotationArray[pos++];
        String LRT_converted_rankscore = functionalAnnotationArray[pos++];
        String LRT_pred = functionalAnnotationArray[pos++];
//        String LRT_score = functionalAnnotationArray[pos++];
//        String M_CAP_pred = functionalAnnotationArray[pos++];
//        String M_CAP_rankscore = functionalAnnotationArray[pos++];
//        String M_CAP_score = functionalAnnotationArray[pos++];

        // N80 MetaLR_pred|MetaLR_rankscore|MetaSVM_pred|MetaSVM_rankscore|MutationAssessor_pred|MutationAssessor_rankscore|
        // N86 MutationTaster_converted_rankscore|MutationTaster_pred|
        // Polyphen2_HVAR_pred|Polyphen2_HVAR_rankscore|REVEL_rankscore|

//        String MetaLR_pred = functionalAnnotationArray[pos++];
//        String MetaLR_rankscore = functionalAnnotationArray[pos++];
//        String MetaLR_score = functionalAnnotationArray[pos++].trim();
//        String MetaSVM_pred = functionalAnnotationArray[pos++].trim();
//        String MetaSVM_rankscore = functionalAnnotationArray[pos++].trim();
//        String MetaSVM_score = functionalAnnotationArray[pos++].trim();
//        String MutPred_AAchange = functionalAnnotationArray[pos++].trim();
//        String MutPred_Top5features = functionalAnnotationArray[pos++].trim();
//        String MutPred_protID = functionalAnnotationArray[pos++].trim();
//        String MutPred_rankscore = functionalAnnotationArray[pos++].trim();
//        String MutPred_score = functionalAnnotationArray[pos++].trim();
//        String MutationAssessor_pred = functionalAnnotationArray[pos++].trim();
//        String MutationAssessor_rankscore = functionalAnnotationArray[pos++].trim();
//        String MutationAssessor_score = functionalAnnotationArray[pos++].trim();
//        String MutationTaster_AAE = functionalAnnotationArray[pos++].trim();
//        String MutationTaster_converted_rankscore = functionalAnnotationArray[pos++].trim();
//        String MutationTaster_model = functionalAnnotationArray[pos++].trim();
//        String MutationTaster_pred = functionalAnnotationArray[pos++].trim();
//        String MutationTaster_score = functionalAnnotationArray[pos++].trim();

//        String PROVEAN_converted_rankscore = functionalAnnotationArray[pos++].trim();
//        String PROVEAN_pred = functionalAnnotationArray[pos++].trim();
//        String PROVEAN_score = functionalAnnotationArray[pos++].trim();
//        String Polyphen2_HDIV = functionalAnnotationArray[pos++].trim();
//        String Polyphen2_HDIV_rankscore = functionalAnnotationArray[pos++].trim();
//        String Polyphen2_HDIV_score = functionalAnnotationArray[pos++];
        String Polyphen2_HVAR_pred = functionalAnnotationArray[pos++].trim();
        String Polyphen2_HVAR_rankscore = functionalAnnotationArray[pos++].trim();
//        String Polyphen2_HVAR_score = functionalAnnotationArray[pos++];
//        String PrimateAI_pred = functionalAnnotationArray[pos++].trim();
//        String PrimateAI_rankscore = functionalAnnotationArray[pos++].trim();
//        String PrimateAI_score = functionalAnnotationArray[pos++].trim();
        String REVEL_rankscore = functionalAnnotationArray[pos++].trim();
//        addStrToJsonObject("REVEL_rankscore", REVEL_rankscore, funcAnnotation, false);
//        String REVEL_score = functionalAnnotationArray[pos++].trim();
//        String Reliability_index = functionalAnnotationArray[pos++].trim();

        //N95 SIFT4G_converted_rankscore|SIFT4G_pred|SIFT_converted_rankscore|SIFT_pred|
        // SIFT_converted_rankscore|SIFT_pred|

//        String SIFT4G_converted_rankscore = functionalAnnotationArray[pos++].trim();
//        String SIFT4G_pred = functionalAnnotationArray[pos++].trim();
//        String SIFT4G_score = functionalAnnotationArray[pos++].trim();
        String SIFT_converted_rankscore = functionalAnnotationArray[pos++].trim();
        String SIFT = functionalAnnotationArray[pos++];
//        String SIFT_score = functionalAnnotationArray[pos++].trim();

        //N99 SiPhy_29way_logOdds_rankscore|UK10K_AC|UK10K_AF|VEST4_rankscore|
        // UK10K_AC|UK10K_AF|

//        String SiPhy_29way_logOdds = functionalAnnotationArray[pos++].trim();
//        String SiPhy_29way_logOdds_rankscore = functionalAnnotationArray[pos++].trim();
//        String SiPhy_29way_pi = functionalAnnotationArray[pos++].trim();
//        String TSL = functionalAnnotationArray[pos++];
//        String TWINSUK_AC = functionalAnnotationArray[pos++];
//        String TWINSUK_AF = functionalAnnotationArray[pos++];
        String UK10K_AC = functionalAnnotationArray[pos++];
        String UK10K_AF = functionalAnnotationArray[pos++];
        addNumberToJsonObject("AC", UK10K_AC , frequencyUk10k, true, 'l');
        addNumberToJsonObject("AF", UK10K_AF , frequencyUk10k, true, 'f');
//        String Uniprot_acc = functionalAnnotationArray[pos++];
//        String Uniprot_entry = functionalAnnotationArray[pos++];
//        String VEP_canonical = functionalAnnotationArray[pos++];
//        String VEST4_rankscore = functionalAnnotationArray[pos++];

//        String VEST4_score = functionalAnnotationArray[pos++];
//        String VindijiaNeandertal = functionalAnnotationArray[pos++];
//        String aaalt = functionalAnnotationArray[pos++];
//        addStrToJsonObject("aaAlt", aaalt, funcAnnotation, false);
//        String aapos = functionalAnnotationArray[pos++];
//        addStrToJsonObject("aaPos", functionalAnnotationArray[pos++] , funcAnnotation, false);
//        String aaref = functionalAnnotationArray[pos++];
//        addStrToJsonObject("aaRef", aaref, funcAnnotation, false);
//        String alt = functionalAnnotationArray[pos++];
//        addStrToJsonObject("alt", functionalAnnotationArray[pos++], funcAnnotation, false);

//        String bStatistic = functionalAnnotationArray[pos++];
//        String bStatistic_rankscore = functionalAnnotationArray[pos++];
//        String cds_strand = functionalAnnotationArray[pos++];
//        String chr = functionalAnnotationArray[pos++];

        //N103 clinvar_MedGen_id|clinvar_OMIM_id|clinvar_Orphanet_id|clinvar_clnsig|clinvar_id|clinvar_review|clinvar_trait|
        // clinvar_MedGen_id|clinvar_OMIM_id|clinvar_Orphanet_id|clinvar_clnsig|clinvar_id|clinvar_trait|

        String clinvar_MedGen_id = functionalAnnotationArray[pos++];
//        addStrToJsonObject("clinvar_MedGen_id", clinvar_MedGen_id, funcAnnotation, false);
        String omim = functionalAnnotationArray[pos++];
        if (addStrToJsonObject("clinvar_OMIM_id", omim, funcAnnotation, false)){
            dbExt[OMIM] = true;
            addStrToListSet(dbExtId, OMIM, omim);
        }
        String orpha = functionalAnnotationArray[pos++];
        if (addStrToJsonObject("clinvar_Orphanet_id", orpha, funcAnnotation, false)){
            dbExt[ORPHANET] = true;
            addStrToListSet(dbExtId, ORPHANET, orpha);
        }
        String clinvarClinsign = functionalAnnotationArray[pos++];
//        addStrToJsonObject("clinvar_clnsig", , funcAnnotation, false);
        if (!clinvarClinsign.isEmpty()) {
            dbExt[CLINVAR_SIG] = true;
            addStrToListSet(dbExtId, CLINVAR_SIG, clinvarClinsign);
        }
//        String clinvar_hgvs = functionalAnnotationArray[pos++];
//        addStrToJsonObject("clinvar_hgvs", , funcAnnotation, false);
        String clinvar = functionalAnnotationArray[pos++];
        if ( !clinvar.isEmpty()){
            dbExt[CLINVAR] = true;
            addStrToListSet(dbExtId, CLINVAR, clinvar);
        }
//        String clinvar_review = functionalAnnotationArray[pos++];
        String clinvarTrait = functionalAnnotationArray[pos++];
        //addStrToJsonObject("clinvar_trait", , funcAnnotation, false);
        if (!clinvarTrait.isEmpty()) {
            dbExt[CLINVAR_TRAIT] = true;
            addStrToListSet(dbExtId, CLINVAR_TRAIT, clinvarTrait);
        }
        /*
        Refactor of the clinvar structure WIP - discarded for now
        String[] clinvarIds = clinvar.split("&");
        String[] omimIds = omim.split("&");
        String[] clinvarClinsignArray = clinvarClinsign.split("&");
        String[] clinvarTraitArray = clinvarTrait.split("&");
        for (int i=0; i<clinvarIds.length; i++) {

            Clinvar clinvar1 = new Clinvar();
            ClinvarInterpretation ci = new ClinvarInterpretation();
            clinvar1.setId(clinvarIds[i]);
            ci.setCondition(clinvarTraitArray[i]);
            ci.setInterpretation(clinvarClinsignArray[i]);
            if (omimIds.length == clinvarIds.length) {
                clinvar1.addOmim(omimIds[i]);
            }
            clinvar1.addInterpretation(ci);
            clinvarSet.add(clinvar1);

        }
        */


        //N110 fathmm-MKL_coding_pred|fathmm-MKL_coding_rankscore|fathmm-XF_coding_pred|fathmm-XF_coding_rankscore|

//        String clinvar_var_source = functionalAnnotationArray[pos++];
//        String codon_degeneracy = functionalAnnotationArray[pos++];
//        String codonpos = functionalAnnotationArray[pos++];
//        String fathmm_MKL_coding_group = functionalAnnotationArray[pos++];
//        String fathmm_MKL_coding_pred = functionalAnnotationArray[pos++];
//        String fathmm_MKL_coding_rankscore = functionalAnnotationArray[pos++];
//        String fathmm_MKL_coding_score = functionalAnnotationArray[pos++];
//        String fathmm_XF_coding_pred = functionalAnnotationArray[pos++];
//        String fathmm_XF_coding_rankscore = functionalAnnotationArray[pos++];
//        String fathmm_XF_coding_score = functionalAnnotationArray[pos++];

        //String geneName = functionalAnnotationArray[pos++];
        if (addStrToJsonObject("geneAffectedSymbol", symbol, funcAnnotation, false)) {
            dbExt[GENES] = true;
            addStrToListSet(dbExtId, GENES, symbol);
        }
        //N114 gnomAD_exomes_controls_AC|gnomAD_exomes_controls_AF|gnomAD_exomes_controls_AN|
        // gnomAD_exomes_AC|gnomAD_exomes_AF|gnomAD_genomes_AC|gnomAD_genomes_AF|

        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');

//        String gnomAD_exomes_AFR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_AFR_nhomalt = functionalAnnotationArray[pos++];

//        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_AMR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_AMR_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_AN = functionalAnnotationArray[pos++];

//        addNumberToJsonObject("ASJ_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("ASJ_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_ASJ_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_ASJ_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_EAS_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_EAS_nhomalt = functionalAnnotationArray[pos++];

//        addNumberToJsonObject("FIN_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("FIN_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_FIN_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_FIN_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_NFE_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_NFE_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("POPMAX_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("POPMAX_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_POPMAX_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_POPMAX_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("SAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'l');
//        addNumberToJsonObject("SAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_SAS_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_SAS_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AC = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("controls_AC", gnomAD_exomes_controls_AC , frequencyGnomadEx, true, 'l');
//        String gnomAD_exomes_controls_AF = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("controls_AF", gnomAD_exomes_controls_AF , frequencyGnomadEx, true, 'f');
//        String gnomAD_exomes_controls_AFR_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AFR_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AFR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AFR_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AMR_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AMR_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AMR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AMR_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_AN = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("controls_AN", gnomAD_exomes_controls_AN , frequencyGnomadEx, true, 'l');
//        String gnomAD_exomes_controls_ASJ_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_ASJ_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_ASJ_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_ASJ_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_EAS_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_EAS_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_EAS_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_EAS_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_FIN_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_FIN_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_FIN_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_FIN_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_NFE_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_NFE_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_NFE_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_NFE_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_POPMAX_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_POPMAX_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_POPMAX_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_POPMAX_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_SAS_AC = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_SAS_AF = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_SAS_AN = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_SAS_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_controls_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_flag = functionalAnnotationArray[pos++];
//        String gnomAD_exomes_nhomalt = functionalAnnotationArray[pos++];
        addNumberToJsonObject("AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'l');
        addNumberToJsonObject("AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'f');
//        addNumberToJsonObject("AFR_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'l');
//        addNumberToJsonObject("AFR_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'f');
//        String gnomAD_genomes_AFR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_AFR_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("AMR_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'l');
//        addNumberToJsonObject("AMR_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'f');
        //N117  gnomAD_genomes_controls_AC|gnomAD_genomes_controls_AF|gnomAD_genomes_controls_AN|
//        String gnomAD_genomes_AMR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_AMR_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_ASJ_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_ASJ_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_ASJ_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_ASJ_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("EAS_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'l');
//        addNumberToJsonObject("EAS_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'f');
//        String gnomAD_genomes_EAS_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_EAS_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_FIN_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_FIN_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_FIN_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_FIN_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("NFE_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'l');
//        addNumberToJsonObject("NFE_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'f');
//        String gnomAD_genomes_NFE_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_NFE_nhomalt = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("POPMAX_AC", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'l');
//        addNumberToJsonObject("POPMAX_AF", functionalAnnotationArray[pos++] , frequencyGnomadGen, true, 'f');
//        String gnomAD_genomes_POPMAX_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_POPMAX_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AF = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("controls_AC", gnomAD_genomes_controls_AC , frequencyGnomadGen, true, 'l');
//        addNumberToJsonObject("controls_AF", gnomAD_genomes_controls_AF , frequencyGnomadGen, true, 'f');
//        String gnomAD_genomes_controls_AFR_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AFR_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AFR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AFR_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AMR_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AMR_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AMR_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AMR_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_AN = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("controls_AN",  gnomAD_genomes_controls_AN, frequencyGnomadGen, true, 'l');
//        String gnomAD_genomes_controls_ASJ_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_ASJ_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_ASJ_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_ASJ_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_EAS_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_EAS_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_EAS_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_EAS_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_FIN_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_FIN_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_FIN_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_FIN_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_NFE_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_NFE_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_NFE_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_NFE_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_POPMAX_AC = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_POPMAX_AF = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_POPMAX_AN = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_POPMAX_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_controls_nhomalt = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_flag = functionalAnnotationArray[pos++];
//        String gnomAD_genomes_nhomalt = functionalAnnotationArray[pos++];

        //N120 integrated_fitCons_rankscore|phastCons100way_vertebrate_rankscore|phastCons17way_primate_rankscore|phastCons30way_mammalian_rankscore|
        // phyloP17way_primate_rankscore|rs_dbSNP151
//        String hg18_chr = functionalAnnotationArray[pos++];
//        String hg18_pos_1based = functionalAnnotationArray[pos++];
//        String hg19_chr = functionalAnnotationArray[pos++];
//        String hg19_pos_1based = functionalAnnotationArray[pos++];
//        String integrated_confidence_value = functionalAnnotationArray[pos++];
//        String integrated_fitCons_rankscore = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("fitCons", integrated_fitCons_rankscore, conservation, true, 'f');
//        String integrated_fitCons_score = functionalAnnotationArray[pos++];
//        String phastCons100way_vertebrate = functionalAnnotationArray[pos++];
//        String phastCons100way_vertebrate_rankscore = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("PhastCons100way", phastCons100way_vertebrate_rankscore, conservation, true, 'f');
//        String phastCons17way_primate = functionalAnnotationArray[pos++];
//        String phastCons17way_primate_rankscore = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("PhastCons17way", phastCons17way_primate_rankscore, conservation, true, 'f');
//        String phastCons30way_mammalian = functionalAnnotationArray[pos++];
//        String phastCons30way_mammalian_rankscore = functionalAnnotationArray[pos++];

        //N124 phyloP100way_vertebrate_rankscore|phyloP17way_primate_rankscore|phyloP30way_mammalian_rankscore|rs_dbSNP151"
//        String phyloP100way_vertebrate = functionalAnnotationArray[pos++];
//        String phyloP100way_vertebrate_rankscore = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("PhyloP100Way", phyloP100way_vertebrate_rankscore, conservation, true, 'f');
//        String phyloP17way_primate = functionalAnnotationArray[pos++];
        String phyloP17way_primate_rankscore = functionalAnnotationArray[pos++];
        addNumberToJsonObject("PhyloP17Way", phyloP17way_primate_rankscore, conservation, true, 'f');
//        String phyloP30way_mammalian = functionalAnnotationArray[pos++];
//        String phyloP30way_mammalian_rankscore = functionalAnnotationArray[pos++];
//        addNumberToJsonObject("PhyloP30Way", phyloP30way_mammalian_rankscore, conservation, true, 'f');
//        String pos_1_based = functionalAnnotationArray[pos++];
//        String ref = functionalAnnotationArray[pos++];
//        addStrToJsonObject("ref", ref, funcAnnotation, false);
//        String refcodon = functionalAnnotationArray[pos++];
//        addStrToJsonObject("ref_codon", functionalAnnotationArray[pos++], funcAnnotation, false);

        String rs_dbSNP151 = functionalAnnotationArray[pos];
        if (!rs_dbSNP151.isEmpty()) {
            dbExt[DBSNP] = true;
            addStrToListSet(dbExtId, DBSNP, rs_dbSNP151);
        }

        // N128
        /* population frequencies Json Obj */
        // Create a map of Map<TranscriptIdFull, deStructuredTranscript>

        if (deStructuredTranscript.isEmpty() && !Ensembl_transcriptid.trim().isEmpty()) {

            String[] ensemblTranscriptIdArray = Ensembl_transcriptid.split("&");
            int length = ensemblTranscriptIdArray.length;
//            String[] ensemblProteinIdArray = Ensembl_proteinid.split("&");
//            String[] aaPosArray = aapos.split("&");
            String[] FATHMM_scoreArray = (FATHMM_converted_rankscore.isEmpty()) ? fillArray(length) : FATHMM_converted_rankscore.split("&");
            String[] FATHMMArray = (FATHMM.isEmpty()) ? fillArray(length) : FATHMM.split("&");
            String[] SIFTArray = SIFT.split("&");
            String[] SIFT_scoreArray = SIFT_converted_rankscore.split("&");
//            String[] MutationAssessor_predArray = MutationAssessor_pred.split("&");
//            String[] MutationAssessor_scoreArray = MutationAssessor_rankscore.split("&");
//            String[] Polyphen2_HDIVArray = Polyphen2_HDIV.split("&");
//            String[] Polyphen2_HDIV_scoreArray = Polyphen2_HDIV_rankscore.split("&");
            String[] Polyphen2_HVAR_predArray = Polyphen2_HVAR_pred.split("&");
            String[] Polyphen2_HVAR_scoreArray = Polyphen2_HVAR_rankscore.split("&");

            //String[] LRT_predArray = LRT_pred.split("&");
            //String[] LRT_scoreArray = LRT_score.split("&");


            for (int i=0; i < ensemblTranscriptIdArray.length; i++) {
                String tid = ensemblTranscriptIdArray[i];
                int scoreLength = SIFT_scoreArray.length;
                int scorePos = (scoreLength > 1) ? i : 0;
                try {
                    /*
                      String ensemblTranscriptId, String aaPosition, String FATHMM, String SIFT, String SIFT_score, String ensemblProteinId,
                      String polyphen2_HVAR_score, String mutationAssessor_score, String polyphen2_HDIV,
                      String polyphen2_HDIV_score, String LRT_pred, String FATHMM_score, String mutationAssessor_pred,
                      String polyphen2_HVAR_pred, String LRT_score,
                      String DANN, String CADD, String REVEL) {
                     */
                    deStructuredTranscript.put(tid, new Transcript(
                            tid, null, FATHMM_PRED_VALUE.get(FATHMMArray[i]), SIFT_PRED_VALUE.get(SIFTArray[i]),
                            SIFT_scoreArray[scorePos], null,
                            Polyphen2_HVAR_scoreArray[scorePos], null, null,
                            null, LRT_PRED_VALUE.get(LRT_pred), FATHMM_scoreArray[scorePos], null,
                            POLYPHEN2_HVAR_PRED_VALUE.get(Polyphen2_HVAR_predArray[i]), LRT_converted_rankscore,
                            DANN_rankscore, CADD_raw_rankscore, REVEL_rankscore));
                } catch (ArrayIndexOutOfBoundsException aioe) {
                    System.err.println("something can be empty");
                    aioe.printStackTrace();
                    System.exit(1);
                }
            }
        }

        if (deStructuredTranscript.containsKey(featureId)) {
            Transcript transcript = deStructuredTranscript.get(featureId);
            //String aaChange = aaref + transcript.getAaPosition() + aaalt;
            String aaChange = null;
            if (Amino_acids !=null && !Amino_acids.isEmpty() && !Amino_acids.startsWith(".")) {
              aaChange = Amino_acids.substring(0,1) + Protein_position + Amino_acids.substring(2);
            }
            //String aaChange = Amino_acids
            addStrToJsonObject("aaChange", aaChange , funcAnnotation, false);
            //addStrToJsonObject("domain", transcript.getDomains() , funcAnnotation, false);

            prediction = transcript.getPrediction();
            if ( (boolean) prediction.get("available") ) {
                prediction.remove("available");
                funcAnnotation.put("predictions", prediction);
            }

        }

        // wrong way of doing it
//        if (!cdnaPos.isEmpty() && !".".startsWith(cdnaPos)) {
//            cdnaChange = ref + cdnaPos + alt;
//            //System.out.println("cdnaChange="+cdnaChange);
//        }
        if (!CDS_position.isEmpty() && !".".startsWith(CDS_position)) {
            cdnaChange = CDS_position + ref + ">" + alt;
            //System.out.println("cdnaChange="+cdnaChange);
        }
        addStrToJsonObject("cdnaChange", cdnaChange , funcAnnotation, false);
        boolean freqAvail = false;
        if ((boolean) frequency1000Gp3.get("available")) {
            frequency1000Gp3.remove("available");
            frequencies.put("1000Gp3", frequency1000Gp3);
            freqAvail = true;
        }
        if ((boolean) frequencyExAc.get("available")) {
            frequencyExAc.remove("available");
            frequencies.put("ExAc", frequencyExAc);
            freqAvail = true;
        }
        if (!frequencyGnomadEx.isNull("available") && (boolean) frequencyGnomadEx.get("available")) {
            frequencyGnomadEx.remove("available");
            frequencies.put("gnomAD_exomes", frequencyGnomadEx);
            freqAvail = true;
        }
        if (!frequencyGnomadGen.isNull("available") && (boolean) frequencyGnomadGen.get("available")) {
            frequencyGnomadGen.remove("available");
            frequencies.put("gnomAD_genomes", frequencyGnomadGen);
            freqAvail = true;
        }
        if (!frequencyUk10k.isNull("available") && (boolean) frequencyUk10k.get("available")) {
            frequencyUk10k.remove("available");
            frequencies.put("Uk10k", frequencyUk10k);
            freqAvail = true;
        }
        if (freqAvail) funcAnnotation.put("frequencies", frequencies);
        funcAnnotation.put("conservationsScores", conservation);

        if (!Gene.isEmpty()) {
            geneSet.add(new Gene(Gene, symbol, biotype));
        }

        return funcAnnotation;

    }

    static boolean addNumberToJsonObject(String name, String var, JSONObject jsonObject, boolean withAvailability, char type) {
        boolean avail = false;
        if (var != null && !var.trim().isEmpty() && !".".equalsIgnoreCase(var)) {
            try {
                if (type == 'l') {
                    jsonObject.put(name, new BigDecimal(var).longValue());
                } else {
                    jsonObject.put(name, new BigDecimal(var).doubleValue());
                }
                avail = true;
                if (withAvailability) jsonObject.put("available", avail);
            } catch (NumberFormatException nfe) {
                nfe.printStackTrace();
                System.out.println("trying to store number name="+ name + " with var="+var);
                toPrint = true;
                if (withAvailability) jsonObject.put("available", avail);
            }
        } else if (withAvailability) {
            jsonObject.put("available", avail);
        }
        return avail;
    }

    static boolean addNumberToJsonObject(String name, BigDecimal var, JSONObject jsonObject, boolean withAvailability, char type) {
        boolean avail = false;
        if (var != null) {
            try {
                if (type == 'l') {
                    jsonObject.put(name, var.longValue());
                } else {
                    jsonObject.put(name, var.doubleValue());
                }
                avail = true;
                if (withAvailability) jsonObject.put("available", avail);
            } catch (NumberFormatException nfe) {
                nfe.printStackTrace();
                System.out.println("trying to store number name="+ name + " with var="+var);
                toPrint = true;
                if (withAvailability) jsonObject.put("available", avail);
            }
        } else if (withAvailability) {
            jsonObject.put("available", avail);

        }
        return avail;
    }


    private static boolean addBooleanToJsonObject(String name, String var, JSONObject jsonObject, boolean withAvailability) {
        boolean avail = false;
        if (var.length()>0 && !".".equalsIgnoreCase(var)) {
            if ("YES".equalsIgnoreCase(var)) {
                jsonObject.put(name, true);
                avail = true;
            } else {
                jsonObject.put(name, false);
                avail = true;
            }
            if (withAvailability) jsonObject.put("available", avail);
        } else if (withAvailability) {
            jsonObject.put("available", avail);
            jsonObject.put(name, false);
        } else {
            jsonObject.put(name, false);
        }
        return avail;
    }

    static boolean addStrToJsonObject(String name, String var, JSONObject jsonObject, boolean withAvailability) {
        boolean avail = false;
        if (var != null && var.length()>0 && !".".equalsIgnoreCase(var)) {
            jsonObject.put(name, var);
            avail = true;
            if (withAvailability) jsonObject.put("available", avail);
        } else if (withAvailability) {
            jsonObject.put("available", avail);
        }
        return avail;
    }

    public static String zygosity(String gt) {
        char[] gtA = gt.toCharArray();
        if (gtA[0] == '.' || gtA[2] == '.') {
            return "UNK";
        }
        if (gtA.length == 3) {
            // Homo
            if ( gtA[0] == gtA[2] ) {
                if ( gtA[0] == '0') {
                    return "HOM REF";
                } else {
                    return "HOM";
                }
            // not equal
            } else {
                return "HET";
            }
        }
        return "UNK";
    }


    public static int countAlternativeAllele(String gt) {
        int count =0;
        char[] gtA = gt.toCharArray();

        if ( gtA[0] == '1') {
            count++;
        }
        if ( gtA[2] == '1') {
            count++;
        }
        return count;
    }
    public static int countNumberOfAllele(String gt) {
        int count =0;
        char[] gtA = gt.toCharArray();
        if ( gtA[0] == '1' || gtA[0] == '0') {
            count++;
        }
        if ( gtA[2] == '1'|| gtA[2] == '0') {
            count++;
        }
        return count;
    }

    // public bc it's used by scala test
    public static int getImpactScore(String impact) {
        if (impact == null) return 0;
        switch (impact.toUpperCase()) {
            case "HIGH" : return 4;
            case "MODERATE" : return 3;
            case "LOW" : return 2;
            case "MODIFIER" : return 1;
            default: return 0;
        }
    }

    public static List<String> getSpecimenList(String filename) {
        List<String> specimenList = new ArrayList<>();

        if (filename != null && !filename.isEmpty()) {
            try (BufferedReader buf = new BufferedReader(new FileReader(filename))) {

                String[] firstLineArray = buf.readLine().split("\t");
                int i=0;
                for (; i<firstLineArray.length;i++) {
                    if (firstLineArray[i].equalsIgnoreCase("FORMAT")) {
                       break;
                    }
                }
                i++;
                for (; i<firstLineArray.length;i++) {
                    if (firstLineArray[i].startsWith("SP")) {
                        specimenList.add(firstLineArray[i]);
                    } else {
                        specimenList.add("SP" + firstLineArray[i]);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return specimenList;
    }

    public static Properties getPropertiesFromFile(String filename) {
        Properties prop = new Properties();

        try (BufferedReader buf = new BufferedReader(new FileReader(filename))) {
            prop.load(buf);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return prop;
    }

    private static void addStrToListSet(List<Set<String>> setList, int position, String value) {

        // Extract the proper set
        Set<String> extractedSet = setList.get(position);
        String[] splitStr = value.split("&");
        for (int i=0; i<splitStr.length;i++) {
            if (splitStr[i].startsWith("_")) {
                splitStr[i] = splitStr[i].replaceFirst("_", "");
            }
        }
        Collections.addAll(extractedSet, splitStr);
    }

    private static String toStringList(Set<String> setStr) {
        return setStr.stream().reduce("", (x,y) -> (x.trim().isEmpty() ? x : x + ",") + y);
    }


    private static void addSetsToArrayToObjs(List<Set<String>> setList, int position, JSONObject jsonObject, String name, JSONObject availObj) {
        JSONArray newIdArray = new JSONArray();
        for (String id: setList.get(position)) {
            newIdArray.put(id);
        }
        if (newIdArray.length() > 0 ) {
            jsonObject.put(name, newIdArray);
            if (availObj != null) availObj.put(name, true);
        }
    }

    static void addGeneSetsToObjs(String ensId, JSONObject jsonObject, JSONObject availObj,
                                  Map<String, Patient> patientMap, List<Pedigree> pedigrees) {

        Set<String> geneSets = getMembersForEnsId(ensId);
        Set<String> aliasSet = new HashSet<>();
        JSONArray hpoGeneSets = new JSONArray();
        JSONArray orphanetGeneSets = new JSONArray();
        JSONArray radboudumc = new JSONArray();
        JSONArray omimSets = new JSONArray();

        geneSets.forEach(member -> {

                if (member.startsWith("symbol:")) {
                    // VEP geneSymbol sometimes put an alias instead
                    String symbol = member.replace("symbol:", "");
                    if (!jsonObject.isNull("geneSymbol")) {
                        String geneSymbol = jsonObject.getString("geneSymbol");
                        if (!symbol.equalsIgnoreCase(geneSymbol)) {
                            aliasSet.add(geneSymbol);
                            jsonObject.put("geneSymbol", symbol);
                            jsonObject.put("oldGeneSymbol", geneSymbol);
                        }
                    } else {
                        jsonObject.put("geneSymbol", symbol);
                    }
                } else if (member.startsWith("HP:")) {
                    String[] hpoTerms = member.split(",");
                    hpoGeneSets.put(hpoTerms[1] + " (" + hpoTerms[0] + ")");
                    for (Pedigree ped : pedigrees) {
                        String phenotype = ped.getPhenotype();

                        // -9 (missing), 0 (missing) ,1 (unaffected)
                        // 2 is affected or could have 1 single hpo term

                        if (!("-9".equalsIgnoreCase(phenotype) || "0".equalsIgnoreCase(phenotype) || "1".equalsIgnoreCase(phenotype))) {
                            if (phenotype.equalsIgnoreCase(hpoTerms[0])) {
                                patientMap.get(ped.getId()).addQtyOfHposTermsFound(1);
                            }
                            for (String hpo : patientMap.get(ped.getId()).getHposTerms()) {
                                if (hpo.equalsIgnoreCase(hpoTerms[0])) {
                                    patientMap.get(ped.getId()).addQtyOfHposTermsFound(1);
                                }
                            }
                        }
                    }
//                if (hpoTermsNeg.contains(hpoTerms[0])) negPosHposTermsFoundArray[0]++;
//                if (hpoTermsPos.contains(hpoTerms[0])) negPosHposTermsFoundArray[1]++;

                } else if (member.startsWith("Orph:")) {
                    String[] orphanetObj = member.split(";");

                    String[] orphTerms = orphanetObj[0].split(",");
                    JSONObject orphanetJsonObj = new JSONObject();
                    orphanetJsonObj.put("dataId", orphanetObj[1].replace("did:", ""));
                    orphanetJsonObj.put("panel", orphTerms[1] + " (" + orphTerms[0] + ")");
                    orphanetGeneSets.put(orphanetJsonObj);
                } else if (member.startsWith("alias:")) {
                    //alias.put(member.replace("alias:", ""));
                    aliasSet.add(member.replace("alias:", ""));
                } else if (member.startsWith("Rad:")) {
                    String[] radTerms = member.split(":");
                    radboudumc.put(radTerms[1] + " (" + member + ")");
                } else if (member.startsWith("geneid:")) {
                    jsonObject.put("geneId", member.replace("geneid:", ""));
                } else if (member.startsWith("map_location:")) {
                    jsonObject.put("location", member.replace("map_location:", ""));
                } else if (member.startsWith("omim:")) {
                    omimSets.put(member.replace("omim:", ""));
                } else if (member.startsWith("hgnc:")) {
                    jsonObject.put("hgnc", member.replace("hgnc:", ""));
                } else if (member.startsWith("name:")) {
                    jsonObject.put("name", member.replace("name:", ""));
                }

            });
        if (hpoGeneSets.length() > 0) {
            jsonObject.put("hpo", hpoGeneSets);
            availObj.put("hpo", true);
        }
        if (orphanetGeneSets.length() > 0) {
            jsonObject.put("orphanet", orphanetGeneSets);
            availObj.put("orphanet", true);
        }
        if (aliasSet.size() > 0) jsonObject.put("alias", aliasSet.toArray());
        if (omimSets.length() > 0) {
            jsonObject.put("omim", omimSets);
        }
        if (radboudumc.length() > 0) {
            jsonObject.put("radboudumc", radboudumc);
            availObj.put("radboudumc", true);
        }

    }

    private static String[] fillArray(int size) {
        String[] array = new String[size];
        Arrays.fill(array, "");
        return array;
    }

    static List<JSONObject> extractGenesFromMutation(JSONObject mutationWithGene, String parentId, boolean andRemove) {
        List<JSONObject> listGenes = new ArrayList<>();

        mutationWithGene.put("mutation_or_gene", "mutation");

        JSONObject joinObj = new JSONObject();
        joinObj.put("name", "gene");
        joinObj.put("parent", parentId);
        if (!mutationWithGene.isNull("genes")) {
            JSONArray geneArray = (JSONArray) (andRemove ? mutationWithGene.remove("genes") : mutationWithGene.get("genes"));
            geneArray.forEach((gene) -> {
                ((JSONObject) gene).put("mutation_or_gene", joinObj);
                ((JSONObject) gene).put("id", parentId);
                listGenes.add((JSONObject) gene);
            });
        }
        return listGenes;
    }

    static Map<String,Family> getFamilyMap(List<String> specimenList, List<Pedigree> pedigrees) {
        Map<String,Family> familyMap = new HashMap<>();

        specimenList.forEach((specimen)-> {
            Pedigree ped = getPedigreeBySpecimenId(specimen, pedigrees);
            if (ped == null) {
                // an error occurs...
                System.err.println("!!! All specimen from the vcf that was found in the patient DB need an entry in the ped file");
                System.exit(1);
            }
            Family family = new Family(ped.getFamilyId());
            familyMap.put(ped.getFamilyId(), family);
        });
        return familyMap;
    }

    public static boolean isAutosomalDominant(List<String> genotypesFamily, List<Pedigree> familyPed, boolean strict) {
        // Gemini
        // All affected must be HET or UNK
        // [affected] No unaffected can be het or homalt (can be unknown)
        // At least 1 affected must have 1 affected parent (or have no parents). <- nope

        // Jannovar
        // at least one affected person has a HET call for this variant,  affected proband is HET
        // no affected person has a REF or HOM call, and  proband cannot be REF or HOM
        // no unaffected person has a HET or HOM call.  parent should not be HET or HOM

        boolean hasParent = false, hasNoParent = false;
        int nbAffected = 0;
        if (familyPed.size() > 1) {
            hasParent = true;
        } else {
            hasNoParent = true;
        }
        int qtyIndividualToTest =1;
        if (strict) {
            qtyIndividualToTest = familyPed.size();
        }
        for (int i=0; i<qtyIndividualToTest; i++) {
            boolean isAffected;
            if (strict) {
                isAffected = (!familyPed.get(i).getPhenotype().equalsIgnoreCase("1"));
            } else {
                isAffected = (i==0 && !familyPed.get(i).getPhenotype().equalsIgnoreCase("1"));
            }
            String zygosity = zygosity(genotypesFamily.get(i));
            if (isAffected) {
                nbAffected++;
                if (! (zygosity.equalsIgnoreCase("HET"))) {
                    if (i>0 && zygosity.equalsIgnoreCase("UNK")) {
                        continue;
                    } else {
                        return false;
                    }
                }
            } else if (zygosity.equalsIgnoreCase("HET") || zygosity.equalsIgnoreCase("HOM") ) {
                return false;
            }
        }
        return (hasNoParent) || (hasParent && nbAffected >= 1);
    }

    public static boolean isAutosomalRecessif(List<String> genotypesFamily, List<Pedigree> familyPed, boolean strict) {
        // Gemini
        // all affecteds must be hom_alt
        // [affected] no unaffected can be hom_alt (can be unknown)
        // [strict] if parents exist they must be unaffected and het for all affected kids

        // Janovar
        // at least one affected person has a HOM call for this variant and
        //  no affected person has a REF or HET call. xx not considered
        // The unaffected parents of affected persons must not be REF or HOM.
        // There is no unaffected person that has a HOM call.

        boolean hasParent = false;
        int nbAffected = 0;

        if (familyPed.size() > 1) {
            hasParent = true;
        }
        int qtyIndividualToTest =1;
        if (strict) {
            qtyIndividualToTest = familyPed.size();
        }

        for (int i=0; i<qtyIndividualToTest; i++) {
            String zygosity = zygosity(genotypesFamily.get(i));
            boolean isAffected;
            if (strict) {
                isAffected = (!familyPed.get(i).getPhenotype().equalsIgnoreCase("1"));
            } else {
                isAffected = (i==0 && !familyPed.get(i).getPhenotype().equalsIgnoreCase("1"));
            }
            if (isAffected) {
                nbAffected++;
                if  (!zygosity.equalsIgnoreCase("HOM"))   {
                    return false;
                }
            } else {
                if (zygosity.startsWith("HOM") ) {
                    return false;
                }
            }
        }
        return !hasParent || nbAffected <= 1;
    }

    public static String isDeNovo(List<String> genotypesFamily, List<Pedigree> familyPed, boolean strict) {

        // Gemini
        // All affected must be HET but they could be ./.
        // [affected] all unaffected must be homref or homalt or unk (could)
        // at least 1 affected kid must have unaffected parents

        boolean hasParent = false;
        boolean isPossible = false;

        int nbAffected = 0;
        if (familyPed.size() > 1) {
            hasParent = true;
        } else {
            return "NO";
        }
        if (familyPed.get(0).getPhenotype().equalsIgnoreCase("1")) {
            return "NO";
        }

        for (int i=0; i<familyPed.size(); i++) {
            boolean isAffected;
            if (strict) {
                isAffected = (!familyPed.get(i).getPhenotype().equalsIgnoreCase("1"));
            } else {
                isAffected = (i==0 && !familyPed.get(i).getPhenotype().equalsIgnoreCase("1"));
            }
            String zygosity = zygosity(genotypesFamily.get(i));

            if (isAffected) {
                nbAffected++;
                if (zygosity.equalsIgnoreCase("UNK")) {
                    isPossible = true;
                    nbAffected--; //unconfirmed affection at this locus
                    continue;
                }
                if (! (zygosity.equalsIgnoreCase("HET"))) {
                    return "NO";
                }
            } else if (zygosity.equalsIgnoreCase("UNK")) {
                isPossible = true;
            } else if (!zygosity.equalsIgnoreCase("HOM REF") ) {
                return "NO";
            }
        }
        if ( (strict && hasParent && nbAffected > 1) ) {
            return "NO";
        }
        if (isPossible) {
            return "DeNovo";
        } else {
            return "DeNovo";
        }
    }

    public static boolean isXRecessif(List<String> genotypesFamily, List<Pedigree> familyPed, boolean strict) {

        // Gemini
        // Affected females must be HOM_ALT
        // Unaffected females are HET or HOM_REF
        // Affected males are not HOM_REF
        // Unaffected males are HOM_REF

        // Janovar (new version)
        // First of all variants must be X-Chromosomal.
        //        Female proband must be HOM
        //        Male proband must be HET or HOM
        //
        //        For the parents of affected females:
            //	•	the father must be affected and
            //	•	the mother cannot have it REF or HOM
        //        For the parents of affected males:
            //	•	the unaffected father cannot have the variant HET or HOM
            //	•	the affected father must be  HET or HOM (can be UNK)
            //	•	the unaffected mother cannot be HOM (can be UNK)
            //	•	the affected mother must be HOM (can be UNK)
        //        There is no unaffected person that has a HOM call.
        //        There is no unaffected male person that has a HET call.

        boolean isGirl=false,isBoy=false;
        boolean xr = false;

        if (familyPed.get(0).sex.equalsIgnoreCase("1")) {
            isBoy=true;
        } else if (familyPed.get(0).sex.equalsIgnoreCase("2")) {
            isGirl=true;
        } else {
            return false;
        }
        if (familyPed.size() == 1) {
            strict = false;
        }
        int qtyIndividualToTest =1;
        if (strict) {
            qtyIndividualToTest = familyPed.size();
        }

        for (int i=0; i<qtyIndividualToTest; i++) {
            String zygosity = zygosity(genotypesFamily.get(i));
            boolean isAffected;
            if (strict) {
                isAffected = isAffected(familyPed.get(i).getPhenotype());
            } else {
                isAffected = (i==0 && isAffected(familyPed.get(i).getPhenotype()));
            }
            if (!strict && i==0) {
                if (isAffected) {

                    if (isGirl) {
                        if (zygosity.equalsIgnoreCase("HOM")) {
                            return true;
                        } else {
                            return false;
                        }
                    } else { //isBoy
                        if (zygosity.equalsIgnoreCase("HOM") || zygosity.equalsIgnoreCase("HET")) {
                            return true;
                        } else {
                            return false;
                        }
                    }
                } else {
                    return false;
                }
            } else {
                // strict mode
                boolean isFth = isBoy(familyPed.get(i).getSex());
                boolean isMth = isGirl(familyPed.get(i).getSex());

                if (i == 0) {
                    if (isGirl) {
                        if (zygosity.equalsIgnoreCase("HOM")) {
                            xr = true;
                        } else {
                            return false;
                        }
                    } else {
                        if (zygosity.equalsIgnoreCase("HOM") || zygosity.equalsIgnoreCase("HET")) {
                            xr = true;
                        } else {
                            return false;
                        }
                    }
                } else { //father or mth
                    if (isGirl) {
                        //For the parents of affected females:
                        //	•	the father must be affected and
                        //	•	the mother cannot have it REF or HOM
                        if (isFth && isAffected) {
                            if (!zygosity.equalsIgnoreCase("HOM REF")) {
                                xr = true;
                            } else {
                                return false;
                            }

                        } else if (isFth && !isAffected) {
                            return false;
                        } else if (isMth && zygosity.startsWith("HOM")) {
                            return false;
                        }
                    } else {
                        // For the parents of affected males:
                        //	•	the unaffected father cannot have the variant HET or HOM
                        //	•	the affected father must be  HET or HOM (can be UNK)
                        //	•	the unaffected mother cannot be HOM (can be UNK)
                        //	•	the affected mother must be HOM (cannot be UNK)
                        if (isFth && isAffected) {
                            if (zygosity.equalsIgnoreCase("HOM")
                                    || zygosity.equalsIgnoreCase("HET")
                                    || zygosity.equalsIgnoreCase("UNK")) {
                                xr = true;
                            } else {
                                return false;
                            }
                         } else if (isFth && !isAffected) {
                            if (zygosity.equalsIgnoreCase("HOM") || zygosity.equalsIgnoreCase("HET")) {
                                return false;
                            } else {
                                xr = true;
                            }
                        } else if (isMth && isAffected) {
                            if (zygosity.equalsIgnoreCase("HOM") || zygosity.equalsIgnoreCase("UNK")) {
                                xr = true;
                            } else {
                                return false;
                            }
                        } else if (isMth && !isAffected) {
                            if (zygosity.equalsIgnoreCase("HOM")) {
                                return false;
                            } else {
                                xr = true;
                            }
                        }
                    }
                }
            }

        }
        return xr;
    }

    public static boolean isXDominant(List<String> genotypesFamily, List<Pedigree> familyPed, boolean strict) {

        // Gemini
            // Affected males are HET or HOM
        // Affected females must be HET
        // Unaffecteds must be HOM_REF
        // girls of affected dad must be affected
        // boys of affected dad must be unaffected
        // mothers of affected males must be het (and affected)
        // at least 1 parent of affected females must be het (and affected)

        // Janovar
//        •	First of all variants must be X-Chromosomal.
//        •	If the pedigree only contains one person then we decide
    //        if * the person is female then the variant call list must contain one HET call.
    //        * else the variant call list must contain a HET or a HOM call.
//        •	If there is more than one person in the pedigree then there must be at least one compatible call, meaning:
//        •	at least one affected male has a HET or HOM call or
//        •	an affected female a HET call for this variant,
//        •	no affected person has a REF call,
//        •	no affected female has a HOM call, and *
//        •	no unaffected person has a HET or HOM call.

        // CHUSJ - no male proband can be XD

        boolean isGirl=false,isBoy=false;
        boolean xd = false;

        if (familyPed.get(0).sex.equalsIgnoreCase("1")) {
            isBoy=true;
            return false;
        } else if (familyPed.get(0).sex.equalsIgnoreCase("2")) {
            isGirl=true;
        } else {
            return false;
        }
        if (familyPed.size() == 1) {
            strict = false;
        }

        for (int i=0; i<familyPed.size(); i++) {

            // at least 1 parent of affected females must be het (and affected)

            String zygosity = zygosity(genotypesFamily.get(i));
            boolean isAffected;
            boolean isFth = isBoy(familyPed.get(i).getSex());
            boolean isMth = isGirl(familyPed.get(i).getSex());
            if (strict) {
                isAffected = isAffected(familyPed.get(i).getPhenotype());

            } else {
                isAffected = (i==0 && isAffected(familyPed.get(i).getPhenotype()));
            }

            if (isAffected) {
                // at least one affected male has a HET or HOM call or a affected female a HET call for this variant
                if (i==0 && !strict && isGirl && zygosity.equalsIgnoreCase("HET")) {
                    return true;
                } else if (i==0 && isGirl && !zygosity.equalsIgnoreCase("HET")) {
                    return false;
                } else if (i==0 && isGirl && zygosity.equalsIgnoreCase("HET")) {
                    xd=true;
                } else if (i==0 && isBoy && !strict && (zygosity.equalsIgnoreCase("HET") || zygosity.equalsIgnoreCase("HOM")) ) {
                    return true;
                } else if (i==0 && isBoy && zygosity.equalsIgnoreCase("HOM REF")) {
                    return false;
                } else if (i==0 && isBoy && (zygosity.equalsIgnoreCase("HET") || zygosity.equalsIgnoreCase("HOM")) ) {
                    xd=true;
                }
                if (strict && i>0) {
                    // girls of affected dad must be affected
                    if (isGirl && isFth) {
                        xd=true;
                    }
                    // at least 1 parent of affected females must be het (and affected)
                    if (isGirl && zygosity.equalsIgnoreCase("HET")) {
                        xd=true;
                    }
                    // mothers of affected males must be het (and affected)
                    if (isMth && isBoy) {
                        if (zygosity.equalsIgnoreCase("HET") ||
                                (zygosity.equalsIgnoreCase("UNK") && zygosity(genotypesFamily.get(0)).equalsIgnoreCase("HET")) ) {
                            xd=true;
                        } else {
                            return false;
                        }
                    }
                    // no affected person has a REF call;
                    if (zygosity.equalsIgnoreCase("HOM REF")) {
                        return false;
                    }
                }

            } else if ( !(zygosity.equalsIgnoreCase("HOM REF") || zygosity.equalsIgnoreCase("UNK")) ) {
                return false;
            // mothers of affected males must be het (and affected)
            } else if (isMth && isBoy ) {
                return false;
            } else {
//                xd = true;
            }
        }
        return xd;
    }

    static boolean isBoy(String sex) {
        if (sex.equalsIgnoreCase("1")) {
            return true;
        } else {
            return false;
        }
    }
    static boolean isGirl(String sex) {
        if (sex.equalsIgnoreCase("2")) {
            return true;
        } else {
            return false;
        }
    }
    static boolean isAffected(String status) {
        if (status.equalsIgnoreCase("1")) {
            return false;
        } else {
            return true;
        }
    }
}
